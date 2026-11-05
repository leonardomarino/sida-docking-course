"""Step 4: redocking validation. Thin glue over installed tools.

    biotite (CCD)   the ligand's identity: topology, stereo and charges come from the PDB component dictionary
    obabel          protonation of that molecule at the configured pH
    RDKit           3D embedding of the protonated molecule, crystal-pose reference, symmetry-aware RMSD
    Meeko           ligand -> PDBQT; docked PDBQT -> RDKit
    Vina (CPU), Uni-Dock (GPU)   the two engines, same receptor, box, seeds and settings

Nothing is tuned here: a failed redock is reported with every change made and its effect, never silently adjusted.
"""
from __future__ import annotations

import re
import shutil
import statistics
from dataclasses import dataclass
from pathlib import Path

import biotite.structure.info as ccd
import numpy as np
from biotite.interface import rdkit as bt_rdkit
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from .boxes import _tool, read_vina_box
from .prep import _run

UNIDOCK_DEFAULT = Path("/opt/miniforge3/envs/sida-docking-gpu/bin/unidock")
EMBED_SEED = 0          # 3D start conformer only; Vina/Uni-Dock randomise the starting pose of every run
MIN_SEEDS = 10          # Step 4 asks for at least 10


# --------------------------------------------------------------------------- the ligand
def ligand_template(resname: str):
    """Heavy-atom molecule (bond orders, stereo) of a PDB component, from biotite's bundled dictionary."""
    try:
        atoms = ccd.residue(resname)
    except Exception as e:                                       # noqa: BLE001 - biotite raises KeyError/others
        raise ValueError("cannot build %r from the PDB component dictionary: %s" % (resname, e)) from e
    mol = bt_rdkit.to_mol(atoms, explicit_hydrogen=True)
    Chem.AssignStereochemistryFrom3D(mol)
    return Chem.RemoveHs(mol)


def protonate(template, ph: float) -> str:
    """Protonation state at `ph` as SMILES (Open Babel). Charges and stereo are returned in the SMILES itself."""
    smi = Chem.MolToSmiles(template)
    import subprocess
    p = subprocess.run([str(_tool("obabel")), "-:" + smi, "-p", "%g" % ph, "-osmi"], capture_output=True, text=True)
    out = p.stdout.split()
    if p.returncode != 0 or not out:
        raise RuntimeError("obabel could not protonate %s at pH %g: %s" % (smi, ph, p.stderr.strip()[-200:]))
    return out[0]


@dataclass(frozen=True)
class LigandPrep:
    pdbqt: str
    smiles: str               # the protonated molecule that was docked
    ph: float | None          # the pH it was protonated at; None when the state was given explicitly
    net_charge: int           # formal charge of that molecule
    n_atoms: int              # ATOM lines in the PDBQT (heavy + polar hydrogens)
    n_heavy: int
    torsdof: int
    partial_charge_sum: float


def ligand_from_smiles(smiles: str, ph: float | None = None, embed_seed: int = EMBED_SEED) -> LigandPrep:
    """A protonated SMILES -> 3D (RDKit; charges survive) -> Meeko PDBQT."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=embed_seed) != 0:
        raise RuntimeError("RDKit could not embed %s in 3D" % smiles)
    AllChem.MMFFOptimizeMolecule(mol)
    setups = MoleculePreparation().prepare(mol)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError("Meeko could not write %s: %s" % (smiles, err))
    atoms = [l for l in pdbqt.splitlines() if l.startswith(("ATOM", "HETATM"))]
    return LigandPrep(
        pdbqt=pdbqt, smiles=smiles, ph=ph, net_charge=Chem.GetFormalCharge(mol), n_atoms=len(atoms),
        n_heavy=sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1),
        torsdof=int(re.search(r"TORSDOF\s+(\d+)", pdbqt).group(1)),
        partial_charge_sum=float(sum(float(l[70:76]) for l in atoms)))


def prepare_ligand(resname: str, ph: float, embed_seed: int = EMBED_SEED) -> LigandPrep:
    """Protonate at `ph`, embed in 3D with RDKit (`obabel -p ... --gen3d` silently dropped the charges), Meeko."""
    return ligand_from_smiles(protonate(ligand_template(resname), ph), ph, embed_seed)


# --------------------------------------------------------------------------- reference pose and RMSD
def crystal_reference(ligand_pdb, template):
    """The deposited ligand as an RDKit molecule: its coordinates, the dictionary's bond orders. Refuses a mismatch."""
    mol = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False, sanitize=False, proximityBonding=True)
    if mol is None:
        raise ValueError("cannot read %s" % ligand_pdb)
    mol = Chem.RemoveHs(mol, sanitize=False)
    try:
        ref = AllChem.AssignBondOrdersFromTemplate(template, mol)
    except Exception as e:                                       # noqa: BLE001
        raise ValueError("%s does not match the template ligand (%d heavy atoms vs %d): %s" % (
            ligand_pdb, mol.GetNumAtoms(), template.GetNumAtoms(), e)) from e
    Chem.AssignStereochemistryFrom3D(ref)
    return ref


def _topology(mol):
    """Heavy-atom graph of elements only (single bonds, no charges): the symmetry group is then every element-preserving
    automorphism, so exchanging the two carboxylate oxygens is free whichever protonation form each side is drawn in."""
    rw = Chem.RWMol(Chem.RemoveHs(mol, sanitize=False))
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0)
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
    rw.UpdatePropertyCache(strict=False)
    return rw.GetMol()


def pose_rmsd(pose, ref) -> float:
    """Symmetry-corrected heavy-atom RMSD, NO re-alignment (both poses are in the receptor frame)."""
    return float(rdMolAlign.CalcRMS(_topology(pose), _topology(ref)))


def read_poses(pdbqt_text: str) -> list:
    """[(score kcal/mol, RDKit heavy-atom molecule with that pose's coordinates)], best score first."""
    scores = [float(x) for x in re.findall(r"REMARK \S+ RESULT:\s+(-?\d+\.?\d*)", pdbqt_text)]
    mols = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(pdbqt_text, skip_typing=True))
    if not mols or mols[0] is None:
        raise ValueError("no poses could be read from the docking output")
    multi = Chem.RemoveHs(mols[0], sanitize=False)
    out = []
    for k, conf in enumerate(multi.GetConformers()):
        m = Chem.Mol(multi)
        m.RemoveAllConformers()
        m.AddConformer(Chem.Conformer(conf), assignId=True)
        out.append((scores[k], m))
    if len(out) != len(scores):
        raise ValueError("%d scores but %d poses in the docking output" % (len(scores), len(out)))
    return sorted(out, key=lambda t: t[0])


def pose_pdb_block(pose, ref, chain: str = "L") -> str:
    """PDB text of `pose` whose atoms carry the names and residue name of the deposited ligand `ref`, so a program that
    wants real residue names (Rosetta) can read it. The atom correspondence is a graph automorphism of the element
    topology; where several exist (the two carboxylate oxygens) they are chemically equivalent, so any one is correct."""
    tp, tr = _topology(pose), _topology(ref)
    match = tp.GetSubstructMatch(tr) if tp.GetNumAtoms() == tr.GetNumAtoms() else ()
    if len(match) != tr.GetNumAtoms():
        raise ValueError("pose (%d atoms) does not match the reference topology (%d atoms)" % (tp.GetNumAtoms(), tr.GetNumAtoms()))
    out = Chem.RWMol(Chem.RemoveHs(pose, sanitize=False))
    for ref_idx, pose_idx in enumerate(match):                      # match[i] = pose atom paired with reference atom i
        info = ref.GetAtomWithIdx(ref_idx).GetPDBResidueInfo()
        new = Chem.AtomPDBResidueInfo()
        new.SetName(info.GetName())
        new.SetResidueName(info.GetResidueName())
        new.SetResidueNumber(1)
        new.SetChainId(chain)
        new.SetIsHeteroAtom(True)
        out.GetAtomWithIdx(pose_idx).SetMonomerInfo(new)
    block = Chem.MolToPDBBlock(out.GetMol())
    return "".join(l for l in block.splitlines(True) if not l.startswith("CONECT"))      # their serials would collide with a receptor file


# --------------------------------------------------------------------------- engines
def _write_ligand(out_dir: Path, pdbqt: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    f = out_dir / "ligand.pdbqt"
    if not f.exists() or f.read_text() != pdbqt:
        f.write_text(pdbqt)
    return f


def run_vina(receptor, ligand_pdbqt: str, box_file, *, seed: int, exhaustiveness: int, num_modes: int, out_dir) -> Path:
    out_dir = Path(out_dir)
    lig = _write_ligand(out_dir, ligand_pdbqt)
    out = out_dir / ("vina_seed%d.pdbqt" % seed)
    _run([_tool("vina"), "--receptor", receptor, "--ligand", lig, "--config", box_file, "--seed", seed,
          "--exhaustiveness", exhaustiveness, "--num_modes", num_modes, "--out", out], out_dir / ("vina_seed%d.log" % seed))
    return out


def run_unidock(receptor, ligand_pdbqt: str, box_file, *, seed: int, exhaustiveness: int, num_modes: int, out_dir,
                binary=UNIDOCK_DEFAULT) -> Path:
    out_dir = Path(out_dir)
    lig = _write_ligand(out_dir, ligand_pdbqt)
    centre, size = read_vina_box(box_file)
    work = out_dir / ("unidock_seed%d" % seed)
    work.mkdir(exist_ok=True)
    cmd = [binary, "--receptor", receptor, "--gpu_batch", lig, "--dir", work, "--seed", seed,
           "--exhaustiveness", exhaustiveness, "--num_modes", num_modes, "--scoring", "vina"]
    for ax, c, s in zip("xyz", centre, size):
        cmd += ["--center_" + ax, "%.3f" % c, "--size_" + ax, "%.3f" % s]
    _run(cmd, out_dir / ("unidock_seed%d.log" % seed))
    produced = sorted(work.glob("*.pdbqt"))
    if len(produced) != 1:
        raise RuntimeError("Uni-Dock wrote %d pdbqt files in %s (expected 1)" % (len(produced), work))
    out = out_dir / ("unidock_seed%d.pdbqt" % seed)
    shutil.move(str(produced[0]), out)
    return out


# --------------------------------------------------------------------------- evaluation
def evaluate_seeds(per_seed: list, accept: float, min_seeds: int | None = None) -> dict:
    """Registered rule (written before any redock was run): PASS if the median over seeds of the top-1 RMSD is
    <= `accept`. The number of seeds within `accept`, the best-in-top-5 and the worst seed are reported beside it."""
    if min_seeds is not None and len(per_seed) < min_seeds:
        raise ValueError("Step 4 needs at least %d seeds, got %d" % (min_seeds, len(per_seed)))
    npos = [s["n_poses"] for s in per_seed if "n_poses" in s]
    t1 = [s["top1_rmsd"] for s in per_seed]
    t5 = [s["best_top5_rmsd"] for s in per_seed]
    return {"n_seeds": len(per_seed), "accept": accept,
            "top1_median": statistics.median(t1), "top1_min": min(t1), "top1_max": max(t1),
            "top1_n_within": sum(x <= accept for x in t1),
            "best_top5_median": statistics.median(t5), "best_top5_max": max(t5),
            "best_top5_n_within": sum(x <= accept for x in t5),
            "top1_score_median": statistics.median(s["top1_score"] for s in per_seed),
            "n_poses_min": min(npos) if npos else None, "n_poses_max": max(npos) if npos else None,
            "passes": bool(statistics.median(t1) <= accept)}


# --------------------------------------------------------------------------- one receptor, one engine, all seeds
def redock_seeds(job_dir: Path, lp: LigandPrep, ref, cfg, engine: str, out_dir: Path, receptor=None, box_file=None) -> list:
    """Dock the prepared ligand into one prepared receptor once per configured seed; RMSD of every pose to the deposited one."""
    runner = {"vina": run_vina, "unidock": run_unidock}[engine]
    rows = []
    for seed in cfg.seeds:
        out = runner(receptor or job_dir / "receptor.pdbqt", lp.pdbqt, box_file or job_dir / "receptor.box.txt", seed=seed,
                     exhaustiveness=cfg.exhaustiveness, num_modes=cfg.num_modes, out_dir=out_dir)
        poses = read_poses(out.read_text())
        rm = [pose_rmsd(m, ref) for _, m in poses]
        rows.append({"seed": seed, "n_poses": len(poses), "scores": [s for s, _ in poses], "rmsds": rm,
                     "top1_rmsd": rm[0], "top1_score": poses[0][0], "best_top5_rmsd": min(rm[:5]),
                     "best_pose_rank": int(np.argmin(rm)) + 1})
    return rows


def benchmark_jobs(cfg, prepared_root: Path) -> list:
    """Prepared receptors of the benchmark structures that kept a deposited ligand, primary first."""
    import json
    keep = []
    for rep in sorted(Path(prepared_root).glob("*/prep_report.json")):
        r = json.loads(rep.read_text())
        job = rep.parent
        if r["pdb_id"] in (cfg.structures.benchmark_primary, cfg.structures.benchmark_secondary) \
                and (job / "ligand_crystal.pdb").exists() and (job / "ligand_crystal.pdb").stat().st_size > 0 \
                and r["nadp_policy"] == cfg.redock_nadp:
            keep.append((0 if r["pdb_id"] == cfg.structures.benchmark_primary else 1, job.name, job))
    return [j for _, _, j in sorted(keep)]


def _engine_version(engine: str) -> str:
    import subprocess
    exe = _tool("vina") if engine == "vina" else UNIDOCK_DEFAULT
    p = subprocess.run([str(exe), "--version"], capture_output=True, text=True)
    return ((p.stdout + p.stderr).strip().splitlines() or ["unknown"])[0]


def redock_all(cfg, prepared_root, out_root, engines=("vina", "unidock")) -> dict:
    import json
    import time
    from .config import write_resolved
    if len(cfg.seeds) < MIN_SEEDS:                                  # fail before any docking is spent
        raise ValueError("Step 4 needs at least %d seeds, config has %d" % (MIN_SEEDS, len(cfg.seeds)))
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    write_resolved(cfg, out_root)
    template = ligand_template(cfg.ligand_resname)
    lp = prepare_ligand(cfg.ligand_resname, cfg.protonation_ph)
    jobs = benchmark_jobs(cfg, Path(prepared_root))
    if not jobs:
        raise FileNotFoundError("no benchmark receptors with a deposited %s under %s" % (cfg.ligand_resname, prepared_root))
    results = []
    for job in jobs:
        ref = crystal_reference(job / "ligand_crystal.pdb", template)
        for eng in engines:
            t0 = time.time()
            rows = redock_seeds(job, lp, ref, cfg, eng, out_root / job.name / eng)
            ev = evaluate_seeds(rows, cfg.redock_rmsd_acceptance, min_seeds=MIN_SEEDS)
            results.append({"job": job.name, "engine": eng, "wall_seconds": round(time.time() - t0, 1), "per_seed": rows, **ev})
    report = {"ligand": {k: getattr(lp, k) for k in ("smiles", "ph", "net_charge", "n_atoms", "n_heavy", "torsdof", "partial_charge_sum")},
              "settings": {"exhaustiveness": cfg.exhaustiveness, "num_modes": cfg.num_modes, "seeds": list(cfg.seeds),
                           "acceptance": cfg.redock_rmsd_acceptance, "box_center_atom": cfg.box_center_atom,
                           "box_size": list(cfg.box_size_angstrom), "protonation_ph": cfg.protonation_ph,
                           "redock_nadp": cfg.redock_nadp},
              "engines": {e: _engine_version(e) for e in engines}, "results": results}
    (out_root / "step4_report.json").write_text(json.dumps(report, indent=1))
    return report


def _fill_pose_counts(r: dict) -> None:
    """Older reports lack n_poses_min/max; the per-seed rows always carry n_poses."""
    if r.get("n_poses_min") is None:
        n = [x["n_poses"] for x in r.get("per_seed", []) if "n_poses" in x]
        if n:
            r["n_poses_min"], r["n_poses_max"] = min(n), max(n)


def render_markdown(report: dict, notes: str = "") -> str:
    s, lg = report["settings"], report["ligand"]
    L = ["# Step 4 - redocking validation", "",
         "Generated by `python -m sidacore.redock`; numbers only. Interpretation is in the notes, if any.", "",
         "- Ligand docked: `%s` (pH %g, net charge %+d, %d atoms in the PDBQT, %d torsions), from the PDB component dictionary."
         % (lg["smiles"], lg["ph"], lg["net_charge"], lg["n_atoms"], lg["torsdof"]),
         "- Box: centre `%s`, %s A. Exhaustiveness %d, %d modes, NADP in the receptor: `%s`. Engines: %s." % (
             s["box_center_atom"], " x ".join("%g" % x for x in s["box_size"]), s["exhaustiveness"], s["num_modes"],
             s["redock_nadp"], "; ".join("%s = %s" % kv for kv in report["engines"].items())),
         "- Seeds: %s." % ", ".join(str(x) for x in s["seeds"]),
         "- RMSD: symmetry-corrected, heavy atoms, deposited pose as reference, no re-alignment. Registered pass rule "
         "(written before any run): the **median over seeds of the top-1 RMSD is <= %g A**; the count of seeds within %g A "
         "is reported beside it." % (s["acceptance"], s["acceptance"]), "",
         "| receptor | engine | top-1 median | min-max (A) | seeds within | best-in-top-5 median | max | top-1 score median | poses returned | verdict |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in report["results"]:
        _fill_pose_counts(r)
        npos = "-" if r.get("n_poses_min") is None else "%d-%d" % (r["n_poses_min"], r["n_poses_max"])
        L.append("| %s | %s | %.2f | %.2f-%.2f | %d/%d | %.2f | %.2f | %.2f | %s | %s |" % (
            r["job"], r["engine"], r["top1_median"], r["top1_min"], r["top1_max"], r["top1_n_within"], r["n_seeds"],
            r["best_top5_median"], r["best_top5_max"], r["top1_score_median"], npos, "**PASS**" if r["passes"] else "**FAIL**"))
    counts = {}
    for r in report["results"]:
        if r.get("n_poses_min") is not None and r["n_poses_min"] < s["num_modes"]:
            lo, hi = counts.get(r["engine"], (r["n_poses_min"], r["n_poses_max"]))
            counts[r["engine"]] = (min(lo, r["n_poses_min"]), max(hi, r["n_poses_max"]))
    short = [(eng, lo, hi) for eng, (lo, hi) in counts.items()]
    for eng, lo, hi in short:
        L.append("")
        L.append("- **%s returned fewer poses than the %d requested** (%d-%d per seed, depending on receptor and seed). Its "
                 "'best-in-top-5' is then over a shorter list (at 5 poses it is the best of all of them), so best-in-list "
                 "numbers are not comparable with an engine that returned %d." % (eng, s["num_modes"], lo, hi, s["num_modes"]))
    L += ["", "## Per-seed top-1 RMSD (A) / score (kcal/mol)", ""]
    for r in report["results"]:
        L.append("- %s, %s: %s" % (r["job"], r["engine"], "; ".join(
            "s%d %.2f / %.2f" % (x["seed"], x["top1_rmsd"], x["top1_score"]) for x in r["per_seed"])))
    if notes.strip():
        L += ["", "## Notes", "", notes.strip()]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    import argparse
    import sys
    import yaml
    from .config import load_config
    ap = argparse.ArgumentParser(description="Step 4: redock the deposited substrate into the benchmark receptors.")
    ap.add_argument("--config")
    ap.add_argument("--prepared")
    ap.add_argument("--out")
    ap.add_argument("--render-only", metavar="OUT_DIR", help="rebuild --report from OUT_DIR/step4_report.json; no docking")
    ap.add_argument("--engine", nargs="+", default=["vina", "unidock"], choices=["vina", "unidock"])
    ap.add_argument("--report")
    ap.add_argument("--notes")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
    a = ap.parse_args(argv)
    if a.render_only:
        import json
        if not a.report:
            ap.error("--render-only needs --report")
        report = json.loads((Path(a.render_only) / "step4_report.json").read_text())
        notes = Path(a.notes).read_text() if a.notes and Path(a.notes).exists() else ""
        Path(a.report).write_text(render_markdown(report, notes))
        return 0
    if not (a.config and a.prepared and a.out):
        ap.error("--config, --prepared and --out are required unless --render-only is given")
    over = {k: yaml.safe_load(v) for k, _, v in (s.partition("=") for s in a.set)}
    report = redock_all(load_config(a.config, overrides=over), a.prepared, a.out, tuple(a.engine))
    if a.report:
        notes = Path(a.notes).read_text() if a.notes and Path(a.notes).exists() else ""
        Path(a.report).write_text(render_markdown(report, notes))
    for r in report["results"]:
        print("%-22s %-8s top-1 median %.2f (%.2f-%.2f) within %d/%d | best-of-5 median %.2f | %s | %.0fs" % (
            r["job"], r["engine"], r["top1_median"], r["top1_min"], r["top1_max"], r["top1_n_within"], r["n_seeds"],
            r["best_top5_median"], "PASS" if r["passes"] else "FAIL", r["wall_seconds"]))
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
