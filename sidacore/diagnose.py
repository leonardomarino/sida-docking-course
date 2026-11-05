"""Step 4 diagnosis: why does a redock fail? Single-factor changes in the order the spec prescribes.

Order: ligand protonation, box centre and size, FAD treatment, receptor protonation. Every variant is run on every
benchmark receptor under the SAME registered pass rule, and every result is reported - including the ones that do not help
and any that break a receptor that worked. A variant that "passes" is a finding, not a fix: it is adopted only with a
physical justification, never because it passed.
"""
from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path

import numpy as np
from rdkit import Chem

from .boxes import inside_box, read_vina_box, write_vina_box
from .config import load_config, write_resolved
from .prep import fetch_structure, prepare_receptor
from .redock import (MIN_SEEDS, _topology, benchmark_jobs, crystal_reference, evaluate_seeds, ligand_from_smiles,
                     ligand_template, prepare_ligand, read_poses, redock_seeds)

ENGINE = "vina"      # the two engines agreed to within 0.05 A in the baseline, so the ladder uses the fast one
FLAVIN, NADP = ("FAD", "FDA"), ("NAP", "NDP")


def _names(mol):
    return [a.GetPDBResidueInfo().GetName().strip() for a in mol.GetAtoms()]


def atom_displacements(pose, ref) -> list:
    """[(atom name, distance in A)] between a pose and the deposited pose, atoms paired by the element-graph automorphism
    that minimises the total displacement (so chemically equivalent atoms, e.g. the carboxylate oxygens, are not penalised).
    Names are the reference's."""
    tp, tr = _topology(pose), _topology(ref)
    matches = tp.GetSubstructMatches(tr, uniquify=False, maxMatches=10000)     # match[ref_idx] = pose_idx
    if not matches:
        raise ValueError("pose does not match the reference topology")
    p, r = pose.GetConformer().GetPositions(), ref.GetConformer().GetPositions()
    best = min(matches, key=lambda m: float(np.sum((p[list(m)] - r) ** 2)))
    d = np.linalg.norm(p[list(best)] - r, axis=1)
    return list(zip(_names(ref), map(float, d)))


# --------------------------------------------------------------------------- rung 1: ligand protonation states
# Written with the same atom order and chirality tag as the baseline so the stereocentre is identical in every state.
LIGAND_STATES = {
    "baseline +1 (pH 7.4)": "[NH3+]CCC[C@H]([NH3+])C(=O)[O-]",
    "delta-amine neutral (net 0)": "NCCC[C@H]([NH3+])C(=O)[O-]",
    "alpha-amine neutral (net 0)": "[NH3+]CCC[C@H](N)C(=O)[O-]",
    "both amines + acid (+2)": "[NH3+]CCC[C@H]([NH3+])C(=O)O",
    "anion (-1)": "NCCC[C@H](N)C(=O)[O-]",
}


# --------------------------------------------------------------------------- rung 3: cofactor treatment (diagnostic edits)
def edit_pdbqt(text: str, drop=(), zero_charge=()) -> str:
    """Receptor PDBQT with whole residues removed and/or the partial charges of residues set to 0. Nothing else changes."""
    out = []
    for line in text.splitlines(True):
        if line.startswith(("ATOM", "HETATM")):
            res = line[17:20].strip()
            if res in drop:
                continue
            if res in zero_charge:
                line = line[:70] + "%+6.3f" % 0.0 + line[76:]
        out.append(line)
    return "".join(out)


# --------------------------------------------------------------------------- rung 2: boxes
def box_variants(base_center, base_size, carried_center, ligand_centroid) -> list:
    """[(label, centre, size)]. One thing changes at a time. The oracle box is centred on the deposited ligand: it is a
    diagnostic upper bound on what a perfectly placed box could do, never a usable protocol."""
    return [("size 18 A", tuple(base_center), (18, 18, 18)), ("size 26 A", tuple(base_center), (26, 26, 26)),
            ("size 30 A", tuple(base_center), (30, 30, 30)),
            ("centre: carried from 6X0I", tuple(carried_center), tuple(base_size)),
            ("centre: deposited-ligand centroid (oracle, diagnostic only)", tuple(ligand_centroid), tuple(base_size))]


def crystal_fraction_inside(ref, center, size) -> float:
    xyz = ref.GetConformer().GetPositions()[[a.GetAtomicNum() > 1 for a in ref.GetAtoms()]]
    return float(inside_box(xyz, center, size).mean())


# --------------------------------------------------------------------------- rung 4: receptor protonation
def receptor_variants() -> list:
    """[(label, pH)]. The ligand stays at the baseline pH; only the receptor's PDB2PQR/PROPKA pH moves."""
    return [("receptor protonation pH 6.0", 6.0), ("receptor protonation pH 8.5", 8.5)]


def _site_states(csv_path: Path) -> dict:
    return {(r["res_name"], r["res_id"]): r["state"] for r in csv.DictReader(open(csv_path))}


def _slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


# --------------------------------------------------------------------------- the ladder
def _steps(cfg, include_baseline: bool, only) -> list:
    steps = [("0 baseline", "baseline", {})] if include_baseline else []
    for label, smi in LIGAND_STATES.items():
        if not label.startswith("baseline"):
            steps.append(("1 ligand protonation", label, {"smiles": smi}))
    for label in ("size 18 A", "size 26 A", "size 30 A", "centre: carried from 6X0I",
                  "centre: deposited-ligand centroid (oracle, diagnostic only)"):
        steps.append(("2 box", label, {"box": label}))
    steps += [("3 FAD treatment", "FAD removed (diagnostic)", {"drop": FLAVIN}),
              ("3 FAD treatment", "FAD charges zeroed (diagnostic)", {"zero": FLAVIN}),
              ("3 FAD treatment", "NADP removed (diagnostic)", {"drop": NADP}),
              ("3 FAD treatment", "FAD and NADP removed (diagnostic)", {"drop": FLAVIN + NADP})]
    steps += [("4 receptor protonation", label, {"ph": ph}) for label, ph in receptor_variants()]
    return [x for x in steps if not only or x[1] in only]


def run_ladder(cfg, prepared_root, boxes_root, out_root, jobs=None, only=None, with_baseline=True, cache_dir="structures_cache") -> dict:
    if len(cfg.seeds) < MIN_SEEDS:
        raise ValueError("Step 4 needs at least %d seeds, config has %d" % (MIN_SEEDS, len(cfg.seeds)))
    prepared_root, boxes_root, out_root = Path(prepared_root), Path(boxes_root), Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    write_resolved(cfg, out_root)
    template = ligand_template(cfg.ligand_resname)
    base_lp = prepare_ligand(cfg.ligand_resname, cfg.protonation_ph)
    size = tuple(cfg.box_size_angstrom)
    results = []
    for job in [j for j in benchmark_jobs(cfg, prepared_root) if not jobs or j.name in jobs]:
        ref = crystal_reference(job / "ligand_crystal.pdb", template)
        centre, _ = read_vina_box(job / "receptor.box.txt")
        carried, _ = read_vina_box(boxes_root / job.name / "box.carried_from_reference.txt")
        lig_centroid = ref.GetConformer().GetPositions().mean(axis=0)
        boxes = {l: (c, s) for l, c, s in box_variants(centre, size, carried, lig_centroid)}
        rep = json.loads((job / "prep_report.json").read_text())
        base_pdbqt = (job / "receptor.pdbqt").read_text()
        for rung, label, spec in _steps(cfg, with_baseline, only):
            t0 = time.time()
            d = out_root / _slug(rung) / _slug(label) / job.name
            d.mkdir(parents=True, exist_ok=True)
            lp, receptor, box_file, extras = base_lp, job / "receptor.pdbqt", job / "receptor.box.txt", {}
            if "smiles" in spec:
                lp = ligand_from_smiles(spec["smiles"])
                extras["ligand"] = {"smiles": lp.smiles, "net_charge": lp.net_charge, "n_atoms": lp.n_atoms}
            if "box" in spec:
                c, s = boxes[spec["box"]]
                box_file = d / "box.txt"
                write_vina_box(box_file, c, s)
                extras["crystal_fraction_inside"] = crystal_fraction_inside(ref, c, s)
                extras["box_centre"], extras["box_size"] = [round(float(x), 3) for x in c], list(s)
            if "drop" in spec or "zero" in spec:
                receptor = d / "receptor.pdbqt"
                receptor.write_text(edit_pdbqt(base_pdbqt, drop=spec.get("drop", ()), zero_charge=spec.get("zero", ())))
            if "ph" in spec:
                cfg_ph = load_config(cfg.source_path, overrides={**dict(cfg.overrides), "protonation_ph": spec["ph"]})
                prep_dir = d / "receptor_prep"
                if not prep_dir.exists():
                    prepare_receptor(fetch_structure(rep["pdb_id"], cache_dir), rep["pdb_id"], rep["chain"], cfg_ph,
                                     rep["nadp_policy"], prep_dir)
                receptor, box_file = prep_dir / "receptor.pdbqt", prep_dir / "receptor.box.txt"
                before, after = _site_states(job / "protonation_states.csv"), _site_states(prep_dir / "protonation_states.csv")
                extras["site_states_changed"] = {"%s%s" % k: "%s -> %s" % (before[k], after[k]) for k in before if k in after and before[k] != after[k]}
            rows = redock_seeds(job, lp, ref, cfg, ENGINE, d, receptor=receptor, box_file=box_file)
            ev = evaluate_seeds(rows, cfg.redock_rmsd_acceptance, min_seeds=MIN_SEEDS)
            results.append({"rung": rung, "label": label, "job": job.name, "engine": ENGINE, "extras": extras,
                            "wall_seconds": round(time.time() - t0, 1), "per_seed": rows, **ev})
    return {"settings": {"engine": ENGINE, "exhaustiveness": cfg.exhaustiveness, "num_modes": cfg.num_modes,
                         "seeds": list(cfg.seeds), "acceptance": cfg.redock_rmsd_acceptance, "protonation_ph": cfg.protonation_ph},
            "variants": results}


def baseline_displacements(cfg, prepared_root, redock_root, seed: int = 1) -> dict:
    """{job: {engine: {"top1_rmsd":..., "atoms": [(name, A)...]}}} for the baseline redock's top-1 pose of one seed."""
    template = ligand_template(cfg.ligand_resname)
    out = {}
    for job in benchmark_jobs(cfg, Path(prepared_root)):
        ref = crystal_reference(job / "ligand_crystal.pdb", template)
        for eng in ("vina", "unidock"):
            f = Path(redock_root) / job.name / eng / ("%s_seed%d.pdbqt" % (eng, seed))
            if f.exists():
                pose = read_poses(f.read_text())[0][1]
                out.setdefault(job.name, {})[eng] = {"atoms": atom_displacements(pose, ref)}
    return out


# --------------------------------------------------------------------------- report
def render_markdown(report: dict, displacements: dict | None = None, notes: str = "") -> str:
    s = report["settings"]
    jobs = sorted({v["job"] for v in report["variants"]})
    short = {j: j.replace("_nadp-keep", "") for j in jobs}
    L = ["# Step 4 - diagnosis of the redocking failures", "",
         "Generated by `python -m sidacore.diagnose`; numbers only. Interpretation is in the notes, if any.", "",
         "- Order, as the spec prescribes: ligand protonation, box centre and size, FAD treatment, receptor protonation. One factor "
         "changes at a time against the baseline; **every** variant is run on **every** benchmark receptor, including those that "
         "pass, so a change that rescues a failure but breaks a pass is visible.",
         "- Engine: %s, exhaustiveness %d, %d modes, seeds %s. Same registered rule as the baseline: **median top-1 RMSD over seeds "
         "<= %g A**. A variant that passes is a finding, not a fix: it is adopted only with a physical justification."
         % (s["engine"], s["exhaustiveness"], s["num_modes"], ", ".join(map(str, s["seeds"])), s["acceptance"]),
         "- Cell: top-1 median RMSD in A (seeds within the acceptance / seeds) and, after a slash, the best-in-top-5 median.", ""]
    for rung in sorted({v["rung"] for v in report["variants"]}):
        L += ["## Rung %s" % rung, "", "| variant | " + " | ".join(short[j] for j in jobs) + " | notes |",
              "|---|" + "---|" * (len(jobs) + 1)]
        for label in dict.fromkeys(v["label"] for v in report["variants"] if v["rung"] == rung):
            cells, notes_l = [], []
            for j in jobs:
                r = next((v for v in report["variants"] if v["rung"] == rung and v["label"] == label and v["job"] == j), None)
                if r is None:
                    cells.append("-")
                    continue
                cells.append("%.2f (%d/%d) / %.2f%s" % (r["top1_median"], r["top1_n_within"], r["n_seeds"], r["best_top5_median"],
                                                      " **PASS**" if r["passes"] else ""))
                e = r["extras"]
                if "crystal_fraction_inside" in e:
                    notes_l.append("%s: %.0f%% of ligand in box" % (short[j], 100 * e["crystal_fraction_inside"]))
                if e.get("site_states_changed"):
                    notes_l.append("%s: %s" % (short[j], ", ".join("%s %s" % kv for kv in e["site_states_changed"].items())))
                if "ligand" in e and j == jobs[0]:
                    notes_l.append("net %+d, %d atoms" % (e["ligand"]["net_charge"], e["ligand"]["n_atoms"]))
            L.append("| %s | %s | %s |" % (label, " | ".join(cells), "; ".join(notes_l) or "-"))
        L.append("")
    if displacements:
        L += ["## Where the baseline top-1 pose differs from the deposited one (seed 1, A per atom)", "",
              "| receptor | engine | " + " | ".join(n for n, _ in next(iter(next(iter(displacements.values())).values()))["atoms"]) + " |",
              "|---|---|" + "---|" * 9]
        for job, per in displacements.items():
            for eng, v in per.items():
                L.append("| %s | %s | %s |" % (short.get(job, job), eng, " | ".join("%.1f" % d for _, d in v["atoms"])))
        L.append("")
    if notes.strip():
        L += ["## Notes", "", notes.strip()]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    import argparse
    import sys
    import yaml
    ap = argparse.ArgumentParser(description="Step 4 diagnosis ladder (single-factor changes, all benchmark receptors).")
    ap.add_argument("--config", required=True)
    ap.add_argument("--prepared", required=True)
    ap.add_argument("--boxes", required=True)
    ap.add_argument("--redock", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default="structures_cache")
    ap.add_argument("--report")
    ap.add_argument("--notes")
    ap.add_argument("--render-only", action="store_true", help="rebuild --report from OUT/step4_diagnosis.json; no docking")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    out = Path(a.out)
    if a.render_only:
        data = json.loads((out / "step4_diagnosis.json").read_text())
        report, disp = data["ladder"], data["displacements"]
    else:
        report = run_ladder(cfg, a.prepared, a.boxes, out, cache_dir=a.cache)
        disp = baseline_displacements(cfg, a.prepared, a.redock)
        (out / "step4_diagnosis.json").write_text(json.dumps({"ladder": report, "displacements": disp}, indent=1))
    if a.report:
        notes = Path(a.notes).read_text() if a.notes and Path(a.notes).exists() else ""
        Path(a.report).write_text(render_markdown(report, disp, notes))
    for v in report["variants"]:
        print("%-24s %-58s %-18s top-1 %.2f (%d/%d) best5 %.2f %s" % (v["rung"], v["label"], v["job"], v["top1_median"],
              v["top1_n_within"], v["n_seeds"], v["best_top5_median"], "PASS" if v["passes"] else "fail"))
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
