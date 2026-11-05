"""Step 2 receptor preparation: thin glue over installed tools.

    biotite            read mmCIF, alternate conformers, atom selection, all distances / gaps
    PDB2PQR + PROPKA   protein hydrogens and pKa at the configured pH
    Meeko              receptor -> PDBQT (cofactor templates, polar H, charges) and the Vina box

FAD is KEPT in the receptor: it is a cofactor, not the ligand, and is never co-docked. The substrate
(config.ligand_resname) is removed and written separately for the redocking RMSD.
Every removed atom is logged with a reason; nothing is dropped silently.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import biotite
import biotite.database.rcsb as rcsb
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbqt as pdbqt
import biotite.structure.io.pdbx as pdbx
import numpy as np

from .config import ConfigError, write_resolved

FLAVIN = ("FAD", "FDA")           # FDA = reduced FAD (6X0J); both carry the isoalloxazine C4X
NADP = ("NAP", "NDP")
IONIZABLE = ("ASP", "GLU", "HIS", "LYS", "ARG", "CYS", "TYR")
TEMPLATES = Path(__file__).with_name("data") / "meeko_cofactor_templates.json"   # FAD, FDA, NAP only; offline
DEFAULT_PDB2PQR = Path("/opt/miniforge3/envs/af2/bin/pdb2pqr30")
FAD_STATEMENT = ("FAD is retained in the receptor: it is a cofactor, not the ligand, and is never co-docked. "
                 "Meeko builds its template from the PDB component dictionary; its protonation state (charge -2) "
                 "does not change with protonation_ph, only the protein's does.")
HIS_NOTE = "PROPKA absolute histidine pKa is no better than assuming 6.5 (Hogues et al. 2025): treat as uncertain"


class SiteAtomError(ValueError):
    """The configured box-centre atom is missing, ambiguous, or is not the intended atom."""


# --------------------------------------------------------------------------- structure in
def fetch_structure(pdb_id: str, cache_dir) -> Path:
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / ("%s.cif" % pdb_id.upper())
    if target.exists():
        return target
    got = Path(rcsb.fetch(pdb_id, "cif", target_path=str(cache_dir)))
    if got != target:
        got.rename(target)
    return target


def read_structure(path, alt_conformers: str):
    """mmCIF -> (AtomArray with bonds, altloc removal log). Alternate conformers resolved by biotite."""
    atoms = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1, altloc="all", extra_fields=["occupancy"])
    if alt_conformers == "highest-occupancy":
        keep = struc.filter_highest_occupancy_altloc(atoms, atoms.altloc_id, atoms.occupancy)
    else:
        keep = struc.filter_first_altloc(atoms, atoms.altloc_id)
    log = [_row(atoms, i, "alternate conformer not selected (policy: %s)" % alt_conformers, altloc=str(atoms.altloc_id[i]))
           for i in np.nonzero(~keep)[0]]
    atoms = atoms[keep]
    atoms.bonds = struc.connect_via_residue_names(atoms)
    return atoms, log


def _row(atoms, i, reason, altloc="", dist=None):
    r = {"chain_id": str(atoms.chain_id[i]), "res_id": int(atoms.res_id[i]), "res_name": str(atoms.res_name[i]),
         "atom_name": str(atoms.atom_name[i]), "altloc": altloc,
         "occupancy": float(atoms.occupancy[i]) if "occupancy" in atoms.get_annotation_categories() else "",
         "reason": reason}
    if dist is not None:
        r["dist_to_site"] = round(float(dist), 2)
    return r


def resolve_site_center(atoms, chain: str, spec: str) -> np.ndarray:
    """'FAD:C4X' -> coordinates. 'FAD' also matches the reduced FDA. Refuses a flavin atom that is not bonded to N5."""
    res, _, name = spec.partition(":")
    names = FLAVIN if res in FLAVIN else NADP if res in NADP else (res,)
    idx = np.nonzero((atoms.chain_id == chain) & np.isin(atoms.res_name, names) & (atoms.atom_name == name))[0]
    if len(idx) != 1:
        raise SiteAtomError("box_center_atom %r: found %d matching atoms in chain %s (need exactly 1)" % (spec, len(idx), chain))
    if res in FLAVIN:
        nbrs = {str(atoms.atom_name[j]) for j in atoms.bonds.get_bonds(int(idx[0]))[0]}
        if "N5" not in nbrs:
            raise SiteAtomError("%r is not the isoalloxazine C4a: it is not bonded to N5 (bonded to %s). The hydroxylation "
                                "site atom is C4X; C4A is the adenine carbon, ~18 A away." % (spec, sorted(nbrs)))
    return atoms.coord[idx[0]].astype(float)


# --------------------------------------------------------------------------- selection
@dataclass(frozen=True)
class Selection:
    kept: struc.AtomArray       # protein + FAD (+NADP, waters, ions per policy)
    ligand: struc.AtomArray     # crystal substrate, set aside (None when absent)
    log: list                   # one row per removed atom, with a reason


def select_atoms(atoms, *, chain, nadp, cfg, center) -> Selection:
    n = atoms.array_length()
    dist = np.linalg.norm(atoms.coord - center, axis=1)
    in_chain = atoms.chain_id == chain
    is_prot = in_chain & struc.filter_canonical_amino_acids(atoms)
    is_flav = in_chain & np.isin(atoms.res_name, FLAVIN)
    is_nadp = in_chain & np.isin(atoms.res_name, NADP)
    is_lig = in_chain & (atoms.res_name == cfg.ligand_resname)
    is_wat = in_chain & struc.filter_solvent(atoms)
    is_ion = in_chain & struc.filter_monoatomic_ions(atoms)
    wat_ok = {"strip": np.zeros(n, bool), "keep": np.ones(n, bool), "keep-within-pocket": dist <= cfg.pocket_radius_angstrom}[cfg.waters]
    keep = (is_prot & (atoms.element != "H")) | is_flav | (is_nadp & (nadp == "keep")) | (is_wat & wat_ok) | (is_ion & (cfg.ions == "keep"))

    reason = np.full(n, "", dtype=object)
    reason[~in_chain] = "other chain"
    for mask, text in ((is_lig, None), (is_nadp & (nadp == "strip"), "NADP stripped (policy)"),
                       (is_wat & ~wat_ok, "water stripped (policy)"), (is_ion & (cfg.ions == "strip"), "ion stripped (policy)"),
                       (is_prot & (atoms.element == "H"), "input hydrogen (re-added by PDB2PQR)")):
        sel = mask & (reason == "") & ~keep
        if text:
            reason[sel] = text
    rest = in_chain & ~keep & ~is_lig & (reason == "")
    reason[rest & atoms.hetero] = "crystallization additive / non-cofactor heterogen"
    reason[rest & ~atoms.hetero] = "non-canonical residue (not one of the 20 standard amino acids)"
    log = [_row(atoms, i, reason[i], dist=dist[i]) for i in np.nonzero(reason != "")[0]]
    ligand = atoms[is_lig] if is_lig.any() else None
    return Selection(kept=atoms[keep], ligand=ligand, log=log)


# --------------------------------------------------------------------------- tools
def _resolve_pdb2pqr(explicit=None) -> Path:
    for cand in (explicit, os.environ.get("SIDACORE_PDB2PQR"), shutil.which("pdb2pqr30"), DEFAULT_PDB2PQR):
        if cand and Path(cand).exists():
            return Path(cand)
    raise FileNotFoundError("pdb2pqr30 not found (tried argument, $SIDACORE_PDB2PQR, PATH, %s)" % DEFAULT_PDB2PQR)


def _meeko_bin() -> Path:
    for cand in (Path(sys.executable).parent / "mk_prepare_receptor.py", shutil.which("mk_prepare_receptor.py")):
        if cand and Path(cand).exists():
            return Path(cand)
    raise FileNotFoundError("mk_prepare_receptor.py not found next to %s or on PATH" % sys.executable)


def _run(cmd, log_path: Path) -> str:
    cmd = [str(c) for c in cmd]
    p = subprocess.run(cmd, capture_output=True, text=True)
    text = "$ %s\n%s%s" % (" ".join(cmd), p.stdout, p.stderr)
    log_path.write_text(text)
    if p.returncode != 0:
        raise RuntimeError("%s failed (rc=%d), see %s:\n%s" % (Path(cmd[0]).name, p.returncode, log_path, "\n".join(text.splitlines()[-12:])))
    return text


def _version(cmd) -> str:
    p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    return ((p.stdout + p.stderr).strip().splitlines() or ["unknown"])[-1]


# --------------------------------------------------------------------------- parsing tool output
def parse_rebuilt_atoms(log: str) -> dict:
    # the same message appears with different prefixes on stdout and in PDB2PQR's own log: de-duplicate
    rows = sorted(set(re.findall(r"Added atom (\S+) to residue (\S+) (\S) +(\d+) at coordinates", log)))
    per = {}
    for atom, rname, chain, rid in rows:
        per.setdefault((rname, chain, int(rid)), []).append(atom)
    return {"n_atoms": len(rows), "n_residues": len(per),
            "residues": [{"res_name": k[0], "chain_id": k[1], "res_id": k[2], "atoms": v}
                         for k, v in sorted(per.items(), key=lambda kv: kv[0][2])]}


def parse_propka(log: str) -> dict:
    """{(res_name, res_id, chain): (pKa, model_pKa)} from the SUMMARY OF THIS PREDICTION block.
    The block also lists N+/C- termini and ends at the next '----' separator or 'Free energy' line."""
    out, on = {}, False
    for line in log.splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            on = True
            continue
        if not on:
            continue
        if line.startswith("-") or line.strip().startswith("Free energy"):
            break
        m = re.match(r"^\s+([A-Z]{3})\s+(\d+)\s+(\S)\s+([-\d.]+)\s+([-\d.]+)", line)
        if m:
            out[(m.group(1), int(m.group(2)), m.group(3))] = (float(m.group(4)), float(m.group(5)))
    return out


def ionization_state(res_name: str, hydrogens: set):
    """(state, net_charge) from which hydrogens PDB2PQR placed."""
    if res_name == "HIS":
        d, e = "HD1" in hydrogens, "HE2" in hydrogens
        return ("HIP", 1) if d and e else ("HID", 0) if d else ("HIE", 0) if e else ("HIS?", 0)
    if res_name == "ASP":
        return ("ASH", 0) if "HD2" in hydrogens else ("ASP", -1)
    if res_name == "GLU":
        return ("GLH", 0) if "HE2" in hydrogens else ("GLU", -1)
    if res_name == "LYS":
        return ("LYS", 1) if len({"HZ1", "HZ2", "HZ3"} & hydrogens) == 3 else ("LYN", 0)
    if res_name == "CYS":
        return ("CYS", 0) if "HG" in hydrogens else ("CYM", -1)
    if res_name == "TYR":
        return ("TYR", 0) if "HH" in hydrogens else ("TYM", -1)
    return ("ARG", 1)


def site_residue_ids(atoms, center, box_size, ligand_xyz=None, contact: float = 4.5) -> dict:
    """The substrate site: {(chain, res_id): {"in_box", "contacts_ligand"}} for every residue of `atoms` (heavy atoms) with an
    atom inside the docking box (an axis-aligned cube about `center`) or within `contact` A of the deposited ligand."""
    xyz = atoms.coord
    in_box = np.all(np.abs(xyz - np.asarray(center, float)) <= np.asarray(box_size, float) / 2, axis=1)
    near = np.zeros(len(xyz), bool)
    if ligand_xyz is not None and len(ligand_xyz):
        hit = struc.CellList(atoms, cell_size=max(contact, 2.0)).get_atoms(np.asarray(ligand_xyz, float), contact)
        near[np.unique(hit[hit >= 0])] = True
    site = {}
    for i in np.nonzero(in_box | near)[0]:
        d = site.setdefault((str(atoms.chain_id[i]), int(atoms.res_id[i])), {"in_box": False, "contacts_ligand": False})
        d["in_box"] = d["in_box"] or bool(in_box[i])
        d["contacts_ligand"] = d["contacts_ligand"] or bool(near[i])
    return site


def protonation_rows(prot_h, site, center, pkas) -> list:
    names = {(str(c), int(r)): str(n) for c, r, n in zip(prot_h.chain_id, prot_h.res_id, prot_h.res_name)}
    keys = sorted((c, r, names[(c, r)]) for (c, r) in site)
    rows = []
    for chain, rid, rname in keys:
        if rname not in IONIZABLE:
            continue
        res = prot_h[(prot_h.chain_id == chain) & (prot_h.res_id == rid)]
        state, q = ionization_state(rname, {str(a) for a in res.atom_name[res.element == "H"]})
        pka, model = pkas.get((rname, rid, chain), ("", ""))
        rows.append({"chain_id": chain, "res_id": rid, "res_name": rname, "state": state, "net_charge": q, "pKa": pka,
                     "model_pKa": model, "dist_to_site": round(float(np.linalg.norm(res.coord - center, axis=1).min()), 2),
                     "in_box": site[(chain, rid)]["in_box"], "contacts_ligand": site[(chain, rid)]["contacts_ligand"],
                     "note": HIS_NOTE if rname == "HIS" else ""})
    return rows


def backbone_gaps(prot, center) -> list:
    out = []
    for i in struc.check_res_id_continuity(prot):
        a, b = prot[i - 1], prot[i]
        da = float(np.linalg.norm(prot.coord[prot.res_id == a.res_id] - center, axis=1).min())
        db = float(np.linalg.norm(prot.coord[prot.res_id == b.res_id] - center, axis=1).min())
        out.append({"after": int(a.res_id), "before": int(b.res_id), "n_missing": int(b.res_id - a.res_id - 1),
                    "flank_dist_to_site": [round(da, 1), round(db, 1)]})
    return out


def pdbqt_cofactors(path: Path) -> dict:
    """Per cofactor: atoms, polar H, summed partial charge, read from Meeko's PDBQT (fixed-width columns)."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            r = line[17:20].strip()
            if r in FLAVIN + NADP:
                d = out.setdefault(r, {"atoms": 0, "polar_h": 0, "charge": 0.0})
                d["atoms"] += 1
                d["polar_h"] += line[77:79].strip() == "HD"
                d["charge"] = round(d["charge"] + float(line[70:76]), 4)
    return out


# --------------------------------------------------------------------------- pipeline
def _write_csv(path: Path, rows: list, cols=None):
    cols = cols or (list(rows[0]) if rows else [])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _residues(a):
    return {(str(c), int(r), str(n)) for c, r, n in zip(a.chain_id, a.res_id, a.res_name)}


def prepare_receptor(cif_path, pdb_id, chain, cfg, nadp, out_dir, *, pdb2pqr=None, templates=TEMPLATES) -> dict:
    t0 = time.time()
    out_dir, cif_path = Path(out_dir), Path(cif_path)
    if out_dir.exists():
        raise ConfigError("output directory already exists, refusing to overwrite: %s" % out_dir)
    work = out_dir / "work"
    work.mkdir(parents=True)
    write_resolved(cfg, out_dir)

    atoms, altloc_log = read_structure(cif_path, cfg.alt_conformers)
    center = resolve_site_center(atoms, chain, cfg.box_center_atom)
    sel = select_atoms(atoms, chain=chain, nadp=nadp, cfg=cfg, center=center)
    log = altloc_log + sel.log
    _write_csv(out_dir / "removal_log.csv", log, ["chain_id", "res_id", "res_name", "atom_name", "altloc", "occupancy", "reason", "dist_to_site"])
    if sel.ligand is not None:
        f = pdb.PDBFile(); pdb.set_structure(f, sel.ligand); f.write(str(out_dir / "ligand_crystal.pdb"))

    # 1. protein (+ kept waters) -> PDB2PQR / PROPKA at the configured pH
    kept = sel.kept
    prot_in = kept[struc.filter_canonical_amino_acids(kept) | struc.filter_solvent(kept)]
    f = pdb.PDBFile(); pdb.set_structure(f, prot_in); f.write(str(work / "protein_in.pdb"))
    p2p = _resolve_pdb2pqr(pdb2pqr)
    cmd1 = [p2p, "--ff=AMBER", "--with-ph=%g" % cfg.protonation_ph, "--titration-state-method=propka", "--keep-chain",
            "--pdb-output", work / "protein_h.pdb", work / "protein_in.pdb", work / "protein.pqr"]
    if not struc.filter_solvent(prot_in).any():
        cmd1.insert(5, "--drop-water")
    log1 = _run(cmd1, work / "pdb2pqr.log")
    plog = work / "protein.log"                      # PDB2PQR writes its full log + the PROPKA summary beside the PQR
    ptext = plog.read_text() if plog.exists() else log1
    prot_h = pdb.PDBFile.read(str(work / "protein_h.pdb")).get_structure(model=1)

    # 2. protein (with H) + cofactors/ions at crystal coordinates -> Meeko
    extra = kept[np.isin(kept.res_name, FLAVIN + NADP) | struc.filter_monoatomic_ions(kept)].copy()
    for cat in ("occupancy", "altloc_id"):
        if cat in extra.get_annotation_categories():
            extra.del_annotation(cat)
    extra.bonds = None
    combo = prot_h + extra
    f = pdb.PDBFile(); pdb.set_structure(f, combo); f.write(str(work / "combined.pdb"))
    box = list(cfg.box_size_angstrom)
    cmd2 = [_meeko_bin(), "--read_pdb", work / "combined.pdb", "-o", out_dir / "receptor", "-p", "-j", "-v",
            "--write_pdb", out_dir / "receptor_prepared.pdb", "--box_center", *["%.3f" % x for x in center],
            "--box_size", *["%g" % x for x in box], "--delete_bad_res_from_box_radius", "%g" % cfg.pocket_radius_angstrom,
            "--add_templates", templates]
    log2 = _run(cmd2, work / "meeko.log")

    # 3. what Meeko kept / deleted, protonation table, gaps, report
    got = pdbqt.PDBQTFile.read(str(out_dir / "receptor.pdbqt")).get_structure(model=1)
    gone = sorted(_residues(combo) - _residues(got), key=lambda k: k[1])
    deleted = [{"chain": c, "res_id": r, "res_name": n,
                "dist_to_site": round(float(np.linalg.norm(combo.coord[(combo.chain_id == c) & (combo.res_id == r)] - center, axis=1).min()), 1)}
               for c, r, n in gone]
    heavy_h = prot_h[prot_h.element != "H"]
    site = site_residue_ids(heavy_h, center, box, None if sel.ligand is None else sel.ligand.coord, cfg.site_contact_angstrom)
    rows = protonation_rows(prot_h, site, center, parse_propka(ptext))
    _write_csv(out_dir / "protonation_states.csv", rows,
               ["chain_id", "res_id", "res_name", "state", "net_charge", "pKa", "model_pKa", "dist_to_site", "in_box", "contacts_ligand", "note"])
    prot_only = kept[struc.filter_canonical_amino_acids(kept)]
    rebuilt = parse_rebuilt_atoms(log1 + "\n" + ptext)
    for r in rebuilt["residues"]:                       # modelled atoms feed the docking: how close to the site are they?
        m = (prot_h.chain_id == r["chain_id"]) & (prot_h.res_id == r["res_id"])
        r["dist_to_site"] = round(float(np.linalg.norm(prot_h.coord[m] - center, axis=1).min()), 1)
    rebuilt["nearest_dist_to_site"] = min((r["dist_to_site"] for r in rebuilt["residues"]), default=None)
    rebuilt["within_pocket_radius"] = [r for r in rebuilt["residues"] if r["dist_to_site"] <= cfg.pocket_radius_angstrom]
    report = {
        "pdb_id": pdb_id, "chain": chain, "nadp_policy": nadp, "input_sha256": hashlib.sha256(cif_path.read_bytes()).hexdigest(),
        "config": cfg.to_plain(), "site_center": [float(x) for x in center], "box_size": box,
        "fad_retained_in_receptor": bool(np.isin(kept.res_name, FLAVIN).any()), "fad_statement": FAD_STATEMENT,
        "removed_by_reason": dict(Counter(r["reason"] for r in log)), "n_atoms_kept": int(kept.array_length()),
        "substrate_set_aside": None if sel.ligand is None else {"res_name": cfg.ligand_resname, "n_atoms": int(sel.ligand.array_length())},
        "site_definition": {"rule": "protein residues with an atom in the docking box, plus residues within site_contact_angstrom of the "
                            "deposited ligand (if the structure has one)", "box_size": box, "contact_angstrom": cfg.site_contact_angstrom,
                            "n_residues": len(site), "n_ionizable": len(rows), "ligand_used": sel.ligand is not None},
        "pdb2pqr_rebuilt_atoms": rebuilt, "meeko_deleted_residues": deleted,
        "cofactors": pdbqt_cofactors(out_dir / "receptor.pdbqt"), "backbone_gaps_not_rebuilt": backbone_gaps(prot_only, center),
        "tools": {"biotite": biotite.__version__, "meeko": _version([_meeko_bin(), "--help"]) and __import__("meeko").__version__,
                  "pdb2pqr": _version([p2p, "--version"]), "propka": _version([p2p.parent / "propka3", "--version"])},
        "commands": [[str(c) for c in cmd1], [str(c) for c in cmd2]], "wall_seconds": round(time.time() - t0, 1),
    }
    (out_dir / "prep_report.json").write_text(json.dumps(report, indent=1, sort_keys=True))
    return report


# --------------------------------------------------------------------------- batch
def job_name(pdb_id, chain, nadp) -> str:
    return "%s_%s_nadp-%s" % (pdb_id, chain, nadp)


def plan_jobs(cfg, flavin_chains: dict) -> list:
    """(pdb_id, chain, nadp) for everything Steps 3-5 need, derived from the config; duplicates collapse."""
    s = cfg.structures
    jobs = [(s.benchmark_primary, cfg.receptor_chain, cfg.redock_nadp)]
    jobs += [(s.benchmark_secondary, c, cfg.redock_nadp) for c in flavin_chains[s.benchmark_secondary]]   # one redock per chain
    jobs += [(a.structure, cfg.receptor_chain, a.nadp) for a in cfg.cross_docking_arms]
    return list(dict.fromkeys(jobs))


def prepare_all(cfg, out_root, cache_dir, **kw) -> dict:
    """Prepare every planned receptor. A failing job is recorded and re-raised at the end; never substituted."""
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    ids = sorted({cfg.structures.benchmark_primary, cfg.structures.benchmark_secondary,
                  cfg.structures.fad_out, cfg.structures.fad_in})
    paths = {i: fetch_structure(i, cache_dir) for i in ids}
    flavin_chains = {}
    for i, pth in paths.items():
        a, _ = read_structure(pth, cfg.alt_conformers)
        flavin_chains[i] = sorted({str(c) for c in a.chain_id[np.isin(a.res_name, FLAVIN)]})
    reports, failures = {}, {}
    for pdb_id, chain, nadp in plan_jobs(cfg, flavin_chains):
        name = job_name(pdb_id, chain, nadp)
        try:
            reports[name] = prepare_receptor(paths[pdb_id], pdb_id, chain, cfg, nadp, out_root / name, **kw)
        except Exception as e:                                   # noqa: BLE001 - recorded, then re-raised below
            failures[name] = "%s: %s" % (type(e).__name__, e)
    if failures:
        (out_root / "prep_failures.json").write_text(json.dumps(failures, indent=1))
        raise RuntimeError("receptor preparation failed for %d job(s): %s" % (len(failures), sorted(failures)))
    return reports


def main(argv=None) -> int:
    import argparse
    from .config import load_config
    ap = argparse.ArgumentParser(description="Prepare every receptor Steps 3-5 need (biotite -> PDB2PQR/PROPKA -> Meeko).")
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True, help="output root; one sub-directory per receptor")
    ap.add_argument("--cache", default="structures_cache")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE", help="override top-level config keys, e.g. protonation_ph=6.5")
    a = ap.parse_args(argv)
    import yaml
    over = {k: yaml.safe_load(v) for k, _, v in (s.partition("=") for s in a.set)}
    reps = prepare_all(load_config(a.config, overrides=over), a.out, a.cache)
    for name, r in reps.items():
        print("%-22s site %s | cofactors %s | meeko-deleted %d | %.0fs" % (name, [round(x, 1) for x in r["site_center"]],
              {k: v["charge"] for k, v in r["cofactors"].items()}, len(r["meeko_deleted_residues"]), r["wall_seconds"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
