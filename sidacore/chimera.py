"""Chimeric-receptor diagnostic: does giving 4B63 the 6X0J side chains of a few residues change its redock?

Design registered before any run (see the generated report): an unmodified round trip through the same write/read/prepare path is the baseline (it must reproduce the
original 4B63 result, or the method itself perturbs it); swapping four residues that differ by only 0.1-0.4 A is the sensitivity control; the test variants swap the residues
the site comparison flagged. Same pass rule as Step 4 (median top-1 <= 2 A over >= 10 seeds, Vina). A variant that "passes" is a finding, not a fix: a chimeric receptor
is never adopted as the benchmark. Only side-chain atoms move: a rotamer transplant onto the target's own backbone (N, CA, C of each residue).
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np

from .config import load_config, write_resolved
from .diagnose import _site_states
from .prep import prepare_receptor, read_structure
from .boxes import protein_chain
from .redock import MIN_SEEDS, benchmark_jobs, crystal_reference, evaluate_seeds, ligand_template, prepare_ligand, redock_seeds
from .sidechain import fit_transform

ENGINE = "vina"
BACKBONE = ("N", "CA", "C", "O")
CLASH_CUTOFF = 2.2                                  # A: closer than this between non-bonded heavy atoms is a steric clash
# (label, donor chain, residues taken from the donor). The first is the baseline, the last the sensitivity control.
VARIANTS = [
    ("round trip, nothing swapped (baseline)", "A", []),
    ("Gln102 only", "A", [102]),
    ("Asp288 only", "A", [288]),
    ("Gln102 + Asp288", "A", [102, 288]),
    ("Gln102 + Asp288, donor chain C", "C", [102, 288]),
    ("all four flagged: Gln102, Arg144, Asp288, Arg328", "A", [102, 144, 288, 328]),
    ("control: four unflagged contacts Lys107, Asn293, Asn323, Ser469", "A", [107, 293, 323, 469]),
]


REVERSE_VARIANTS = [
    ("round trip, nothing swapped (baseline)", "A", []),
    ("4B63 Gln102", "A", [102]),
    ("4B63 Asp288", "A", [288]),
    ("4B63 Gln102 + Asp288", "A", [102, 288]),
    ("all four that differed: 4B63 Gln102, Arg144, Asp288, Arg328", "A", [102, 144, 288, 328]),
    ("control: four unflagged contacts from 4B63 (Lys107, Asn293, Asn323, Ser469)", "A", [107, 293, 323, 469]),
]
REVERSE_TARGETS = (("6X0J", "A"), ("6X0J", "B"), ("6X0J", "D"))       # the chains that pass the baseline redock


def job_for(cfg, prepared_root, target):
    """The prepared benchmark job of one structure AND chain (a name prefix is not enough: 6X0J has four chains)."""
    want = "%s_%s_" % target
    for j in benchmark_jobs(cfg, Path(prepared_root)):
        if j.name.startswith(want):
            return j
    raise ValueError("no prepared benchmark job for %s chain %s" % target)


def swap_sidechains(target, donor, res_ids):
    """(new target atoms, log): a rotamer transplant. For each residue the donor residue is superimposed on the TARGET residue's own backbone atoms (N, CA, C) and its
    side-chain atoms replace the target's, so the donor's side-chain torsions sit on the target's backbone and the CA-CB junction keeps the target's geometry. (Copying the donor's
    coordinates after a global fit stretched CA-CB to 2.05 A, because the two backbones differ by ~0.5 A there.) Residue types must match and every side-chain atom must exist in the donor."""
    out, log = target.copy(), []
    for rid in res_ids:
        t_res, d_res = target[target.res_id == rid], donor[donor.res_id == rid]
        if len(t_res) == 0 or len(d_res) == 0:
            raise ValueError("residue %d is missing from the target or the donor" % rid)
        if str(t_res.res_name[0]) != str(d_res.res_name[0]):
            raise ValueError("residue %d: target is %s but donor is %s" % (rid, t_res.res_name[0], d_res.res_name[0]))
        tb = t_res[np.isin(t_res.atom_name, ("N", "CA", "C"))]
        db = d_res[np.isin(d_res.atom_name, ("N", "CA", "C"))]
        if len(tb) != 3 or len(db) != 3 or list(tb.atom_name) != list(db.atom_name):
            raise ValueError("residue %d needs N, CA and C in both structures" % rid)
        _, tf = struc.superimpose(tb, db)
        moved = d_res.copy()
        moved.coord = tf.apply(d_res.coord)
        t_idx = np.nonzero((target.res_id == rid) & ~np.isin(target.atom_name, BACKBONE))[0]
        shifts = []
        for i in t_idx:
            hit = np.nonzero(moved.atom_name == target.atom_name[i])[0]
            if len(hit) != 1:
                raise ValueError("donor residue %d has %d atoms named %s" % (rid, len(hit), target.atom_name[i]))
            shifts.append(float(np.linalg.norm(moved.coord[hit[0]] - target.coord[i])))
            out.coord[i] = moved.coord[hit[0]]
        log.append({"res_id": int(rid), "res_name": str(t_res.res_name[0]), "n_atoms": int(len(t_idx)), "max_shift": max(shifts) if shifts else 0.0,
                    "backbone_rmsd": float(struc.rmsd(tb, tf.apply(db)))})
    return out, log


def closest_contact(atoms, res_ids, cutoff: float = CLASH_CUTOFF):
    """The closest approach (below `cutoff`) of any side-chain atom of `res_ids` to an atom of ANOTHER residue, or None."""
    heavy = atoms[atoms.element != "H"]
    cells = struc.CellList(heavy, cell_size=cutoff)
    best = None
    for rid in res_ids:
        for i in np.nonzero((heavy.res_id == rid) & ~np.isin(heavy.atom_name, BACKBONE))[0]:
            for j in cells.get_atoms(heavy.coord[i], cutoff):
                if int(heavy.res_id[j]) == rid and str(heavy.chain_id[j]) == str(heavy.chain_id[i]):
                    continue
                d = float(np.linalg.norm(heavy.coord[i] - heavy.coord[j]))
                if d < cutoff and (best is None or d < best["distance"]):
                    best = {"res_id": int(rid), "atom": str(heavy.atom_name[i]), "other_res_id": int(heavy.res_id[j]),
                            "other_res_name": str(heavy.res_name[j]), "other_atom": str(heavy.atom_name[j]), "distance": d}
    return best


def build_chimera_cif(cache_dir, out_path, target=("4B63", "A"), donor=("6X0J", "A"), res_ids=()) -> dict:
    """Write `target` with the side chains of `res_ids` taken from `donor`. Returns what was swapped and the closest non-bonded contact the swap created."""
    cache_dir = Path(cache_dir)
    t_all, _ = read_structure(cache_dir / (target[0] + ".cif"), "highest-occupancy")
    d_all, _ = read_structure(cache_dir / (donor[0] + ".cif"), "highest-occupancy")
    t_mask = (t_all.chain_id == target[1]) & struc.filter_canonical_amino_acids(t_all)
    d_prot = protein_chain(d_all, donor[1])
    new_prot, log = swap_sidechains(t_all[t_mask], d_prot, list(res_ids))
    out = t_all.copy()
    out.coord[t_mask] = new_prot.coord
    contact = closest_contact(out[~struc.filter_solvent(out)], list(res_ids)) if res_ids else None
    out.bonds = None
    f = pdbx.CIFFile()
    pdbx.set_structure(f, out, data_block="chimera")
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    f.write(str(out_path))
    return {"path": str(out_path), "swapped": log, "closest_contact": contact, "donor": "%s:%s" % donor}


def _slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def run_chimeras(cfg, prepared_root, out_root, cache_dir="structures_cache", target=("4B63", "A"), variants=VARIANTS, only=None, checkpoint=None, donor_pid="6X0J") -> dict:
    if len(cfg.seeds) < MIN_SEEDS:
        raise ValueError("Step 4 needs at least %d seeds, config has %d" % (MIN_SEEDS, len(cfg.seeds)))
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    write_resolved(cfg, out_root)
    base_job = job_for(cfg, prepared_root, target)
    ref = crystal_reference(base_job / "ligand_crystal.pdb", ligand_template(cfg.ligand_resname))
    lp = prepare_ligand(cfg.ligand_resname, cfg.protonation_ph)
    base_states = _site_states(base_job / "protonation_states.csv")
    results = []
    for label, donor_chain, res_ids in variants:
        if only and label not in only:
            continue
        t0, slug = time.time(), _slug(label)
        built = build_chimera_cif(cache_dir, out_root / "cif" / (slug + ".cif"), target=target, donor=(donor_pid, donor_chain), res_ids=res_ids)
        prep_dir = out_root / "receptors" / slug
        try:
            if not prep_dir.exists():
                prepare_receptor(built["path"], target[0], target[1], cfg, "keep", prep_dir)
        except RuntimeError as e:                    # e.g. a swap that overlaps a neighbour makes Meeko infer a spurious bond and refuse: record it, do not drop it
            results.append({"label": label, "donor_chain": donor_chain, "res_ids": list(res_ids), "swapped": built["swapped"], "closest_contact": built["closest_contact"],
                            "prepare_error": str(e)})
            if checkpoint:
                Path(checkpoint).write_text(json.dumps({"settings": {"partial": True}, "variants": results}, indent=1))
            continue
        states = _site_states(prep_dir / "protonation_states.csv")
        changed = {"%s%s" % k: "%s -> %s" % (base_states[k], states[k]) for k in base_states if k in states and base_states[k] != states[k]}
        rows = redock_seeds(prep_dir, lp, ref, cfg, ENGINE, out_root / "redock" / slug, receptor=prep_dir / "receptor.pdbqt", box_file=prep_dir / "receptor.box.txt")
        ev = evaluate_seeds(rows, cfg.redock_rmsd_acceptance, min_seeds=MIN_SEEDS)
        results.append({"label": label, "donor_chain": donor_chain, "res_ids": list(res_ids), "swapped": built["swapped"], "closest_contact": built["closest_contact"],
                        "site_states_changed": changed, "wall_seconds": round(time.time() - t0, 1), "per_seed": rows, **ev})
        if checkpoint:
            Path(checkpoint).write_text(json.dumps({"settings": {"partial": True}, "variants": results}, indent=1))
    return {"settings": {"engine": ENGINE, "exhaustiveness": cfg.exhaustiveness, "num_modes": cfg.num_modes, "seeds": list(cfg.seeds),
                         "acceptance": cfg.redock_rmsd_acceptance, "protonation_ph": cfg.protonation_ph, "target": "%s:%s" % target,
                         "direction": "forward" if target[0] == "4B63" else "reverse", "donor": donor_pid}, "variants": results}


def run_reverse(cfg, prepared_root, out_root, cache_dir="structures_cache", targets=REVERSE_TARGETS, donor=("4B63", "A"), checkpoint=None) -> dict:
    """The reverse chimera: each passing 6X0J chain receives the donor's (4B63's) side-chain conformations. Same variants, rule and controls as the forward run."""
    out_root, merged, settings = Path(out_root), [], None
    for tgt in targets:
        rep = run_chimeras(cfg, prepared_root, out_root / ("%s_%s" % tgt), cache_dir=cache_dir, target=tgt, variants=REVERSE_VARIANTS, donor_pid=donor[0])
        settings = rep["settings"]
        merged += [{"target": "%s:%s" % tgt, **v} for v in rep["variants"]]
        if checkpoint:
            Path(checkpoint).write_text(json.dumps({"settings": {"partial": True}, "variants": merged}, indent=1))
    settings = {**settings, "target": ", ".join("%s:%s" % t for t in targets), "direction": "reverse", "donor": "%s:%s" % donor}
    return {"settings": settings, "variants": merged}


def render_markdown(report: dict, notes: str = "") -> str:
    s = report["settings"]
    reverse = s.get("direction") == "reverse"
    multi = any("target" in v for v in report["variants"])
    if reverse:
        title = "# Step 4: reverse chimeras - are 4B63's Gln102 and Asp288 necessary for a 6X0J chain to pass?"
        what = ("Targets %s; each receives the %s side-chain conformations of the named residues as a rotamer transplant (the donor residue is superimposed on the target residue's own N, CA and C "
                "and only its side-chain atoms are copied)." % (s["target"], s.get("donor", "4B63")))
        reading = ("(5) Reading, fixed in advance: a swap that makes a passing chain fail while the control does not means those residues are necessary for the pass; if every chain still "
                   "passes, they are not. The baseline of each chain must reproduce its original pass.")
    else:
        title = "# Step 4: chimeric receptors - do Gln102 and Asp288 explain the 4B63 failure?"
        what = ("Target %s; donor residues come from the named 6X0J chain as a rotamer transplant: each donor residue is superimposed on the 4B63 residue's own N, CA and C and only its "
                "side-chain atoms are copied." % s["target"])
        reading = ""
    L = [title, "", "Generated by `python -m sidacore.chimera`; numbers only. " + what, "",
         "**Design registered before the run.** (1) The round trip with nothing swapped is the baseline: it goes through the same write, read and prepare path and must reproduce the original "
         "result, otherwise the method itself changes the answer. (2) The control swaps four residues that differ by only 0.1 to 0.5 A: if it changes the result as much as the real swap, "
         "docking is simply sensitive to sub-angstrom changes and the real swap means little. (3) Same pass rule as Step 4: median top-1 RMSD over %d seeds <= %g A (%s, exhaustiveness %d, %d modes). "
         "(4) A variant that passes is a finding, never a fix; a chimeric receptor is not adopted as the benchmark. %s" % (
             len(s["seeds"]), s["acceptance"], s["engine"], s["exhaustiveness"], s["num_modes"], reading), "",
         ("| Target " if multi else "") + "| Variant | residues swapped | closest new contact (A) | top-1 median RMSD (within / seeds) | best-in-top-5 median | verdict | protonation changes in the site |",
         ("|---" if multi else "") + "|---|---|---|---|---|---|---|"]
    for v in report["variants"]:
        c = v["closest_contact"]
        head = ("| %s " % v["target"]) if multi else ""
        swapped = ", ".join("%s%d" % (x["res_name"], x["res_id"]) for x in v["swapped"]) or "none"
        contact = "none under 2.2" if not c else "%s%d %s - %s%d %s: %.2f" % (next(x["res_name"] for x in v["swapped"] if x["res_id"] == c["res_id"]), c["res_id"], c["atom"],
                                                                           c["other_res_name"], c["other_res_id"], c["other_atom"], c["distance"])
        if "prepare_error" in v:
            why = "; ".join([l for l in v["prepare_error"].splitlines() if "excess inter-residue" in l] or v["prepare_error"].splitlines()[:1])
            L.append("%s| %s | %s | %s | not prepared: %s | - | - | - |" % (head, v["label"], swapped, contact, why))
            continue
        L.append("%s| %s | %s | %s | %.2f (%d/%d) | %.2f | %s | %s |" % (
            head, v["label"], swapped, contact, v["top1_median"], v["top1_n_within"], v["n_seeds"], v["best_top5_median"], "**PASS**" if v["passes"] else "fail",
            "; ".join("%s: %s" % kv for kv in v["site_states_changed"].items()) or "none"))
    return "\n".join(L) + "\n" + ("\n" + notes if notes else "")


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Chimeric-receptor diagnostic for the 4B63 redock failure.")
    ap.add_argument("--config", required=True)
    ap.add_argument("--prepared", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default="structures_cache")
    ap.add_argument("--report")
    ap.add_argument("--notes")
    ap.add_argument("--direction", choices=("forward", "reverse"), default="forward", help="forward: 6X0J side chains into 4B63; reverse: 4B63 side chains into the passing 6X0J chains")
    ap.add_argument("--render-only", action="store_true", help="rebuild --report from OUT/step4_chimera.json; no docking")
    a = ap.parse_args(argv)
    cfg, out = load_config(a.config), Path(a.out)
    if a.render_only:
        report = json.loads((out / "step4_chimera.json").read_text())
    else:
        part = out / "step4_chimera.partial.json"
        out.mkdir(parents=True, exist_ok=True)
        report = (run_reverse(cfg, a.prepared, out, cache_dir=a.cache, checkpoint=part) if a.direction == "reverse"
                  else run_chimeras(cfg, a.prepared, out, cache_dir=a.cache, checkpoint=part))
        (out / "step4_chimera.json").write_text(json.dumps(report, indent=1))
    if a.report:
        notes = Path(a.notes).read_text() if a.notes and Path(a.notes).exists() else ""
        Path(a.report).write_text(render_markdown(report, notes))
    for v in report["variants"]:
        if "prepare_error" in v:
            print("%-8s %-72s NOT PREPARED" % (v.get("target", ""), v["label"][:72]))
            continue
        print("%-8s %-72s top-1 %.2f (%d/%d) best5 %.2f %s" % (v.get("target", ""), v["label"][:72], v["top1_median"], v["top1_n_within"], v["n_seeds"], v["best_top5_median"], "PASS" if v["passes"] else "fail"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
