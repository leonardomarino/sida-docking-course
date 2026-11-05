"""Inputs for the Rosetta rescoring (confirm/rosetta_rescore.py), made with the licence-free sidacore tools.

Writes, per benchmark receptor: receptor.pdb (protein + cofactors, heavy atoms, original residue names), poses/crystal.pdb
(the deposited ligand) and poses/<engine>_s<seed>_m<rank>.pdb for every docked mode of one seed, all carrying the deposited
atom names on chain L. Returns the manifest rows (also written as manifest.csv): engine, rank, Vina-style score and the
heavy-atom symmetry-corrected RMSD to the deposited pose.
"""
from __future__ import annotations

import csv
from pathlib import Path

from sidacore.redock import benchmark_jobs, crystal_reference, ligand_template, pose_pdb_block, pose_rmsd, read_poses


def _heavy_receptor(job_dir: Path) -> str:
    keep = []
    for l in (job_dir / "work" / "combined.pdb").read_text().splitlines(True):
        if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() != "H":
            keep.append(l)
    return "".join(keep)


def _hetatm_only(block: str) -> str:
    return "".join(l for l in block.splitlines(True) if l.startswith("HETATM")) + "END\n"


def make_inputs(cfg, prepared_root, redock_root, out_root, engines=("vina", "unidock"), seed: int = 1, only=None) -> list:
    prepared_root, redock_root, out_root = Path(prepared_root), Path(redock_root), Path(out_root)
    template = ligand_template(cfg.ligand_resname)
    rows = []
    for job in benchmark_jobs(cfg, prepared_root):
        if only and job.name not in only:
            continue
        ref = crystal_reference(job / "ligand_crystal.pdb", template)
        d = out_root / job.name
        (d / "poses").mkdir(parents=True, exist_ok=True)
        (d / "receptor.pdb").write_text(_heavy_receptor(job))
        (d / "poses" / "crystal.pdb").write_text(_hetatm_only(pose_pdb_block(ref, ref)))
        rows.append({"job": job.name, "engine": "crystal", "rank": 0, "score": "", "rmsd_to_crystal": 0.0, "pose": "crystal"})
        for eng in engines:
            src = redock_root / job.name / eng / ("%s_seed%d.pdbqt" % (eng, seed))
            for rank, (score, mol) in enumerate(read_poses(src.read_text()), 1):
                name = "%s_s%d_m%02d" % (eng, seed, rank)
                (d / "poses" / (name + ".pdb")).write_text(_hetatm_only(pose_pdb_block(mol, ref)))
                rows.append({"job": job.name, "engine": eng, "rank": rank, "score": score,
                             "rmsd_to_crystal": pose_rmsd(mol, ref), "pose": name})
    out_root.mkdir(parents=True, exist_ok=True)
    with open(out_root / "manifest.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["job", "engine", "rank", "score", "rmsd_to_crystal", "pose"])
        w.writeheader()
        w.writerows(rows)
    return rows


def main(argv=None) -> int:
    import argparse
    import sys
    from sidacore.config import load_config
    ap = argparse.ArgumentParser(description="Make Rosetta rescoring inputs from the Step 2 receptors and the Step 4 poses.")
    ap.add_argument("--config", required=True)
    ap.add_argument("--prepared", required=True)
    ap.add_argument("--redock", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    rows = make_inputs(cfg, a.prepared, a.redock, a.out, seed=a.seed)
    starts = make_start_ligands(cfg, a.prepared, a.out, cfg.seeds)
    print("%d poses over %d receptors, %d start ligands -> %s" % (len(rows), len({r["job"] for r in rows}), len(starts), a.out))
    return 0



def make_start_ligands(cfg, prepared_root, out_root, seeds, only=None) -> list:
    """One start ligand per seed for the Rosetta redock: an RDKit conformer (embedded with that seed) of the protonated
    ligand, atom-named like the deposited ligand, its heavy-atom centroid placed at the centre of the C4X box. The crystal
    pose is never used; the rows record each start's RMSD to it so a leak would be visible."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from sidacore.boxes import read_vina_box
    from sidacore.redock import prepare_ligand
    smiles = prepare_ligand(cfg.ligand_resname, cfg.protonation_ph).smiles
    template = ligand_template(cfg.ligand_resname)
    rows = []
    for job in benchmark_jobs(cfg, Path(prepared_root)):
        if only and job.name not in only:
            continue
        ref = crystal_reference(job / "ligand_crystal.pdb", template)
        centre, _ = read_vina_box(job / "receptor.box.txt")
        d = Path(out_root) / job.name / "starts"
        d.mkdir(parents=True, exist_ok=True)
        for seed in seeds:
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            if AllChem.EmbedMolecule(mol, randomSeed=int(seed)) != 0:
                raise RuntimeError("could not embed %s with seed %d" % (smiles, seed))
            AllChem.MMFFOptimizeMolecule(mol)
            conf = mol.GetConformer()
            xyz = conf.GetPositions()
            heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
            shift = centre - xyz[heavy].mean(axis=0)
            for i in range(mol.GetNumAtoms()):
                p = xyz[i] + shift
                conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))
            (d / ("start_s%d.pdb" % seed)).write_text(_hetatm_only(pose_pdb_block(mol, ref)))
            rows.append({"job": job.name, "seed": int(seed), "rmsd_to_crystal": pose_rmsd(mol, ref)})
    return rows


if __name__ == "__main__":
    import sys
    sys.exit(main())
