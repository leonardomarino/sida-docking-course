"""Rosetta (PyRosetta) rescoring of docked ligand poses: an independent check on the Vina / Uni-Dock results.

LICENCE: PyRosetta is free for academic use and needs a licence for commercial use. This module is deliberately outside
`sidacore` (the shared layer, which must stay free of non-commercial components) and `sidacore` never imports it.

Protocol, fixed before any docked pose was scored and identical for the deposited pose and every docked pose:
  1. receptor (protein + FAD/FDA + NADP, heavy atoms) and ligand are read together; residue types for the cofactors and
     the zwitterionic ornithine come from the PDB component dictionary; Rosetta adds its own hydrogens;
  2. ref2015 minimisation: the ligand's rigid-body jump and torsions, and the side chains of protein residues within
     `shell` angstrom of it. Backbone and cofactors stay fixed;
  3. InterfaceAnalyzer, ligand against everything else, with the partners separated and repacked (dG_separated).
Lower dG_separated is a better interface. It is a second, independent scoring function, not ground truth.
"""
from __future__ import annotations

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pyrosetta
from pyrosetta import rosetta as R

OPTIONS = "-load_PDB_components true -ignore_unrecognized_res false -mute all -ex1 -ex2 -use_input_sc -constant_seed true -jran 1"
LIGAND_CHAIN = "L"
_ready = False


def init_rosetta():
    global _ready
    if not _ready:
        pyrosetta.init(OPTIONS, silent=True)
        _ready = True


def _heavy_xyz(res):
    return np.array([[res.xyz(a).x, res.xyz(a).y, res.xyz(a).z] for a in range(1, res.natoms() + 1)
                     if not res.atom_is_hydrogen(a)])


def _ligand_jump(pose, lig):
    ft = pose.fold_tree()
    for j in range(1, pose.num_jump() + 1):
        if lig in (ft.downstream_jump_residue(j), ft.upstream_jump_residue(j)):
            return j
    raise RuntimeError("the ligand residue %d has no jump in the fold tree" % lig)


def _receptor_chains(text: str) -> str:
    """Chain letters of the receptor file (protein and cofactors share one). Never assumed to be 'A'."""
    chains = sorted({l[21] for l in text.splitlines() if l.startswith(("ATOM", "HETATM"))} - {LIGAND_CHAIN, " "})
    if not chains:
        raise ValueError("the receptor file has no chain letters other than the ligand's %r" % LIGAND_CHAIN)
    return "".join(chains)


def rescore(receptor_pdb, ligand_pdb, shell: float = 8.0, out_ligand_pdb=None) -> dict:
    init_rosetta()
    receptor_text = "".join(l for l in Path(receptor_pdb).read_text().splitlines(True) if not l.startswith(("END", "TER")))
    partners = "%s_%s" % (LIGAND_CHAIN, _receptor_chains(receptor_text))
    text = receptor_text + Path(ligand_pdb).read_text()
    pose = pyrosetta.Pose()
    R.core.import_pose.pose_from_pdbstring(pose, text)
    lig = next(i for i in range(1, pose.total_residue() + 1) if pose.pdb_info().chain(i) == LIGAND_CHAIN)
    sf = pyrosetta.create_score_function("ref2015")
    start, xyz0 = sf(pose), _heavy_xyz(pose.residue(lig))

    sel = R.core.select.residue_selector
    near = sel.NeighborhoodResidueSelector(sel.ChainSelector(LIGAND_CHAIN), shell, False).apply(pose)
    mm = R.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    mm.set_jump(_ligand_jump(pose, lig), True)
    mm.set_bb(lig, True)
    mm.set_chi(lig, True)
    for i in range(1, pose.total_residue() + 1):
        if near[i] and pose.residue(i).is_protein():
            mm.set_chi(i, True)
    R.protocols.minimization_packing.MinMover(mm, sf, "lbfgs_armijo_nonmonotone", 0.001, True).apply(pose)
    minimised = sf(pose)

    ia = R.protocols.analysis.InterfaceAnalyzerMover(
        R.core.pose.DockingPartners.docking_partners_from_string(partners), False, sf, False, False, True, False, False)
    ia.apply(pose)
    data = ia.get_all_data()
    xyz1 = _heavy_xyz(pose.residue(lig))
    if out_ligand_pdb:
        sub = pyrosetta.Pose()
        sub.append_residue_by_jump(pose.residue(lig), 1)
        sub.dump_pdb(str(out_ligand_pdb))
    return {"score_start": start, "score_min": minimised,
            "dG_separated": ia.get_separated_interface_energy(), "dSASA": ia.get_interface_delta_sasa(),
            "sc_value": data.sc_value, "unsat_hbonds": ia.get_interface_delta_hbond_unsat(),
            "n_interface_residues": ia.get_num_interface_residues(),
            "ligand_shift": float(np.sqrt(np.mean(np.sum((xyz1 - xyz0) ** 2, axis=1))))}


def _one(args):
    receptor, ligand, out_dir = args
    name = Path(ligand).stem
    row = rescore(receptor, ligand, out_ligand_pdb=Path(out_dir) / ("%s.min.pdb" % name))
    return {"pose": name, **row}


def rescore_many(receptor, ligands, out_dir, jobs: int = 8) -> list:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work = [(str(receptor), str(l), str(out_dir)) for l in ligands]
    with ProcessPoolExecutor(max_workers=jobs, initializer=init_rosetta) as ex:
        rows = list(ex.map(_one, work))
    if not any(r["n_interface_residues"] > 0 for r in rows):
        raise RuntimeError("no pose has any interface residue with the receptor: wrong chain letters or ligand far from the site?")
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Rescore docked ligand poses with Rosetta ref2015 (academic licence).")
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--ligands", nargs="+", required=True, help="ligand PDB files (chain L)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args(argv)
    rows = rescore_many(a.receptor, a.ligands, a.out_dir, a.jobs)
    with open(a.csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("%d poses rescored -> %s" % (len(rows), a.csv))
    return 0


if __name__ == "__main__":
    sys.exit(main())
