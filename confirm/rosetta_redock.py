"""RosettaLigand redock of the deposited substrate: an independent search AND scoring function for the Step 4 question.

LICENCE: PyRosetta is free for academic use and needs a licence for commercial use. This module is outside `sidacore`
(which must stay free of non-commercial components) and `sidacore` never imports it.

Protocol, fixed before any decoy was run:
  1. relax the apo receptor (protein + FAD/FDA + NADP) with FastRelax restrained to the crystal coordinates, once per receptor;
  2. for each seed, a start ligand (one RDKit conformer per seed, centred on the C4X box, never the crystal pose) is moved by
     Transform (rigid-body Monte Carlo on a scoring grid, sphere of `box_size` A), then HighResDocker (Monte Carlo repack +
     minimise), then FinalMinimizer, then InterfaceScoreCalculator;
  3. `mode="native"` runs only steps HighResDocker/FinalMinimizer from the deposited pose: the refinement score of the crystal
     pose under the same protocol, to separate a search failure from a scoring failure.
The interface energy is evaluated with the coordinate-constraint term switched off: left in, it returned -2.8 million.
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


def protocol_xml(mode: str = "dock", box_size: float = 11.0, transform_cycles: int = 500, transform_repeats: int = 3,
                 hrd_cycles: int = 6) -> str:
    if mode not in ("dock", "native"):
        raise ValueError("mode must be 'dock' or 'native', got %r" % mode)
    steps = ('      <Add mover_name="transform"/>\n' if mode == "dock" else "") + \
            '      <Add mover_name="high_res_docker"/>\n      <Add mover_name="final"/>\n      <Add mover_name="add_scores"/>\n'
    hard = '''      <Reweight scoretype="fa_intra_rep" weight="0.004"/>
      <Reweight scoretype="fa_elec" weight="0.42"/>
      <Reweight scoretype="hbond_bb_sc" weight="1.3"/>
      <Reweight scoretype="hbond_sc" weight="1.3"/>
      <Reweight scoretype="rama" weight="0.2"/>
'''
    return '''<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="ligand_soft_rep" weights="ligand_soft_rep">
      <Reweight scoretype="fa_elec" weight="0.42"/>
      <Reweight scoretype="hbond_bb_sc" weight="1.3"/>
      <Reweight scoretype="hbond_sc" weight="1.3"/>
      <Reweight scoretype="rama" weight="0.2"/>
    </ScoreFunction>
    <ScoreFunction name="hard_rep" weights="ligand">
%s    </ScoreFunction>
    <ScoreFunction name="iface" weights="ligand">
%s      <Reweight scoretype="coordinate_constraint" weight="0"/>
    </ScoreFunction>
  </SCOREFXNS>
  <SCORINGGRIDS ligand_chain="%s" width="25.0">
    <ClassicGrid grid_name="classic" weight="1.0"/>
  </SCORINGGRIDS>
  <LIGAND_AREAS>
    <LigandArea name="docking_sidechain" chain="%s" cutoff="6.0" add_nbr_radius="true" all_atom_mode="false"/>
    <LigandArea name="final_sidechain" chain="%s" cutoff="6.0" add_nbr_radius="true" all_atom_mode="false"/>
    <LigandArea name="final_backbone" chain="%s" cutoff="7.0" add_nbr_radius="false" all_atom_mode="true" Calpha_restraints="0.3"/>
  </LIGAND_AREAS>
  <INTERFACE_BUILDERS>
    <InterfaceBuilder name="side_chain_for_docking" ligand_areas="docking_sidechain"/>
    <InterfaceBuilder name="side_chain_for_final" ligand_areas="final_sidechain"/>
    <InterfaceBuilder name="backbone" ligand_areas="final_backbone" extension_window="3"/>
  </INTERFACE_BUILDERS>
  <MOVEMAP_BUILDERS>
    <MoveMapBuilder name="docking" sc_interface="side_chain_for_docking"/>
    <MoveMapBuilder name="final" sc_interface="side_chain_for_final" bb_interface="backbone"/>
  </MOVEMAP_BUILDERS>
  <MOVERS>
    <Transform name="transform" chain="%s" box_size="%s" move_distance="0.2" angle="10" cycles="%d" repeats="%d" temperature="5"/>
    <HighResDocker name="high_res_docker" cycles="%d" repack_every_Nth="3" scorefxn="ligand_soft_rep" movemap_builder="docking"/>
    <FinalMinimizer name="final" scorefxn="hard_rep" movemap_builder="final"/>
    <InterfaceScoreCalculator name="add_scores" chains="%s" scorefxn="iface"/>
    <ParsedProtocol name="run">
%s    </ParsedProtocol>
  </MOVERS>
  <PROTOCOLS>
    <Add mover_name="run"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>
''' % (hard, hard, LIGAND_CHAIN, LIGAND_CHAIN, LIGAND_CHAIN, LIGAND_CHAIN, LIGAND_CHAIN, "%.1f" % box_size, transform_cycles,
       transform_repeats, hrd_cycles, LIGAND_CHAIN, steps)


def _complex(receptor_pdb, ligand_pdb):
    rec = "".join(l for l in Path(receptor_pdb).read_text().splitlines(True) if not l.startswith(("END", "TER")))
    pose = pyrosetta.Pose()
    R.core.import_pose.pose_from_pdbstring(pose, rec + Path(ligand_pdb).read_text())
    return pose


def _ligand_residue(pose) -> int:
    return next(i for i in range(1, pose.total_residue() + 1) if pose.pdb_info().chain(i) == LIGAND_CHAIN)


def relax_receptor(receptor_pdb, out_pdb, repeats: int = 3):
    """FastRelax of the apo receptor, restrained to the starting coordinates (ref2015 + coordinate constraints)."""
    init_rosetta()
    pose = pyrosetta.Pose()
    R.core.import_pose.pose_from_pdbstring(pose, "".join(l for l in Path(receptor_pdb).read_text().splitlines(True)
                                                          if not l.startswith(("END", "TER"))))
    sf = pyrosetta.create_score_function("ref2015_cst")
    fr = R.protocols.relax.FastRelax(sf, repeats)
    fr.constrain_relax_to_start_coords(True)
    fr.coord_constrain_sidechains(True)
    fr.ramp_down_constraints(False)
    before = sf(pose)
    fr.apply(pose)
    pose.dump_pdb(str(out_pdb))
    return {"score_before": before, "score_after": sf(pose)}


def dock_decoy(receptor_pdb, start_ligand_pdb, seed: int, out_ligand_pdb=None, mode: str = "dock", box_size: float = 11.0,
               transform_cycles: int = 500, transform_repeats: int = 3, hrd_cycles: int = 6) -> dict:
    init_rosetta()
    pose = _complex(receptor_pdb, start_ligand_pdb)
    objs = R.protocols.rosetta_scripts.XmlObjects.create_from_string(
        protocol_xml(mode, box_size, transform_cycles, transform_repeats, hrd_cycles))
    R.numeric.random.rg().set_seed(int(seed))
    objs.get_mover("run").apply(pose)
    sc = dict(pose.scores)
    if out_ligand_pdb:
        sub = pyrosetta.Pose()
        sub.append_residue_by_jump(pose.residue(_ligand_residue(pose)), 1)
        sub.dump_pdb(str(out_ligand_pdb))
    return {"total_score": float(sc["total_score"]), "interface_delta": float(sc["interface_delta_%s" % LIGAND_CHAIN]),
            "grid_score": float(sc.get("Grid_score", float("nan"))), "touching": float(sc.get("ligand_is_touching_%s" % LIGAND_CHAIN, float("nan")))}


def _one(task):
    receptor, start, seed, decoy, mode, out_dir, kw = task
    name = "%s_s%d_d%02d" % (mode, seed, decoy)
    try:
        row = dock_decoy(receptor, start, seed * 1000 + decoy, Path(out_dir) / (name + ".pdb"), mode=mode, **kw)
        return {"decoy": name, "mode": mode, "seed": seed, "index": decoy, "error": "", **row}
    except Exception as e:                                       # noqa: BLE001 - recorded in the CSV, never dropped silently
        return {"decoy": name, "mode": mode, "seed": seed, "index": decoy, "error": str(e).strip().splitlines()[-1][:120],
                "total_score": float("nan"), "interface_delta": float("nan"), "grid_score": float("nan"), "touching": float("nan")}


def run_decoys(receptor, starts: dict, out_dir, decoys_per_seed: int, mode: str = "dock", jobs: int = 8, **kw) -> list:
    """starts: {seed: start_ligand_pdb}. One decoy per (seed, index); the RNG seed is seed*1000+index."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tasks = [(str(receptor), str(p), int(s), d, mode, str(out_dir), kw) for s, p in sorted(starts.items())
             for d in range(1, decoys_per_seed + 1)]
    with ProcessPoolExecutor(max_workers=jobs, initializer=init_rosetta) as ex:
        rows = list(ex.map(_one, tasks))
    if all(r["error"] for r in rows):
        raise RuntimeError("every decoy failed; first error: %s" % rows[0]["error"])
    return rows


def _write_csv(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="RosettaLigand redock (academic licence).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("relax", help="restrained FastRelax of the apo receptor")
    r.add_argument("--receptor", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--repeats", type=int, default=3)
    d = sub.add_parser("dock", help="decoys from per-seed start ligands, or native-start refinement")
    d.add_argument("--receptor", required=True, help="the relaxed receptor")
    d.add_argument("--starts", nargs="+", required=True, help="start ligand PDBs named start_s<seed>.pdb (dock) or any file (native)")
    d.add_argument("--mode", choices=["dock", "native"], default="dock")
    d.add_argument("--decoys-per-seed", type=int, default=10)
    d.add_argument("--out-dir", required=True)
    d.add_argument("--csv", required=True)
    d.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args(argv)
    if a.cmd == "relax":
        print(relax_receptor(a.receptor, a.out, a.repeats))
        return 0
    starts = {}
    for p in a.starts:
        stem = Path(p).stem
        starts[int(stem.rsplit("_s", 1)[1]) if "_s" in stem else 0] = p
    rows = run_decoys(a.receptor, starts, a.out_dir, a.decoys_per_seed, a.mode, a.jobs)
    _write_csv(rows, a.csv)
    print("%d decoys (%d failed) -> %s" % (len(rows), sum(1 for x in rows if x["error"]), a.csv))
    return 0


if __name__ == "__main__":
    sys.exit(main())
