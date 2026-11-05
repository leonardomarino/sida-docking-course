"""Step 4: redocking validation. The state of the ligand that is actually docked is checked, not assumed."""
import shutil
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from sidacore.config import load_config
from sidacore.redock import (UNIDOCK_DEFAULT, benchmark_jobs, crystal_reference, evaluate_seeds, ligand_template,
                             main as redock_main, pose_pdb_block, pose_rmsd, prepare_ligand, protonate, read_poses, redock_all, redock_seeds, render_markdown,
                             run_unidock, run_vina)

ROOT = Path(__file__).resolve().parent.parent
PREPARED = ROOT / "prepared"
JOB = PREPARED / "4B63_A_nadp-keep"
DEFAULT = ROOT / "config" / "default.yaml"
need_prepared = pytest.mark.skipif(not (JOB / "receptor.pdbqt").exists(), reason="run `python -m sidacore.prep` first")
need_vina = pytest.mark.skipif(shutil.which("vina") is None, reason="vina not on PATH")
need_unidock = pytest.mark.skipif(not Path(UNIDOCK_DEFAULT).exists(), reason="Uni-Dock env missing")


@pytest.fixture(scope="module")
def cfg():
    return load_config(DEFAULT)


# ------------------------------------------------------------------ ligand identity comes from the PDB dictionary
def test_template_is_l_ornithine_from_the_component_dictionary():
    t = ligand_template("ORN")
    assert t.GetNumAtoms() == 9 and Chem.MolToSmiles(t) == Chem.MolToSmiles(Chem.MolFromSmiles("N[C@@H](CCCN)C(=O)O"))


def test_unknown_residue_name_is_an_error_not_an_empty_molecule():
    with pytest.raises(ValueError, match="XXXX"):
        ligand_template("XXXX")


# ------------------------------------------------------------------ protonation: the first thing to suspect when a redock fails
@pytest.mark.parametrize("ph, charge", [(7.4, 1), (1.0, 2), (12.0, -1)])
def test_protonation_follows_the_ph(ph, charge):
    """alpha-NH3+ and delta-NH3+ at 7.4 with a carboxylate (net +1); both amines and the acid at pH 1; free amines at 12."""
    smi = protonate(ligand_template("ORN"), ph)
    assert Chem.GetFormalCharge(Chem.MolFromSmiles(smi)) == charge


def test_the_docked_pdbqt_carries_the_protonation_state_not_just_the_smiles():
    """Regression: obabel `-p 7.4 --gen3d -> Meeko` silently delivered the NEUTRAL molecule (14 atoms) to Vina."""
    lp = prepare_ligand("ORN", 7.4)
    assert lp.net_charge == 1 and lp.n_atoms == 15 and lp.n_heavy == 9
    assert sum(1 for l in lp.pdbqt.splitlines() if l.startswith("ATOM")) == 15
    assert lp.partial_charge_sum == pytest.approx(1.0, abs=0.05)


def test_ph_change_is_a_visible_difference_in_the_prepared_ligand():
    a, b = prepare_ligand("ORN", 7.4), prepare_ligand("ORN", 12.0)
    assert (a.net_charge, a.n_atoms) != (b.net_charge, b.n_atoms) and a.smiles != b.smiles


def test_ligand_preparation_is_reproducible():
    assert prepare_ligand("ORN", 7.4).pdbqt == prepare_ligand("ORN", 7.4).pdbqt


# ------------------------------------------------------------------ the crystal reference and the RMSD
@need_prepared
def test_crystal_reference_has_the_deposited_coordinates_and_the_dictionary_bonds():
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    assert ref.GetNumAtoms() == 9
    first = ref.GetConformer().GetAtomPosition(0)
    assert (first.x, first.y, first.z) == pytest.approx((-24.797, -40.544, 11.806), abs=1e-3)
    assert Chem.MolToSmiles(ref) == Chem.MolToSmiles(ligand_template("ORN"))


@need_prepared
def test_crystal_reference_refuses_a_ligand_that_does_not_match_the_template(tmp_path):
    bad = tmp_path / "bad.pdb"
    bad.write_text("\n".join(l for l in (JOB / "ligand_crystal.pdb").read_text().splitlines() if " OXT " not in l) + "\n")
    with pytest.raises(ValueError):
        crystal_reference(bad, ligand_template("ORN"))


@need_prepared
def test_rmsd_of_the_reference_with_itself_is_zero():
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    assert pose_rmsd(ref, ref) == pytest.approx(0.0, abs=1e-6)


@need_prepared
def test_swapping_the_two_carboxylate_oxygens_costs_nothing():
    """Symmetry-corrected: O and OXT are the same chemical group, so their exchange is not a pose error."""
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    swapped = Chem.Mol(ref)
    names = [a.GetPDBResidueInfo().GetName().strip() for a in ref.GetAtoms()]
    i, j = names.index("O"), names.index("OXT")
    conf = swapped.GetConformer()
    pi, pj = conf.GetAtomPosition(i), conf.GetAtomPosition(j)
    conf.SetAtomPosition(i, pj)
    conf.SetAtomPosition(j, pi)
    naive = np.sqrt(np.mean(np.sum((swapped.GetConformer().GetPositions() - ref.GetConformer().GetPositions()) ** 2, axis=1)))
    assert naive > 0.2                                          # without symmetry the swap would be counted as an error
    assert pose_rmsd(swapped, ref) == pytest.approx(0.0, abs=1e-6)


@need_prepared
def test_a_displaced_pose_has_the_displacement_as_rmsd():
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    moved = Chem.Mol(ref)
    c = moved.GetConformer()
    for k in range(moved.GetNumAtoms()):
        p = c.GetAtomPosition(k)
        c.SetAtomPosition(k, (p.x + 3.0, p.y, p.z))
    assert pose_rmsd(moved, ref) == pytest.approx(3.0, abs=1e-4)


# ------------------------------------------------------------------ aggregation (pure)
def _seed(top1, best5, score=-5.0):
    return {"top1_rmsd": top1, "best_top5_rmsd": best5, "top1_score": score}


def test_evaluate_seeds_reports_median_fraction_and_the_registered_verdict():
    r = evaluate_seeds([_seed(1.0, 0.8), _seed(1.5, 1.0), _seed(2.5, 1.2), _seed(0.9, 0.9)], accept=2.0)
    assert r["n_seeds"] == 4 and r["top1_median"] == pytest.approx(1.25) and r["top1_n_within"] == 3
    assert r["passes"] is True and r["best_top5_median"] == pytest.approx(0.95)


def test_evaluate_seeds_fails_when_the_median_exceeds_the_acceptance():
    r = evaluate_seeds([_seed(3.0, 1.0), _seed(2.5, 0.9), _seed(1.0, 0.5)], accept=2.0)
    assert r["passes"] is False and r["top1_n_within"] == 1


def test_evaluate_seeds_rejects_fewer_than_ten_seeds_when_asked_for_the_registered_minimum():
    with pytest.raises(ValueError, match="at least 10"):
        evaluate_seeds([_seed(1.0, 1.0)] * 3, accept=2.0, min_seeds=10)


# ------------------------------------------------------------------ real engines
@need_prepared
@need_vina
def test_vina_redock_returns_ranked_poses_and_a_finite_rmsd(tmp_path, cfg):
    lp = prepare_ligand("ORN", cfg.protonation_ph)
    out = run_vina(JOB / "receptor.pdbqt", lp.pdbqt, JOB / "receptor.box.txt", seed=1, exhaustiveness=8, num_modes=5,
                   out_dir=tmp_path)
    poses = read_poses(out.read_text())
    assert 1 <= len(poses) <= 5 and all(poses[i][0] <= poses[i + 1][0] + 1e-6 for i in range(len(poses) - 1))
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    assert all(np.isfinite(pose_rmsd(m, ref)) for _, m in poses)


@need_prepared
@need_unidock
def test_unidock_redock_returns_poses_for_the_same_inputs(tmp_path, cfg):
    lp = prepare_ligand("ORN", cfg.protonation_ph)
    out = run_unidock(JOB / "receptor.pdbqt", lp.pdbqt, JOB / "receptor.box.txt", seed=1, exhaustiveness=8, num_modes=5,
                      out_dir=tmp_path)
    poses = read_poses(out.read_text())
    assert len(poses) >= 1 and poses[0][0] < 0


# ------------------------------------------------------------------ orchestration
@need_prepared
def test_benchmark_jobs_are_the_primary_then_the_secondary_chains_and_nothing_else(cfg):
    names = [j.name for j in benchmark_jobs(cfg, PREPARED)]
    assert names[0].startswith("4B63") and [n[:6] for n in names[1:]] == ["6X0J_A", "6X0J_B", "6X0J_C", "6X0J_D"]
    assert not any(n.startswith(("6X0H", "6X0I")) for n in names)


@need_prepared
def test_fewer_than_ten_seeds_is_refused_before_any_docking(cfg, tmp_path):
    few = load_config(DEFAULT, overrides={"seeds": [1, 2, 3]})
    with pytest.raises(ValueError, match="at least 10"):
        redock_all(few, PREPARED, tmp_path / "out")
    assert not list((tmp_path / "out").glob("*/")) if (tmp_path / "out").exists() else True


@need_prepared
@need_vina
def test_redock_seeds_gives_one_row_per_configured_seed(tmp_path):
    cfg2 = load_config(DEFAULT, overrides={"seeds": [3, 4], "exhaustiveness": 8})
    lp = prepare_ligand("ORN", cfg2.protonation_ph)
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    rows = redock_seeds(JOB, lp, ref, cfg2, "vina", tmp_path)
    assert [r["seed"] for r in rows] == [3, 4]
    assert all(r["top1_rmsd"] == r["rmsds"][0] and r["best_top5_rmsd"] <= r["top1_rmsd"] + 1e-9 for r in rows)


def test_markdown_lists_the_seeds_and_states_the_registered_rule():
    per = [{"seed": i, "top1_rmsd": 1.0 + i / 10, "top1_score": -5.0, "best_top5_rmsd": 0.9, "n_poses": 20} for i in range(1, 11)]
    ev = evaluate_seeds(per, 2.0, min_seeds=10)
    rep = {"ligand": {"smiles": "X", "ph": 7.4, "net_charge": 1, "n_atoms": 15, "n_heavy": 9, "torsdof": 7, "partial_charge_sum": 1.0},
           "settings": {"exhaustiveness": 32, "num_modes": 20, "seeds": list(range(1, 11)), "acceptance": 2.0,
                        "box_center_atom": "FAD:C4X", "box_size": [22, 22, 22], "protonation_ph": 7.4, "redock_nadp": "keep"},
           "engines": {"vina": "v"}, "results": [{"job": "J", "engine": "vina", "wall_seconds": 1.0, "per_seed": per, **ev}]}
    md = render_markdown(rep, notes="n")
    assert "median over seeds of the top-1 RMSD is <= 2 A" in md and "**PASS**" in md and "s10 2.00" in md and md.rstrip().endswith("n")


# ------------------------------------------------------------------ handing a pose to another program
@need_prepared
def test_pose_pdb_block_carries_the_deposited_atom_names_and_keeps_the_coordinates():
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    block = pose_pdb_block(ref, ref, chain="L")
    lines = [l for l in block.splitlines() if l.startswith("HETATM")]
    assert sorted(l[12:16].strip() for l in lines) == sorted(["N", "CA", "CB", "CG", "CD", "NE", "C", "O", "OXT"])
    assert {l[17:20] for l in lines} == {"ORN"} and {l[21] for l in lines} == {"L"}
    assert not [l for l in block.splitlines() if l.startswith("CONECT")]     # serials would collide with a receptor file
    xyz = sorted((float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in lines)
    want = sorted(map(tuple, np.round(ref.GetConformer().GetPositions(), 3)))
    np.testing.assert_allclose(xyz, want, atol=1e-3)


@need_prepared
@need_vina
def test_a_docked_pose_gets_the_same_names_as_the_deposited_ligand(tmp_path, cfg):
    lp = prepare_ligand("ORN", cfg.protonation_ph)
    out = run_vina(JOB / "receptor.pdbqt", lp.pdbqt, JOB / "receptor.box.txt", seed=1, exhaustiveness=8, num_modes=3, out_dir=tmp_path)
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    pose = read_poses(out.read_text())[0][1]
    block = pose_pdb_block(pose, ref)
    back = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False, proximityBonding=True)
    assert back.GetNumAtoms() == 9
    # naming must not change the geometry: RMSD of the re-read pose to the deposited one equals the original number
    assert pose_rmsd(back, ref) == pytest.approx(pose_rmsd(pose, ref), abs=1e-2)


def test_render_only_rebuilds_the_report_with_notes_and_does_not_dock(tmp_path):
    import json
    per = [{"seed": i, "top1_rmsd": 1.0, "top1_score": -5.0, "best_top5_rmsd": 0.9, "n_poses": 20} for i in range(1, 11)]
    ev = evaluate_seeds(per, 2.0, min_seeds=10)
    rep = {"ligand": {"smiles": "X", "ph": 7.4, "net_charge": 1, "n_atoms": 15, "n_heavy": 9, "torsdof": 7, "partial_charge_sum": 1.0},
           "settings": {"exhaustiveness": 32, "num_modes": 20, "seeds": list(range(1, 11)), "acceptance": 2.0,
                        "box_center_atom": "FAD:C4X", "box_size": [22, 22, 22], "protonation_ph": 7.4, "redock_nadp": "keep"},
           "engines": {"vina": "v"}, "results": [{"job": "J", "engine": "vina", "wall_seconds": 1.0, "per_seed": per, **ev}]}
    (tmp_path / "step4_report.json").write_text(json.dumps(rep))
    (tmp_path / "notes.md").write_text("hand-written note")
    assert redock_main(["--render-only", str(tmp_path), "--report", str(tmp_path / "r.md"), "--notes", str(tmp_path / "notes.md")]) == 0
    assert "hand-written note" in (tmp_path / "r.md").read_text()
