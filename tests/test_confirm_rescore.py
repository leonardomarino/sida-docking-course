"""confirm/rescore.py: the registered Step 4b rescoring rule (reports/step4_rescoring_rule.md), pure logic plus command construction."""
import os
import stat
from pathlib import Path

import numpy as np
import pytest

from confirm.rescore import (main, _read_box_center, logged_interaction_energies, outcome_statement, render_markdown, AD4_LIGAND_TYPES, AD4_RECEPTOR_TYPES, docked_scores, g0_check, headline_pose, is_admissible,
                             outcome_key, outside_map, parse_estimated_energy, pdbqt_coords, run_vina_score, seed_result,
                             split_models, summarise_job, write_gpf)

TWO_MODELS = """MODEL 1
REMARK VINA RESULT:    -4.922      0.000      0.000
REMARK SMILES [NH3+]CCC
ROOT
ATOM      1  C   UNL     1     -23.561 -39.042  11.668  1.00  0.00    +0.068 C
ENDROOT
BRANCH   1   2
ATOM      2  HD  UNL     1     -24.156 -39.546  12.984  1.00  0.00    +0.259 HD
ENDBRANCH   1   2
TORSDOF 1
ENDMDL
MODEL 2
REMARK VINA RESULT:    -4.761      1.606      1.766
ROOT
ATOM      1  C   UNL     1     -20.000 -30.000  10.000  1.00  0.00    +0.068 C
ENDROOT
TORSDOF 0
ENDMDL
"""

VINA_OUT = """Scoring function : vina
Estimated Free Energy of Binding   : -4.924 (kcal/mol) [=(1)+(2)+(3)-(4)]
(1) Final Intermolecular Energy    : -6.363 (kcal/mol)
"""


# ----------------------------------------------------------------------------- pose files
def test_split_models_gives_single_model_files_without_model_tags():
    parts = split_models(TWO_MODELS)
    assert len(parts) == 2
    assert all("MODEL" not in p for p in parts)                      # Vina rejects a MODEL tag in a ligand file
    assert "-4.922" in parts[0] and "-4.761" in parts[1] and "-20.000" in parts[1]
    assert "-20.000" not in parts[0]


def test_docked_scores_in_file_order():
    assert docked_scores(TWO_MODELS) == [-4.922, -4.761]


def test_pdbqt_coords_include_polar_hydrogens_and_skip_remarks():
    xyz = pdbqt_coords(TWO_MODELS)
    assert xyz.shape == (3, 3)
    assert xyz[1].tolist() == [-24.156, -39.546, 12.984]


# ----------------------------------------------------------------------------- scores
def test_parse_estimated_energy():
    assert parse_estimated_energy(VINA_OUT) == -4.924


def test_parse_estimated_energy_refuses_output_without_the_line():
    with pytest.raises(ValueError):
        parse_estimated_energy("PDBQT parsing error: Unexpected multi-MODEL tag")


def _fake_vina(tmp_path, body="echo 'Estimated Free Energy of Binding   : -1.500 (kcal/mol)'"):
    exe = tmp_path / "vina"
    exe.write_text("#!/bin/bash\necho \"$@\" > %s/args.txt\n%s\n" % (tmp_path, body))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    return exe


def test_vina_and_vinardo_use_the_receptor_and_autobox(tmp_path):
    e = run_vina_score(_fake_vina(tmp_path), "vinardo", "lig.pdbqt", receptor="rec.pdbqt")
    assert e == -1.5
    args = (tmp_path / "args.txt").read_text().split()
    assert args[args.index("--scoring") + 1] == "vinardo" and "--score_only" in args and "--autobox" in args
    assert args[args.index("--receptor") + 1] == "rec.pdbqt" and "--maps" not in args


def test_ad4_uses_maps_and_never_a_receptor(tmp_path):
    run_vina_score(_fake_vina(tmp_path), "ad4", "lig.pdbqt", maps="m/rec")
    args = (tmp_path / "args.txt").read_text().split()
    assert args[args.index("--maps") + 1] == "m/rec" and "--receptor" not in args           # Vina refuses --receptor with ad4
    assert args[args.index("--scoring") + 1] == "ad4"


def test_a_vina_failure_is_an_error_not_a_score(tmp_path):
    with pytest.raises(RuntimeError):
        run_vina_score(_fake_vina(tmp_path, "echo 'PDBQT parsing error' >&2; exit 1"), "vina", "l.pdbqt", receptor="r.pdbqt")


def test_scoring_function_is_validated(tmp_path):
    with pytest.raises(ValueError):
        run_vina_score(_fake_vina(tmp_path), "gnina", "l.pdbqt", receptor="r.pdbqt")


# ----------------------------------------------------------------------------- AutoGrid input
def test_box_centre_is_read_from_the_prepared_box_file(tmp_path):
    f = tmp_path / "receptor.box.txt"
    f.write_text("center_x = -25.372\ncenter_y = -31.449\ncenter_z = 14.541\nsize_x = 22.000\nsize_y = 22.000\nsize_z = 22.000\n")
    assert _read_box_center(f) == (-25.372, -31.449, 14.541)


def test_gpf_carries_exactly_the_registered_values():
    t = write_gpf("rec.pdbqt", "rec", (-25.372, -31.449, 14.541), npts=60)
    lines = [l.split("#")[0].strip() for l in t.splitlines()]
    for want in ("npts 60 60 60", "spacing 0.375", "gridcenter -25.372 -31.449 14.541", "smooth 0.5", "dielectric -0.1465",
                 "receptor rec.pdbqt", "gridfld rec.maps.fld", "receptor_types A C HD N NA OA P SA", "ligand_types C HD N OA",
                 "elecmap rec.e.map", "dsolvmap rec.d.map", "map rec.C.map", "map rec.HD.map", "map rec.N.map", "map rec.OA.map"):
        assert want in lines, want
    assert AD4_RECEPTOR_TYPES == ("A", "C", "HD", "N", "NA", "OA", "P", "SA") and AD4_LIGAND_TYPES == ("C", "HD", "N", "OA")


def test_gpf_fallback_grid_is_80_points():
    assert "npts 80 80 80" in write_gpf("r.pdbqt", "r", (0, 0, 0), npts=80)


def test_out_of_map_rule_uses_the_grid_extent():
    c = (10.0, 20.0, 30.0)
    inside = np.array([[10, 20, 30], [21.25, 20, 30], [10, 8.75, 41.25]], float)           # 60 * 0.375 / 2 = 11.25
    assert outside_map(inside, c, npts=60) is False
    assert outside_map(np.array([[21.26, 20, 30]]), c, npts=60) is True
    assert outside_map(np.array([[21.26, 20, 30]]), c, npts=80) is False                # 80 points reach 15 A


# ----------------------------------------------------------------------------- ranking
def test_clear_minimum_is_the_headline_pose():
    assert headline_pose([-4.0, -5.0, -4.5]) == (1, False)


def test_a_tie_within_0_05_goes_to_the_pose_vina_ranked_higher():
    assert headline_pose([-4.00, -5.00, -5.04]) == (1, True)                                  # pose 2 is better by 0.04: tie, lower index wins
    assert headline_pose([-5.04, -4.00, -5.00]) == (0, True)
    assert headline_pose([-5.00, -4.00, -5.06]) == (2, False)                                 # 0.06 apart: a real difference


def test_exactly_0_05_apart_counts_as_a_tie():
    assert headline_pose([-5.00, -5.05]) == (0, True)


def test_seed_result_reports_rank_of_best_pose_and_ceiling():
    r = seed_result(rescored=[-5.0, -4.2, -4.8, -4.0], rmsds=[4.8, 1.5, 3.0, 6.0], accept=2.0)
    assert r["top1_idx"] == 0 and r["top1_rmsd"] == 4.8 and r["tied"] is False
    assert r["best_idx"] == 1 and r["best_rank"] == 3 and r["best_rmsd"] == 1.5           # two poses score better than the best one
    assert r["has_near_native"] is True


def test_seed_result_without_any_near_native_pose():
    r = seed_result(rescored=[-5.0, -4.0], rmsds=[4.8, 3.0], accept=2.0)
    assert r["has_near_native"] is False


def test_summarise_job_applies_the_step_4_rule():
    mk = lambda rm: seed_result([-5.0, -4.0], [rm, 9.0], accept=2.0)
    s = summarise_job([mk(1.0), mk(1.5), mk(5.0)], accept=2.0)
    assert s["top1_median"] == 1.5 and s["passes"] is True and s["n_within"] == 2 and s["n_seeds"] == 3
    s = summarise_job([mk(1.0), mk(4.0), mk(5.0)], accept=2.0)
    assert s["top1_median"] == 4.0 and s["passes"] is False
    assert s["ceiling_n"] == 1 and s["n_tied"] == 0                  # only the 1.0 A seed has any pose within 2 A


# ----------------------------------------------------------------------------- gates and outcomes
def test_g0_passes_when_every_pose_reproduces_the_log():
    ok = g0_check([-4.922, -4.761, -4.5], [-4.924, -4.760, -4.5])
    assert ok["ok"] is True and ok["n_bad"] == 0 and ok["max_abs_diff"] == pytest.approx(0.002)


def test_g0_fails_on_one_pose_beyond_0_05():
    bad = g0_check([-4.9, -4.8], [-4.9, -4.7])
    assert bad["ok"] is False and bad["n_bad"] == 1 and bad["max_abs_diff"] == pytest.approx(0.1)


def test_g0_refuses_mismatched_lengths():
    with pytest.raises(ValueError):
        g0_check([-4.9], [-4.9, -4.8])


def test_function_is_admissible_only_if_all_three_controls_pass():
    ok = {"6X0J_A": True, "6X0J_B": True, "6X0J_D": True}
    assert is_admissible(ok)["admissible"] is True
    broken = dict(ok, **{"6X0J_B": False})
    r = is_admissible(broken)
    assert r["admissible"] is False and r["broken"] == ["6X0J_B"]


def test_admissibility_demands_all_three_controls_be_present():
    with pytest.raises(KeyError):
        is_admissible({"6X0J_A": True, "6X0J_B": True})


@pytest.mark.parametrize("admissible,passes,key", [(True, True, "rescues"), (True, False, "no_rescue"),
                                                    (False, True, "cannot_arbitrate"), (False, False, "cannot_arbitrate")])
def test_outcome_table_is_the_registered_one(admissible, passes, key):
    assert outcome_key(admissible, passes) == key



# ----------------------------------------------------------------------------- Amendment 1: interaction energy
REMARKS = """MODEL 1
REMARK VINA RESULT:    -4.922      0.000      0.000
REMARK INTER + INTRA:          -6.761
REMARK INTER:                  -6.361
REMARK INTRA:                  -0.399
REMARK UNBOUND:                -0.399
ENDMDL
MODEL 2
REMARK VINA RESULT:    -4.761      1.606      1.766
REMARK INTER + INTRA:          -6.552
REMARK INTER:                  -6.447
REMARK INTRA:                  -0.105
REMARK UNBOUND:                -0.399
ENDMDL
"""


def test_logged_interaction_energy_drops_the_internal_term():
    e = logged_interaction_energies(REMARKS)
    assert e[0] == pytest.approx(-4.922, abs=0.002)                 # pose 1: INTRA == UNBOUND, so equal to the logged score
    assert e[1] == pytest.approx(-4.989, abs=0.002)                 # pose 2: the log says -4.761, interaction only is -4.989 (what --score_only returns)
    assert docked_scores(REMARKS)[1] == -4.761


def test_logged_interaction_energy_does_not_confuse_the_combined_line_with_inter_or_intra():
    assert logged_interaction_energies(REMARKS.replace("-6.761", "-9.999"))[0] == pytest.approx(-4.922, abs=0.002)


def test_logged_interaction_energy_refuses_a_file_whose_first_pose_is_not_its_own_reference():
    with pytest.raises(ValueError):
        logged_interaction_energies(REMARKS.replace("REMARK UNBOUND:                -0.399\nENDMDL\nMODEL 2", "REMARK UNBOUND:                -0.250\nENDMDL\nMODEL 2", 1))


def test_logged_interaction_energy_refuses_a_file_without_remarks():
    with pytest.raises(ValueError):
        logged_interaction_energies(TWO_MODELS)


def test_baseline_qualification_is_added_to_f1_and_f2_only_when_vina_interaction_only_already_passes():
    plain = outcome_statement("ad4", "rescues", moot=False)
    moot = outcome_statement("ad4", "rescues", moot=True)
    assert "internal-strain" not in plain and "accounts for it" in moot
    assert "internal-strain" not in outcome_statement("vina", "rescues", moot=True)          # the baseline itself is not qualified


def test_report_states_the_gate_the_baseline_and_every_outcome():
    summ = {"4B63_A": {"top1_median": 4.8, "top1_min": 4.7, "top1_max": 4.9, "n_within": 0, "n_seeds": 10, "passes": False, "n_tied": 1,
                       "ceiling_n": 10, "best_rank_median": 4, "accept": 2.0}}
    fn = {"jobs": summ, "admissible": {"admissible": False, "broken": ["6X0J_B"]}, "outcome": "cannot_arbitrate",
          "statement": "AutoDock4 scoring cannot arbitrate this site: it breaks control(s) 6X0J_B."}
    res = {"g0": {"ok": True, "n": 1000, "n_bad": 0, "max_abs_diff": 0.002, "tol": 0.05},
           "grid": {"npts": 60, "n_poses_outside_60": {"4B63_A": 0}, "n_poses_outside_final": {"4B63_A": 0}},
           "functions": {"vina": fn, "ad4": fn, "vinardo": fn}}
    md = render_markdown(res, notes="NOTE")
    assert "Vina, interaction energy only" in md and "AutoDock4" in md and "Vinardo" in md and "cannot arbitrate" in md
    assert "0.002" in md and "| 4B63_A |" in md and md.rstrip().endswith("NOTE")

def test_render_only_rewrites_the_report_from_saved_results_without_scoring(tmp_path):
    import json
    fn = {"jobs": {"4B63_A": {"top1_median": 4.9, "top1_min": 4.8, "top1_max": 5.0, "n_within": 0, "n_seeds": 10, "passes": False, "n_tied": 0,
                              "ceiling_n": 10, "best_rank_median": 3, "accept": 2.0}},
          "admissible": {"admissible": True, "broken": []}, "outcome": "no_rescue", "statement": "X does not rescue 4B63 either."}
    (tmp_path / "result.json").write_text(json.dumps({"g0": {"ok": True, "n": 1, "n_bad": 0, "max_abs_diff": 0.0, "tol": 0.05},
                                                      "functions": {"vina": fn, "ad4": fn, "vinardo": fn}}))
    (tmp_path / "n.md").write_text("HAND-WRITTEN\n")
    assert main(["--out", str(tmp_path), "--render-only", "--report", str(tmp_path / "r.md"), "--notes", str(tmp_path / "n.md")]) == 0
    assert (tmp_path / "r.md").read_text().rstrip().endswith("HAND-WRITTEN")


# ----------------------------------------------------------------------------- real data (skipped where absent)
COURSE = Path(__file__).resolve().parent.parent
POSES = COURSE / "redock/4B63_A_nadp-keep/vina/vina_seed1.pdbqt"
RECEPTOR = COURSE / "prepared/4B63_A_nadp-keep/receptor.pdbqt"


def _vina():
    from sidacore.boxes import _tool
    try:
        return _tool("vina")
    except FileNotFoundError:
        return None


@pytest.mark.skipif(not (POSES.exists() and RECEPTOR.exists()) or _vina() is None, reason="needs the 4B63 redock and vina")
def test_vina_function_reproduces_the_logged_interaction_energy_on_real_poses(tmp_path):
    text = POSES.read_text()
    logged = logged_interaction_energies(text)[:3]
    got = []
    for i, part in enumerate(split_models(text)[:3]):
        f = tmp_path / ("p%d.pdbqt" % i)
        f.write_text(part)
        got.append(run_vina_score(_vina(), "vina", f, receptor=RECEPTOR))
    assert g0_check(logged, got)["ok"] is True
