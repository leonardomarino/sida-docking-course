"""confirm/summarise.py: what a second scoring function says about the docked poses (pure logic on the manifest + scores)."""
import pytest

from confirm.summarise import render_markdown, summarise_job


def _m(pose, engine, rank, score, rmsd):
    return {"job": "J", "engine": engine, "rank": rank, "score": score, "rmsd_to_crystal": rmsd, "pose": pose}


MANIFEST = [_m("crystal", "crystal", 0, "", 0.0),
            _m("v1", "vina", 1, -5.0, 4.8), _m("v2", "vina", 2, -4.9, 1.2), _m("v3", "vina", 3, -4.7, 3.0), _m("v4", "vina", 4, -4.5, 6.0)]
ROSETTA = {"crystal": -13.0, "v1": -5.0, "v2": -12.0, "v3": -8.0, "v4": 0.0}      # dG_separated, lower is better


def test_rosetta_can_prefer_a_pose_the_docking_score_ranked_second():
    s = summarise_job(MANIFEST, ROSETTA)["vina"]
    assert s["engine_top1_rmsd"] == 4.8 and s["rosetta_top1_rank_in_engine"] == 2 and s["rosetta_top1_rmsd"] == 1.2


def test_crystal_pose_versus_the_engines_top_pose_separates_search_from_scoring():
    s = summarise_job(MANIFEST, ROSETTA)["vina"]
    assert s["crystal_dG"] == -13.0 and s["engine_top1_dG"] == -5.0
    assert s["crystal_better_than_engine_top1"] is True and s["crystal_rank_among_all"] == 1


def test_best_rmsd_pose_and_rank_correlation():
    s = summarise_job(MANIFEST, ROSETTA)["vina"]
    assert s["best_rmsd_pose_rosetta_rank"] == 1 and s["best_rmsd"] == 1.2
    assert s["spearman_score_vs_dG"] == pytest.approx(0.4, abs=0.01)     # checked independently with scipy: 0.4


def test_refuses_a_pose_without_a_rosetta_score():
    with pytest.raises(KeyError):
        summarise_job(MANIFEST, {"crystal": -13.0})


def test_markdown_states_the_protocol_and_every_receptor():
    md = render_markdown({"J": summarise_job(MANIFEST, ROSETTA)}, notes="n")
    assert "ref2015" in md and "| J |" in md and md.rstrip().endswith("n")


def test_a_pose_with_no_interface_cannot_win_the_rosetta_ranking():
    """dG is 0 when nothing touches the receptor, which would beat any real binder whose dG is positive."""
    man = MANIFEST + [_m("v5", "vina", 5, -4.0, 15.0)]
    dG = {"crystal": 1.0, "v1": 13.0, "v2": 12.0, "v3": 14.0, "v4": 15.0, "v5": 0.0}
    contact = {"crystal": 17, "v1": 21, "v2": 20, "v3": 19, "v4": 18, "v5": 0}
    s = summarise_job(man, dG, contact=contact)["vina"]
    assert s["n_no_contact"] == 1 and s["rosetta_top1_pose"] == "v2" and s["rosetta_top1_rmsd"] == 1.2


def test_total_energy_view_is_reported_separately_and_compares_with_the_crystal_pose():
    total = {"crystal": -337.0, "v1": -325.0, "v2": -334.0, "v3": -320.0, "v4": -300.0}
    s = summarise_job(MANIFEST, ROSETTA, total=total)["vina"]["total"]
    assert s["crystal_better_than_engine_top1"] is True and s["engine_top1_minus_crystal"] == pytest.approx(12.0)
    assert s["rosetta_top1_pose"] == "v2" and s["rosetta_top1_rmsd"] == 1.2


def test_markdown_shows_the_pose_as_docked_and_after_minimisation_and_never_a_rank_zero():
    s = summarise_job(MANIFEST, ROSETTA, min_rmsd={"v2": 2.1})
    s["vina"]["best_rmsd_pose_rosetta_rank"] = None            # the best-RMSD pose was excluded for having no contact
    md = render_markdown({"J": s})
    assert "after minimisation" in md and "1.20 / 2.10" in md and "[0]" not in md
