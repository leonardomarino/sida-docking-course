"""confirm/summarise_redock.py: the registered rule applied to the RosettaLigand decoys, and the native-start comparison."""
import pytest

from confirm.summarise_redock import render_markdown, summarise_job


def _row(seed, idx, total, iface, err=""):
    return {"decoy": "dock_s%d_d%02d" % (seed, idx), "seed": seed, "index": idx, "total_score": total,
            "interface_delta": iface, "error": err}


def _data(good_first=True):
    """10 seeds x 3 decoys. In every seed the lowest total score is the near-native decoy (d01, 1.0 A) when good_first."""
    rows, rmsd = [], {}
    for s in range(1, 11):
        for d, (tot, rm) in enumerate([(-100.0, 1.0), (-90.0, 6.0), (-80.0, 7.0)], 1):
            if not good_first:
                tot, rm = (-100.0 - d, 6.0 + d) if d == 1 else (tot, 1.0 if d == 3 else rm)
            r = _row(s, d, tot, tot / 10)
            rows.append(r)
            rmsd[r["decoy"]] = rm
    return rows, rmsd


NATIVE = [{"decoy": "native_s0_d%02d" % i, "seed": 0, "index": i, "total_score": -110.0 - i, "interface_delta": -12.0, "error": ""} for i in range(1, 4)]


def test_the_registered_rule_passes_when_the_lowest_scoring_decoy_of_most_seeds_is_near_native():
    rows, rmsd = _data(True)
    s = summarise_job(rows, NATIVE, rmsd, accept=2.0)
    assert s["by_total_score"]["passes"] is True and s["by_total_score"]["top1_median"] == 1.0 and s["by_total_score"]["top1_n_within"] == 10


def test_it_fails_when_the_lowest_scoring_decoy_is_wrong_even_though_a_near_native_decoy_exists():
    rows, rmsd = _data(False)
    s = summarise_job(rows, NATIVE, rmsd, accept=2.0)
    assert s["by_total_score"]["passes"] is False and s["by_total_score"]["top1_median"] > 2.0
    assert s["near_native_fraction"] == pytest.approx(1 / 3)            # the search did find it, the score did not pick it


def test_native_start_comparison_separates_search_from_scoring():
    rows, rmsd = _data(False)
    s = summarise_job(rows, NATIVE, rmsd, accept=2.0)
    assert s["native"]["best_total_score"] == -113.0 and s["native"]["lower_than_best_docked_total"] is True
    assert s["native"]["margin_to_best_docked"] == pytest.approx(113.0 - 101.0)       # best docked is -101, best native -113


def test_failed_decoys_are_counted_and_excluded_never_ranked():
    rows, rmsd = _data(True)
    rows[0] = _row(1, 1, float("nan"), float("nan"), err="NAN occurred in H-bonding")
    s = summarise_job(rows, NATIVE, rmsd, accept=2.0)
    assert s["n_failed"] == 1 and s["by_total_score"]["n_seeds"] == 10


def test_a_seed_whose_decoys_all_failed_is_an_error():
    rows, rmsd = _data(True)
    for r in rows:
        if r["seed"] == 3:
            r.update(error="boom", total_score=float("nan"), interface_delta=float("nan"))
    with pytest.raises(ValueError, match="seed 3"):
        summarise_job(rows, NATIVE, rmsd, accept=2.0)


def test_markdown_states_the_registered_rule_and_the_verdict():
    rows, rmsd = _data(True)
    md = render_markdown({"J": summarise_job(rows, NATIVE, rmsd, accept=2.0)}, notes="n")
    assert "lowest total score" in md and "**PASS**" in md and "| J |" in md and md.rstrip().endswith("n")
