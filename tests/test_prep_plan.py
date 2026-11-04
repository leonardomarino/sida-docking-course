"""Which receptors Steps 3-5 need, derived from the config (no hardcoded list)."""
from pathlib import Path

from sidacore.config import load_config
from sidacore.prep import job_name, plan_jobs

DEFAULT = Path(__file__).resolve().parent.parent / "config" / "default.yaml"
FLAVIN_CHAINS = {"4B63": ["A"], "6X0J": ["A", "B", "C", "D"], "6X0H": ["A", "B", "C", "D"], "6X0I": ["A", "B", "C", "D"]}


def test_default_plan_covers_primary_secondary_chains_and_the_three_arms():
    jobs = plan_jobs(load_config(DEFAULT), FLAVIN_CHAINS)
    assert jobs == [("4B63", "A", "keep"),
                    ("6X0J", "A", "keep"), ("6X0J", "B", "keep"), ("6X0J", "C", "keep"), ("6X0J", "D", "keep"),
                    ("6X0H", "A", "keep"), ("6X0I", "A", "keep"), ("6X0I", "A", "strip")]


def test_duplicates_collapse_and_nadp_policy_follows_the_config():
    cfg = load_config(DEFAULT, overrides={"redock_nadp": "strip"})
    jobs = plan_jobs(cfg, FLAVIN_CHAINS)
    assert ("4B63", "A", "strip") in jobs and len(jobs) == len(set(jobs))


def test_job_name_is_stable_and_readable():
    assert job_name("6X0I", "A", "strip") == "6X0I_A_nadp-strip"
