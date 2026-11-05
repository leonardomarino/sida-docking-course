"""confirm/make_inputs.py: the licence-free half of the Rosetta check (RDKit and sidacore only; no PyRosetta)."""
import csv
import json
from pathlib import Path

import pytest

from confirm.make_inputs import make_inputs
from sidacore.config import load_config

ROOT = Path(__file__).resolve().parent.parent
PREPARED, REDOCK = ROOT / "prepared", ROOT / "redock"
JOB = "4B63_A_nadp-keep"
pytestmark = pytest.mark.skipif(not (REDOCK / JOB / "vina" / "vina_seed1.pdbqt").exists(), reason="run `python -m sidacore.redock` first")


@pytest.fixture(scope="module")
def made(tmp_path_factory):
    out = tmp_path_factory.mktemp("rin")
    cfg = load_config(ROOT / "config" / "default.yaml")
    return out, make_inputs(cfg, PREPARED, REDOCK, out, engines=("vina",), seed=1, only=[JOB])


def test_receptor_is_heavy_atoms_with_original_cofactor_names(made):
    out, _ = made
    lines = (out / JOB / "receptor.pdb").read_text().splitlines()
    assert lines and all(l[76:78].strip() != "H" for l in lines if l.startswith(("ATOM", "HETATM")))
    assert {l[17:20] for l in lines if l.startswith("HETATM")} >= {"FAD", "NAP"}


def test_every_ligand_file_is_nine_named_atoms_on_chain_L_with_no_conect(made):
    out, rows = made
    for r in rows:
        txt = (out / JOB / "poses" / r["pose"]).with_suffix(".pdb").read_text()
        het = [l for l in txt.splitlines() if l.startswith("HETATM")]
        assert len(het) == 9 and {l[21] for l in het} == {"L"} and "CONECT" not in txt


def test_manifest_has_the_crystal_pose_at_zero_and_the_docked_modes_ranked(made):
    _, rows = made
    crystal = [r for r in rows if r["pose"] == "crystal"]
    assert len(crystal) == 1 and crystal[0]["rmsd_to_crystal"] == pytest.approx(0.0, abs=1e-6)
    vina = [r for r in rows if r["engine"] == "vina"]
    assert [r["rank"] for r in vina] == list(range(1, len(vina) + 1))
    assert all(vina[i]["score"] <= vina[i + 1]["score"] + 1e-6 for i in range(len(vina) - 1))


def test_manifest_top1_rmsd_agrees_with_the_redock_report(made):
    _, rows = made
    rep = json.loads((REDOCK / "step4_report.json").read_text())
    seed1 = next(r for r in rep["results"] if r["job"] == JOB and r["engine"] == "vina")["per_seed"][0]
    top1 = next(r for r in rows if r["engine"] == "vina" and r["rank"] == 1)
    assert top1["rmsd_to_crystal"] == pytest.approx(seed1["top1_rmsd"], abs=1e-3)
