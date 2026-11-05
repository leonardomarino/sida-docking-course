"""confirm/make_inputs.py: the licence-free half of the Rosetta check (RDKit and sidacore only; no PyRosetta)."""
import csv
import json
from pathlib import Path

import pytest

from confirm.make_inputs import make_inputs, make_start_ligands
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


# ------------------------------------------------------------------ start ligands for the Rosetta redock
@pytest.fixture(scope="module")
def starts(tmp_path_factory):
    out = tmp_path_factory.mktemp("starts")
    cfg = load_config(ROOT / "config" / "default.yaml")
    return out, make_start_ligands(cfg, PREPARED, out, seeds=(1, 2, 3), only=[JOB])


def _xyz(path):
    return [(l[12:16].strip(), float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in Path(path).read_text().splitlines() if l.startswith("HETATM")]


def test_one_named_start_ligand_per_seed_centred_on_the_box(starts):
    out, rows = starts
    box = [float(x.split("=")[1]) for x in (PREPARED / JOB / "receptor.box.txt").read_text().splitlines()[:3]]
    assert [r["seed"] for r in rows] == [1, 2, 3]
    for r in rows:
        atoms = _xyz(out / JOB / "starts" / ("start_s%d.pdb" % r["seed"]))
        assert sorted(a[0] for a in atoms) == sorted(["N", "CA", "CB", "CG", "CD", "NE", "C", "O", "OXT"])
        c = [sum(a[i] for a in atoms) / 9 for i in (1, 2, 3)]
        assert max(abs(c[i] - box[i]) for i in range(3)) < 0.01


def test_no_start_ligand_is_near_the_deposited_pose(starts):
    """The start must not leak the answer: every start is far (symmetry-corrected RMSD) from the crystal pose."""
    _, rows = starts
    assert all(r["rmsd_to_crystal"] > 3.0 for r in rows)


def test_the_seeds_give_more_than_one_internal_geometry(starts):
    """Compared by interatomic distances, which do not depend on orientation, so this tests the conformers themselves."""
    import numpy as np
    out, _ = starts
    def dm(seed):
        xyz = np.array([[x, y, z] for _, x, y, z in sorted(_xyz(out / JOB / "starts" / ("start_s%d.pdb" % seed)))])
        return np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
    d = [dm(s) for s in (1, 2, 3)]
    assert max(np.abs(d[i] - d[j]).max() for i in range(3) for j in range(i + 1, 3)) > 0.1
