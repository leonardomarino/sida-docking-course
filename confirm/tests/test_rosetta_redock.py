"""Rosetta (RosettaLigand) redock. Runs in the PyRosetta env (academic licence); skipped elsewhere. Not part of sidacore."""
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pyrosetta")

from confirm.rosetta_redock import dock_decoy, protocol_xml  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RIN = ROOT / "rosetta_in" / "4B63_A_nadp-keep"
pytestmark = pytest.mark.skipif(not (RIN / "receptor.pdb").exists(), reason="run `python -m confirm.make_inputs` first")
FAST = dict(box_size=11.0, transform_cycles=50, transform_repeats=1, hrd_cycles=2)


def _heavy(pdb_text):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in pdb_text.splitlines()
                     if l.startswith("HETATM") and l[76:78].strip() != "H"])


def _by_name(pdb_text):
    """{atom name: xyz} of the heavy atoms. Rosetta writes the residue's atoms in its own order, so compare by name."""
    return {l[12:16].strip(): np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) for l in pdb_text.splitlines()
            if l.startswith("HETATM") and l[76:78].strip() != "H"}


def _shifted_crystal(tmp_path, dx):
    out = []
    for l in (RIN / "poses" / "crystal.pdb").read_text().splitlines():
        if l.startswith("HETATM"):
            l = l[:30] + "%8.3f" % (float(l[30:38]) + dx) + l[38:]
            out.append(l)
    p = tmp_path / ("start_%g.pdb" % dx)
    p.write_text("\n".join(out) + "\nEND\n")
    return p


def test_protocol_xml_carries_the_registered_settings_and_native_mode_has_no_transform():
    dock, native = protocol_xml(mode="dock", **FAST), protocol_xml(mode="native", **FAST)
    for needed in ("<Transform", "<HighResDocker", "<FinalMinimizer", "<InterfaceScoreCalculator", "<SCORINGGRIDS", 'box_size="11.0"'):
        assert needed in dock
    assert '<Add mover_name="transform"/>' in dock and '<Add mover_name="transform"/>' not in native and "<HighResDocker" in native
    # the interface energy is evaluated with the coordinate-constraint term switched off (it returned -2.8 million otherwise)
    assert 'coordinate_constraint" weight="0"' in dock


def test_a_decoy_has_finite_sane_scores_and_a_ligand_file(tmp_path):
    r = dock_decoy(RIN / "receptor.pdb", _shifted_crystal(tmp_path, 4.0), seed=1, out_ligand_pdb=tmp_path / "d.pdb", **FAST)
    assert np.isfinite([r["total_score"], r["interface_delta"]]).all()
    assert abs(r["interface_delta"]) < 1000                      # regression: a constraint term once made it -2,771,178
    xyz = _heavy((tmp_path / "d.pdb").read_text())
    assert len(xyz) == 9


def test_the_same_seed_gives_the_same_decoy_and_a_different_seed_a_different_one(tmp_path):
    s = _shifted_crystal(tmp_path, 4.0)
    a = dock_decoy(RIN / "receptor.pdb", s, seed=7, out_ligand_pdb=tmp_path / "a.pdb", **FAST)
    b = dock_decoy(RIN / "receptor.pdb", s, seed=7, out_ligand_pdb=tmp_path / "b.pdb", **FAST)
    c = dock_decoy(RIN / "receptor.pdb", s, seed=8, out_ligand_pdb=tmp_path / "c.pdb", **FAST)
    assert a["total_score"] == pytest.approx(b["total_score"], abs=1e-6)
    assert not np.allclose(_heavy((tmp_path / "a.pdb").read_text()), _heavy((tmp_path / "c.pdb").read_text()), atol=1e-3)


def test_native_mode_refines_in_place_and_does_not_search(tmp_path):
    start = RIN / "poses" / "crystal.pdb"
    r = dock_decoy(RIN / "receptor.pdb", start, seed=1, out_ligand_pdb=tmp_path / "n.pdb", mode="native", **FAST)
    a, b = _by_name(start.read_text()), _by_name((tmp_path / "n.pdb").read_text())
    assert set(a) == set(b) and len(b) == 9
    assert np.sqrt(np.mean([np.sum((a[k] - b[k]) ** 2) for k in a])) < 2.0
    assert np.isfinite(r["total_score"])
