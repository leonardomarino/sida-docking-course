"""Rosetta confirmation of the docking poses. Runs in the PyRosetta env (academic licence); skipped elsewhere.

Not part of `sidacore` and never imported by it.
"""
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pyrosetta")

from confirm.rosetta_rescore import rescore, rescore_many  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
JOB = ROOT / "prepared" / "4B63_A_nadp-keep"
pytestmark = pytest.mark.skipif(not (JOB / "work" / "combined.pdb").exists(), reason="run `python -m sidacore.prep` first")


def _receptor(tmp_path):
    """Protein + FAD + NADP heavy atoms with their original residue names; Rosetta places its own hydrogens."""
    rec = tmp_path / "receptor.pdb"
    rec.write_text("".join(l for l in (JOB / "work" / "combined.pdb").read_text().splitlines(True)
                           if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() != "H"))
    return rec


def _ligand(tmp_path, shift=(0.0, 0.0, 0.0), name="lig.pdb"):
    out = []
    for l in (JOB / "ligand_crystal.pdb").read_text().splitlines():
        if l.startswith("HETATM"):
            x, y, z = (float(l[30:38]) + shift[0], float(l[38:46]) + shift[1], float(l[46:54]) + shift[2])
            l = l[:21] + "L" + l[22:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + l[54:]
            out.append(l)
    p = tmp_path / name
    p.write_text("\n".join(out) + "\nEND\n")
    return p


def test_crystal_pose_is_scored_and_minimisation_never_raises_the_energy(tmp_path):
    r = rescore(_receptor(tmp_path), _ligand(tmp_path), out_ligand_pdb=tmp_path / "min.pdb")
    assert r["score_min"] <= r["score_start"] + 1e-6
    assert r["n_interface_residues"] > 10 and r["dSASA"] > 100
    assert np.isfinite([r["dG_separated"], r["sc_value"]]).all()
    assert (tmp_path / "min.pdb").read_text().count("HETATM") >= 9      # the minimised ligand is written for the RMSD stage


def test_rescoring_is_deterministic(tmp_path):
    a = rescore(_receptor(tmp_path), _ligand(tmp_path))
    b = rescore(_receptor(tmp_path), _ligand(tmp_path, name="again.pdb"))
    assert a["dG_separated"] == pytest.approx(b["dG_separated"], abs=1e-6)


def test_a_ligand_moved_out_of_the_site_has_no_interface(tmp_path):
    far = rescore(_receptor(tmp_path), _ligand(tmp_path, shift=(60.0, 60.0, 60.0), name="far.pdb"))
    assert far["n_interface_residues"] == 0 and abs(far["dG_separated"]) < 1.0


def test_the_minimiser_moves_the_ligand_only_a_little_from_a_good_pose(tmp_path):
    r = rescore(_receptor(tmp_path), _ligand(tmp_path))
    assert 0.0 <= r["ligand_shift"] < 2.0


def _rechain(src, new, dst):
    """Same receptor with its chain letter changed (6X0J's chains are B, C and D, not A)."""
    dst.write_text("".join(l[:21] + new + l[22:] if l.startswith(("ATOM", "HETATM")) else l for l in src.read_text().splitlines(True)))
    return dst


def test_the_receptor_chain_letter_does_not_matter(tmp_path):
    """Regression: the interface was hard-coded as ligand-vs-chain-A, so chains B/C/D silently scored dG = 0."""
    a = rescore(_receptor(tmp_path), _ligand(tmp_path))
    rec_b = _rechain(_receptor(tmp_path), "B", tmp_path / "rec_b.pdb")
    b = rescore(rec_b, _ligand(tmp_path, name="lig_b.pdb"))
    assert b["n_interface_residues"] == a["n_interface_residues"] > 10
    assert b["dG_separated"] == pytest.approx(a["dG_separated"], abs=1e-3)


def test_a_batch_in_which_no_pose_touches_the_receptor_fails_loudly(tmp_path):
    far = _ligand(tmp_path, shift=(60.0, 60.0, 60.0), name="far.pdb")
    with pytest.raises(RuntimeError, match="no pose"):
        rescore_many(_receptor(tmp_path), [far], tmp_path / "o", jobs=1)
