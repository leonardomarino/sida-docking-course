"""All-residue comparison of the substrate site, 4B63 against the four 6X0J chains (biotite fit; existing site definition)."""
from pathlib import Path

import numpy as np
import pytest
import biotite.structure as struc

from sidacore import sitecompare as SC


def _res(name, rows, res_id=1, chain="A"):
    a = struc.AtomArray(len(rows))
    for i, (an, el, x, y, z) in enumerate(rows):
        a.atom_name[i], a.element[i], a.res_name[i], a.res_id[i], a.chain_id[i] = an, el, name, res_id, chain
        a.coord[i] = (x, y, z)
    return a


LYS = [("N", "N", 0, 0, 0), ("CA", "C", 1, 0, 0), ("C", "C", 2, 0, 0), ("O", "O", 3, 0, 0),
       ("CB", "C", 1, 1, 0), ("CG", "C", 1, 2, 0), ("CD", "C", 1, 3, 0), ("CE", "C", 1, 4, 0), ("NZ", "N", 1, 5, 0)]


def test_identical_side_chains_have_zero_rmsd_and_backbone_is_ignored():
    a = _res("LYS", LYS)
    moved_backbone = a.copy()
    moved_backbone.coord[:4] += 5.0                                             # N, CA, C, O moved: not a side-chain atom
    assert SC.sidechain_rmsd(a, moved_backbone) == pytest.approx(0.0)


def test_one_displaced_atom_gives_the_expected_rmsd():
    a = _res("LYS", LYS)
    b = a.copy()
    b.coord[b.atom_name == "NZ"] += (3.0, 0.0, 0.0)
    assert SC.sidechain_rmsd(a, b) == pytest.approx(3.0 / np.sqrt(5), rel=1e-6)     # 5 side-chain atoms: CB CG CD CE NZ


def test_a_symmetric_carboxylate_swap_is_not_a_difference():
    asp = [("N", "N", 0, 0, 0), ("CA", "C", 1, 0, 0), ("C", "C", 2, 0, 0), ("O", "O", 3, 0, 0), ("CB", "C", 1, 1, 0),
           ("CG", "C", 1, 2, 0), ("OD1", "O", 0.2, 3, 0), ("OD2", "O", 2.2, 3, 0)]
    a = _res("ASP", asp)
    b = a.copy()
    b.atom_name[b.atom_name == "OD1"], b.atom_name[b.atom_name == "OD2"] = "OD2", "OD1"      # same atoms, labels exchanged
    assert SC.sidechain_rmsd(a, b) == pytest.approx(0.0, abs=1e-9)


def test_glycine_has_no_side_chain_and_a_name_mismatch_is_an_error():
    gly = _res("GLY", LYS[:4])
    assert SC.sidechain_rmsd(gly, gly) is None
    with pytest.raises(ValueError):
        SC.sidechain_rmsd(_res("LYS", LYS), _res("ARG", LYS))


def test_outlier_rule_needs_to_beat_the_whole_reference_cloud_and_the_floor():
    assert SC.is_outlier([1.5, 1.8], [0.3, 0.4, 0.9], floor=1.0) is True
    assert SC.is_outlier([1.5, 1.8], [0.3, 1.6], floor=1.0) is False           # a 6X0J pair differs as much: not specific to 4B63
    assert SC.is_outlier([0.8, 0.9], [0.1], floor=1.0) is False                # below the floor
    assert SC.is_outlier([1.5, 0.2], [0.1], floor=1.0) is False                # one 6X0J chain matches 4B63: not outside the cloud


DATA = Path(__file__).resolve().parents[1] / "structures_cache" / "4B63.cif"


@pytest.mark.skipif(not DATA.exists(), reason="structure cache not present")
def test_real_comparison_covers_the_site_and_both_cofactors():
    res = SC.compare(DATA.parent)
    ids = {r["res_id"] for r in res["residues"]}
    assert {107, 293, 323, 469} <= ids                                           # the four ornithine contacts
    for r in res["residues"]:
        assert len(r["vs_ref"]) == 4                                             # one value per 6X0J chain
        if r["resname"] != "GLY":
            assert all(v is not None and v >= 0.0 for v in r["vs_ref"])
    assert {c["name"] for c in res["cofactors"]} >= {"flavin", "NADP"}
    assert "LYS107" in SC.render_markdown(res)
