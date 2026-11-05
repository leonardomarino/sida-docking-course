"""Side-chain comparison of one residue across structures (biotite superposition and dihedrals; no hand-written geometry)."""
from pathlib import Path

import numpy as np
import pytest
import biotite.structure as struc

from sidacore import sidechain as S


def _lys(origin=(0.0, 0.0, 0.0), chain="A", res_id=107):
    """An all-trans lysine in the xy plane: every chi is 180 degrees."""
    names = ["N", "CA", "CB", "CG", "CD", "CE", "NZ"]
    xy = [(0.0, 0.0), (1.2, 0.7), (2.4, 0.0), (3.6, 0.7), (4.8, 0.0), (6.0, 0.7), (7.2, 0.0)]
    a = struc.AtomArray(len(names))
    for i, (n, (x, y)) in enumerate(zip(names, xy)):
        a.atom_name[i], a.res_name[i], a.res_id[i], a.chain_id[i] = n, "LYS", res_id, chain
        a.element[i] = "N" if n in ("N", "NZ") else "C"
        a.coord[i] = (x + origin[0], y + origin[1], origin[2])
    return a


def test_chi_angles_of_an_all_trans_chain_are_180():
    chis = S.chi_angles(_lys())
    assert len(chis) == 4
    assert all(abs(abs(c) - 180.0) < 1e-3 for c in chis)


def test_chi_angles_change_when_a_torsion_changes():
    a = _lys()
    a.coord[3:] = struc.rotate_about_axis(a.coord[3:], a.coord[2] - a.coord[1], 60.0, support=a.coord[2])   # twist about CA-CB
    chis = S.chi_angles(a)
    assert abs(abs(chis[0]) - 180.0) > 30.0                                    # chi1 moved
    assert abs(abs(chis[1]) - 180.0) < 1e-3                                     # the others did not


def test_a_rigidly_moved_copy_has_zero_displacement_after_fitting():
    ref = _lys()
    moved = struc.rotate(struc.translate(ref, (5.0, -3.0, 2.0)), (0.3, 0.2, 0.1))
    got = S.displacement_after_fit(ref, moved, atom_name="NZ")
    assert got["displacement"] < 1e-3
    assert got["n_anchors"] == 7


def test_a_displaced_side_chain_is_measured_against_a_common_backbone():
    ref = _lys()
    bent = ref.copy()
    bent.coord[5:] += (0.0, 0.0, 3.0)                                           # CE and NZ lifted 3 A, backbone untouched
    got = S.displacement_after_fit(ref, bent, atom_name="NZ", anchors=("N", "CA", "CB"))
    assert got["displacement"] == pytest.approx(3.0, abs=1e-2)


def test_missing_atom_raises():
    ref = _lys()
    with pytest.raises(ValueError):
        S.displacement_after_fit(ref, ref[ref.atom_name != "NZ"], atom_name="NZ")


CACHE = Path(__file__).resolve().parents[1] / "structures_cache" / "4B63.cif"


@pytest.mark.skipif(not CACHE.exists(), reason="structure cache not present")
def test_real_structures_report_lys107_for_every_chain():
    rows = S.compare_structures(CACHE.parent, res_id=107, res_name="LYS", reference=("4B63", "A"),
                                members=[("4B63", "A"), ("6X0J", "A"), ("6X0J", "C")], ligand_name="ORN")
    assert [r["label"] for r in rows] == ["4B63:A", "6X0J:A", "6X0J:C"]
    assert rows[0]["displacement"] == pytest.approx(0.0, abs=1e-6)             # the reference against itself
    assert all(len(r["chi"]) == 4 and 0.0 < r["nz_to_carboxylate"] < 4.0 for r in rows)
    md = S.render_markdown(rows, res_label="Lys107")
    assert "6X0J:C" in md and "chi1" in md
