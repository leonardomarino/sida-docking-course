"""Polar-contact retention of docked poses against the deposited ligand (biotite cell list; no hand-written distance loops)."""
from pathlib import Path

import numpy as np
import pytest
import biotite.structure as struc

from sidacore import contacts as C


def _atoms(rows):
    """rows: (res_name, chain, res_id, atom_name, element, x, y, z)"""
    a = struc.AtomArray(len(rows))
    for i, (rn, ch, ri, an, el, x, y, z) in enumerate(rows):
        a.res_name[i], a.chain_id[i], a.res_id[i], a.atom_name[i], a.element[i] = rn, ch, ri, an, el
        a.coord[i] = (x, y, z)
    return a


REC = _atoms([
    ("LYS", "A", 107, "NZ", "N", 0.0, 0.0, 0.0),
    ("LYS", "A", 107, "CE", "C", 1.0, 0.0, 0.0),     # carbon: never a polar contact
    ("SER", "A", 469, "OG", "O", 10.0, 0.0, 0.0),
    ("ASN", "A", 293, "ND2", "N", 0.0, 10.0, 0.0),
    ("HOH", "A", 900, "H1", "H", 0.5, 0.0, 0.0),     # hydrogen: ignored
])
LIG = _atoms([
    ("ORN", "L", 1, "O", "O", 2.8, 0.0, 0.0),       # 2.8 from Lys NZ
    ("ORN", "L", 1, "OXT", "O", 10.0, 3.0, 0.0),    # 3.0 from Ser OG
    ("ORN", "L", 1, "N", "N", 0.0, 14.0, 0.0),      # 4.0 from Asn ND2: beyond a 3.5 cutoff
    ("ORN", "L", 1, "CA", "C", 0.5, 0.0, 0.0),
])


def test_polar_contacts_keep_only_heavy_n_o_pairs_within_the_cutoff():
    got = C.polar_contacts(REC, LIG, cutoff=3.5)
    assert {(c.res_name, c.res_id, c.rec_atom, c.lig_atom) for c in got} == {
        ("LYS", 107, "NZ", "O"), ("SER", 469, "OG", "OXT")}
    assert all(c.distance <= 3.5 for c in got)


def test_the_cutoff_is_respected_in_both_directions():
    assert len(C.polar_contacts(REC, LIG, cutoff=2.0)) == 0
    assert len(C.polar_contacts(REC, LIG, cutoff=4.5)) == 3


def test_contact_residues_are_identified_by_name_chain_and_number():
    res = C.contact_residues(C.polar_contacts(REC, LIG, cutoff=3.5))
    assert res == {("LYS", "A", 107), ("SER", "A", 469)}


def test_retention_separates_kept_from_lost_residues():
    native = {("LYS", "A", 107), ("SER", "A", 469), ("ASN", "A", 293)}
    pose = {("LYS", "A", 107), ("THR", "A", 5)}
    kept, lost = C.retention(native, pose)
    assert kept == {("LYS", "A", 107)}
    assert lost == {("SER", "A", 469), ("ASN", "A", 293)}      # a new contact (THR 5) does not rescue a lost one


def test_protein_only_drops_cofactor_residues():
    res = {("LYS", "A", 107), ("NAP", "A", 601), ("FAD", "A", 600)}
    assert C.protein_only(res) == {("LYS", "A", 107)}


def test_min_distance_uses_only_the_named_atoms():
    d = C.min_distance(REC, LIG, res_id=107, rec_atom="NZ", lig_atoms=("O", "OXT"))
    assert d == pytest.approx(2.8)
    assert C.min_distance(REC, LIG, res_id=107, rec_atom="NZ", lig_atoms=("OXT",)) > 9.0


def test_min_distance_raises_for_a_missing_atom_instead_of_returning_nan():
    with pytest.raises(ValueError):
        C.min_distance(REC, LIG, res_id=999, rec_atom="NZ", lig_atoms=("O",))


DATA = Path(__file__).resolve().parents[1] / "rosetta_in" / "4B63_A_nadp-keep"


@pytest.mark.skipif(not DATA.exists(), reason="Rosetta inputs not generated")
def test_deposited_4b63_has_the_published_lys107_salt_bridge_and_retains_itself():
    rec, lig = C.read_job(DATA, "crystal")
    native = C.contact_residues(C.polar_contacts(rec, lig))
    assert {r[2] for r in native} >= {107, 293, 323, 469}                    # Robinson 2015: Lys107, Asn293, Asn323, Ser469
    assert C.min_distance(rec, lig, 107, "NZ", ("O", "OXT")) < 3.2
    kept, lost = C.retention(native, native)
    assert not lost


@pytest.mark.skipif(not DATA.exists(), reason="Rosetta inputs not generated")
def test_analysis_and_report_on_the_real_inputs():
    res = C.analyse(DATA.parent)
    assert len(res) == 5
    ok = res["6X0J_A_nadp-keep"]["engines"]["vina"]["top1"]                  # a registered PASS keeps every contact
    assert ok["lost"] == [] and ok["lys107"] < 3.5
    bad = res["4B63_A_nadp-keep"]["engines"]["vina"]["top1"]                 # a registered FAIL loses the Lys107 salt bridge
    assert 107 in bad["lost"] and bad["lys107"] > 5.0
    md = C.render_markdown(res)
    assert "| 4B63_A_nadp-keep | vina |" in md and "every protein contact" in md
