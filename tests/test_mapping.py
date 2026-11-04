"""Construct check: explicit residue mapping between two structures. Never renumbers anything."""
import csv

import pytest

from sidacore.mapping import map_residues, write_mapping_csv


def res(seq, start=1, skip=()):
    """[(num, icode, aa)] for a sequence, author numbers from `start`, omitting numbers in `skip`."""
    return [(start + i, "", a) for i, a in enumerate(seq) if (start + i) not in skip]


SEQ = "MKTAYIAKQRQISFVKSHFSRQ"


def test_identical_numbering_with_different_coverage_is_shared():
    m = map_residues(res(SEQ, 10), res(SEQ, 10, skip={14, 15}))
    assert m.numbering_shared and m.n_mismatch == 0 and m.n_different_number == 0
    assert m.n_ref == len(SEQ) and m.n_qry == len(SEQ) - 2
    assert m.n_ref_only == 2 and m.n_qry_only == 0          # residues resolved in one structure only


def test_point_mutation_is_reported_not_hidden():
    q = res(SEQ, 10); q[5] = (15, "", "W")
    m = map_residues(res(SEQ, 10), q)
    assert m.n_mismatch == 1 and m.numbering_shared
    bad = [r for r in m.rows if not r["identical"] and r["ref_aa"] and r["qry_aa"]]
    assert [(r["ref_num"], r["ref_aa"], r["qry_aa"]) for r in bad] == [("15", "I", "W")]


def test_shifted_numbering_is_detected_with_its_offset():
    m = map_residues(res(SEQ, 10), res(SEQ, 12))
    assert not m.numbering_shared
    assert m.offsets == {2: len(SEQ)}                       # every aligned pair is off by exactly +2


def test_insertion_in_query_is_a_gap_row():
    q = res(SEQ[:8] + "GGG" + SEQ[8:], 10)
    m = map_residues(res(SEQ, 10), q)
    assert m.n_qry_only == 3 and m.n_mismatch == 0


def test_every_input_residue_appears_exactly_once_in_rows():
    r, q = res(SEQ, 10), res(SEQ, 10, skip={20})
    m = map_residues(r, q)
    assert sum(1 for x in m.rows if x["ref_aa"]) == len(r) and sum(1 for x in m.rows if x["qry_aa"]) == len(q)


def test_csv_round_trip(tmp_path):
    m = map_residues(res(SEQ, 10), res(SEQ, 10, skip={14}))
    out = write_mapping_csv(m, tmp_path / "map.csv")
    rows = list(csv.DictReader(out.open()))
    assert len(rows) == len(m.rows)
    assert set(rows[0]) == {"ref_num", "ref_aa", "qry_num", "qry_aa", "identical", "same_number"}


def test_empty_input_is_an_error_not_a_silent_empty_mapping():
    with pytest.raises(ValueError):
        map_residues([], res(SEQ))


# ---- the real construct check (Amendment 1): 4B63 must share numbering with the 6X0* batch -------------------
from pathlib import Path

from sidacore.mapping import load_polymer_residues

CACHE = Path(__file__).resolve().parent.parent / "structures_cache"
needs_cache = pytest.mark.skipif(not all((CACHE / ("%s.cif" % i)).exists() for i in ("4B63", "6X0H", "6X0I", "6X0J")),
                                 reason="structures_cache incomplete")


@needs_cache
@pytest.mark.parametrize("qry,identical", [("6X0H", 448), ("6X0I", 443), ("6X0J", 447)])
def test_4b63_shares_numbering_and_residues_with_the_6x0_batch(qry, identical):
    m = map_residues(load_polymer_residues(CACHE / "4B63.cif", "A"), load_polymer_residues(CACHE / ("%s.cif" % qry), "A"))
    assert m.numbering_shared and m.n_mismatch == 0 and m.offsets == {}
    assert m.n_aligned == identical                        # every aligned residue is identical: a point mutation would show here


@needs_cache
def test_loading_a_missing_chain_is_an_error():
    with pytest.raises(Exception):
        load_polymer_residues(CACHE / "4B63.cif", "Z")
