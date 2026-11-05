"""Chimeric receptors: side chains of chosen residues taken from a donor structure, everything else left as deposited."""
from pathlib import Path

import numpy as np
import pytest
import biotite.structure as struc

from sidacore import chimera as CH


def _protein(n=5, lys_at=3, shift=(0.0, 0.0, 0.0), nz_extra=(0.0, 0.0, 0.0)):
    """n residues on a gentle helix-like path; residue `lys_at` is a lysine with a side chain, the others are alanines."""
    atoms = []
    for i in range(1, n + 1):
        base = np.array([3.8 * i, 2.0 * np.sin(i), 2.0 * np.cos(i)]) + np.array(shift)
        name = "LYS" if i == lys_at else "ALA"
        rows = [("N", "N", base + (-1.2, 0.3, 0.0)), ("CA", "C", base), ("C", "C", base + (1.2, 0.4, 0.0)), ("O", "O", base + (1.4, 1.5, 0.0)),
                ("CB", "C", base + (0.0, -1.0, 1.2))]
        if name == "LYS":
            for k, (an, el) in enumerate([("CG", "C"), ("CD", "C"), ("CE", "C"), ("NZ", "N")]):
                extra = np.array(nz_extra) if an == "NZ" else 0.0
                rows.append((an, el, base + (0.0, -1.5 - 1.2 * k, 2.2 + 0.4 * k) + extra))
        for an, el, xyz in rows:
            atoms.append((i, name, an, el, xyz))
    a = struc.AtomArray(len(atoms))
    for j, (rid, rn, an, el, xyz) in enumerate(atoms):
        a.res_id[j], a.res_name[j], a.atom_name[j], a.element[j], a.chain_id[j] = rid, rn, an, el, "A"
        a.coord[j] = xyz
    return a


def test_swap_takes_the_donor_side_chain_and_keeps_everything_else():
    target = _protein()
    donor = struc.rotate(struc.translate(_protein(nz_extra=(0.0, 0.0, 2.0)), (5.0, -3.0, 1.0)), (0.2, 0.1, 0.3))   # moved rigidly, NZ also displaced 2 A
    out, log = CH.swap_sidechains(target, donor, res_ids=[3])
    nz = lambda a: a[(a.res_id == 3) & (a.atom_name == "NZ")].coord[0]
    assert np.linalg.norm(nz(out) - (nz(target) + (0.0, 0.0, 2.0))) < 1e-3          # donor NZ, expressed in the target frame
    keep = ~((target.res_id == 3) & ~np.isin(target.atom_name, ["N", "CA", "C", "O"]))
    assert np.allclose(out.coord[keep], target.coord[keep])                           # backbone and other residues untouched
    assert log[0]["res_id"] == 3 and log[0]["n_atoms"] == 5
    assert log[0]["backbone_rmsd"] < 1e-3                                           # the donor residue sat on the target backbone exactly


def test_the_junction_with_the_backbone_keeps_the_target_geometry_even_if_the_donor_backbone_is_offset():
    """The bug this guards: copying donor coordinates after a global fit stretched CA-CB to 2.05 A when the two backbones differ by ~0.5 A."""
    target = _protein()
    donor = target.copy()
    res = donor.res_id == 3
    donor.coord[res] += (0.5, 0.4, 0.0)                                             # this residue sits 0.64 A away in the donor's frame
    out, _ = CH.swap_sidechains(target, donor, res_ids=[3])
    ca = lambda a: a[(a.res_id == 3) & (a.atom_name == "CA")].coord[0]
    cb = lambda a: a[(a.res_id == 3) & (a.atom_name == "CB")].coord[0]
    assert np.linalg.norm(cb(out) - ca(out)) == pytest.approx(np.linalg.norm(cb(target) - ca(target)), abs=1e-6)
    assert np.allclose(out.coord, target.coord, atol=1e-6)                           # same rotamer on the same backbone: nothing moves


def test_a_residue_type_mismatch_is_an_error():
    target, donor = _protein(), _protein()
    donor.res_name[donor.res_id == 3] = "ARG"
    with pytest.raises(ValueError):
        CH.swap_sidechains(target, donor, res_ids=[3])


def test_a_missing_donor_atom_is_an_error():
    target, donor = _protein(), _protein()
    donor = donor[~((donor.res_id == 3) & (donor.atom_name == "NZ"))]
    with pytest.raises(ValueError):
        CH.swap_sidechains(target, donor, res_ids=[3])


def test_clash_check_finds_a_planted_overlap_and_ignores_covalent_neighbours():
    a = _protein()
    assert CH.closest_contact(a, res_ids=[3], cutoff=2.2) is None                   # nothing closer than a bond-length excluded
    b = a.copy()
    b.coord[(b.res_id == 3) & (b.atom_name == "NZ")] = b.coord[(b.res_id == 5) & (b.atom_name == "CA")] + (0.0, 0.0, 1.0)
    hit = CH.closest_contact(b, res_ids=[3], cutoff=2.2)
    assert hit is not None and hit["distance"] == pytest.approx(1.0, abs=1e-3) and hit["other_res_id"] == 5


CACHE = Path(__file__).resolve().parents[1] / "structures_cache"


@pytest.mark.skipif(not (CACHE / "4B63.cif").exists(), reason="structure cache not present")
def test_round_trip_through_a_cif_changes_nothing_and_a_real_swap_changes_only_the_named_residue(tmp_path):
    from sidacore.prep import read_structure
    same = CH.build_chimera_cif(CACHE, tmp_path / "rt.cif", target=("4B63", "A"), donor=("6X0J", "A"), res_ids=[])
    base, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    rt, _ = read_structure(tmp_path / "rt.cif", "highest-occupancy")
    assert len(rt) == len(base) and np.abs(rt.coord - base.coord).max() < 1e-3
    assert same["swapped"] == []
    CH.build_chimera_cif(CACHE, tmp_path / "q.cif", target=("4B63", "A"), donor=("6X0J", "A"), res_ids=[102])
    q, _ = read_structure(tmp_path / "q.cif", "highest-occupancy")
    moved = np.linalg.norm(q.coord - base.coord, axis=1) > 1e-3
    assert set(q.res_id[moved]) == {102} and q.chain_id[moved].tolist() == ["A"] * int(moved.sum())
    assert not np.isin(q.atom_name[moved], ["N", "CA", "C", "O"]).any()               # side chain only
    ca_cb = lambda a: np.linalg.norm(a[(a.res_id == 102) & (a.atom_name == "CB")].coord[0] - a[(a.res_id == 102) & (a.atom_name == "CA")].coord[0])
    assert abs(ca_cb(q) - ca_cb(base)) < 0.05                                          # the swap must not stretch the CA-CB bond (it did: 2.05 A)


@pytest.mark.skipif(not (CACHE / "4B63.cif").exists() or not (CACHE.parent / "prepared").exists(), reason="inputs not present")
def test_a_variant_that_cannot_be_prepared_is_recorded_not_dropped_and_the_run_continues(tmp_path, monkeypatch):
    from sidacore.config import load_config
    cfg = load_config(CACHE.parent / "config" / "default.yaml")

    def boom(*a, **k):
        raise RuntimeError("mk_prepare_receptor.py failed (rc=1)\nmatched with excess inter-residue bond(s): A:102")

    monkeypatch.setattr(CH, "prepare_receptor", boom)
    rep = CH.run_chimeras(cfg, CACHE.parent / "prepared", tmp_path, cache_dir=CACHE, variants=[("two bad ones", "A", [102]), ("another", "A", [288])])
    assert [v["label"] for v in rep["variants"]] == ["two bad ones", "another"]          # both recorded, none dropped, the loop did not stop
    assert all("excess inter-residue" in v["prepare_error"] and "top1_median" not in v for v in rep["variants"])
    md = CH.render_markdown(rep)
    assert "not prepared" in md and "two bad ones" in md


@pytest.mark.skipif(not (CACHE / "4B63.cif").exists() or not (CACHE.parent / "prepared").exists(), reason="inputs not present")
def test_the_job_for_a_target_is_chosen_by_structure_and_chain_not_by_prefix():
    from sidacore.config import load_config
    cfg = load_config(CACHE.parent / "config" / "default.yaml")
    assert CH.job_for(cfg, CACHE.parent / "prepared", ("6X0J", "C")).name.startswith("6X0J_C_")
    assert CH.job_for(cfg, CACHE.parent / "prepared", ("6X0J", "D")).name.startswith("6X0J_D_")
    with pytest.raises(ValueError):
        CH.job_for(cfg, CACHE.parent / "prepared", ("6X0J", "Z"))


@pytest.mark.skipif(not (CACHE / "4B63.cif").exists(), reason="structure cache not present")
def test_the_reverse_direction_gives_a_6x0j_chain_the_4b63_residue_and_nothing_else(tmp_path):
    from sidacore.prep import read_structure
    CH.build_chimera_cif(CACHE, tmp_path / "r.cif", target=("6X0J", "B"), donor=("4B63", "A"), res_ids=[102])
    base, _ = read_structure(CACHE / "6X0J.cif", "highest-occupancy")
    q, _ = read_structure(tmp_path / "r.cif", "highest-occupancy")
    moved = np.linalg.norm(q.coord - base.coord, axis=1) > 1e-3
    assert set(zip(q.chain_id[moved], q.res_id[moved])) == {("B", 102)}                  # one residue of one chain
    assert not np.isin(q.atom_name[moved], ["N", "CA", "C", "O"]).any()
    ca_cb = lambda a: np.linalg.norm(a[(a.chain_id == "B") & (a.res_id == 102) & (a.atom_name == "CB")].coord[0] - a[(a.chain_id == "B") & (a.res_id == 102) & (a.atom_name == "CA")].coord[0])
    assert abs(ca_cb(q) - ca_cb(base)) < 0.05


def test_a_multi_target_report_names_the_target_and_the_direction():
    rep = {"settings": {"engine": "vina", "exhaustiveness": 32, "num_modes": 20, "seeds": list(range(1, 11)), "acceptance": 2.0, "protonation_ph": 7.4,
                        "target": "6X0J:A, 6X0J:B", "direction": "reverse", "donor": "4B63:A"},
           "variants": [{"target": "6X0J:A", "label": "4B63 Gln102", "donor_chain": "A", "res_ids": [102], "swapped": [{"res_id": 102, "res_name": "GLN"}], "closest_contact": None,
                         "site_states_changed": {}, "top1_median": 1.0, "top1_n_within": 10, "n_seeds": 10, "best_top5_median": 1.0, "passes": True}]}
    md = CH.render_markdown(rep)
    assert "reverse" in md.lower()
    assert "| Target |" in md                                                           # a Target column heads the table
    assert "| 6X0J:A | 4B63 Gln102 |" in md                                              # and the row starts with the target
