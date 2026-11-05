"""Step 2 glue: biotite selection -> PDB2PQR/PROPKA -> Meeko. Integration tests on the real structures."""
import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from sidacore.config import ConfigError, load_config
from sidacore.prep import site_residue_ids
from sidacore.prep import (SiteAtomError, fetch_structure, prepare_receptor, read_structure,
                           resolve_site_center, select_atoms)

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "structures_cache"
DEFAULT = ROOT / "config" / "default.yaml"
PDB2PQR = Path("/opt/miniforge3/envs/af2/bin/pdb2pqr30")

need = pytest.mark.skipif(not all((CACHE / ("%s.cif" % i)).exists() for i in ("4B63", "6X0I")),
                          reason="structures_cache incomplete")
need_tools = pytest.mark.skipif(not PDB2PQR.exists() or shutil.which("obabel") is None, reason="pdb2pqr30 not available")


@pytest.fixture(scope="module")
def cfg():
    return load_config(DEFAULT)


@need
def test_fetch_uses_existing_cache_without_network():
    p = fetch_structure("4B63", CACHE)
    assert p == CACHE / "4B63.cif" and len(p.read_bytes()) > 1000


@need
def test_site_center_is_c4x_and_matches_measured_coordinates(cfg):
    atoms, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    c = resolve_site_center(atoms, "A", "FAD:C4X")
    assert c == pytest.approx([-25.372, -31.449, 14.541], abs=1e-3)


@need
def test_adenine_c4a_is_refused_not_silently_used(cfg):
    # C4A exists in FAD but is the ADENINE carbon, 18 A from the hydroxylation site
    atoms, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    with pytest.raises(SiteAtomError, match="N5"):
        resolve_site_center(atoms, "A", "FAD:C4A")


@need
def test_missing_site_atom_or_cofactor_is_an_error(cfg):
    atoms, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    with pytest.raises(SiteAtomError):
        resolve_site_center(atoms, "A", "FAD:NOPE")
    with pytest.raises(SiteAtomError):
        resolve_site_center(atoms, "Z", "FAD:C4X")


@need
def test_altloc_removals_are_logged_per_atom():
    atoms, log = read_structure(CACHE / "6X0I.cif", "highest-occupancy")
    assert any(r["res_name"] == "MET" and r["res_id"] == 101 and r["chain_id"] == "A" for r in log)
    assert all(r["reason"].startswith("alternate conformer") for r in log) and len(log) == 131


@need
def test_selection_logs_every_removed_atom_with_a_reason(cfg):
    atoms, _ = read_structure(CACHE / "6X0I.cif", "highest-occupancy")
    c = resolve_site_center(atoms, "A", "FAD:C4X")
    sel = select_atoms(atoms, chain="A", nadp="strip", cfg=cfg, center=c)
    assert sel.kept.array_length() + len(sel.log) + (sel.ligand.array_length() if sel.ligand is not None else 0) == atoms.array_length()
    reasons = {r["reason"] for r in sel.log}
    assert {"other chain", "NADP stripped (policy)", "water stripped (policy)", "ion stripped (policy)",
            "crystallization additive / non-cofactor heterogen"} <= reasons
    assert not any(r["reason"] == "" for r in sel.log)
    assert (sel.kept.res_name == "FAD").sum() == 53                    # the cofactor is KEPT
    assert not np.isin(sel.kept.res_name, ["NAP", "ACT", "CA", "HOH"]).any()
    nadp = [r for r in sel.log if r["reason"] == "NADP stripped (policy)"]
    assert len(nadp) == 48 and all("dist_to_site" in r for r in nadp)


@need
def test_nadp_keep_retains_it_and_substrate_is_set_aside_not_in_receptor(cfg):
    atoms, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    c = resolve_site_center(atoms, "A", "FAD:C4X")
    sel = select_atoms(atoms, chain="A", nadp="keep", cfg=cfg, center=c)
    assert (sel.kept.res_name == "NAP").sum() == 48 and (sel.kept.res_name == "FAD").sum() == 53
    assert not (sel.kept.res_name == "ORN").any()                      # the answer must not be in the receptor
    assert sel.ligand is not None and sel.ligand.array_length() == 9 and set(sel.ligand.res_name) == {"ORN"}


@need
def test_waters_keep_within_pocket_only_keeps_those_near_the_site():
    cfg = load_config(DEFAULT, overrides={"waters": "keep-within-pocket"})
    atoms, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    c = resolve_site_center(atoms, "A", "FAD:C4X")
    sel = select_atoms(atoms, chain="A", nadp="keep", cfg=cfg, center=c)
    w = sel.kept[sel.kept.res_name == "HOH"]
    assert w.array_length() > 0
    assert (np.linalg.norm(w.coord - c, axis=1) <= cfg.pocket_radius_angstrom).all()
    assert any(r["reason"] == "water stripped (policy)" for r in sel.log)


@need
@need_tools
def test_prepare_4b63_end_to_end(cfg, tmp_path):
    rep = prepare_receptor(CACHE / "4B63.cif", "4B63", "A", cfg, "keep", tmp_path / "4B63_A_nadp-keep")
    d = tmp_path / "4B63_A_nadp-keep"
    for f in ("receptor.pdbqt", "receptor.box.txt", "receptor_prepared.pdb", "ligand_crystal.pdb", "removal_log.csv",
              "protonation_states.csv", "prep_report.json", "resolved_config.yaml"):
        assert (d / f).exists(), f
    pdbqt = (d / "receptor.pdbqt").read_text()
    assert " FAD " in pdbqt and " NAP " in pdbqt and " ORN " not in pdbqt      # cofactors in, substrate OUT
    assert rep["cofactors"]["FAD"]["charge"] == pytest.approx(-2.0, abs=0.01)
    assert rep["cofactors"]["NAP"]["charge"] == pytest.approx(-3.0, abs=0.01)
    assert rep["fad_retained_in_receptor"] is True and "cofactor" in rep["fad_statement"]
    assert rep["site_center"] == pytest.approx([-25.372, -31.449, 14.541], abs=1e-3)
    assert rep["meeko_deleted_residues"] == [] or all(x["dist_to_site"] > 8 + 11 for x in rep["meeko_deleted_residues"])
    assert rep["pdb2pqr_rebuilt_atoms"]["n_atoms"] > 0                          # 4B63 has truncated surface side chains
    rows = list(csv.DictReader((d / "protonation_states.csv").open()))
    his91 = [r for r in rows if r["res_name"] == "HIS" and r["res_id"] == "91"]
    assert len(his91) == 1 and float(his91[0]["pKa"]) == pytest.approx(4.5, abs=0.2) and "6.5" in his91[0]["note"]
    assert rep["tools"]["meeko"] and rep["tools"]["pdb2pqr"] and rep["tools"]["biotite"]
    assert all(isinstance(c, str) for cmd in rep["commands"] for c in cmd)
    assert len(rep["commands"]) == 2 and "--add_templates" in rep["commands"][1]


@need
@need_tools
def test_prepare_refuses_to_overwrite_an_existing_run(cfg, tmp_path):
    out = tmp_path / "run"
    prepare_receptor(CACHE / "4B63.cif", "4B63", "A", cfg, "keep", out)
    with pytest.raises(ConfigError, match="exists"):
        prepare_receptor(CACHE / "4B63.cif", "4B63", "A", cfg, "keep", out)


@need
@need_tools
def test_prepare_6x0i_nadp_strip_has_flavin_but_no_nadp(cfg, tmp_path):
    rep = prepare_receptor(CACHE / "6X0I.cif", "6X0I", "A", cfg, "strip", tmp_path / "x")
    pdbqt = (tmp_path / "x" / "receptor.pdbqt").read_text()
    assert " FAD " in pdbqt and " NAP " not in pdbqt
    assert "NAP" not in rep["cofactors"] and rep["nadp_policy"] == "strip"


@need
@need_tools
def test_ph_change_produces_a_visible_diff_in_protonation_table(tmp_path):
    lo = load_config(DEFAULT, overrides={"protonation_ph": 4.0})
    hi = load_config(DEFAULT)
    prepare_receptor(CACHE / "4B63.cif", "4B63", "A", lo, "keep", tmp_path / "lo")
    prepare_receptor(CACHE / "4B63.cif", "4B63", "A", hi, "keep", tmp_path / "hi")
    st = lambda d: {(r["res_name"], r["res_id"]): r["state"] for r in csv.DictReader((tmp_path / d / "protonation_states.csv").open())}
    a, b = st("lo"), st("hi")
    assert a != b                                                            # pH 4.0 vs 7.4 must change some site residue


# ------------------------------------------------------------------ the substrate site (review finding): not a sphere about the cofactor
def _res(chain, rid, name, xyz):
    import biotite.structure as struc
    a = struc.AtomArray(1)
    a.coord = np.array([xyz], float)
    a.chain_id[:], a.res_id[:], a.res_name[:], a.atom_name[:], a.element[:] = chain, rid, name, "CA", "C"
    return a


def test_site_is_every_residue_in_the_docking_box_or_touching_the_ligand():
    atoms = _res("A", 1, "LYS", (5, 0, 0)) + _res("A", 2, "ASP", (30, 0, 0)) + _res("A", 3, "GLU", (60, 0, 0))
    lig = np.array([[33.0, 0.0, 0.0]])                          # 3 A from residue 2, which is outside the box
    site = site_residue_ids(atoms, center=(0, 0, 0), box_size=(22, 22, 22), ligand_xyz=lig, contact=4.5)
    assert site == {("A", 1): {"in_box": True, "contacts_ligand": False}, ("A", 2): {"in_box": False, "contacts_ligand": True}}


def test_without_a_ligand_the_site_is_the_box_alone():
    atoms = _res("A", 1, "LYS", (5, 0, 0)) + _res("A", 2, "ASP", (30, 0, 0))
    assert set(site_residue_ids(atoms, (0, 0, 0), (22, 22, 22))) == {("A", 1)}


def test_the_box_is_a_cube_not_a_sphere():
    atoms = _res("A", 1, "LYS", (10.5, 10.5, 10.5))                 # 18.2 A from the centre but inside the 22 A cube
    assert ("A", 1) in site_residue_ids(atoms, (0, 0, 0), (22, 22, 22))


@need
@need_tools
def test_the_protonation_table_now_lists_the_ionizable_residues_that_touch_the_substrate(cfg, tmp_path):
    """Regression: Lys107 (and Asp288 in 6X0J) contact the ornithine but lie 9-12 A from C4X, outside the old 8 A sphere."""
    prepare_receptor(CACHE / "4B63.cif", "4B63", "A", cfg, "keep", tmp_path / "r")
    rows = {(r["res_name"], r["res_id"]): r for r in csv.DictReader((tmp_path / "r" / "protonation_states.csv").open())}
    assert ("LYS", "107") in rows and rows[("LYS", "107")]["contacts_ligand"] == "True" and rows[("LYS", "107")]["in_box"] == "True"
    assert {("HIS", "91"), ("LYS", "100"), ("ARG", "144")} <= set(rows)          # the old three are still there
    assert all(r["state"] for r in rows.values())
