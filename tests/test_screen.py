"""Generic target x ligand screening: clash measurement, target preparation, ligand preparation, the control gate.

Pure-logic tests use hand-built geometry. Real-data tests use the repo's cached structures and prepared receptors and
skip when those are absent. Tests marked `slow` run PDB2PQR, Meeko and Vina (minutes); run `pytest -m "not slow"` for the rest."""
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from sidacore import screen as S

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "structures_cache"
PREPARED = REPO / "prepared"
needs_cache = pytest.mark.skipif(not (CACHE / "4B63.cif").exists(), reason="structures_cache/4B63.cif not present")
needs_prepared = pytest.mark.skipif(not (PREPARED / "4B63_A_nadp-keep" / "receptor_prepared.pdb").exists(),
                                    reason="validated prepared receptors not present")
def _have(name):
    return bool(shutil.which(name)) or (Path(sys.executable).parent / name).exists()


needs_tools = pytest.mark.skipif(not (_have("vina") and (_have("pdb2pqr30") or Path("/opt/miniforge3/envs/af2/bin/pdb2pqr30").exists())),
                                 reason="vina / pdb2pqr30 not available")


def atoms(rows):
    """rows: (res_name, res_id, atom_name, element, x, y, z) -> Atoms in chain A."""
    return S.Atoms(coord=np.array([r[4:7] for r in rows], float), element=np.array([r[3].upper() for r in rows]),
                   res_name=np.array([r[0] for r in rows]), res_id=np.array([r[1] for r in rows]),
                   chain=np.array(["A"] * len(rows)), atom_name=np.array([r[2] for r in rows]))


def pose_from_xyz(smiles, xyz):
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, tuple(float(x) for x in p))
    mol.AddConformer(conf)
    return mol


# --------------------------------------------------------------------------- clash arithmetic
LIG_C = atoms([("LIG", 1, "C1", "C", 0.0, 0.0, 0.0)])


def test_a_carbon_carbon_pair_clashes_when_overlap_exceeds_the_tolerance():
    rec = atoms([("LEU", 5, "CD1", "C", 2.9, 0.0, 0.0)])                  # sum 3.40, d 2.9 -> overlap 0.50 > 0.4
    got = S.find_clashes(LIG_C, rec)
    assert len(got) == 1
    assert (got[0].res_name, got[0].res_id, got[0].rec_atom) == ("LEU", 5, "CD1")
    assert got[0].overlap == pytest.approx(0.50)


def test_a_pair_inside_the_tolerance_is_not_a_clash():
    rec = atoms([("LEU", 5, "CD1", "C", 3.1, 0.0, 0.0)])                  # overlap 0.30 <= 0.4
    assert S.find_clashes(LIG_C, rec) == []


def test_the_tolerance_is_a_parameter():
    rec = atoms([("LEU", 5, "CD1", "C", 3.1, 0.0, 0.0)])
    assert len(S.find_clashes(LIG_C, rec, tolerance=0.1)) == 1
    assert S.find_clashes(LIG_C, rec, tolerance=0.5) == []


def test_polar_pairs_are_judged_by_a_minimum_distance_not_by_vdw_overlap():
    lig = atoms([("LIG", 1, "N1", "N", 0.0, 0.0, 0.0)])
    ok = atoms([("SER", 7, "OG", "O", 2.6, 0.0, 0.0)])                    # an ordinary hydrogen bond: overlap 0.47 but allowed
    bad = atoms([("SER", 7, "OG", "O", 2.2, 0.0, 0.0)])
    assert S.find_clashes(lig, ok) == []
    assert len(S.find_clashes(lig, bad)) == 1
    assert len(S.find_clashes(lig, bad, polar_min=2.0)) == 0


def test_metals_and_receptor_hydrogens_never_clash():
    rec = atoms([("ZN", 900, "ZN", "ZN", 2.0, 0.0, 0.0), ("LEU", 5, "HD11", "H", 1.0, 0.0, 0.0)])
    assert S.find_clashes(LIG_C, rec) == []


def test_clashes_are_sorted_worst_first_and_name_the_receptor_residue():
    rec = atoms([("LEU", 5, "CD1", "C", 2.9, 0.0, 0.0), ("PHE", 9, "CZ", "C", 0.0, 2.5, 0.0)])
    got = S.find_clashes(LIG_C, rec)
    assert [c.res_id for c in got] == [9, 5]
    assert "PHE9" in got[0].label() and "CZ" in got[0].label()


def test_a_pose_far_from_the_receptor_has_no_clashes():
    assert S.find_clashes(atoms([("LIG", 1, "C1", "C", 30.0, 0.0, 0.0)]), atoms([("LEU", 5, "CD1", "C", 0.0, 0.0, 0.0)])) == []


def test_pose_atoms_reads_heavy_atoms_and_elements_from_an_rdkit_pose():
    pose = pose_from_xyz("CCO", [(0, 0, 0), (1.5, 0, 0), (2.4, 1.0, 0)])
    a = S.pose_atoms(pose)
    assert list(a.element) == ["C", "C", "O"]
    assert np.allclose(a.coord[2], (2.4, 1.0, 0))


def test_pose_pdb_text_round_trips_an_aromatic_pose_without_rdkit_kekulisation(tmp_path):
    ring = [(1.4 * np.cos(a), 1.4 * np.sin(a), 0.0) for a in np.linspace(0, 2 * np.pi, 7)[:-1]]
    pose = pose_from_xyz("c1ccccc1", ring)
    f = tmp_path / "p.pdb"
    f.write_text(S.pose_pdb_text(pose))
    back = S.read_pdb_heavy(f)
    assert list(back.element) == ["C"] * 6 and np.allclose(back.coord, np.array(ring), atol=1e-3)
    assert sum(l.startswith("CONECT") for l in f.read_text().splitlines()) == 6        # every ring atom lists its two neighbours


# --------------------------------------------------------------------------- site metrics
REC_SITE = atoms([("AAA", 1, "CA", "C", 0.0, 0.0, 0.0), ("BBB", 2, "CA", "C", 4.0, 0.0, 0.0), ("CCC", 3, "CA", "C", 40.0, 0.0, 0.0)])


def test_a_pose_on_top_of_the_crystal_ligand_has_zero_displacement_and_full_site_contact():
    crystal = atoms([("LIG", 1, "C1", "C", 2.0, 0.0, 0.0)])
    m = S.site_metrics(crystal, crystal, REC_SITE)
    assert m["centroid_dist"] == pytest.approx(0.0)
    assert m["site_contact_fraction"] == pytest.approx(1.0)


def test_a_pose_in_another_part_of_the_box_misses_the_site():
    crystal = atoms([("LIG", 1, "C1", "C", 2.0, 0.0, 0.0)])
    far = atoms([("LIG", 1, "C1", "C", 40.0, 3.0, 0.0)])                # touches only residue 3
    m = S.site_metrics(far, crystal, REC_SITE)
    assert m["centroid_dist"] > 30
    assert m["site_contact_fraction"] == pytest.approx(0.0)


def test_the_site_fraction_is_the_share_of_the_deposited_contacts_that_the_pose_also_touches():
    crystal = atoms([("LIG", 1, "C1", "C", 2.0, 0.0, 0.0)])                                   # touches residues 1 and 2
    sprawling = atoms([("LIG", 1, "C1", "C", 2.0, 0.0, 0.0), ("LIG", 1, "C2", "C", 40.0, 3.0, 0.0)])   # those two plus residue 3
    half = atoms([("LIG", 1, "C1", "C", -2.0, 0.0, 0.0)])                                     # touches residue 1 only
    assert S.site_metrics(sprawling, crystal, REC_SITE)["site_contact_fraction"] == pytest.approx(1.0)
    assert S.site_metrics(half, crystal, REC_SITE)["site_contact_fraction"] == pytest.approx(0.5)


# --------------------------------------------------------------------------- the plain-language reading
SET = S.Settings()


def row(**kw):
    base = {"n_clashes": 0, "n_clash_free": 60, "n_poses": 60, "site_contact_fraction": 0.9}
    base.update(kw)
    return base


def test_reading_says_clashes_in_every_pose_when_no_pose_is_clash_free():
    assert S.reading(row(n_clashes=8, n_clash_free=0), "validated", SET).startswith("clashes in every pose")


def test_reading_says_the_best_pose_clashes_when_a_clash_free_pose_exists():
    assert S.reading(row(n_clashes=2, n_clash_free=5), "validated", SET).startswith("best pose clashes")


def test_reading_calls_a_clash_free_pose_at_the_rim_outside_the_site_not_a_fit():
    assert S.reading(row(site_contact_fraction=0.33), "validated", SET).startswith("clash-free but outside the site")


def test_reading_says_fits_for_a_clash_free_pose_in_the_site():
    assert S.reading(row(), "validated", SET) == "fits: clash-free, in the site"


def test_reading_is_prefixed_when_the_target_failed_its_control():
    assert S.reading(row(), "FAILED", SET).startswith("(target failed its control) ")


def test_reading_without_a_site_ligand_does_not_judge_site_occupancy():
    r = S.reading(row(site_contact_fraction=float("nan")), "no control (no site ligand)", SET)
    assert "not assessed" in r


# --------------------------------------------------------------------------- ligands in
def test_ligands_load_from_a_dict_a_csv_and_a_smi_file(tmp_path):
    assert S.load_ligands({"a": "CCO"}) == {"a": "CCO"}
    (tmp_path / "x.csv").write_text("name,smiles\nethanol,CCO\nglycine,NCC(=O)O\n")
    assert S.load_ligands(tmp_path / "x.csv") == {"ethanol": "CCO", "glycine": "NCC(=O)O"}
    (tmp_path / "y.smi").write_text("CCO ethanol\nNCC(=O)O glycine\n")
    assert S.load_ligands(tmp_path / "y.smi") == {"ethanol": "CCO", "glycine": "NCC(=O)O"}


def test_duplicate_ligand_names_are_refused(tmp_path):
    (tmp_path / "x.csv").write_text("name,smiles\na,CCO\na,CCN\n")
    with pytest.raises(ValueError, match="duplicate"):
        S.load_ligands(tmp_path / "x.csv")


def test_an_unreadable_smiles_is_recorded_as_a_failure_not_raised():
    got = S.prepare_ligands({"bad": "not-a-smiles", "ok": "CCO"}, ph=None)
    by = {e.name: e for e in got}
    assert by["bad"].error and by["bad"].pdbqt is None
    assert by["ok"].error is None and by["ok"].pdbqt and by["ok"].n_heavy == 3


def test_ligand_warnings_cover_undefined_stereo_and_high_flexibility():
    got = {e.name: e for e in S.prepare_ligands({"ala": "CC(N)C(=O)O", "chain": "CCCCCCCCCCCCCCCC(=O)O", "ok": "C[C@H](N)C(=O)O"}, ph=None)}
    assert any("stereo" in w for w in got["ala"].warnings)
    assert any("rotatable" in w for w in got["chain"].warnings)
    assert got["ok"].warnings == ()


def test_protonation_at_a_ph_changes_the_charge_state_of_an_amino_acid():
    zw = S.prepare_ligands({"gly": "NCC(=O)O"}, ph=7.4)[0]
    neutral = S.prepare_ligands({"gly": "NCC(=O)O"}, ph=None)[0]
    assert zw.net_charge == 0 and "[O-]" in zw.smiles          # zwitterion
    assert neutral.smiles == "NCC(=O)O"


# --------------------------------------------------------------------------- target selection (real structure, no tools)
@needs_cache
def test_heterogens_of_4b63_chain_a_are_listed_with_waters_counted_apart():
    from sidacore.prep import read_structure
    a, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    het = {h["res_name"]: h for h in S.list_heterogens(a, "A")}
    assert {"FAD", "NAP", "ORN", "SO4", "GOL"} <= set(het)
    assert "HOH" not in het and S.count_waters(a, "A") > 100
    assert het["GOL"]["n_residues"] == 4 and het["SO4"]["n_residues"] == 1


@needs_cache
def test_auto_keeps_the_flavin_and_nadp_and_strips_additives_and_sets_the_site_ligand_aside():
    from sidacore.prep import read_structure
    a, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    sel = S.select_target_atoms(a, S.TargetSpec("4B63", "A", site_ligand="ORN"))
    assert sorted(sel.kept_cofactors) == ["FAD", "NAP"]
    assert sel.ligand is not None and set(sel.ligand.res_name) == {"ORN"}
    reasons = {r["res_name"]: r["reason"] for r in sel.removed}
    assert "additive" in reasons["SO4"] and "additive" in reasons["GOL"] and "water" in reasons["HOH"]
    assert not np.isin(sel.kept.res_name, ["SO4", "GOL", "HOH", "ORN"]).any()


@needs_cache
def test_an_explicit_keep_list_overrides_auto():
    from sidacore.prep import read_structure
    a, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    sel = S.select_target_atoms(a, S.TargetSpec("4B63", "A", site_ligand="ORN", keep=["FAD"]))
    assert sel.kept_cofactors == ["FAD"]
    assert any(r["res_name"] == "NAP" for r in sel.removed)


@needs_cache
def test_the_box_is_centred_on_the_site_ligand_unless_a_centre_is_given():
    from sidacore.prep import read_structure
    a, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    sel = S.select_target_atoms(a, S.TargetSpec("4B63", "A", site_ligand="ORN"))
    assert np.allclose(sel.center, sel.ligand.coord[sel.ligand.element != "H"].mean(axis=0))
    given = S.select_target_atoms(a, S.TargetSpec("4B63", "A", site_ligand="ORN", center=(1.0, 2.0, 3.0)))
    assert np.allclose(given.center, (1.0, 2.0, 3.0))


@needs_cache
def test_a_missing_site_ligand_and_no_centre_is_an_error():
    from sidacore.prep import read_structure
    a, _ = read_structure(CACHE / "6X0H.cif", "highest-occupancy")        # apo of ornithine
    with pytest.raises(ValueError, match="site_ligand|centre|center"):
        S.select_target_atoms(a, S.TargetSpec("6X0H", "A", site_ligand="ORN"))


@needs_cache
def test_the_default_box_is_site_sized_and_a_box_too_small_for_the_site_ligand_is_noted():
    from sidacore.prep import read_structure
    assert S.TargetSpec("4B63", "A").box_size == 12.0
    a, _ = read_structure(CACHE / "4B63.cif", "highest-occupancy")
    assert S.select_target_atoms(a, S.TargetSpec("4B63", "A", site_ligand="ORN")).notes == []
    tiny = S.select_target_atoms(a, S.TargetSpec("4B63", "A", site_ligand="ORN", box_size=4.0))
    assert any("extends outside the box" in n for n in tiny.notes)


@needs_cache
def test_preview_target_reports_kept_and_removed_without_running_any_tool():
    sel = S.preview_target(S.TargetSpec("4B63", "A", site_ligand="ORN"), CACHE)
    assert sel.kept_cofactors == ["FAD", "NAP"] and {r["res_name"] for r in sel.removed} >= {"SO4", "GOL", "HOH", "ORN"}


def test_check_tools_names_what_is_missing_and_finds_what_is_installed():
    got = S.check_tools()
    assert set(got) == {"vina", "obabel", "pdb2pqr30", "mk_prepare_receptor.py"}
    assert (Path(sys.executable).parent / "vina").exists() == (got["vina"] is not None)


def test_a_cofactor_without_a_meeko_template_is_refused_early():
    with pytest.raises(ValueError, match="template"):
        S.check_keep(["HEM"])


# --------------------------------------------------------------------------- calibration on crystal poses
@needs_prepared
@pytest.mark.parametrize("job", ["4B63_A_nadp-keep", "6X0J_A_nadp-keep", "6X0J_B_nadp-keep", "6X0J_C_nadp-keep", "6X0J_D_nadp-keep"])
def test_the_deposited_ornithine_has_no_clashes_at_the_default_tolerance(job):
    d = PREPARED / job
    rec = S.read_pdb_heavy(d / "receptor_prepared.pdb")
    lig = S.read_pdb_heavy(d / "ligand_crystal.pdb")
    assert S.find_clashes(lig, rec) == []


@needs_prepared
def test_a_crystal_ligand_pushed_into_the_protein_is_flagged():
    d = PREPARED / "4B63_A_nadp-keep"
    rec = S.read_pdb_heavy(d / "receptor_prepared.pdb")
    lig = S.read_pdb_heavy(d / "ligand_crystal.pdb")
    # a deliberately wrong pose: move the ligand 1.5 A toward its closest receptor atom
    dist = np.linalg.norm(rec.coord[:, None, :] - lig.coord[None, :, :], axis=2)
    j, i = np.unravel_index(dist.argmin(), dist.shape)
    v = rec.coord[j] - lig.coord[i]
    moved = S.Atoms(lig.coord + 1.5 * v / np.linalg.norm(v), lig.element, lig.res_name, lig.res_id, lig.chain, lig.atom_name)
    assert len(S.find_clashes(moved, rec)) >= 1


def test_reading_meeko_style_pdb_drops_hydrogens_and_reads_fixed_columns(tmp_path):
    f = tmp_path / "r.pdb"
    f.write_text("ATOM      1  C   PRO A  29     -33.197 -32.028  42.546                       C \n"
                 "ATOM      2  HA  PRO A  29     -33.504 -33.206  44.267                       H \n"
                 "ATOM      3 HO4' FAD A1492     -22.353 -28.204  20.289                       H \n"
                 "ATOM      4  O4' FAD A1492     -22.325 -28.707  19.437                       O \n")
    a = S.read_pdb_heavy(f)
    assert list(a.atom_name) == ["C", "O4'"] and list(a.res_name) == ["PRO", "FAD"] and list(a.res_id) == [29, 1492]
    assert a.coord[0] == pytest.approx((-33.197, -32.028, 42.546))


# --------------------------------------------------------------------------- slow: the tools
@pytest.mark.slow
@needs_cache
@needs_prepared
@needs_tools
def test_generic_preparation_of_4b63_reproduces_the_validated_receptor(tmp_path):
    t = S.prepare_target(S.TargetSpec("4B63", "A", site_ligand="ORN"), cache_dir=CACHE, work_dir=tmp_path)
    ref = (PREPARED / "4B63_A_nadp-keep" / "receptor.pdbqt").read_text().splitlines()
    got = t.receptor_pdbqt.read_text().splitlines()
    atom = lambda L: [l for l in L if l.startswith(("ATOM", "HETATM"))]
    assert len(atom(got)) == len(atom(ref))
    assert [l[12:30] for l in atom(got)] == [l[12:30] for l in atom(ref)]        # same atoms, residues, numbering
    assert t.kept_cofactors == ["FAD", "NAP"]
    assert {r["res_name"] for r in t.removed} >= {"SO4", "GOL", "HOH"}
    assert t.crystal_ligand_pdb.stat().st_size > 0


@pytest.mark.slow
@needs_cache
@needs_tools
def test_preparing_the_same_target_twice_reuses_it_and_a_changed_spec_is_refused(tmp_path):
    spec = S.TargetSpec("4B63", "A", site_ligand="ORN")
    first = S.prepare_target(spec, cache_dir=CACHE, work_dir=tmp_path)
    stamp = first.receptor_pdbqt.stat().st_mtime_ns
    again = S.prepare_target(spec, cache_dir=CACHE, work_dir=tmp_path)
    assert again.receptor_pdbqt.stat().st_mtime_ns == stamp
    with pytest.raises(ValueError, match="different settings"):
        S.prepare_target(S.TargetSpec("4B63", "A", site_ligand="ORN", keep=["FAD"]), cache_dir=CACHE, work_dir=tmp_path)


@pytest.mark.slow
@needs_cache
@needs_tools
def test_a_small_screen_runs_end_to_end_and_the_oversized_ligand_clashes_where_the_substrate_does_not(tmp_path):
    settings = S.Settings(seeds=(1, 2), exhaustiveness=8, num_modes=10, control_seeds=(1, 2, 3))
    spec = S.TargetSpec("4B63", "A", site_ligand="ORN")
    res = S.run_screen([spec], {"ornithine": "NCCC[C@H](N)C(=O)O",
                                "tetraphenylmethane": "C(c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1", "bad": "zzz"},
                       settings=settings, cache_dir=CACHE, work_dir=tmp_path)
    r = res.results.set_index("ligand")
    assert r.loc["bad", "status"].startswith("failed")
    assert r.loc["ornithine", "status"] == "ok" and r.loc["tetraphenylmethane", "status"] == "ok"
    assert r.loc["ornithine", "n_clashes"] == 0
    assert r.loc["tetraphenylmethane", "n_clashes"] > 0 and r.loc["tetraphenylmethane", "n_clash_free"] == 0
    assert "clashes in every pose" in r.loc["tetraphenylmethane", "reading"]
    assert res.controls.loc[0, "target"] == "4B63_A" and "median_top1_rmsd" in res.controls.columns
    assert (tmp_path / "results.csv").exists() and (tmp_path / "run.json").exists()
    assert json.loads((tmp_path / "run.json").read_text())["settings"]["seeds"] == [1, 2]
    assert (tmp_path / "poses" / "4B63_A__tetraphenylmethane.pdb").stat().st_size > 0
