"""Step 3: box by C4X, the same box carried into other structures by superposition, and an fpocket cross-check."""
import shutil
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from sidacore.boxes import (Pocket, analyze_all, assess_pockets, box_overlap_fraction, heavy_atom_pdb, ligand_pocket_contact, main as boxes_main,
                            ligand_vs_box, parse_fpocket_info, pocket_lining, protein_chain, read_vina_box, render_markdown, run_fpocket,
                            transfer_site, write_vina_box)
from sidacore.config import load_config
from sidacore.prep import read_structure, resolve_site_center

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "structures_cache"
DEFAULT = ROOT / "config" / "default.yaml"
need = pytest.mark.skipif(not all((CACHE / ("%s.cif" % i)).exists() for i in ("4B63", "6X0H", "6X0I", "6X0J")),
                          reason="structures_cache incomplete")
FPOCKET = Path("/opt/miniforge3/envs/sida-docking/bin/fpocket")
need_fpocket = pytest.mark.skipif(not FPOCKET.exists(), reason="fpocket not installed")


# ------------------------------------------------------------------ box geometry (pure)
def test_overlap_is_one_for_identical_boxes_and_zero_for_disjoint():
    c = np.zeros(3)
    assert box_overlap_fraction(c, c, (22, 22, 22)) == pytest.approx(1.0)
    assert box_overlap_fraction(c, c + [30, 0, 0], (22, 22, 22)) == 0.0


def test_overlap_of_boxes_shifted_by_half_a_width_is_one_half():
    assert box_overlap_fraction(np.zeros(3), np.array([11.0, 0, 0]), (22, 22, 22)) == pytest.approx(0.5)


def test_overlap_shifted_on_two_axes_multiplies():
    assert box_overlap_fraction(np.zeros(3), np.array([11.0, 11.0, 0]), (22, 22, 22)) == pytest.approx(0.25)


# ------------------------------------------------------------------ superposition on a synthetic chain
def _chain(n=40, seed=0):
    rng = np.random.default_rng(seed)
    names = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
             "SER", "THR", "TRP", "TYR", "VAL"]
    atoms = struc.AtomArray(n)
    atoms.coord = np.cumsum(rng.normal(0, 2.0, (n, 3)), axis=0)
    atoms.atom_name[:] = "CA"
    atoms.element[:] = "C"
    atoms.chain_id[:] = "A"
    atoms.res_id[:] = np.arange(1, n + 1)
    atoms.res_name[:] = [names[i] for i in rng.integers(0, 20, n)]
    atoms.hetero[:] = False
    return atoms


def test_transfer_site_reproduces_a_known_rigid_motion():
    ref = _chain()
    rot = Rotation.from_euler("xyz", [31, -47, 112], degrees=True)
    shift = np.array([10.0, -20.0, 5.0])
    qry = ref.copy()
    # 0.1 A of noise: an exactly rigid copy has a zero interquartile range, so float dust would be trimmed as outliers
    qry.coord = rot.apply(ref.coord) + shift + np.random.default_rng(1).normal(0, 0.1, ref.coord.shape)
    site = ref.coord[10] + [1.0, 0.5, -0.3]                      # a point near, not on, a residue
    got = transfer_site(ref, qry, site, radius=8.0, outlier_threshold=1.5)
    np.testing.assert_allclose(got.center_in_query, rot.apply(site) + shift, atol=0.2)
    assert got.rmsd_anchors < 0.2
    assert got.n_aligned == 40 and got.n_anchors >= 30


def test_a_residue_that_moved_is_reported_in_site_rmsd_even_if_trimmed_as_an_outlier():
    """Trimming must not hide a conformational change: the site RMSD is over ALL aligned pairs near the site."""
    ref = _chain()
    qry = ref.copy()
    qry.coord = ref.coord + [3.0, 0, 0] + np.random.default_rng(2).normal(0, 0.1, ref.coord.shape)
    qry.coord[10] += [0.0, 9.0, 0.0]                              # one residue next to the site moves 9 A
    site = ref.coord[10].copy()
    got = transfer_site(ref, qry, site, radius=6.0, outlier_threshold=1.5)
    assert got.n_trimmed >= 1                                     # it was dropped from the anchors ...
    assert got.site_max_displacement > 8.0                        # ... but the site numbers still show it
    assert got.n_site >= 1


def test_the_reference_carried_onto_itself_is_the_identity_with_nothing_trimmed():
    ref = _chain()
    got = transfer_site(ref, ref.copy(), ref.coord[3] + 0.5, radius=8.0, outlier_threshold=1.5)
    np.testing.assert_allclose(got.center_in_query, ref.coord[3] + 0.5, atol=1e-6)
    assert got.n_trimmed == 0 and got.n_anchors == 40 and got.rmsd_anchors == 0.0 and got.site_rmsd == 0.0


# ------------------------------------------------------------------ pocket agreement (pure)
def _pocket(rank, centre, n=30, spread=1.0, seed=0):
    rng = np.random.default_rng(seed + rank)
    sp = np.asarray(centre, float) + rng.normal(0, spread, (n, 3))
    return Pocket(rank=rank, score=1.0 / rank, druggability=0.5, volume=100.0 * n, spheres=sp)


def test_parse_fpocket_info_reads_every_pocket_in_order():
    text = ("Pocket 1 :\n\tScore : \t0.927\n\tDruggability Score : \t0.971\n\tNumber of Alpha Spheres : \t171\n"
            "\tVolume : \t1603.625\n\nPocket 2 :\n\tScore : \t0.5\n\tDruggability Score : \t0.1\n"
            "\tNumber of Alpha Spheres : \t20\n\tVolume : \t300.0\n")
    rows = parse_fpocket_info(text)
    assert [r["rank"] for r in rows] == [1, 2]
    assert rows[0] == {"rank": 1, "score": 0.927, "druggability": 0.971, "n_spheres": 171, "volume": 1603.625}


def test_agreement_when_the_top_pocket_holds_the_box_and_is_centred():
    pockets = [_pocket(1, [0, 0, 0]), _pocket(2, [40, 0, 0])]
    a = assess_pockets(pockets, np.zeros(3), (22, 22, 22), max_offset=5.0)
    assert a.site_rank == 1 and a.agrees and a.reasons == []


def test_disagreement_when_the_box_sits_on_a_lower_ranked_pocket():
    pockets = [_pocket(1, [40, 0, 0]), _pocket(2, [0, 0, 0])]
    a = assess_pockets(pockets, np.zeros(3), (22, 22, 22), max_offset=5.0)
    assert a.site_rank == 2 and not a.agrees
    assert any("rank 2" in r for r in a.reasons)


def test_disagreement_when_the_pocket_centroid_is_off_the_box_centre():
    pockets = [_pocket(1, [8, 0, 0])]
    a = assess_pockets(pockets, np.zeros(3), (22, 22, 22), max_offset=5.0)
    assert a.site_rank == 1 and not a.agrees
    assert a.centroid_offset == pytest.approx(8.0, abs=1.0)
    assert any("offset" in r for r in a.reasons)


def test_no_pocket_in_the_box_is_a_disagreement_not_a_crash():
    a = assess_pockets([_pocket(1, [60, 0, 0])], np.zeros(3), (22, 22, 22), max_offset=5.0)
    assert a.site_rank is None and not a.agrees and a.reasons


# ------------------------------------------------------------------ ligand vs box (pure)
def test_ligand_vs_box_counts_atoms_inside_and_centroid_distance():
    lig = struc.AtomArray(4)
    lig.coord = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [30, 0, 0]], float)
    lig.element[:] = "C"
    r = ligand_vs_box(lig, np.zeros(3), (22, 22, 22))
    assert r["n_heavy"] == 4 and r["n_inside"] == 3
    assert r["centroid_distance"] == pytest.approx(np.linalg.norm(lig.coord.mean(0)))


def test_ligand_pocket_contact_names_the_pocket_that_owns_the_nearest_spheres():
    pockets = [_pocket(1, [30, 0, 0]), _pocket(2, [0, 0, 0])]
    lig = struc.AtomArray(3)
    lig.coord = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    lig.element[:] = "C"
    r = ligand_pocket_contact(lig, pockets)
    assert r["nearest_sphere_owner"] == {2: 3}
    assert r["nearest_sphere_distance_max"] < 3.0


def test_ligand_vs_box_ignores_hydrogens():
    lig = struc.AtomArray(2)
    lig.coord = np.array([[0, 0, 0], [50, 0, 0]], float)
    lig.element[:] = ["C", "H"]
    assert ligand_vs_box(lig, np.zeros(3), (22, 22, 22))["n_heavy"] == 1


# ------------------------------------------------------------------ real structures
@pytest.fixture(scope="module")
def cfg():
    return load_config(DEFAULT)


@pytest.fixture(scope="module")
def sites():
    out = {}
    for pid in ("4B63", "6X0H", "6X0I", "6X0J"):
        atoms, _ = read_structure(CACHE / ("%s.cif" % pid), "highest-occupancy")
        out[pid] = (protein_chain(atoms, "A"), resolve_site_center(atoms, "A", "FAD:C4X"))
    return out


@need
def test_protein_chain_excludes_ornithine_and_cofactors(sites):
    prot, _ = sites["4B63"]
    assert "ORN" not in set(prot.res_name) and "FAD" not in set(prot.res_name)
    assert set(prot.atom_name) >= {"CA"}


@need
def test_fad_in_structures_agree_with_each_other_expectation_registered_before_measuring(sites, cfg):
    """Pre-registered: carrying 6X0I's C4X into 6X0J (also FAD-in) lands within 2 A of 6X0J's own C4X."""
    ref, ref_c = sites["6X0I"]
    qry, qry_c = sites["6X0J"]
    got = transfer_site(ref, qry, ref_c, cfg.pocket_radius_angstrom, cfg.superposition_outlier_threshold)
    assert np.linalg.norm(got.center_in_query - qry_c) < 2.0


@need
def test_transfer_across_the_conformational_pair_is_well_anchored(sites, cfg):
    ref, ref_c = sites["6X0I"]
    qry, _ = sites["6X0H"]
    got = transfer_site(ref, qry, ref_c, cfg.pocket_radius_angstrom, cfg.superposition_outlier_threshold)
    n_ca = min((ref.atom_name == "CA").sum(), (qry.atom_name == "CA").sum())
    assert got.n_aligned >= 0.9 * n_ca                      # the alignment covers the chain, not a fragment
    assert got.n_anchors >= 0.5 * got.n_aligned and got.rmsd_anchors < 3.0
    assert np.isfinite(got.center_in_query).all()


# ------------------------------------------------------------------ fpocket and the whole step, on the real receptors
PREPARED = ROOT / "prepared"
need_prepared = pytest.mark.skipif(not (PREPARED / "6X0I_A_nadp-keep" / "prep_report.json").exists(),
                                   reason="run `python -m sidacore.prep` first")


def test_vina_box_file_roundtrip(tmp_path):
    write_vina_box(tmp_path / "b.txt", [1.5, -2.25, 3.0], [22, 22, 20])
    c, sz = read_vina_box(tmp_path / "b.txt")
    np.testing.assert_allclose(c, [1.5, -2.25, 3.0])
    np.testing.assert_allclose(sz, [22, 22, 20])


@need_prepared
@need_fpocket
def test_fpocket_on_the_prepared_receptor_returns_ranked_pockets_with_their_spheres(tmp_path):
    job = PREPARED / "6X0I_A_nadp-keep"
    heavy, _ = heavy_atom_pdb(job / "work" / "combined.pdb", tmp_path / "heavy.pdb")
    atoms = pdb.PDBFile.read(str(heavy)).get_structure(model=1)
    assert (atoms.element != "H").all() and (atoms.element == "").sum() == 0 and "FAD" in set(atoms.res_name)
    pockets = run_fpocket(heavy, tmp_path / "fp.log")
    assert [p.rank for p in pockets] == list(range(1, len(pockets) + 1))
    assert pockets[0].volume > 100 and len(pockets[0].spheres) > 15
    lining = pocket_lining(heavy, 1)                                # fpocket's own list of lining atoms
    assert lining["n_protein_residues"] > 5 and lining["cofactor_atoms"].get("FAD", 0) > 0


@need_prepared
@need_fpocket
def test_fpocket_sees_the_reduced_flavin_of_6X0J(tmp_path):
    """fpocket silently drops hetero residues it does not know by name; it knows FAD but not FDA (reduced FAD).
    Without the relabel, 6X0J's pockets were computed with the flavin missing and FAD never lined the site pocket."""
    job = PREPARED / "6X0J_A_nadp-keep"
    heavy, rec = heavy_atom_pdb(job / "work" / "combined.pdb", tmp_path / "heavy.pdb")
    assert rec["relabelled_atoms"] == {"FDA->FAD": 53} and rec["hetero_unverified_in_fpocket"] == []
    pockets = run_fpocket(heavy, tmp_path / "fp.log")
    assert pocket_lining(heavy, pockets[0].rank)["cofactor_atoms"].get("FAD", 0) > 0


@need
@need_prepared
@need_fpocket
def test_analyze_all_end_to_end(cfg, tmp_path):
    rep = analyze_all(cfg, PREPARED, CACHE, tmp_path / "boxes")
    by = {j["job"]: j for j in rep["jobs"]}
    assert len(by) == len(list(PREPARED.glob("*/prep_report.json")))
    ref = by["6X0I_A_nadp-keep"]                                   # the reference carried onto itself must not move
    assert ref["offset_c4x_vs_carried"] < 1e-3 and ref["box_overlap_c4x_vs_carried"] == pytest.approx(1.0)
    assert ref["superposition"]["n_trimmed"] == 0
    assert (tmp_path / "boxes" / "resolved_config.yaml").exists() and (tmp_path / "boxes" / "step3_report.json").exists()
    assert by["4B63_A_nadp-keep"]["ligand"] is not None and by["6X0H_A_nadp-keep"]["ligand"] is None
    md = render_markdown(rep, notes="a note")
    assert "6X0H_A_nadp-keep" in md and "carried" in md and md.rstrip().endswith("a note")


@need
@need_prepared
@need_fpocket
def test_a_step2_box_that_disagrees_with_the_recomputed_c4x_is_refused(cfg, tmp_path):
    import shutil as sh
    root = tmp_path / "prep"
    sh.copytree(PREPARED / "6X0I_A_nadp-keep", root / "6X0I_A_nadp-keep")
    box = root / "6X0I_A_nadp-keep" / "receptor.box.txt"
    box.write_text(box.read_text().replace("center_x = ", "center_x = 9"))
    with pytest.raises(RuntimeError, match="disagrees with the C4X"):
        analyze_all(cfg, root, CACHE, tmp_path / "out")


def test_render_only_rebuilds_the_step3_report_from_the_retained_json_without_rerunning_fpocket(tmp_path):
    import json
    src = ROOT / "boxes" / "step3_report.json"
    if not src.exists():
        pytest.skip("run `python -m sidacore.boxes` first")
    out = tmp_path / "boxes"
    out.mkdir()
    (out / "step3_report.json").write_text(src.read_text())
    from sidacore.boxes import main as bmain
    assert bmain(["--render-only", str(out), "--report", str(tmp_path / "r.md")]) == 0
    assert "## Boxes" in (tmp_path / "r.md").read_text()


def test_site_rmsd_can_be_restricted_to_named_reference_residues():
    ref = _chain()
    qry = ref.copy()
    qry.coord = ref.coord + [3.0, 0, 0] + np.random.default_rng(3).normal(0, 0.1, ref.coord.shape)
    qry.coord[20] += [0.0, 5.0, 0.0]                              # residue 21 moves 5 A
    site = {("A", 21), ("A", 22)}                                 # res_id = index + 1
    got = transfer_site(ref, qry, ref.coord[20], radius=1.0, outlier_threshold=1.5, site_residues=site)
    assert got.n_site == 2 and got.site_max_displacement > 4.0
    far = transfer_site(ref, qry, ref.coord[20], radius=1.0, outlier_threshold=1.5, site_residues={("A", 5)})
    assert far.n_site == 1 and far.site_max_displacement < 1.0
