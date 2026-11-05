"""Step 4 diagnosis: what is different about a failing pose, and the registered ladder of single-factor changes."""
from pathlib import Path

import pytest
from rdkit import Chem

from sidacore.config import load_config
from sidacore.diagnose import (LIGAND_STATES, atom_displacements, box_variants, crystal_fraction_inside, edit_pdbqt,
                               receptor_variants, render_markdown, run_ladder)
from sidacore.redock import crystal_reference, ligand_from_smiles, ligand_template, prepare_ligand, protonate

ROOT = Path(__file__).resolve().parent.parent
PREPARED = ROOT / "prepared"
JOB = PREPARED / "4B63_A_nadp-keep"
need_prepared = pytest.mark.skipif(not (JOB / "receptor.pdbqt").exists(), reason="run `python -m sidacore.prep` first")


def _names(mol):
    return [a.GetPDBResidueInfo().GetName().strip() for a in mol.GetAtoms()]


def _moved(mol, name, dx):
    out = Chem.Mol(mol)
    i = _names(out).index(name)
    p = out.GetConformer().GetAtomPosition(i)
    out.GetConformer().SetAtomPosition(i, (p.x + dx, p.y, p.z))
    return out


@need_prepared
def test_displacements_of_the_reference_with_itself_are_all_zero():
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    d = dict(atom_displacements(ref, ref))
    assert set(d) == {"N", "CA", "CB", "CG", "CD", "NE", "C", "O", "OXT"} and max(d.values()) < 1e-6


@need_prepared
def test_one_moved_atom_is_the_only_nonzero_displacement():
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    d = dict(atom_displacements(_moved(ref, "NE", 3.0), ref))
    assert d["NE"] == pytest.approx(3.0, abs=1e-3) and all(v < 1e-3 for k, v in d.items() if k != "NE")


@need_prepared
def test_the_equivalent_carboxylate_oxygens_are_not_reported_as_displaced():
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    sw = Chem.Mol(ref)
    n = _names(sw)
    c = sw.GetConformer()
    pi, pj = c.GetAtomPosition(n.index("O")), c.GetAtomPosition(n.index("OXT"))
    c.SetAtomPosition(n.index("O"), pj)
    c.SetAtomPosition(n.index("OXT"), pi)
    assert max(v for _, v in atom_displacements(sw, ref)) < 1e-3


# ------------------------------------------------------------------ rung 1: ligand protonation states
def test_the_baseline_state_is_exactly_what_the_pipeline_docks_at_the_configured_ph():
    base = LIGAND_STATES["baseline +1 (pH 7.4)"]
    assert Chem.CanonSmiles(base) == Chem.CanonSmiles(protonate(ligand_template("ORN"), 7.4))


@pytest.mark.parametrize("label, charge", [("baseline +1 (pH 7.4)", 1), ("delta-amine neutral (net 0)", 0),
                                           ("alpha-amine neutral (net 0)", 0), ("both amines + acid (+2)", 2), ("anion (-1)", -1)])
def test_every_state_is_L_ornithine_with_the_stated_charge(label, charge):
    m = Chem.MolFromSmiles(LIGAND_STATES[label])
    assert Chem.GetFormalCharge(m) == charge and sum(a.GetAtomicNum() > 1 for a in m.GetAtoms()) == 9
    assert Chem.FindMolChiralCenters(m, useLegacyImplementation=False) == [(4, "S")]


def test_a_ligand_can_be_prepared_from_an_explicit_state_and_the_pdbqt_reflects_it():
    lp = ligand_from_smiles(LIGAND_STATES["delta-amine neutral (net 0)"])
    assert lp.net_charge == 0 and lp.n_atoms == 14 and lp.ph is None and lp.n_heavy == 9
    assert lp.partial_charge_sum == pytest.approx(0.0, abs=0.05)


# ------------------------------------------------------------------ rung 3: FAD treatment edits of the receptor PDBQT
@need_prepared
def test_dropping_a_cofactor_removes_exactly_its_atoms():
    txt = (JOB / "receptor.pdbqt").read_text()
    n = lambda t, r: sum(1 for l in t.splitlines() if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() in r)
    out = edit_pdbqt(txt, drop=("FAD", "FDA"))
    assert n(txt, ("FAD", "FDA")) == 61 and n(out, ("FAD", "FDA")) == 0
    assert n(out, ("NAP",)) == n(txt, ("NAP",))


@need_prepared
def test_zeroing_a_cofactors_charges_changes_nothing_else():
    txt = (JOB / "receptor.pdbqt").read_text()
    out = edit_pdbqt(txt, zero_charge=("FAD",))
    a, b = txt.splitlines(), out.splitlines()
    assert len(a) == len(b)
    for x, y in zip(a, b):
        if x.startswith(("ATOM", "HETATM")) and x[17:20].strip() == "FAD":
            assert float(y[70:76]) == 0.0 and x[:70] == y[:70] and x[76:] == y[76:]
        else:
            assert x == y


# ------------------------------------------------------------------ rung 2: boxes
def test_box_variants_vary_one_thing_at_a_time_and_include_the_oracle_marked_as_such():
    vs = box_variants(base_center=(0, 0, 0), base_size=(22, 22, 22), carried_center=(0.2, 0.1, 0), ligand_centroid=(5, 5, 5))
    labels = [v[0] for v in vs]
    assert "size 18 A" in labels and "size 30 A" in labels and any("oracle" in l for l in labels) and any("carried" in l for l in labels)
    by = {v[0]: v for v in vs}
    assert by["size 18 A"][2] == (18, 18, 18) and by["size 18 A"][1] == (0, 0, 0)
    assert by["centre: carried from 6X0I"][2] == (22, 22, 22)


@need_prepared
def test_crystal_fraction_inside_reports_how_much_of_the_ligand_a_box_contains():
    ref = crystal_reference(JOB / "ligand_crystal.pdb", ligand_template("ORN"))
    c = ref.GetConformer().GetPositions().mean(axis=0)
    assert crystal_fraction_inside(ref, c, (22, 22, 22)) == 1.0
    assert crystal_fraction_inside(ref, c + 50.0, (22, 22, 22)) == 0.0


def test_receptor_variants_are_the_two_registered_ph_values():
    assert [v[1] for v in receptor_variants()] == [6.0, 8.5]


# ------------------------------------------------------------------ the ladder end to end (small)
@need_prepared
def test_ladder_runs_a_chosen_variant_on_a_chosen_receptor_and_records_the_pass_rule(tmp_path):
    cfg = load_config(ROOT / "config" / "default.yaml", overrides={"exhaustiveness": 8})
    rep = run_ladder(cfg, PREPARED, ROOT / "boxes", tmp_path, jobs=["4B63_A_nadp-keep"], only=["anion (-1)"], with_baseline=False)
    r = rep["variants"][0]
    assert r["rung"] == "1 ligand protonation" and r["label"] == "anion (-1)" and r["job"] == "4B63_A_nadp-keep"
    assert r["n_seeds"] == 10 and r["accept"] == 2.0 and "passes" in r and len(r["per_seed"]) == 10


@need_prepared
def test_a_receptor_ph_variant_reruns_step_2_and_reports_which_site_states_changed(tmp_path):
    cfg = load_config(ROOT / "config" / "default.yaml", overrides={"exhaustiveness": 8})
    rep = run_ladder(cfg, PREPARED, ROOT / "boxes", tmp_path, jobs=["4B63_A_nadp-keep"], only=["receptor protonation pH 6.0"],
                     with_baseline=False)
    r = rep["variants"][0]
    assert r["rung"] == "4 receptor protonation" and isinstance(r["extras"]["site_states_changed"], dict)
    assert (tmp_path / "4_receptor_protonation" / "receptor_protonation_pH_6_0" / "4B63_A_nadp-keep" / "receptor_prep" / "receptor.pdbqt").exists()


def test_fewer_than_ten_seeds_is_refused_by_the_ladder_too(tmp_path):
    with pytest.raises(ValueError, match="at least 10"):
        run_ladder(load_config(ROOT / "config" / "default.yaml", overrides={"seeds": [1, 2]}), PREPARED, ROOT / "boxes", tmp_path)


def test_markdown_has_one_section_per_rung_one_column_per_receptor_and_flags_a_pass():
    per = [{"seed": i, "top1_rmsd": 1.0, "top1_score": -5.0, "best_top5_rmsd": 0.9, "n_poses": 20} for i in range(1, 11)]
    from sidacore.redock import evaluate_seeds
    ev = evaluate_seeds(per, 2.0, min_seeds=10)
    v = lambda job, rung, label: {"rung": rung, "label": label, "job": job, "engine": "vina", "extras": {}, "wall_seconds": 1.0,
                                  "per_seed": per, **ev}
    rep = {"settings": {"engine": "vina", "exhaustiveness": 32, "num_modes": 20, "seeds": list(range(1, 11)), "acceptance": 2.0,
                        "protonation_ph": 7.4},
           "variants": [v("A_nadp-keep", "1 ligand protonation", "anion (-1)"), v("B_nadp-keep", "1 ligand protonation", "anion (-1)")]}
    md = render_markdown(rep, notes="hello")
    assert "## Rung 1 ligand protonation" in md and "| variant | A | B | notes |" in md
    assert "**PASS**" in md and md.rstrip().endswith("hello")


def _variant(job, label, frac=None, passes=False):
    per = [{"seed": i, "top1_rmsd": 4.8, "top1_score": -5.0, "best_top5_rmsd": 4.8, "n_poses": 20} for i in range(1, 11)]
    from sidacore.redock import evaluate_seeds
    ev = evaluate_seeds(per, 2.0, min_seeds=10)
    extras = {} if frac is None else {"crystal_fraction_inside": frac}
    return {"rung": "2 box", "label": label, "job": job, "engine": "vina", "extras": extras, "wall_seconds": 1.0, "per_seed": per, **ev}


_SETTINGS = {"engine": "vina", "exhaustiveness": 32, "num_modes": 20, "seeds": list(range(1, 11)), "acceptance": 2.0, "protonation_ph": 7.4}


def test_a_box_that_cannot_contain_the_deposited_ligand_is_excluded_not_reported_as_a_failure():
    md = render_markdown({"settings": _SETTINGS, "variants": [_variant("A_nadp-keep", "size 18 A", frac=0.56)]})
    assert "excluded" in md and "56%" in md and "4.80 (0/10)" not in md


def test_the_charge_blind_control_is_labelled_as_a_control():
    v = _variant("A_nadp-keep", "FAD charges zeroed (diagnostic)")
    v["rung"] = "3 FAD treatment"
    md = render_markdown({"settings": _SETTINGS, "variants": [v]})
    assert "control" in md and "no electrostatic" in md
