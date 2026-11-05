"""Where does one residue's side chain sit in each structure, relative to a common backbone? (biotite superposition and dihedrals.)

Residues are matched by number and atom name. That is the published residue mapping for this series: the construct check (reports/step3, `mapping`) found
4B63 and 6X0H/I/J numbered identically with 0 mismatches over the aligned residues, so no renumbering happens here.
"""
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np

CHI_ATOMS = {"LYS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "CE"), ("CG", "CD", "CE", "NZ")]}
CARBOXYLATE = ("O", "OXT")
POLAR_PARTNER_CUTOFF = 3.5


def _atom(res, name):
    sel = res[res.atom_name == name]
    if len(sel) == 0:
        raise ValueError("atom %s missing from residue %s %s" % (name, res.res_name[0], res.res_id[0]))
    return sel.coord[0]


def chi_angles(res) -> list:
    """Side-chain torsions in degrees (chi1..chi4 for lysine)."""
    table = CHI_ATOMS[str(res.res_name[0])]
    return [float(np.degrees(struc.dihedral(*[_atom(res, n) for n in quad]))) for quad in table]


def _matched(ref, qry, anchors):
    """Index arrays into ref and qry of atoms present in both (same residue number and atom name), limited to `anchors` atom names if given."""
    pos = {}
    for i, (r, n) in enumerate(zip(ref.res_id, ref.atom_name)):
        if anchors is None or n in anchors:
            pos[(int(r), str(n))] = i
    ri, qi = [], []
    for j, (r, n) in enumerate(zip(qry.res_id, qry.atom_name)):
        i = pos.get((int(r), str(n)))
        if i is not None:
            ri.append(i)
            qi.append(j)
    return np.array(ri, dtype=int), np.array(qi, dtype=int)


def fit_transform(ref, qry, anchors=None, within=None):
    """Rigid fit of `qry` onto `ref` over matched anchor atoms. Returns (transform, n_anchors, fit_rmsd).
    `anchors`: atom names to fit on (None = every matched atom). `within`: (centre, radius) keeps only anchors whose reference position lies in that sphere."""
    ri, qi = _matched(ref, qry, anchors)
    if within is not None:
        keep = struc.distance(ref.coord[ri], np.asarray(within[0])) <= within[1]
        ri, qi = ri[keep], qi[keep]
    if len(ri) < 3:
        raise ValueError("only %d anchor atoms to fit on" % len(ri))
    ref_a, qry_a = ref[ri], qry[qi]
    _, tf = struc.superimpose(ref_a, qry_a)
    return tf, int(len(ri)), float(struc.rmsd(ref_a, tf.apply(qry_a)))


def displacement_after_fit(ref, qry, atom_name: str, anchors=None, within=None, res_id=None) -> dict:
    """Fit `qry` onto `ref` (see fit_transform), then return how far `atom_name` of the same residue is from its reference position.
    `res_id`: which residue's `atom_name` to measure (None = it must be unique)."""
    tf, n, rmsd = fit_transform(ref, qry, anchors, within)
    pick = (lambda a: a[(a.atom_name == atom_name) & ((a.res_id == res_id) if res_id is not None else True)])
    ref_atom, qry_atom = pick(ref), pick(qry)
    if len(ref_atom) != 1 or len(qry_atom) != 1:
        raise ValueError("atom %s must occur exactly once in each structure's residue (found %d and %d)" % (atom_name, len(ref_atom), len(qry_atom)))
    moved = tf.apply(qry_atom.coord)[0]
    return {"displacement": float(np.linalg.norm(moved - ref_atom.coord[0])), "n_anchors": n, "fit_rmsd": rmsd}


# ---------------------------------------------------------------------------------------------- real structures
def _load(cache_dir, pid):
    f = pdbx.CIFFile.read(str(Path(cache_dir) / (pid + ".cif")))
    return pdbx.get_structure(f, model=1, altloc="occupancy", extra_fields=["b_factor", "occupancy"], include_bonds=False)


def compare_structures(cache_dir, res_id: int, res_name: str, reference, members, ligand_name: str, local_radius: float = 12.0) -> list:
    """One row per (structure, chain): torsions, B-factor, distance to the deposited ligand's carboxylate, and the NZ-like end atom's displacement from the
    reference after fitting (a) every protein CA of the chain and (b) only the CAs within `local_radius` of the reference ligand's centroid."""
    cache, end_atom = {}, {"LYS": "NZ"}[res_name]

    def parts(pid, chain):
        if pid not in cache:
            cache[pid] = _load(cache_dir, pid)
        a = cache[pid]
        chain_atoms = a[a.chain_id == chain]
        protein = chain_atoms[struc.filter_canonical_amino_acids(chain_atoms)]
        lig = chain_atoms[chain_atoms.res_name == ligand_name]
        if len(lig) == 0:
            raise ValueError("no %s in %s chain %s" % (ligand_name, pid, chain))
        return protein, lig, chain_atoms

    ref_prot, ref_lig, _ = parts(*reference)
    centre = ref_lig.coord.mean(axis=0)
    rows = []
    for pid, chain in members:
        prot, lig, chain_atoms = parts(pid, chain)
        res = prot[prot.res_id == res_id]
        if len(res) == 0 or str(res.res_name[0]) != res_name:
            raise ValueError("%s chain %s has no %s %d" % (pid, chain, res_name, res_id))
        end = res[res.atom_name == end_atom]
        carb = lig[np.isin(lig.atom_name, CARBOXYLATE)]
        near = chain_atoms[(chain_atoms.res_name != ligand_name) & np.isin(chain_atoms.element, ["N", "O"])]
        d_near = struc.distance(near.coord, end.coord[0])
        partners = sorted("%s%d %s" % (str(near.res_name[k]).title(), int(near.res_id[k]), str(near.atom_name[k])) for k in np.where(d_near <= POLAR_PARTNER_CUTOFF)[0]
                          if str(near.res_name[k]) != res_name or int(near.res_id[k]) != res_id)
        glob = displacement_after_fit(ref_prot, prot, end_atom, anchors=("CA",), res_id=res_id)
        loc = displacement_after_fit(ref_prot, prot, end_atom, anchors=("CA",), within=(centre, local_radius), res_id=res_id)
        rows.append({"label": "%s:%s" % (pid, chain), "chi": chi_angles(res), "nz_b": float(end.b_factor[0]), "nz_occupancy": float(end.occupancy[0]),
                     "nz_to_carboxylate": float(min(struc.distance(end.coord[0], c) for c in carb.coord)), "partners": partners,
                     "displacement": glob["displacement"], "fit_rmsd": glob["fit_rmsd"], "n_anchors": glob["n_anchors"],
                     "displacement_local": loc["displacement"], "fit_rmsd_local": loc["fit_rmsd"], "n_anchors_local": loc["n_anchors"]})
    return rows


def render_markdown(rows: list, res_label: str) -> str:
    L = ["# Step 4: where does %s sit in each structure?" % res_label, "",
         "Generated by `python -m sidacore.sidechain`; do not edit. Deposited coordinates (highest-occupancy alternate), chain-by-chain; residues matched by number (the "
         "construct check found identical numbering). Displacement = distance of the end atom from its 4B63 position after fitting CA atoms of the whole chain "
         "(`all CA`) or only those within 12 A of the deposited ligand (`site CA`); the fit RMSD in brackets is the fit's own quality.", "",
         "| Structure:chain | chi1 | chi2 | chi3 | chi4 | end atom B (A^2) | to carboxylate (A) | displacement, all CA (A) | displacement, site CA (A) | other polar neighbours within 3.5 A |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append("| %s | %s | %.1f | %.2f | %.2f [%.2f] | %.2f [%.2f] | %s |" % (
            r["label"], " | ".join("%.0f" % c for c in r["chi"]), r["nz_b"], r["nz_to_carboxylate"], r["displacement"], r["fit_rmsd"],
            r["displacement_local"], r["fit_rmsd_local"], ", ".join(r["partners"]) or "none"))
    return "\n".join(L) + "\n"


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache", default="structures_cache")
    ap.add_argument("--out", default="reports/step4_lys107.md")
    a = ap.parse_args(argv)
    members = [("4B63", "A")] + [("6X0J", c) for c in "ABCD"]
    rows = compare_structures(a.cache, 107, "LYS", ("4B63", "A"), members, "ORN")
    Path(a.out).write_text(render_markdown(rows, "Lys107"))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
