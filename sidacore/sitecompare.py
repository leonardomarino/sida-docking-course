"""Which residues of the substrate site differ between 4B63 and the four 6X0J chains?

Registered before the first run: a residue is flagged "4B63 differs" only if it lies OUTSIDE the cloud of the four 6X0J chains, i.e. the smallest side-chain RMSD from 4B63
to any 6X0J chain is larger than the largest RMSD between any two 6X0J chains AND larger than a 1.0 A floor. The site is the one Steps 2 and 3 use (docking box plus
residues within site_contact_angstrom of the deposited ligand, in the 4B63 frame). Descriptive: nothing here is a pass rule, and a flagged residue is a lead, not a cause.

All structures are fitted onto 4B63 over the CA atoms within 12 A of its deposited ornithine; residues are matched by number (the construct check found identical numbering).
"""
from itertools import combinations
from pathlib import Path

import biotite.structure as struc
import numpy as np

from .prep import FLAVIN, NADP, read_structure, resolve_site_center, site_residue_ids
from .boxes import protein_chain
from .sidechain import fit_transform

BACKBONE = ("N", "CA", "C", "O")
FLOOR = 1.0                                                                       # A
RING = ("N1", "C2", "O2", "N3", "C4", "O4", "C4X", "N5", "C5X", "C6", "C7", "C7M", "C8", "C8M", "C9", "C9A", "N10", "C10")
_SWAP = {"ASP": {"OD1": "OD2", "OD2": "OD1"}, "GLU": {"OE1": "OE2", "OE2": "OE1"}, "ARG": {"NH1": "NH2", "NH2": "NH1"},
         "PHE": {"CD1": "CD2", "CD2": "CD1", "CE1": "CE2", "CE2": "CE1"}, "TYR": {"CD1": "CD2", "CD2": "CD1", "CE1": "CE2", "CE2": "CE1"}}


def _rmsd_by_name(a, b, names) -> float:
    da = {str(n): c for n, c in zip(a.atom_name, a.coord)}
    db = {str(n): c for n, c in zip(b.atom_name, b.coord)}
    common = [n for n in names if n in da and n in db]
    if not common:
        return None
    return float(np.sqrt(np.mean([np.sum((da[n] - db[n]) ** 2) for n in common])))


def sidechain_rmsd(a, b):
    """Heavy-atom side-chain RMSD between two same-type residues already in one frame (no refit). Symmetric groups (carboxylate, guanidinium, rings) are compared
    with their labels both ways and the lower value kept. None for a residue with no side-chain atoms (glycine)."""
    name = str(a.res_name[0])
    if name != str(b.res_name[0]):
        raise ValueError("residue types differ: %s and %s" % (name, b.res_name[0]))
    names = [str(n) for n in a.atom_name if str(n) not in BACKBONE and not str(n).startswith("H")]
    if not names:
        return None
    best = _rmsd_by_name(a, b, names)
    swap = _SWAP.get(name)
    if swap and best is not None:
        b2 = b.copy()
        b2.atom_name = np.array([swap.get(str(n), str(n)) for n in b.atom_name])
        alt = _rmsd_by_name(a, b2, names)
        best = min(best, alt) if alt is not None else best
    return best


def is_outlier(vs_ref: list, within: list, floor: float = FLOOR) -> bool:
    """The registered rule: outside the 6X0J cloud and above the floor."""
    if not vs_ref or any(v is None for v in vs_ref):
        return False
    lo = min(vs_ref)
    return bool(lo > floor and lo > (max(within) if within else 0.0))


def _moved(atoms, tf):
    out = atoms.copy()
    out.coord = tf.apply(out.coord)
    return out


def _cofactor(atoms, chain, names):
    return atoms[(atoms.chain_id == chain) & np.isin(atoms.res_name, names) & (atoms.element != "H")]


def compare(cache_dir, reference=("4B63", "A"), members=(("6X0J", "A"), ("6X0J", "B"), ("6X0J", "C"), ("6X0J", "D")), ligand_name="ORN",
            center_spec="FAD:C4X", box: float = 22.0, contact: float = 4.5, local_radius: float = 12.0, floor: float = FLOOR) -> dict:
    cache_dir, loaded = Path(cache_dir), {}

    def load(pid):
        if pid not in loaded:
            loaded[pid] = read_structure(cache_dir / (pid + ".cif"), "highest-occupancy")[0]
        return loaded[pid]

    rpid, rchain = reference
    ref_all = load(rpid)
    ref_prot = protein_chain(ref_all, rchain)
    ref_prot = ref_prot[ref_prot.element != "H"]
    lig = ref_all[(ref_all.chain_id == rchain) & (ref_all.res_name == ligand_name)]
    if len(lig) == 0:
        raise ValueError("no %s in %s chain %s" % (ligand_name, rpid, rchain))
    centre = resolve_site_center(ref_all, rchain, center_spec)
    site = site_residue_ids(ref_prot, centre, box, lig.coord, contact)

    frames = {}                                                                    # label -> (protein, flavin, NADP) in the 4B63 frame
    for pid, chain in members:
        all_atoms = load(pid)
        prot = protein_chain(all_atoms, chain)
        prot = prot[prot.element != "H"]
        tf, n, rmsd = fit_transform(ref_prot, prot, anchors=("CA",), within=(lig.coord.mean(axis=0), local_radius))
        frames[pid + ":" + chain] = {"prot": _moved(prot, tf), "flavin": _moved(_cofactor(all_atoms, chain, FLAVIN), tf),
                                     "nadp": _moved(_cofactor(all_atoms, chain, NADP), tf), "fit_rmsd": rmsd, "n_anchors": n}
    labels = list(frames)

    rows = []
    for (chain, res_id), flags in sorted(site.items(), key=lambda kv: kv[0][1]):
        ref_res = ref_prot[ref_prot.res_id == res_id]
        if len(ref_res) == 0:
            continue

        def sc(label):
            r = frames[label]["prot"]
            r = r[r.res_id == res_id]
            return sidechain_rmsd(ref_res, r) if len(r) and str(r.res_name[0]) == str(ref_res.res_name[0]) else None

        def pair(l1, l2):
            r1, r2 = frames[l1]["prot"], frames[l2]["prot"]
            r1, r2 = r1[r1.res_id == res_id], r2[r2.res_id == res_id]
            return sidechain_rmsd(r1, r2) if len(r1) and len(r2) and str(r1.res_name[0]) == str(r2.res_name[0]) else None

        vs_ref = [sc(l) for l in labels]
        within = [v for v in (pair(a, b) for a, b in combinations(labels, 2)) if v is not None]
        rows.append({"chain": chain, "res_id": int(res_id), "resname": str(ref_res.res_name[0]), "in_box": flags["in_box"], "contacts_ligand": flags["contacts_ligand"],
                     "vs_ref": vs_ref, "within_max": max(within) if within else None, "outlier": is_outlier(vs_ref, within, floor)})

    cofactors = []
    ref_cof = {"flavin": _cofactor(ref_all, rchain, FLAVIN), "NADP": _cofactor(ref_all, rchain, NADP)}
    for name, key, names_of in (("flavin", "flavin", None), ("isoalloxazine ring", "flavin", RING), ("NADP", "nadp", None)):
        base = ref_cof["NADP" if name == "NADP" else "flavin"]
        use = names_of or sorted({str(n) for n in base.atom_name})
        vs = [_rmsd_by_name(base, frames[l][key], use) for l in labels]
        w = [v for v in (_rmsd_by_name(frames[a][key], frames[b][key], use) for a, b in combinations(labels, 2)) if v is not None]
        cofactors.append({"name": name, "vs_ref": vs, "within_max": max(w) if w else None,
                          "outlier": is_outlier(vs, w, floor)})
    return {"reference": "%s:%s" % reference, "members": labels, "floor": floor, "contact": contact, "box": box, "local_radius": local_radius,
            "fits": {l: {"n_anchors": f["n_anchors"], "fit_rmsd": f["fit_rmsd"]} for l, f in frames.items()}, "residues": rows, "cofactors": cofactors}


def _fmt(v):
    return "-" if v is None else "%.2f" % v


def render_markdown(res: dict) -> str:
    chains = [m.split(":")[1] for m in res["members"]]
    flagged = [r for r in res["residues"] if r["outlier"]]
    L = ["# Step 4: which site residues differ between 4B63 and 6X0J?", "",
         "Generated by `python -m sidacore.sitecompare`; do not edit. Reference %s. The 6X0J chains are fitted onto it over the CA atoms within %.0f A of its deposited "
         "ornithine (fit RMSD: %s). Side-chain heavy-atom RMSD per residue, symmetric groups compared both ways." % (
             res["reference"], res["local_radius"], ", ".join("%s %.2f" % (l, f["fit_rmsd"]) for l, f in res["fits"].items())), "",
         "**Rule registered before the run:** a residue is flagged if the smallest 4B63-to-6X0J RMSD exceeds both the largest RMSD between two 6X0J chains and %.1f A. "
         "The site (%d residues) is the docking box plus residues within %.1f A of the deposited ornithine." % (res["floor"], len(res["residues"]), res["contact"]), "",
         "**%d of %d site residues are flagged.**" % (len(flagged), len(res["residues"])), "",
         "| Residue | in box | contacts ornithine | " + " | ".join("4B63 vs 6X0J:%s" % c for c in chains) + " | largest 6X0J-6X0J | flagged |",
         "|---|---|---|" + "---|" * len(chains) + "---|---|"]
    shown = [r for r in res["residues"] if r["outlier"] or r["contacts_ligand"]]
    for r in sorted(shown, key=lambda r: (not r["outlier"], not r["contacts_ligand"], r["res_id"])):
        L.append("| %s%d | %s | %s | %s | %s | %s |" % (r["resname"], r["res_id"], "yes" if r["in_box"] else "no", "yes" if r["contacts_ligand"] else "no",
                                                        " | ".join(_fmt(v) for v in r["vs_ref"]), _fmt(r["within_max"]), "**yes**" if r["outlier"] else "no"))
    L += ["", "Only flagged residues and residues that contact the deposited ornithine are listed; the other %d are in the `--json` output." % (
        len(res["residues"]) - len(shown)), "",
          "| Cofactor (matched heavy atoms) | " + " | ".join("4B63 vs 6X0J:%s" % c for c in chains) + " | largest 6X0J-6X0J | flagged |", "|---|" + "---|" * len(chains) + "---|---|"]
    for c in res["cofactors"]:
        L.append("| %s | %s | %s | %s |" % (c["name"], " | ".join(_fmt(v) for v in c["vs_ref"]), _fmt(c["within_max"]), "**yes**" if c["outlier"] else "no"))
    return "\n".join(L) + "\n"


def main(argv=None):
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache", default="structures_cache")
    ap.add_argument("--out", default="reports/step4_site_comparison.md")
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    res = compare(a.cache)
    Path(a.out).write_text(render_markdown(res))
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
