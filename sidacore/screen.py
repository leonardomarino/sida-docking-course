"""Generic screening: any PDB target x any ligand list -> Vina poses, a redock control per target, and a direct clash count.

What this answers: "can this molecule be placed in this site without overlapping the protein, and does the docking
protocol recover the deposited ligand in this target at all?". What it does NOT answer: binding affinity. Vina scores are
reported because they come for free; on SidA they mis-rank (4B63: the deposited pose scores 0.4-0.5 kcal/mol below a decoy).

Pieces, all thin glue over the validated Step 2/4 code (prep.py, redock.py), which is not modified:

    target    biotite (select) -> PDB2PQR/PROPKA (hydrogens at a pH) -> Meeko (PDBQT, box). FAD/FDA/NADP kept, waters and
              crystallization additives stripped, every removal logged. The box is centred on the deposited site ligand.
    ligand    SMILES -> Open Babel protonation at a pH -> RDKit 3D -> Meeko PDBQT.
    control   the deposited ligand is redocked into each target (same box, same settings): a target whose median top-1 RMSD
              exceeds `accept` is labelled FAILED and every result against it must be read with that in mind.
    clashes   heavy-atom pairs closer than the sum of Bondi vdW radii minus `tolerance`; N/O..N/O pairs (hydrogen bonds) are
              judged by a minimum distance instead; metals and hydrogens never count.
"""
from __future__ import annotations

import csv
import json
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from rdkit import Chem, RDLogger
from scipy.spatial import cKDTree

import biotite.structure as struc
import biotite.structure.io.pdb as bpdb

from .boxes import _tool
from .prep import TEMPLATES, _meeko_bin, _resolve_pdb2pqr, _run, fetch_structure, read_structure
from .redock import (crystal_reference, evaluate_seeds, ligand_from_smiles, ligand_template, protonate,
                     read_poses, redock_seeds, run_vina)

# Bondi (1964) van der Waals radii, angstrom. Unknown elements get 1.8.
BONDI = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98}
DEFAULT_RADIUS = 1.8
POLAR = frozenset({"N", "O"})
METALS = frozenset({"LI", "NA", "MG", "K", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG"})
MAX_PAIR = 2 * 1.98 + 0.1          # larger than any Bondi sum above: the neighbour search radius

# Calibrated on the five deposited ornithine poses in the SidA receptors (4B63, 6X0J A-D): their worst non-polar overlap
# is 0.16 A and their closest N/O..N/O hydrogen bond is 2.50 A. 0.4 clears the first by 0.24 A, 2.4 sits under the second.
# The calibration set is small (one 9-atom ligand); change these for other targets only on evidence from crystal poses.
DEFAULT_TOLERANCE = 0.4
DEFAULT_POLAR_MIN = 2.4

SUPPORTED_COFACTORS = ("FAD", "FDA", "NAP")      # residues Meeko has a shipped template for (data/meeko_cofactor_templates.json)
KNOWN_NO_TEMPLATE = {"NDP": "reduced NADPH", "NAD": "NAD+", "NAI": "reduced NAD", "FMN": "FMN", "HEM": "heme"}
SITE_CONTACT = 4.5                  # a receptor residue "contacts" a ligand within this many angstrom (as in prep.py)
MAX_TORSDOF = 10                    # beyond this a ligand is flagged as highly flexible
POCKET_DELETE_RADIUS = 8.0          # Meeko deletes residues it cannot parameterise within this radius of the box centre (prep.py value)


# --------------------------------------------------------------------------- atoms
@dataclass(frozen=True)
class Atoms:
    coord: np.ndarray
    element: np.ndarray       # upper-case
    res_name: np.ndarray
    res_id: np.ndarray
    chain: np.ndarray
    atom_name: np.ndarray

    def __len__(self):
        return len(self.coord)


def read_pdb_heavy(path) -> Atoms:
    """Heavy atoms of a PDB file by fixed columns. Works for Meeko's receptor_prepared.pdb (blank occupancy, all ATOM
    records) and for ordinary PDB files; hydrogens/deuterium are dropped."""
    rows = [l for l in Path(path).read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    rows = [l for l in rows if l[76:78].strip().upper() not in ("H", "D", "")]
    return Atoms(coord=np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in rows], float).reshape(-1, 3),
                 element=np.array([l[76:78].strip().upper() for l in rows]),
                 res_name=np.array([l[17:20].strip() for l in rows]), res_id=np.array([int(l[22:26]) for l in rows]),
                 chain=np.array([l[21] for l in rows]), atom_name=np.array([l[12:16].strip() for l in rows]))


def pose_atoms(mol) -> Atoms:
    """Heavy atoms of an RDKit pose (first conformer)."""
    xyz = mol.GetConformer().GetPositions()
    keep = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    sym = [mol.GetAtomWithIdx(i).GetSymbol().upper() for i in keep]
    return Atoms(coord=xyz[keep], element=np.array(sym), res_name=np.array(["LIG"] * len(keep)),
                 res_id=np.ones(len(keep), int), chain=np.array(["L"] * len(keep)),
                 atom_name=np.array(["%s%d" % (s.capitalize(), i + 1) for s, i in zip(sym, keep)]))


def pose_pdb_text(mol) -> str:
    """PDB text of a pose from its coordinates and bond graph (heavy atoms). RDKit's own PDB writer cannot kekulise the
    unsanitised aromatic poses Meeko returns, so the bonds are written as CONECT records instead."""
    xyz = mol.GetConformer().GetPositions()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    serial = {i: k + 1 for k, i in enumerate(heavy)}
    lines = []
    for i in heavy:
        sym = mol.GetAtomWithIdx(i).GetSymbol().upper()
        lines.append("HETATM%5d %-4s LIG L   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                     % (serial[i], ("%s%d" % (sym, i + 1))[:4], xyz[i][0], xyz[i][1], xyz[i][2], sym))
    for i in heavy:
        nbrs = [serial[b.GetOtherAtomIdx(i)] for b in mol.GetAtomWithIdx(i).GetBonds() if b.GetOtherAtomIdx(i) in serial]
        if nbrs:
            lines.append("CONECT%5d" % serial[i] + "".join("%5d" % n for n in nbrs))
    return "\n".join(lines + ["END", ""])


def _biotite_atoms(a) -> Atoms:
    a = a[a.element != "H"]
    return Atoms(a.coord.astype(float), np.char.upper(a.element.astype(str)), a.res_name.astype(str), a.res_id.astype(int),
                 a.chain_id.astype(str), a.atom_name.astype(str))


# --------------------------------------------------------------------------- clashes
@dataclass(frozen=True)
class Clash:
    lig_atom: str
    lig_element: str
    res_name: str
    res_id: int
    chain: str
    rec_atom: str
    distance: float
    overlap: float        # r_lig + r_rec - distance (Bondi); for N/O..N/O pairs the same number, though they are judged by distance

    def label(self) -> str:
        return "%s%d:%s - %s, %.2f A (overlap %.2f)" % (self.res_name, self.res_id, self.rec_atom, self.lig_atom, self.distance, self.overlap)


def find_clashes(lig: Atoms, rec: Atoms, tolerance: float = DEFAULT_TOLERANCE, polar_min: float = DEFAULT_POLAR_MIN) -> list:
    """Ligand-receptor heavy-atom clashes, worst overlap first. Receptor hydrogens are not considered (pass heavy atoms)."""
    if not len(lig) or not len(rec):
        return []
    tree = cKDTree(rec.coord)
    out = []
    for i, nbrs in enumerate(tree.query_ball_point(lig.coord, MAX_PAIR)):
        ei = lig.element[i]
        if ei in METALS:
            continue
        for j in nbrs:
            ej = rec.element[j]
            if ej in METALS or ej == "H" or ei == "H":
                continue
            d = float(np.linalg.norm(lig.coord[i] - rec.coord[j]))
            overlap = BONDI.get(ei, DEFAULT_RADIUS) + BONDI.get(ej, DEFAULT_RADIUS) - d
            bad = d < polar_min if (ei in POLAR and ej in POLAR) else overlap > tolerance
            if bad:
                out.append(Clash(str(lig.atom_name[i]), ei, str(rec.res_name[j]), int(rec.res_id[j]), str(rec.chain[j]),
                                 str(rec.atom_name[j]), d, overlap))
    return sorted(out, key=lambda c: -c.overlap)


def _contact_residues(centres: np.ndarray, rec: Atoms, cutoff: float) -> set:
    tree = cKDTree(rec.coord)
    hit = {j for nb in tree.query_ball_point(centres, cutoff) for j in nb}
    return {(str(rec.chain[j]), int(rec.res_id[j])) for j in hit}


def site_metrics(pose: Atoms, crystal: Atoms, rec: Atoms, contact: float = SITE_CONTACT) -> dict:
    """Is the pose in the deposited ligand's site? Box volume includes surface and solvent, so a clash-free pose is not
    automatically a pose in the pocket. centroid_dist: pose centroid to deposited-ligand centroid (A).
    site_contact_fraction: share of the receptor residues touching the deposited ligand that also touch the pose."""
    site = _contact_residues(crystal.coord, rec, contact)
    mine = _contact_residues(pose.coord, rec, contact)
    return {"centroid_dist": float(np.linalg.norm(pose.coord.mean(axis=0) - crystal.coord.mean(axis=0))),
            "site_contact_fraction": len(site & mine) / len(site) if site else float("nan")}


# --------------------------------------------------------------------------- targets
@dataclass(frozen=True)
class TargetSpec:
    pdb_id: str
    chain: str
    site_ligand: str | None = None          # residue name of the deposited ligand that defines the site (and the control)
    center: tuple | None = None             # box centre (x, y, z) when there is no site ligand, or to override its centroid
    keep: object = "auto"                   # "auto" = FAD/FDA/NAP present in the chain; or an explicit list of residue names
    box_size: float = 12.0                  # cubic box edge, angstrom. 12 = site-sized: a ligand that does not fit is forced to clash.
                                            # 22 (the Step 3 box) leaves room at the rim, where an oversized ligand docks clash-free.
    ph: float = 7.4                         # protein protonation (PROPKA); ligand pH is a separate Settings value
    label: str = ""

    @property
    def name(self) -> str:
        return "%s_%s%s" % (self.pdb_id.upper(), self.chain, "_" + self.label if self.label else "")


@dataclass(frozen=True)
class Selection:
    kept: object                # biotite AtomArray: protein heavy atoms + kept cofactors
    ligand: object              # biotite AtomArray of the site ligand (one copy), or None
    center: np.ndarray
    kept_cofactors: list
    removed: list               # one row per removed residue type: res_name, n_residues, n_atoms, reason, min_dist_to_site
    notes: list


def check_keep(keep) -> None:
    bad = [k for k in keep if k not in SUPPORTED_COFACTORS]
    if bad:
        raise ValueError("cannot keep %s: Meeko has no template for it (shipped: %s); add a template or strip it" % (bad, list(SUPPORTED_COFACTORS)))


def _hetero_kind(res_name: str) -> str:
    return "water" if res_name in ("HOH", "DOD", "WAT") else "heterogen"


def list_heterogens(atoms, chain: str) -> list:
    """Non-polymer residues of one chain (waters excluded; see count_waters): res_name, n_residues, n_atoms."""
    m = (atoms.chain_id == chain) & atoms.hetero & ~struc.filter_solvent(atoms)
    out = {}
    for rn, rid in sorted({(str(n), int(r)) for n, r in zip(atoms.res_name[m], atoms.res_id[m])}):
        out.setdefault(rn, {"res_name": rn, "n_residues": 0, "n_atoms": 0})["n_residues"] += 1
    for rn in out:
        out[rn]["n_atoms"] = int((m & (atoms.res_name == rn)).sum())
    return list(out.values())


def count_waters(atoms, chain: str) -> int:
    return len({int(r) for r in atoms.res_id[(atoms.chain_id == chain) & struc.filter_solvent(atoms)]})


def select_target_atoms(atoms, spec: TargetSpec) -> Selection:
    chain = spec.chain
    in_chain = atoms.chain_id == chain
    if not in_chain.any():
        raise ValueError("%s has no chain %r (chains: %s)" % (spec.pdb_id, chain, sorted({str(c) for c in atoms.chain_id})))
    notes = []
    ligand = None
    if spec.site_ligand:
        ids = sorted({int(r) for r in atoms.res_id[in_chain & (atoms.res_name == spec.site_ligand)]})
        if ids:
            ligand = atoms[in_chain & (atoms.res_name == spec.site_ligand) & (atoms.res_id == ids[0]) & (atoms.element != "H")]
            if len(ids) > 1:
                notes.append("%d copies of %s in chain %s (residue numbers %s); the first (%d) defines the site" % (len(ids), spec.site_ligand, chain, ids, ids[0]))
    if spec.center is not None:
        center = np.asarray(spec.center, float)
    elif ligand is not None:
        center = ligand.coord.mean(axis=0)
    else:
        raise ValueError("%s chain %s: no deposited %r to centre the box on; give site_ligand present in the structure or an explicit center=(x, y, z)"
                         % (spec.pdb_id, chain, spec.site_ligand))
    if ligand is not None and spec.center is None:
        outside = int((np.abs(ligand.coord - center) > spec.box_size / 2).any(axis=1).sum())
        if outside:
            notes.append("the site ligand %s extends outside the box (%d of %d heavy atoms beyond %.0f A): enlarge box_size, or the control and every fit result are unreliable"
                         % (spec.site_ligand, outside, len(ligand), spec.box_size))
    present = {str(n) for n in atoms.res_name[in_chain & atoms.hetero]}
    keep = [c for c in SUPPORTED_COFACTORS if c in present] if spec.keep == "auto" else list(spec.keep)
    check_keep(keep)

    is_prot = in_chain & struc.filter_canonical_amino_acids(atoms) & (atoms.element != "H")
    is_keep = in_chain & np.isin(atoms.res_name, keep) & (atoms.element != "H")
    kept = atoms[is_prot | is_keep]

    dist = np.linalg.norm(atoms.coord - center, axis=1)
    removed = []
    rest = in_chain & ~is_prot & ~is_keep & ~((atoms.res_name == spec.site_ligand) if spec.site_ligand else np.zeros(len(atoms), bool)) & (atoms.element != "H")
    for rn in sorted({str(n) for n in atoms.res_name[rest]}):
        m = rest & (atoms.res_name == rn)
        if _hetero_kind(rn) == "water":
            reason = "water stripped"
        elif rn in KNOWN_NO_TEMPLATE:
            reason = "cofactor (%s) stripped: no Meeko template" % KNOWN_NO_TEMPLATE[rn]
            notes.append("%s present in chain %s but stripped: Meeko has no template for it" % (rn, chain))
        elif rn in keep:
            continue
        elif not bool(atoms.hetero[m].any()):
            reason = "non-canonical residue stripped"
        else:
            reason = "crystallization additive / other heterogen stripped"
        removed.append({"res_name": rn, "n_residues": len({int(r) for r in atoms.res_id[m]}), "n_atoms": int(m.sum()),
                        "reason": reason, "min_dist_to_site": round(float(dist[m].min()), 1)})
    if spec.site_ligand and ligand is not None:
        copies = len({int(r) for r in atoms.res_id[in_chain & (atoms.res_name == spec.site_ligand)]})
        removed.append({"res_name": spec.site_ligand, "n_residues": copies, "n_atoms": int(len(ligand)),
                        "reason": "site ligand set aside (redock control reference)" + ("; the other copies are dropped" if copies > 1 else ""),
                        "min_dist_to_site": 0.0})
    return Selection(kept=kept, ligand=ligand, center=center, kept_cofactors=keep, removed=removed, notes=notes)


def preview_target(spec: TargetSpec, cache_dir) -> Selection:
    """What prepare_target would keep and strip, without running PDB2PQR or Meeko (downloads the structure if needed)."""
    atoms, _ = read_structure(fetch_structure(spec.pdb_id, cache_dir), "highest-occupancy")
    return select_target_atoms(atoms, spec)


def check_tools() -> dict:
    """{tool: path or None} for every external program the screen calls."""
    out = {}
    for name in ("vina", "obabel", "mk_prepare_receptor.py"):
        try:
            out[name] = str(_tool(name))
        except FileNotFoundError:
            out[name] = None
    exe = Path(sys.executable).parent / "pdb2pqr30"
    try:
        out["pdb2pqr30"] = str(_resolve_pdb2pqr(exe if exe.exists() else None))
    except FileNotFoundError:
        out["pdb2pqr30"] = None
    return out


@dataclass(frozen=True)
class PreparedTarget:
    spec: TargetSpec
    dir: Path
    receptor_pdbqt: Path
    receptor_pdb: Path
    box_file: Path
    center: np.ndarray
    box_size: tuple
    crystal_ligand_pdb: Path | None
    kept_cofactors: list
    removed: list
    meeko_deleted: list
    notes: list

    @property
    def name(self) -> str:
        return self.spec.name


def _spec_json(spec: TargetSpec) -> dict:
    d = asdict(spec)
    d["center"] = list(d["center"]) if d["center"] is not None else None
    d["keep"] = list(d["keep"]) if not isinstance(d["keep"], str) else d["keep"]
    return d


def _load_prepared(spec: TargetSpec, out: Path) -> PreparedTarget:
    r = json.loads((out / "target.json").read_text())
    return PreparedTarget(spec=spec, dir=out, receptor_pdbqt=out / "receptor.pdbqt", receptor_pdb=out / "receptor_prepared.pdb",
                          box_file=out / "receptor.box.txt", center=np.array(r["center"]), box_size=tuple(r["box_size"]),
                          crystal_ligand_pdb=(out / "ligand_crystal.pdb") if r["has_site_ligand"] else None,
                          kept_cofactors=r["kept_cofactors"], removed=r["removed"], meeko_deleted=r["meeko_deleted"], notes=r["notes"])


def prepare_target(spec: TargetSpec, cache_dir, work_dir) -> PreparedTarget:
    """Download (cached), select, protonate and convert one target. Reuses an earlier preparation of the identical spec;
    refuses to reuse a directory prepared with different settings."""
    out = Path(work_dir) / "targets" / spec.name
    if (out / "target.json").exists():
        prev = json.loads((out / "target.json").read_text())["spec"]
        if prev != _spec_json(spec):
            raise ValueError("%s was already prepared with different settings (%s); change `label` or remove %s" % (spec.name, prev, out))
        return _load_prepared(spec, out)

    atoms, _ = read_structure(fetch_structure(spec.pdb_id, cache_dir), "highest-occupancy")
    sel = select_target_atoms(atoms, spec)
    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)

    if sel.ligand is not None:
        f = bpdb.PDBFile(); bpdb.set_structure(f, sel.ligand); f.write(str(out / "ligand_crystal.pdb"))
    prot = sel.kept[struc.filter_canonical_amino_acids(sel.kept)]
    f = bpdb.PDBFile(); bpdb.set_structure(f, prot); f.write(str(work / "protein_in.pdb"))
    exe = Path(sys.executable).parent / "pdb2pqr30"
    p2p = _resolve_pdb2pqr(exe if exe.exists() else None)
    _run([p2p, "--ff=AMBER", "--with-ph=%g" % spec.ph, "--titration-state-method=propka", "--keep-chain", "--drop-water",
          "--pdb-output", work / "protein_h.pdb", work / "protein_in.pdb", work / "protein.pqr"], work / "pdb2pqr.log")
    prot_h = bpdb.PDBFile.read(str(work / "protein_h.pdb")).get_structure(model=1)

    extra = sel.kept[np.isin(sel.kept.res_name, sel.kept_cofactors)].copy()
    for cat in ("occupancy", "altloc_id"):
        if cat in extra.get_annotation_categories():
            extra.del_annotation(cat)
    extra.bonds = None
    f = bpdb.PDBFile(); bpdb.set_structure(f, prot_h + extra); f.write(str(work / "combined.pdb"))
    box = [float(spec.box_size)] * 3
    _run([_meeko_bin(), "--read_pdb", work / "combined.pdb", "-o", out / "receptor", "-p", "-j", "-v",
          "--write_pdb", out / "receptor_prepared.pdb", "--box_center", *["%.3f" % x for x in sel.center],
          "--box_size", *["%g" % x for x in box], "--delete_bad_res_from_box_radius", "%g" % POCKET_DELETE_RADIUS,
          "--add_templates", TEMPLATES], work / "meeko.log")

    import biotite.structure.io.pdbqt as bpdbqt
    combo = prot_h + extra
    got = bpdbqt.PDBQTFile.read(str(out / "receptor.pdbqt")).get_structure(model=1)
    res = lambda a: {(str(c), int(r), str(n)) for c, r, n in zip(a.chain_id, a.res_id, a.res_name)}
    deleted = [{"chain": c, "res_id": r, "res_name": n,
                "dist_to_site": round(float(np.linalg.norm(combo.coord[(combo.chain_id == c) & (combo.res_id == r)] - sel.center, axis=1).min()), 1)}
               for c, r, n in sorted(res(combo) - res(got), key=lambda k: k[1])]
    notes = list(sel.notes)
    if deleted:
        notes.append("Meeko deleted %d residue(s) it could not parameterise: %s" % (len(deleted), ["%s%d" % (d["res_name"], d["res_id"]) for d in deleted]))
    record = {"spec": _spec_json(spec), "center": [float(x) for x in sel.center], "box_size": box, "has_site_ligand": sel.ligand is not None,
              "kept_cofactors": sel.kept_cofactors, "removed": sel.removed, "meeko_deleted": deleted, "notes": notes}
    (out / "target.json").write_text(json.dumps(record, indent=1))
    return _load_prepared(spec, out)


# --------------------------------------------------------------------------- ligands
@dataclass(frozen=True)
class LigandEntry:
    name: str
    input_smiles: str
    smiles: str | None = None
    pdbqt: str | None = None
    net_charge: int | None = None
    n_heavy: int | None = None
    torsdof: int | None = None
    warnings: tuple = ()
    error: str | None = None


def load_ligands(source) -> dict:
    """{name: SMILES} from a dict, a CSV with columns name,smiles, or a .smi/.txt file of `SMILES name` lines."""
    if isinstance(source, dict):
        return dict(source)
    path = Path(source)
    pairs = []
    if path.suffix.lower() == ".csv":
        with open(path, newline="") as fh:
            pairs = [(r["name"].strip(), r["smiles"].strip()) for r in csv.DictReader(fh)]
    elif path.suffix.lower() in (".smi", ".txt"):
        for line in path.read_text().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                pairs.append((parts[1], parts[0]))
    else:
        raise ValueError("ligand file %s: use .csv (name,smiles) or .smi/.txt (SMILES name)" % path)
    names = [n for n, _ in pairs]
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        raise ValueError("duplicate ligand name(s) in %s: %s" % (path, dup))
    return dict(pairs)


def prepare_ligands(ligands: dict, ph: float | None) -> list:
    """Protonate at `ph` (None = take the SMILES as written), embed, convert. A ligand that cannot be prepared is returned
    with `error` set; it is never dropped silently and never stops the others."""
    out = []
    for name, smi in ligands.items():
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                raise ValueError("RDKit cannot read this SMILES")
            warn = []
            centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
            if any(c[1] == "?" for c in centres):
                warn.append("%d undefined stereocentre(s): the 3D start geometry picks one arbitrarily" % sum(c[1] == "?" for c in centres))
            use = protonate(mol, ph) if ph is not None else smi
            lp = ligand_from_smiles(use, ph)
            if lp.torsdof > MAX_TORSDOF:
                warn.append("%d rotatable bonds: sampling is harder, treat the poses with caution" % lp.torsdof)
            out.append(LigandEntry(name, smi, lp.smiles, lp.pdbqt, lp.net_charge, lp.n_heavy, lp.torsdof, tuple(warn)))
        except Exception as e:                                   # noqa: BLE001 - recorded, not hidden
            out.append(LigandEntry(name, smi, error="%s: %s" % (type(e).__name__, e)))
    return out


# --------------------------------------------------------------------------- docking, control, analysis
@dataclass(frozen=True)
class Settings:
    ligand_ph: float | None = 7.4           # None = SMILES taken as written
    seeds: tuple = (1, 2, 3)
    exhaustiveness: int = 16
    num_modes: int = 20
    control_seeds: tuple = (1, 2, 3, 4, 5)
    accept: float = 2.0                     # control passes when the median top-1 RMSD over seeds is within this (Step 4 rule)
    tolerance: float = DEFAULT_TOLERANCE
    polar_min: float = DEFAULT_POLAR_MIN
    site_min: float = 0.6                   # a pose is 'in the site' when it touches at least this share of the residues touching the deposited ligand
                                            # (ornithine poses 0.83-1.0; an oversized ligand at the rim 0.33-0.55 in the SidA probe)


def dock(target: PreparedTarget, entry: LigandEntry, settings: Settings, out_dir) -> list:
    """[(seed, rank, score, RDKit pose)] over all seeds, one list."""
    rows = []
    for seed in settings.seeds:
        f = run_vina(target.receptor_pdbqt, entry.pdbqt, target.box_file, seed=seed, exhaustiveness=settings.exhaustiveness,
                     num_modes=settings.num_modes, out_dir=out_dir)
        rows += [(seed, k + 1, s, m) for k, (s, m) in enumerate(read_poses(f.read_text()))]
    return rows


def analyse(poses: list, rec: Atoms, crystal: Atoms | None, settings: Settings):
    """One result row for a ligand in a target, plus the clash list of the reported pose. The reported pose is the
    lowest-scoring pose over all seeds; the clash-free column asks whether ANY docked pose avoids the protein."""
    per = []
    for seed, rank, score, mol in poses:
        a = pose_atoms(mol)
        per.append((score, seed, rank, mol, a, find_clashes(a, rec, settings.tolerance, settings.polar_min)))
    per.sort(key=lambda t: t[0])
    score, seed, rank, mol, a, cl = per[0]
    free = [p for p in per if not p[5]]
    row = {"best_score": score, "best_seed": seed, "best_rank": rank, "n_clashes": len(cl),
           "worst_overlap": cl[0].overlap if cl else 0.0, "worst_clash": cl[0].label() if cl else "",
           "best_clash_free_score": free[0][0] if free else float("nan"),
           "n_clash_free": len(free), "n_poses": len(per)}
    if crystal is not None:
        m = site_metrics(a, crystal, rec)
        row.update(centroid_dist=m["centroid_dist"], site_contact_fraction=m["site_contact_fraction"])
        if free:
            row["clash_free_site_contact_fraction"] = site_metrics(free[0][4], crystal, rec)["site_contact_fraction"]
    else:
        row.update(centroid_dist=float("nan"), site_contact_fraction=float("nan"))
    return row, cl, mol


def control(target: PreparedTarget, settings: Settings, out_dir) -> dict:
    """Redock the deposited site ligand into its own target (same box, same engine). One row for the controls table."""
    base = {"target": target.name, "site_ligand": target.spec.site_ligand}
    if target.crystal_ligand_pdb is None:
        return {**base, "status": "no control (no site ligand)", "passes": None}
    try:
        template = ligand_template(target.spec.site_ligand)
        lp = ligand_from_smiles(protonate(template, settings.ligand_ph), settings.ligand_ph) if settings.ligand_ph is not None \
            else ligand_from_smiles(Chem.MolToSmiles(template))
        RDLogger.DisableLog("rdApp.warning")         # "More than one matching pattern found" = the symmetric carboxylate oxygens; harmless
        try:
            ref = crystal_reference(target.crystal_ligand_pdb, template)
        finally:
            RDLogger.EnableLog("rdApp.warning")
        cfg = SimpleNamespace(seeds=settings.control_seeds, exhaustiveness=settings.exhaustiveness, num_modes=settings.num_modes)
        per = redock_seeds(target.dir, lp, ref, cfg, "vina", Path(out_dir), receptor=target.receptor_pdbqt, box_file=target.box_file)
        ev = evaluate_seeds(per, settings.accept)
    except Exception as e:                                       # noqa: BLE001 - recorded, not hidden
        return {**base, "status": "control failed to run: %s: %s" % (type(e).__name__, e), "passes": None}
    rec = read_pdb_heavy(target.receptor_pdb)
    crystal = read_pdb_heavy(target.crystal_ligand_pdb)
    return {**base, "status": "validated" if ev["passes"] else "FAILED", "passes": ev["passes"],
            "median_top1_rmsd": ev["top1_median"], "top1_min": ev["top1_min"], "top1_max": ev["top1_max"],
            "seeds_within_accept": "%d/%d" % (ev["top1_n_within"], ev["n_seeds"]), "best_top5_median": ev["best_top5_median"],
            "crystal_ligand_clashes": len(find_clashes(crystal, rec, settings.tolerance, settings.polar_min))}


def reading(row: dict, control_status: str, settings: Settings) -> str:
    """The table in words, by fixed rules (in this order): no clash-free pose -> 'clashes in every pose'; best pose clashes
    -> 'best pose clashes'; clash-free but touching under `site_min` of the site residues -> 'outside the site'; else 'fits'.
    A target that failed its control gets a prefix: the reading is then about the model, not the molecule."""
    frac = row["site_contact_fraction"]
    if row["n_clash_free"] == 0:
        text = "clashes in every pose (all %d)" % row["n_poses"]
    elif row["n_clashes"] > 0:
        text = "best pose clashes (%d clash-free poses of %d)" % (row["n_clash_free"], row["n_poses"])
    elif frac != frac:
        text = "clash-free (site occupancy not assessed: no site ligand)"
    elif frac < settings.site_min:
        text = "clash-free but outside the site (touches %.0f%% of the site residues)" % (100 * frac)
    else:
        text = "fits: clash-free, in the site"
    return ("(target failed its control) " if control_status == "FAILED" else "") + text


# --------------------------------------------------------------------------- the screen
@dataclass
class ScreenResult:
    results: object                 # pandas DataFrame, one row per target x ligand
    controls: object                # pandas DataFrame, one row per target
    targets: dict                   # name -> PreparedTarget
    ligands: list                   # [LigandEntry]
    clashes: dict                   # (target, ligand) -> [Clash] of the reported pose
    poses: dict                     # (target, ligand) -> Path of the reported pose (PDB)
    work_dir: Path


def _slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_") or "x"


def tool_versions() -> dict:
    """Versions of the tools a run depends on (recorded in run.json and shown in the notebook's run record)."""
    import subprocess
    import biotite, meeko, rdkit
    v = subprocess.run([str(_tool("vina")), "--version"], capture_output=True, text=True)
    return {"vina": ((v.stdout + v.stderr).strip().splitlines() or ["unknown"])[0], "rdkit": rdkit.__version__,
            "meeko": meeko.__version__, "biotite": biotite.__version__, "python": sys.version.split()[0]}


def run_screen(targets: list, ligands, *, settings: Settings = Settings(), cache_dir, work_dir) -> ScreenResult:
    import pandas as pd
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    entries = prepare_ligands(load_ligands(ligands), settings.ligand_ph)
    prepared, results, controls, clashes, poses = {}, [], [], {}, {}
    for spec in targets:
        t = prepare_target(spec, cache_dir, work_dir)
        prepared[t.name] = t
        rec = read_pdb_heavy(t.receptor_pdb)
        crystal = read_pdb_heavy(t.crystal_ligand_pdb) if t.crystal_ligand_pdb else None
        ctl = control(t, settings, work_dir / "control" / t.name)
        controls.append(ctl)
        for e in entries:
            base = {"target": t.name, "ligand": e.name, "box": t.box_size[0], "target_control": ctl["status"], "net_charge": e.net_charge,
                    "torsdof": e.torsdof, "warnings": "; ".join(e.warnings)}
            if e.error:
                results.append({**base, "status": "failed: " + e.error})
                continue
            try:
                d = work_dir / "dock" / t.name / _slug(e.name)
                row, cl, mol = analyse(dock(t, e, settings, d), rec, crystal, settings)
            except Exception as ex:                              # noqa: BLE001 - recorded, not hidden
                results.append({**base, "status": "failed: %s: %s" % (type(ex).__name__, ex)})
                continue
            pose_path = work_dir / "poses" / ("%s__%s.pdb" % (t.name, _slug(e.name)))
            pose_path.parent.mkdir(parents=True, exist_ok=True)
            pose_path.write_text(pose_pdb_text(mol))
            clashes[(t.name, e.name)], poses[(t.name, e.name)] = cl, pose_path
            results.append({**base, "status": "ok", **row, "reading": reading(row, ctl["status"], settings)})
    res = pd.DataFrame(results)
    ctl = pd.DataFrame(controls)
    res.to_csv(work_dir / "results.csv", index=False)
    ctl.to_csv(work_dir / "controls.csv", index=False)
    (work_dir / "run.json").write_text(json.dumps({
        "time": time.strftime("%Y-%m-%d %H:%M:%S"), "settings": {k: (list(v) if isinstance(v, tuple) else v) for k, v in asdict(settings).items()},
        "targets": [_spec_json(s) for s in targets], "ligands": {e.name: {"input": e.input_smiles, "docked": e.smiles, "error": e.error} for e in entries},
        "versions": tool_versions()}, indent=1))
    return ScreenResult(res, ctl, prepared, entries, clashes, poses, work_dir)


# --------------------------------------------------------------------------- view
def show_pose(target: PreparedTarget, pose_pdb, clash_list: list, shell: float = 6.0, width: int = 640, height: int = 460):
    """py3Dmol view: receptor cartoon, residues within `shell` of the pose, the pose (green), the deposited ligand (grey),
    and the receptor atoms that clash (red spheres)."""
    import py3Dmol
    rec = read_pdb_heavy(target.receptor_pdb)
    pose = read_pdb_heavy(pose_pdb)
    near = sorted({(c, r) for c, r in _contact_residues(pose.coord, rec, shell)})
    v = py3Dmol.view(width=width, height=height)
    v.addModel(Path(target.receptor_pdb).read_text(), "pdb")
    v.setStyle({"model": 0}, {"cartoon": {"color": "#c8d0d8", "opacity": 0.55}})
    for chain, rid in near:
        v.addStyle({"model": 0, "chain": chain, "resi": rid}, {"stick": {"radius": 0.12, "colorscheme": "lightGrayCarbon"}})
    v.addModel(Path(pose_pdb).read_text(), "pdb")
    v.setStyle({"model": 1}, {"stick": {"radius": 0.22, "colorscheme": "greenCarbon"}})
    if target.crystal_ligand_pdb is not None:
        v.addModel(Path(target.crystal_ligand_pdb).read_text(), "pdb")
        v.setStyle({"model": 2}, {"stick": {"radius": 0.1, "color": "#777777"}})
    for c in clash_list:
        v.addStyle({"model": 0, "chain": c.chain, "resi": c.res_id, "atom": c.rec_atom},
                   {"sphere": {"radius": 0.45, "color": "red", "opacity": 0.8}})
    v.zoomTo({"model": 1})
    return v
