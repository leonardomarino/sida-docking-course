"""Step 3: the substrate-site box, the same box carried into other structures, and an fpocket cross-check.

    biotite   superimpose_homologs (sequence alignment, CA anchors, outlier trimming), rmsd, distance, sphere files
    fpocket   independent pocket detection on the prepared receptor (heavy atoms, cofactors included)

The box is centred on the isoalloxazine C4a (PDB atom C4X). To compare conformers the FAD-in box is carried into
every other structure by superposition, never by eye. Where the C4X box and the transferred box differ, or where
fpocket disagrees with a box, both numbers are reported; neither is silently preferred.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path

import biotite
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import numpy as np

from .config import write_resolved
from .prep import _run, read_structure, resolve_site_center


# --------------------------------------------------------------------------- geometry
def protein_chain(atoms, chain: str):
    """Canonical amino acids of one chain. NOT filter_amino_acids: that one counts L-ornithine as protein."""
    return atoms[(atoms.chain_id == chain) & struc.filter_canonical_amino_acids(atoms)]


def inside_box(coords, center, size):
    return np.all(np.abs(np.asarray(coords, float) - np.asarray(center, float)) <= np.asarray(size, float) / 2, axis=1)


def box_overlap_fraction(a, b, size) -> float:
    """Fraction of one axis-aligned box (same size) that the other covers."""
    s = np.asarray(size, float)
    d = np.abs(np.asarray(a, float) - np.asarray(b, float))
    return float(np.prod(np.clip(s - d, 0.0, None) / s))


def ligand_vs_box(ligand, center, size) -> dict:
    heavy = ligand.coord[ligand.element != "H"]
    return {"n_heavy": int(len(heavy)), "n_inside": int(inside_box(heavy, center, size).sum()),
            "centroid_distance": float(np.linalg.norm(heavy.mean(axis=0) - np.asarray(center, float)))}


# --------------------------------------------------------------------------- superposition
@dataclass(frozen=True, eq=False)
class SiteTransfer:
    center_in_query: np.ndarray   # the reference site point, expressed in the query's own coordinate frame
    n_aligned: int                # CA pairs from the sequence alignment
    n_anchors: int                # pairs kept after outlier trimming; the fit uses only these
    n_trimmed: int
    rmsd_anchors: float
    n_site: int                   # aligned pairs whose reference CA lies within `radius` of the site
    n_site_trimmed: int           # of those, how many were trimmed as outliers (i.e. moved)
    site_rmsd: float | None       # over ALL aligned pairs near the site, trimmed or not
    site_max_displacement: float | None


def transfer_site(ref, qry, ref_site, radius: float, outlier_threshold: float) -> SiteTransfer:
    """Carry `ref_site` (a point in the reference frame) into the query frame.

    The fit is trimmed of conformational outliers, so a moved loop cannot drag it. That would also hide the very
    change we want to see, so the site numbers are computed over all aligned pairs, after the trimmed fit."""
    ref_site = np.asarray(ref_site, float)
    if ref.array_length() == qry.array_length() and np.array_equal(ref.coord, qry.coord):
        # the reference carried onto itself: the identity. A fit would have a zero interquartile range and trim
        # floating-point dust as "outliers", reporting a meaningless anchor count.
        ca = np.nonzero(ref.atom_name == "CA")[0]
        fix_i = mob_i = fix_all = mob_all = ca
        move = lambda c: np.asarray(c, float)                      # noqa: E731
    else:
        _, tf, fix_i, mob_i = struc.superimpose_homologs(qry, ref, outlier_threshold=outlier_threshold)
        _, _, fix_all, mob_all = struc.superimpose_homologs(qry, ref, max_iterations=1)
        move = tf.apply
    disp = struc.distance(qry.coord[fix_all], move(ref.coord[mob_all]))
    near = np.linalg.norm(ref.coord[mob_all] - ref_site, axis=1) <= radius
    trimmed = ~np.isin(mob_all, mob_i)
    return SiteTransfer(
        center_in_query=move(ref_site[None, :])[0],
        n_aligned=int(len(mob_all)), n_anchors=int(len(mob_i)), n_trimmed=int(trimmed.sum()),
        rmsd_anchors=float(struc.rmsd(qry.coord[fix_i], move(ref.coord[mob_i]))),
        n_site=int(near.sum()), n_site_trimmed=int((near & trimmed).sum()),
        site_rmsd=float(np.sqrt(np.mean(disp[near] ** 2))) if near.any() else None,
        site_max_displacement=float(disp[near].max()) if near.any() else None)


# --------------------------------------------------------------------------- fpocket
@dataclass(frozen=True, eq=False)
class Pocket:
    rank: int
    score: float
    druggability: float
    volume: float
    spheres: np.ndarray           # alpha-sphere centres, (n, 3)

    @property
    def centroid(self):
        return self.spheres.mean(axis=0)


_INFO_FIELDS = (("score", r"Score\s*:", float), ("druggability", r"Druggability Score\s*:", float),
                ("n_spheres", r"Number of Alpha Spheres\s*:", int), ("volume", r"Volume\s*:", float))


def parse_fpocket_info(text: str) -> list:
    parts = re.split(r"(?m)^Pocket\s+(\d+)\s*:\s*$", text)      # [preamble, rank, body, rank, body, ...]
    out = []
    for rank, body in zip(parts[1::2], parts[2::2]):
        row = {"rank": int(rank)}
        for key, pat, cast in _INFO_FIELDS:
            m = re.search(r"(?m)^\s*" + pat + r"\s*(-?[\d.]+)", body)
            if not m:
                raise ValueError("fpocket info: pocket %s has no %r" % (rank, key))
            row[key] = cast(m.group(1))
        out.append(row)
    return out


def _fpocket_banner() -> str:
    """The banner fpocket prints (ANSI codes stripped) and the binary used; the exact package build is in the env lock."""
    p = subprocess.run([str(_tool("fpocket")), "-h"], capture_output=True, text=True)
    m = re.search(r"fpocket\s+[\d.]+", re.sub(r"\x1b\[[0-9;]*m", "", p.stdout + p.stderr))
    return "%s (banner) at %s" % (m.group(0) if m else "unknown", _tool("fpocket"))


def _tool(name: str) -> Path:
    for cand in (Path(sys.executable).parent / name, shutil.which(name)):
        if cand and Path(cand).exists():
            return Path(cand)
    raise FileNotFoundError("%s not found next to %s or on PATH" % (name, sys.executable))


# fpocket silently drops hetero residues whose names it does not know (no warning in its log). Checked on 6X0J A:
# with the reduced flavin left as FDA its rank-1 site pocket had 158 alpha spheres and no FAD in its lining; relabelled
# FAD it has 257 and FAD lines it. FAD and NAP are kept (NAP: the keep/strip receptors score differently). The flavin's
# heavy atoms are identical in FAD and FDA (the difference is a hydrogen), so the relabel changes only the lookup.
FPOCKET_ALIASES = {"FDA": "FAD", "NDP": "NAP"}
FPOCKET_KNOWN_HETERO = {"FAD", "NAP", "HOH"}          # verified to be kept (HOH is dropped by design)


def heavy_atom_pdb(combined_pdb: Path, out_path: Path):
    """Step 2's work/combined.pdb (protein + cofactors, exactly what Meeko was given) minus hydrogens, for fpocket.

    Returns (path, record). The record lists what was relabelled and which hetero residues fpocket was given whose
    treatment is unverified (e.g. crystallisation ions), so that nothing about its input is silent."""
    atoms = pdb.PDBFile.read(str(combined_pdb)).get_structure(model=1)
    if (atoms.element == "").any():
        raise ValueError("%s has atoms without an element; cannot tell hydrogens from heavy atoms" % combined_pdb)
    atoms = atoms[atoms.element != "H"]
    relabelled = {}
    for old, new in FPOCKET_ALIASES.items():
        hit = atoms.res_name == old
        if hit.any():
            relabelled["%s->%s" % (old, new)] = int(hit.sum())
            atoms.res_name[hit] = new
    hetero = sorted(set(map(str, atoms.res_name[atoms.hetero])))
    f = pdb.PDBFile()
    f.set_structure(atoms)
    f.write(str(out_path))
    return out_path, {"relabelled_atoms": relabelled, "hetero_residues": hetero,
                      "hetero_unverified_in_fpocket": [h for h in hetero if h not in FPOCKET_KNOWN_HETERO]}


def run_fpocket(heavy_pdb: Path, log: Path) -> list:
    out = heavy_pdb.with_name(heavy_pdb.stem + "_out")
    shutil.rmtree(out, ignore_errors=True)                       # our own previous output; fpocket would add _out1
    _run([_tool("fpocket"), "-f", heavy_pdb], log)
    info = out / (heavy_pdb.stem + "_info.txt")
    if not info.exists():
        raise RuntimeError("fpocket wrote no %s (see %s)" % (info, log))
    pockets = []
    for r in parse_fpocket_info(info.read_text()):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")                      # biotite guesses elements of the STP pseudo-atoms
            sph = pdb.PDBFile.read(str(out / "pockets" / ("pocket%d_vert.pqr" % r["rank"]))).get_structure(model=1).coord
        if len(sph) != r["n_spheres"]:
            raise RuntimeError("pocket %d: info says %d alpha spheres, file has %d" % (r["rank"], r["n_spheres"], len(sph)))
        pockets.append(Pocket(rank=r["rank"], score=r["score"], druggability=r["druggability"], volume=r["volume"], spheres=sph))
    return pockets


@dataclass(frozen=True)
class PocketAgreement:
    site_rank: int | None
    n_spheres_in_box: int
    fraction_of_pocket_in_box: float | None
    centroid_offset: float | None          # whole-pocket alpha-sphere centroid to box centre
    in_box_centroid_offset: float | None
    agrees: bool
    reasons: list
    top: list = field(default_factory=list)   # rank, score, druggability, volume, n_spheres, n_in_box, centroid_offset


def assess_pockets(pockets, center, size, max_offset: float, n_top: int = 3) -> PocketAgreement:
    """The pocket holding most of the box is the site pocket. The box and fpocket agree only if it is rank 1 AND its
    centroid is within `max_offset` of the box centre. Anything else is listed as a reason, never smoothed over."""
    center = np.asarray(center, float)
    counts = {p.rank: int(inside_box(p.spheres, center, size).sum()) for p in pockets}
    top = [{"rank": p.rank, "score": p.score, "druggability": p.druggability, "volume": p.volume,
            "n_spheres": len(p.spheres), "n_in_box": counts[p.rank],
            "centroid_offset": float(np.linalg.norm(p.centroid - center))}
           for p in sorted(pockets, key=lambda q: q.rank)]
    best = max(counts.values(), default=0)
    if best == 0:
        return PocketAgreement(None, 0, None, None, None, False,
                               ["no fpocket pocket has an alpha sphere inside the box"], top[:n_top])
    site = min((p for p in pockets if counts[p.rank] == best), key=lambda q: q.rank)
    off = float(np.linalg.norm(site.centroid - center))
    inbox = site.spheres[inside_box(site.spheres, center, size)]
    reasons = []
    if site.rank != 1:
        r1 = next((t for t in top if t["rank"] == 1), None)
        reasons.append("the pocket holding the box is fpocket rank %d, not 1%s" % (
            site.rank, " (rank 1 is %.1f A from the box centre)" % r1["centroid_offset"] if r1 else ""))
    if off > max_offset:
        reasons.append("alpha-sphere centroid offset %.1f A exceeds %.1f A from the box centre" % (off, max_offset))
    shown = top[:n_top] + [t for t in top[n_top:] if t["rank"] == site.rank]
    return PocketAgreement(site.rank, best, best / len(site.spheres), off,
                           float(np.linalg.norm(inbox.mean(axis=0) - center)), not reasons, reasons, shown)


def ligand_pocket_contact(ligand, pockets) -> dict:
    """Which fpocket pocket owns the alpha sphere nearest to each ligand heavy atom, and how near that sphere is.
    Threshold-free: it asks where the deposited substrate is, not whether it passes a cut-off."""
    heavy = ligand.coord[ligand.element != "H"]
    owner = np.concatenate([np.full(len(p.spheres), p.rank) for p in pockets])
    spheres = np.concatenate([p.spheres for p in pockets])
    d = np.linalg.norm(heavy[:, None, :] - spheres[None, :, :], axis=2)
    nearest = d.argmin(axis=1)
    ranks, counts = np.unique(owner[nearest], return_counts=True)
    dist = d.min(axis=1)
    return {"nearest_sphere_owner": {int(r): int(c) for r, c in zip(ranks, counts)},
            "nearest_sphere_distance_median": float(np.median(dist)), "nearest_sphere_distance_max": float(dist.max())}


def pocket_lining(heavy_pdb: Path, rank: int) -> dict:
    """fpocket's own list of the atoms that line a pocket (pocketN_atm.pdb): protein residues and cofactor atoms."""
    f = heavy_pdb.with_name(heavy_pdb.stem + "_out") / "pockets" / ("pocket%d_atm.pdb" % rank)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = pdb.PDBFile.read(str(f)).get_structure(model=1)
    prot = struc.filter_canonical_amino_acids(a)
    cof = {}
    for name in np.unique(a.res_name[~prot]):
        cof[str(name)] = int((a.res_name == name).sum())
    return {"n_protein_residues": int(len(set(zip(a.chain_id[prot], a.res_id[prot])))),
            "cofactor_atoms": cof,
            "protein_residues": sorted({"%s%d" % (n, i) for n, i in zip(a.res_name[prot], a.res_id[prot])},
                                       key=lambda x: int(x[3:]))}


# --------------------------------------------------------------------------- per-structure analysis
def write_vina_box(path: Path, center, size):
    lines = ["center_%s = %.3f" % (a, c) for a, c in zip("xyz", center)] + ["size_%s = %.3f" % (a, s) for a, s in zip("xyz", size)]
    Path(path).write_text("\n".join(lines) + "\n")


def read_vina_box(path: Path):
    kv = dict(re.findall(r"(\w+)\s*=\s*(-?[\d.]+)", Path(path).read_text()))
    return np.array([float(kv["center_" + a]) for a in "xyz"]), np.array([float(kv["size_" + a]) for a in "xyz"])


def _asdict_plain(obj):
    d = asdict(obj)
    return {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in d.items()}


def analyze_job(job_dir: Path, cfg, ref_protein, ref_center, cif_cache: dict, out_dir: Path) -> dict:
    rep = json.loads((job_dir / "prep_report.json").read_text())
    pid, chain = rep["pdb_id"], rep["chain"]
    atoms = cif_cache[pid]
    size = np.array(cfg.box_size_angstrom, float)
    direct = resolve_site_center(atoms, chain, cfg.box_center_atom)
    prep_center, prep_size = read_vina_box(job_dir / "receptor.box.txt")
    if np.abs(prep_center - direct).max() > 2e-3 or np.abs(prep_size - size).max() > 1e-6:
        raise RuntimeError("%s: Step 2's receptor.box.txt (%s) disagrees with the C4X recomputed here (%s)" % (
            job_dir.name, prep_center.round(3), direct.round(3)))
    prot = protein_chain(atoms, chain)
    tr = transfer_site(ref_protein, prot, ref_center, cfg.pocket_radius_angstrom, cfg.superposition_outlier_threshold)
    carried = np.asarray(tr.center_in_query, float)
    job_out = out_dir / job_dir.name
    job_out.mkdir(parents=True, exist_ok=True)
    write_vina_box(job_out / "box.c4x.txt", direct, size)
    write_vina_box(job_out / "box.carried_from_reference.txt", carried, size)
    heavy, fp_input = heavy_atom_pdb(job_dir / "work" / "combined.pdb", job_out / "receptor_heavy.pdb")
    pockets = run_fpocket(heavy, job_out / "fpocket.log")
    res = {"job": job_dir.name, "pdb_id": pid, "chain": chain, "nadp": rep["nadp_policy"],
           "c4x_center": direct.round(3).tolist(), "carried_center": carried.round(3).tolist(),
           "offset_c4x_vs_carried": float(np.linalg.norm(direct - carried)),
           "box_overlap_c4x_vs_carried": box_overlap_fraction(direct, carried, size),
           "superposition": _asdict_plain(tr), "n_pockets": len(pockets), "fpocket_input": fp_input,
           "fpocket": {"c4x_box": asdict(assess_pockets(pockets, direct, size, cfg.pocket_agreement_offset_angstrom)),
                       "carried_box": asdict(assess_pockets(pockets, carried, size, cfg.pocket_agreement_offset_angstrom))},
           "site_pocket_lining": None, "ligand": None}
    site_rank = res["fpocket"]["c4x_box"]["site_rank"]
    if site_rank is not None:
        res["site_pocket_lining"] = pocket_lining(heavy, site_rank)
    lig_path = job_dir / "ligand_crystal.pdb"
    if lig_path.exists() and lig_path.stat().st_size > 0:
        lig = pdb.PDBFile.read(str(lig_path)).get_structure(model=1)
        res["ligand"] = {"c4x_box": ligand_vs_box(lig, direct, size), "carried_box": ligand_vs_box(lig, carried, size),
                         "fpocket": ligand_pocket_contact(lig, pockets)}
    return res


def analyze_all(cfg, prepared_root, cache_dir, out_root) -> dict:
    prepared_root, out_root, cache_dir = Path(prepared_root), Path(out_root), Path(cache_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    jobs = sorted(p.parent for p in prepared_root.glob("*/prep_report.json"))
    if not jobs:
        raise FileNotFoundError("no prepared receptors under %s; run `python -m sidacore.prep` first" % prepared_root)
    ids = {json.loads((j / "prep_report.json").read_text())["pdb_id"] for j in jobs} | {cfg.structures.fad_in}
    cif = {}
    for i in sorted(ids):
        cif[i], _ = read_structure(cache_dir / ("%s.cif" % i), cfg.alt_conformers)
    ref_atoms = cif[cfg.structures.fad_in]
    ref_protein = protein_chain(ref_atoms, cfg.receptor_chain)
    ref_center = resolve_site_center(ref_atoms, cfg.receptor_chain, cfg.box_center_atom)
    results = [analyze_job(j, cfg, ref_protein, ref_center, cif, out_root) for j in jobs]
    report = {"reference": {"structure": cfg.structures.fad_in, "chain": cfg.receptor_chain,
                            "c4x_center": ref_center.round(3).tolist()},
              "box_size": list(cfg.box_size_angstrom),
              "settings": {"box_center_atom": cfg.box_center_atom, "pocket_radius_angstrom": cfg.pocket_radius_angstrom,
                           "superposition_outlier_threshold": cfg.superposition_outlier_threshold,
                           "pocket_agreement_offset_angstrom": cfg.pocket_agreement_offset_angstrom},
              "tools": {"biotite": biotite.__version__, "fpocket": _fpocket_banner()},
              "jobs": results}
    write_resolved(cfg, out_root)
    (out_root / "step3_report.json").write_text(json.dumps(report, indent=1))
    return report


# --------------------------------------------------------------------------- report
def _f(x, nd=1):
    return "-" if x is None else ("%.*f" % (nd, x))


def render_markdown(report: dict, notes: str = "") -> str:
    ref, s = report["reference"], report["settings"]
    L = ["# Step 3 - box definition", "",
         "Generated by `python -m sidacore.boxes`; numbers only. Observations are in the hand-written notes below, if any.", "",
         "- Box centre atom: `%s` (the isoalloxazine C4a). Box size: %s A. Reference frame: %s chain %s, C4X at %s."
         % (s["box_center_atom"], " x ".join("%g" % x for x in report["box_size"]), ref["structure"], ref["chain"], ref["c4x_center"]),
         "- 'C4X box' is centred on that structure's own FAD C4X. 'Carried box' is the reference C4X moved into the structure by CA "
         "superposition (biotite `superimpose_homologs`, outlier threshold %g). Where FAD sits in a different place the two differ; "
         "both are reported." % s["superposition_outlier_threshold"], "",
         "## Boxes", "", "| receptor | C4X box centre | carried box centre | offset (A) | box overlap |", "|---|---|---|---|---|"]
    for j in report["jobs"]:
        L.append("| %s | %s | %s | %s | %.0f%% |" % (j["job"], j["c4x_center"], j["carried_center"], _f(j["offset_c4x_vs_carried"]),
                                                    100 * j["box_overlap_c4x_vs_carried"]))
    L += ["", "## Superposition onto the reference", "",
          "| receptor | aligned CA | trimmed | anchor RMSD (A) | site CA (<=%g A) | site trimmed | site RMSD (A) | site max (A) |" % s["pocket_radius_angstrom"],
          "|---|---|---|---|---|---|---|---|"]
    for j in report["jobs"]:
        t = j["superposition"]
        L.append("| %s | %d | %d | %s | %d | %d | %s | %s |" % (j["job"], t["n_aligned"], t["n_trimmed"], _f(t["rmsd_anchors"], 2),
                 t["n_site"], t["n_site_trimmed"], _f(t["site_rmsd"], 2), _f(t["site_max_displacement"], 2)))
    L += ["", "## fpocket against each box", "",
          "Agreement means the pocket holding the box is fpocket rank 1 and its centroid is within %g A of the box centre." % s["pocket_agreement_offset_angstrom"], "",
          "| receptor | box | site-pocket rank | spheres in box | centroid offset (A) | agrees | reasons |", "|---|---|---|---|---|---|---|"]
    for j in report["jobs"]:
        for key, label in (("c4x_box", "C4X"), ("carried_box", "carried")):
            a = j["fpocket"][key]
            L.append("| %s | %s | %s | %d | %s | %s | %s |" % (j["job"], label, a["site_rank"] or "-", a["n_spheres_in_box"],
                     _f(a["centroid_offset"]), "yes" if a["agrees"] else "**no**", "; ".join(a["reasons"]) or "-"))
    L += ["", "### Top fpocket pockets per receptor (C4X box)", "",
          "| receptor | rank | score | druggability | volume (A^3) | spheres | in box | centroid offset (A) |", "|---|---|---|---|---|---|---|---|"]
    for j in report["jobs"]:
        for t in j["fpocket"]["c4x_box"]["top"]:
            L.append("| %s | %d | %.3f | %.3f | %.0f | %d | %d | %s |" % (j["job"], t["rank"], t["score"], t["druggability"],
                     t["volume"], t["n_spheres"], t["n_in_box"], _f(t["centroid_offset"])))
    L += ["", "### Site pocket lining (fpocket's own atom list, C4X box)", "",
          "| receptor | site pocket rank | protein residues | cofactor atoms lining it |", "|---|---|---|---|"]
    for j in report["jobs"]:
        g = j["site_pocket_lining"]
        if g:
            L.append("| %s | %s | %d | %s |" % (j["job"], j["fpocket"]["c4x_box"]["site_rank"], g["n_protein_residues"],
                     ", ".join("%s %d" % kv for kv in sorted(g["cofactor_atoms"].items())) or "none"))
    withlig = [j for j in report["jobs"] if j["ligand"]]
    if withlig:
        L += ["", "## Deposited substrate against each box", "",
              "| receptor | box | heavy atoms inside | ligand centroid to box centre (A) |", "|---|---|---|---|"]
        for j in withlig:
            for key, label in (("c4x_box", "C4X"), ("carried_box", "carried")):
                g = j["ligand"][key]
                L.append("| %s | %s | %d / %d | %s |" % (j["job"], label, g["n_inside"], g["n_heavy"], _f(g["centroid_distance"])))
        L += ["", "Which fpocket pocket owns the alpha sphere nearest each deposited-ORN heavy atom (no threshold):", "",
              "| receptor | atoms per owning pocket rank | nearest-sphere distance, median / max (A) |", "|---|---|---|"]
        for j in withlig:
            g = j["ligand"]["fpocket"]
            L.append("| %s | %s | %s / %s |" % (j["job"], ", ".join("rank %s: %d" % (k, v) for k, v in sorted(g["nearest_sphere_owner"].items(), key=lambda kv: int(kv[0]))),
                     _f(g["nearest_sphere_distance_median"]), _f(g["nearest_sphere_distance_max"])))
    rel = {k: v for j in report["jobs"] for k, v in j["fpocket_input"]["relabelled_atoms"].items()}
    unv = sorted({h for j in report["jobs"] for h in j["fpocket_input"]["hetero_unverified_in_fpocket"]})
    L += ["", "## fpocket input", "",
          "fpocket silently drops hetero residues it does not know by name. The reduced flavin `FDA` (6X0J) was relabelled `FAD` "
          "for fpocket only (relabelled: %s); the heavy atoms are identical. Hetero residues given to fpocket whose handling is "
          "unverified: %s." % (rel or "none", unv or "none")]
    L += ["", "## Tools", "", "- " + "\n- ".join("%s %s" % kv for kv in report["tools"].items())]
    if notes.strip():
        L += ["", "## Notes", "", notes.strip()]
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    import argparse
    import yaml
    from .config import load_config
    ap = argparse.ArgumentParser(description="Step 3: C4X box, carried box by superposition, fpocket cross-check.")
    ap.add_argument("--config")
    ap.add_argument("--prepared", help="Step 2 output root")
    ap.add_argument("--cache", default="structures_cache")
    ap.add_argument("--out")
    ap.add_argument("--render-only", metavar="OUT_DIR", help="rebuild --report from OUT_DIR/step3_report.json; no fpocket run")
    ap.add_argument("--report", help="write a markdown report here")
    ap.add_argument("--notes", help="hand-written notes file appended to the report verbatim")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
    a = ap.parse_args(argv)
    if a.render_only:
        if not a.report:
            ap.error("--render-only needs --report")
        notes = Path(a.notes).read_text() if a.notes and Path(a.notes).exists() else ""
        Path(a.report).write_text(render_markdown(json.loads((Path(a.render_only) / "step3_report.json").read_text()), notes))
        return 0
    if not (a.config and a.prepared and a.out):
        ap.error("--config, --prepared and --out are required unless --render-only is given")
    over = {k: yaml.safe_load(v) for k, _, v in (s.partition("=") for s in a.set)}
    report = analyze_all(load_config(a.config, overrides=over), a.prepared, a.cache, a.out)
    if a.report:
        notes = Path(a.notes).read_text() if a.notes and Path(a.notes).exists() else ""
        Path(a.report).write_text(render_markdown(report, notes))
    for j in report["jobs"]:
        print("%-22s C4X %s | carried %s | offset %.1f A | site RMSD %s | fpocket agrees: C4X %s, carried %s" % (
            j["job"], j["c4x_center"], j["carried_center"], j["offset_c4x_vs_carried"], _f(j["superposition"]["site_rmsd"], 2),
            j["fpocket"]["c4x_box"]["agrees"], j["fpocket"]["carried_box"]["agrees"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
