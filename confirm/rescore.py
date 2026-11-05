"""Step 4b: rescore the existing Vina poses with other scoring functions. The rule is `reports/step4_rescoring_rule.md`
(registered in commit 92d2d0b before any pose was rescored); this module implements it and decides nothing of its own.

Functions (exactly these): Vina (reproducibility gate G0 only), AutoDock4 (PRIMARY), Vinardo (secondary).

LICENCE: AutoGrid 4.2.9 (GPL-2.0-only) is run as a separate program, from the `sida-docking-ad4` env, and only the map files it
writes are used; `vina` (Apache-2.0) comes from `sida-docking`. `sidacore` imports nothing from here and nothing here is linked to
anything GPL. Delete `confirm/` and nothing else changes.

Poses are scored as docked: no minimisation, no re-docking, no flexibility.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- registered constants
AD4_RECEPTOR_TYPES = ("A", "C", "HD", "N", "NA", "OA", "P", "SA")
AD4_LIGAND_TYPES = ("C", "HD", "N", "OA")
SPACING = 0.375
NPTS_DEFAULT, NPTS_FALLBACK = 60, 80
SMOOTH, DIELECTRIC = 0.5, -0.1465
TIE_TOL = 0.05                  # kcal/mol: poses this close to the best are tied
G0_TOL = 0.05                   # kcal/mol: F0 must reproduce the docking log this closely
ACCEPT = 2.0                    # angstrom, the Step 4 pass rule
N_SEEDS = 10
LIGAND_RESNAME = "ORN"
PRIMARY = "4B63_A"
CONTROL_LABELS = ("6X0J_A", "6X0J_B", "6X0J_D")
REPORTED_LABELS = (PRIMARY,) + CONTROL_LABELS + ("6X0J_C",)       # chain C: reported, gates nothing
FUNCTIONS = ("vina", "ad4", "vinardo")
AUTOGRID_DEFAULT = Path("/opt/miniforge3/envs/sida-docking-ad4/bin/autogrid4")


# --------------------------------------------------------------------------- pose files
def split_models(pdbqt_text: str) -> list:
    """Single-model PDBQT texts, in file order. Vina rejects a MODEL tag in a ligand file ("Unexpected multi-MODEL tag")."""
    parts, cur = [], None
    for ln in pdbqt_text.splitlines():
        if ln.startswith("MODEL"):
            cur = []
            parts.append(cur)
        elif ln.startswith("ENDMDL"):
            cur = None
        elif cur is not None:
            cur.append(ln)
    if not parts:
        return [pdbqt_text]
    return ["\n".join(p) + "\n" for p in parts]


def docked_scores(pdbqt_text: str) -> list:
    return [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)", pdbqt_text)]


def logged_interaction_energies(pdbqt_text: str, tol: float = 0.01) -> list:
    """Amendment 1: the interaction-only score `vina --score_only` returns, read from the docking output's REMARK lines.
    Logged score_k = (INTER_k + INTRA_k - UNBOUND) / f with UNBOUND the first pose's internal energy; --score_only takes each input as
    its own reference, so it returns INTER_k / f. f = INTER_1 / logged_1 is valid only if INTRA_1 == UNBOUND, which is checked."""
    def grab(key):
        return [float(x) for x in re.findall(r"(?m)^REMARK %s:\s+(-?\d+\.?\d*)" % key, pdbqt_text)]
    res, inter, intra, unb = docked_scores(pdbqt_text), grab("INTER"), grab("INTRA"), grab("UNBOUND")
    if not res or not (len(res) == len(inter) == len(intra) == len(unb)):
        raise ValueError("the docking output lacks the INTER / INTRA / UNBOUND remarks for every pose (%d scores, %d INTER)" % (len(res), len(inter)))
    if abs(intra[0] - unb[0]) > tol:
        raise ValueError("the first pose is not its own unbound reference (INTRA %.3f, UNBOUND %.3f)" % (intra[0], unb[0]))
    f = inter[0] / res[0]
    return [x / f for x in inter]


def pdbqt_coords(pdbqt_text: str) -> np.ndarray:
    """Coordinates of every ATOM/HETATM line, polar hydrogens included (they carry AD4 maps)."""
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in pdbqt_text.splitlines() if l.startswith(("ATOM", "HETATM"))], dtype=float).reshape(-1, 3)


# --------------------------------------------------------------------------- scoring
def parse_estimated_energy(vina_output: str) -> float:
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.?\d*)", vina_output)
    if not m:
        raise ValueError("no 'Estimated Free Energy of Binding' line in the vina output")
    return float(m.group(1))


def run_vina_score(vina, scoring: str, ligand, receptor=None, maps=None) -> float:
    """One `--score_only` evaluation. vina/vinardo take the receptor PDBQT; ad4 takes AutoGrid maps and refuses --receptor."""
    if scoring not in FUNCTIONS:
        raise ValueError("scoring must be one of %s, got %r" % (FUNCTIONS, scoring))
    cmd = [str(vina), "--ligand", str(ligand), "--scoring", scoring, "--score_only", "--cpu", "1"]
    if scoring == "ad4":
        if maps is None:
            raise ValueError("ad4 scoring needs maps")
        cmd += ["--maps", str(maps)]
    else:
        if receptor is None:
            raise ValueError("%s scoring needs a receptor" % scoring)
        cmd += ["--receptor", str(receptor), "--autobox"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = p.stdout + p.stderr
    if p.returncode != 0:
        raise RuntimeError("vina failed (%d) on %s: %s" % (p.returncode, ligand, out.strip()[-300:]))
    try:
        return parse_estimated_energy(out)
    except ValueError as e:
        raise RuntimeError("%s on %s: %s" % (e, ligand, out.strip()[-300:])) from e


# --------------------------------------------------------------------------- AutoGrid
def write_gpf(receptor_file: str, prefix: str, center, npts: int = NPTS_DEFAULT) -> str:
    """The grid parameter file, with exactly the registered values and no others."""
    L = ["npts %d %d %d" % (npts, npts, npts),
         "gridfld %s.maps.fld" % prefix,
         "spacing %g" % SPACING,
         "receptor_types %s" % " ".join(AD4_RECEPTOR_TYPES),
         "ligand_types %s" % " ".join(AD4_LIGAND_TYPES),
         "receptor %s" % receptor_file,
         "gridcenter %.3f %.3f %.3f" % tuple(center),
         "smooth %g" % SMOOTH]
    L += ["map %s.%s.map" % (prefix, t) for t in AD4_LIGAND_TYPES]
    L += ["elecmap %s.e.map" % prefix, "dsolvmap %s.d.map" % prefix, "dielectric %g" % DIELECTRIC]
    return "\n".join(L) + "\n"


def outside_map(coords, center, npts: int = NPTS_DEFAULT) -> bool:
    """True if any atom lies outside the cubic map of `npts` points about `center` (the registered out-of-map rule)."""
    half = npts * SPACING / 2
    return bool(np.any(np.abs(np.asarray(coords, float) - np.asarray(center, float)) > half + 1e-9))


def run_autogrid(autogrid, receptor_pdbqt: Path, center, npts: int, work: Path) -> str:
    """Maps for one receptor; returns the map prefix to pass to `vina --maps`."""
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy(receptor_pdbqt, work / "receptor.pdbqt")
    (work / "receptor.gpf").write_text(write_gpf("receptor.pdbqt", "receptor", center, npts))
    p = subprocess.run([str(autogrid), "-p", "receptor.gpf", "-l", "receptor.glg"], cwd=work, capture_output=True, text=True)
    needed = ["receptor.maps.fld", "receptor.e.map", "receptor.d.map"] + ["receptor.%s.map" % t for t in AD4_LIGAND_TYPES]
    missing = [n for n in needed if not (work / n).exists()]
    if p.returncode != 0 or missing:
        raise RuntimeError("autogrid4 failed in %s (exit %d, missing %s): %s" % (work, p.returncode, missing, (p.stdout + p.stderr)[-300:]))
    return str(work / "receptor")


# --------------------------------------------------------------------------- ranking and the registered rule
def headline_pose(scores, tol: float = TIE_TOL):
    """(index of the headline top-1 pose, tied?). Lowest score wins; poses within `tol` of the best are tied and the one Vina
    ranked highest (lowest index, since the pose files are in Vina order) is the headline."""
    best = min(scores)
    within = [i for i, s in enumerate(scores) if s <= best + tol + 1e-9]
    return within[0], len(within) > 1


def seed_result(rescored, rmsds, accept: float = ACCEPT) -> dict:
    if len(rescored) != len(rmsds):
        raise ValueError("%d scores but %d RMSDs" % (len(rescored), len(rmsds)))
    top, tied = headline_pose(rescored)
    best = int(np.argmin(rmsds))
    return {"top1_idx": top, "top1_rmsd": float(rmsds[top]), "tied": tied,
            "best_idx": best, "best_rmsd": float(rmsds[best]),
            "best_rank": 1 + sum(s < rescored[best] for s in rescored),        # rank of the best-RMSD pose under this function
            "has_near_native": bool(min(rmsds) <= accept)}


def summarise_job(seeds: list, accept: float = ACCEPT) -> dict:
    t1 = [s["top1_rmsd"] for s in seeds]
    return {"n_seeds": len(seeds), "accept": accept, "top1_median": statistics.median(t1), "top1_min": min(t1), "top1_max": max(t1),
            "n_within": sum(x <= accept for x in t1), "passes": bool(statistics.median(t1) <= accept),
            "n_tied": sum(s["tied"] for s in seeds), "ceiling_n": sum(s["has_near_native"] for s in seeds),
            "best_rank_median": statistics.median(s["best_rank"] for s in seeds)}


def g0_check(logged, rescored, tol: float = G0_TOL) -> dict:
    """Gate G0: the pipeline, scoring with Vina's own function, must reproduce the docking log."""
    if len(logged) != len(rescored):
        raise ValueError("%d logged scores but %d rescored" % (len(logged), len(rescored)))
    d = np.abs(np.asarray(logged, float) - np.asarray(rescored, float))
    n_bad = int((d > tol + 1e-9).sum())
    return {"ok": n_bad == 0, "n": len(d), "n_bad": n_bad, "max_abs_diff": float(d.max()), "tol": tol}


def is_admissible(control_passes: dict) -> dict:
    """Gate G1: all three Vina-passing controls must still pass. A missing control is an error, never a silent pass."""
    broken = [c for c in CONTROL_LABELS if not control_passes[c]]
    return {"admissible": not broken, "broken": broken}


def outcome_key(admissible: bool, passes: bool) -> str:
    if not admissible:
        return "cannot_arbitrate"
    return "rescues" if passes else "no_rescue"


FUNCTION_NAMES = {"ad4": "AutoDock4", "vinardo": "Vinardo", "vina": "Vina"}


BASELINE_SENTENCE = (" Vina without its internal-strain term already does the same, so dropping that term, not the function, accounts for it.")


def outcome_statement(fn: str, key: str, broken=(), moot: bool = False) -> str:
    """The registered statements (Amendment 1 adds the baseline sentence to F1 / F2 when Vina interaction-only already passes)."""
    n = FUNCTION_NAMES[fn]
    text = {"rescues": "%s scoring ranks a near-native pose first on 4B63 among Vina's poses and keeps the controls. A finding about "
                       "scoring functions, not a fix: 4B63 stays primary and its Vina failure stays in the report." % n,
            "no_rescue": "%s scoring does not rescue 4B63 either." % n,
            "cannot_arbitrate": "%s scoring cannot arbitrate this site: it breaks control(s) %s." % (n, ", ".join(broken))}[key]
    return text + (BASELINE_SENTENCE if moot and fn != "vina" else "")


# --------------------------------------------------------------------------- the run
def _label(job_dir: Path) -> str:
    return job_dir.name.replace("_nadp-keep", "")


def _read_box_center(box_file: Path):
    v = {k.strip(): x.strip() for k, x in (l.split("=") for l in box_file.read_text().split("\n") if "=" in l)}
    return tuple(float(v["center_" + a]) for a in "xyz")


def load_job(course: Path, label: str, work: Path) -> dict:
    """Poses of every seed (split to single-model files), the deposited-ligand reference and the box centre."""
    from sidacore.redock import crystal_reference, ligand_template, pose_rmsd, read_poses   # noqa: F401 - sidacore is the only import
    job = "%s_nadp-keep" % label
    prep, dock = course / "prepared" / job, course / "redock" / job / "vina"
    seeds = {}
    files = sorted(dock.glob("vina_seed*.pdbqt"), key=lambda f: int(re.search(r"seed(\d+)", f.name).group(1)))
    if len(files) != N_SEEDS:
        raise ValueError("%s: the rule needs %d seeds, found %d" % (label, N_SEEDS, len(files)))
    for f in files:
        seed = int(re.search(r"seed(\d+)", f.name).group(1))
        text = f.read_text()
        parts, logged = split_models(text), docked_scores(text)
        if len(parts) != len(logged):
            raise ValueError("%s seed %d: %d models but %d scores" % (label, seed, len(parts), len(logged)))
        paths = []
        for i, part in enumerate(parts):
            p = work / label / ("seed%d" % seed) / ("pose%02d.pdbqt" % (i + 1))
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(part)
            paths.append(p)
        mols = read_poses(text)
        if [s for s, _ in mols] != logged:
            raise ValueError("%s seed %d: pose file is not in Vina's score order, RMSDs would be misaligned" % (label, seed))
        seeds[seed] = {"paths": paths, "logged": logged, "inter": logged_interaction_energies(text), "mols": [m for _, m in mols],
                       "coords": [pdbqt_coords(t) for t in parts]}
    ref = crystal_reference(str(prep / "ligand_crystal.pdb"), ligand_template(LIGAND_RESNAME))
    return {"label": label, "receptor": prep / "receptor.pdbqt", "center": _read_box_center(prep / "receptor.box.txt"),
            "ref": ref, "seeds": seeds}


def score_all(vina, jobs: dict, fn: str, n_jobs: int, maps: dict | None = None) -> dict:
    """{label: {seed: [score per pose]}} for one function."""
    tasks = [(lab, seed, i, p) for lab, j in jobs.items() for seed, s in j["seeds"].items() for i, p in enumerate(s["paths"])]

    def one(t):
        lab, seed, i, p = t
        if fn == "ad4":
            return run_vina_score(vina, "ad4", p, maps=maps[lab])
        return run_vina_score(vina, fn, p, receptor=jobs[lab]["receptor"])

    with ThreadPoolExecutor(n_jobs) as ex:
        vals = list(ex.map(one, tasks))
    out = {lab: {seed: [None] * len(s["paths"]) for seed, s in j["seeds"].items()} for lab, j in jobs.items()}
    for (lab, seed, i, _), v in zip(tasks, vals):
        out[lab][seed][i] = v
    return out


def run(course: Path, out: Path, vina, autogrid, n_jobs: int, g0_only: bool = False) -> dict:
    from sidacore.redock import pose_rmsd
    out.mkdir(parents=True, exist_ok=True)
    work = out / "poses"
    jobs = {lab: load_job(course, lab, work) for lab in REPORTED_LABELS}

    # ---- G0: nothing else is computed or read unless Vina's own function reproduces the docking log
    f0 = score_all(vina, jobs, "vina", n_jobs)
    logged = [s for lab in jobs for seed in jobs[lab]["seeds"] for s in jobs[lab]["seeds"][seed]["inter"]]       # Amendment 1
    got = [s for lab in jobs for seed in jobs[lab]["seeds"] for s in f0[lab][seed]]
    g0 = g0_check(logged, got)
    (out / "g0.json").write_text(json.dumps(g0, indent=1))
    print("G0:", g0, flush=True)
    result = {"g0": g0, "functions": {}}
    if not g0["ok"] or g0_only:
        (out / "result.json").write_text(json.dumps(result, indent=1))
        return result

    # ---- maps: out-of-map rule decided once for all jobs
    npts = NPTS_DEFAULT
    outside = {lab: sum(outside_map(c, j["center"], NPTS_DEFAULT) for s in j["seeds"].values() for c in s["coords"]) for lab, j in jobs.items()}
    if any(outside.values()):
        npts = NPTS_FALLBACK
    result["grid"] = {"npts": npts, "n_poses_outside_60": outside,
                      "n_poses_outside_final": {lab: sum(outside_map(c, j["center"], npts) for s in j["seeds"].values() for c in s["coords"])
                                                for lab, j in jobs.items()}}
    print("grid:", result["grid"], flush=True)
    maps = {lab: run_autogrid(autogrid, j["receptor"], j["center"], npts, out / "maps" / lab) for lab, j in jobs.items()}

    scores = {"vina": f0, "ad4": score_all(vina, jobs, "ad4", n_jobs, maps), "vinardo": score_all(vina, jobs, "vinardo", n_jobs)}
    rows, per_fn = [], {fn: {} for fn in FUNCTIONS}
    for lab, j in jobs.items():
        for seed, s in j["seeds"].items():
            rm = [pose_rmsd(m, j["ref"]) for m in s["mols"]]
            vina_top = s["mols"][0]
            for fn in FUNCTIONS:
                r = seed_result(scores[fn][lab][seed], rm)
                r["top1_is_vina_top1"] = bool(pose_rmsd(s["mols"][r["top1_idx"]], vina_top) <= ACCEPT)
                per_fn[fn].setdefault(lab, []).append(r)
            for i in range(len(rm)):
                rows.append({"job": lab, "seed": seed, "vina_rank": i + 1, "logged": s["logged"][i], "vina": scores["vina"][lab][seed][i],
                             "ad4": scores["ad4"][lab][seed][i], "vinardo": scores["vinardo"][lab][seed][i], "rmsd": round(rm[i], 3)})
    with (out / "poses.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    moot = False
    for fn in FUNCTIONS:                                  # FUNCTIONS starts with vina, so the baseline is known before F1 / F2
        summ = {lab: summarise_job(per_fn[fn][lab]) for lab in jobs}
        adm = is_admissible({c: summ[c]["passes"] for c in CONTROL_LABELS})
        key = outcome_key(adm["admissible"], summ[PRIMARY]["passes"])
        if fn == "vina":
            moot = key == "rescues"
        result["functions"][fn] = {"jobs": summ, "admissible": adm, "outcome": key,
                                   "statement": outcome_statement(fn, key, adm["broken"], moot=moot), "seeds": per_fn[fn]}
    (out / "result.json").write_text(json.dumps(result, indent=1, default=float))
    return result


FN_TITLES = {"vina": "Vina, interaction energy only (baseline)", "ad4": "AutoDock4 (PRIMARY)", "vinardo": "Vinardo (secondary)"}


def render_markdown(result: dict, notes: str = "") -> str:
    g0, L = result["g0"], []
    L += ["# Step 4b: rescoring the Vina poses (generated by `confirm.rescore`; rule `step4_rescoring_rule.md` + Amendment 1)", "",
          "Poses scored as docked, interaction energy only. Pass = median over %d seeds of the rescored top-1 RMSD <= %g A (the Step 4 rule). "
          "Controls: %s. Chain 6X0J_C is reported and gates nothing." % (N_SEEDS, ACCEPT, ", ".join(CONTROL_LABELS)), "",
          "**G0** (Vina reproduces INTER/f of the docking log on every pose): %s; %d poses, %d beyond %g kcal/mol, maximum difference %.3f." % (
              "PASSED" if g0["ok"] else "FAILED", g0["n"], g0["n_bad"], g0["tol"], g0["max_abs_diff"]), ""]
    if "grid" in result:
        gr = result["grid"]
        L += ["**AutoGrid map:** %d points per axis (%.1f A). Poses outside the 60-point map, per job: %s; outside the final map: %s." % (
            gr["npts"], gr["npts"] * SPACING, gr["n_poses_outside_60"], gr["n_poses_outside_final"]), ""]
    for fn in FUNCTIONS:
        if fn not in result.get("functions", {}):
            continue
        r = result["functions"][fn]
        adm = r["admissible"]
        L += ["## %s" % FN_TITLES[fn], "",
              "| job | top-1 median (A) | min - max | seeds within %g A | pass | seeds with any pose within %g A | tied seeds | median rank of best pose |" % (ACCEPT, ACCEPT),
              "|---|---|---|---|---|---|---|---|"]
        for lab in REPORTED_LABELS:
            if lab in r["jobs"]:
                j = r["jobs"][lab]
                L.append("| %s | %.2f | %.2f - %.2f | %d/%d | %s | %d/%d | %d | %g |" % (
                    lab, j["top1_median"], j["top1_min"], j["top1_max"], j["n_within"], j["n_seeds"], "**PASS**" if j["passes"] else "fail",
                    j["ceiling_n"], j["n_seeds"], j["n_tied"], j["best_rank_median"]))
        L += ["", "Admissible (all three controls pass): **%s**%s." % ("yes" if adm["admissible"] else "no",
                                                                          "" if adm["admissible"] else "; broken: " + ", ".join(adm["broken"])),
              "", "> %s" % r["statement"], ""]
    return "\n".join(L) + ("\n" + notes if notes else "")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Step 4b: rescore the Vina poses (rule: reports/step4_rescoring_rule.md).")
    ap.add_argument("--course", type=Path, default=Path("."))
    ap.add_argument("--out", type=Path, default=Path("rescoring"))
    ap.add_argument("--vina", type=Path, default=None)
    ap.add_argument("--autogrid", type=Path, default=AUTOGRID_DEFAULT)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--report", type=Path, default=None, help="write the generated Markdown report here")
    ap.add_argument("--notes", type=Path, default=None, help="hand-written notes appended to the report")
    ap.add_argument("--render-only", action="store_true", help="rewrite the report from <out>/result.json; scores nothing")
    ap.add_argument("--g0-only", action="store_true", help="run gate G0 and stop, computing and reading nothing else")
    a = ap.parse_args(argv)
    if a.render_only:
        if not a.report:
            ap.error("--render-only needs --report")
        a.report.write_text(render_markdown(json.loads((a.out / "result.json").read_text()), a.notes.read_text() if a.notes else ""))
        return 0
    from sidacore.boxes import _tool
    r = run(a.course.resolve(), a.out, a.vina or _tool("vina"), a.autogrid, a.jobs, a.g0_only)
    if not r["g0"]["ok"]:
        print("G0 FAILED: the pipeline does not reproduce the docking log; nothing else was computed.", file=sys.stderr)
        return 2
    if a.report and "functions" in r:
        a.report.write_text(render_markdown(r, a.notes.read_text() if a.notes else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
