# Screening notebook: dock literature ligands into a target and look for clashes

`screen.ipynb` takes crystal structures (PDB ID + chain) and ligands (SMILES, e.g. from a paper), docks every ligand into every
target with AutoDock Vina, and reports for each pair whether the best pose **clashes** with the protein, whether it sits **in the site**,
and whether the protocol could be trusted on that target (a **redock control** of the deposited ligand). The worked example is SidA with
L-ornithine. Code: `sidacore/screen.py`; tests: `tests/test_screen.py`, `tests/test_screen_notebook.py`.

It answers *can this molecule be placed here without overlapping the protein*. It does **not** rank affinity: on SidA, Vina scores mis-rank
(see the notebook's "How to read the results" and `reports/step4_*`).

## Set up (Linux or macOS, once)

```bash
git clone git@github.com:leonardomarino/sida-docking-course.git && cd sida-docking-course
git checkout feat/sida-docking-prototype
conda env create -f notebooks/environment.yml     # name: sida-screen (conda-forge: vina, meeko, rdkit, biotite, openbabel; pip: pdb2pqr, propka)
conda activate sida-screen
pip install -e .                                   # makes `sidacore` importable
cd notebooks && jupyter lab screen.ipynb
```

Run all cells (about 75 s for the example on a 24-core Linux machine; expect several minutes on a laptop). Structures are downloaded from
RCSB into `notebooks/structures_cache/` (internet needed the first time). Outputs go to `notebooks/screen_runs/<RUN_NAME>/`:
`results.csv`, `controls.csv`, `run.json` (settings and tool versions), `poses/` (best pose per pair, PDB), `targets/` (prepared receptors).

**macOS:** the same pins resolve for osx-arm64 (checked with a dry-run solve on Linux; no Mac has run the notebook yet). The last notebook cell compares your
table with `reference_results.csv`, produced on Linux with CPU Vina 1.2.7. Vina is deterministic per seed on one machine; small score differences across
machines are possible, a different `reading` or a different number of clashes is worth a look.

## Your own run

Edit the "parameters" cell (section 2): `TARGETS` (one `S.TargetSpec` per structure), `LIGANDS` (dict, or the path of a `.csv` with columns `name,smiles`, or a `.smi` file),
`SETTINGS`. Section 3 shows what will be kept and stripped from each structure before anything is docked.

* `box_size` is the setting that matters most: **12 Å = site-sized** (a ligand that does not fit is forced to clash); **22 Å = roomy** (an oversized ligand docks clash-free at the rim, which is not a fit). Compare both before concluding.
* New targets: `site_ligand` is the residue name of the ligand deposited in the site you care about. For a structure without one, give `center=(x, y, z)` and there is no control.
* Only FAD, FDA and NADP can be kept as cofactors (the only ones with a shipped Meeko template); NADPH (`NDP`), heme, FMN, etc. are stripped, and section 3 says so.

## Tests

```bash
pytest -m "not slow" tests/test_screen.py     # 39 tests, ~2 s: clash arithmetic, ligand/target handling, no external programs
pytest -m slow tests/test_screen.py tests/test_screen_notebook.py   # run PDB2PQR, Meeko, Vina and the whole notebook (~2 min)
```

## What was measured while building it (SidA, ligand-centred box, mj-ai 2026-11-05)

* The generic preparation reproduces the validated Step 2 receptor for 4B63 atom for atom (test `test_generic_preparation_of_4b63_reproduces_the_validated_receptor`).
* Redock control, 5 seeds: 4B63 **fails** at every box size 12-22 Å (median top-1 RMSD 4.7-4.9 Å); 6X0J chain A passes at every size (about 1 Å);
  6X0J chain C **fails at 22 Å (5.96 Å, as in Step 4) and passes at 12 Å (1.06 Å, 4 of 5 seeds)**: the control gets easier as the box tightens around the deposited ligand.
* Oversized ligand (tetraphenylmethane): at 12 Å it is forced into the pocket (best-pose score +0.4 to +9 kcal/mol, 8-15 clashes, no clash-free pose); at 16-22 Å it docks at the rim with 0 clashes and a score better than the substrate's.
* Clash thresholds (0.4 Å overlap, 2.4 Å N/O) come from five deposited ornithine poses: a small calibration set.
