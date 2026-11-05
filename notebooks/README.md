# Screening notebook: dock literature ligands into a target and look for clashes

`screen.ipynb` takes crystal structures (PDB ID + chain) and ligands (SMILES, e.g. from a paper), docks every ligand into every
target with AutoDock Vina, and reports for each pair whether the best pose **clashes** with the protein, whether it sits **in the site**,
and whether the protocol could be trusted on that target (a **redock control** of the deposited ligand). The worked example is SidA with
L-ornithine. Code: `sidacore/screen.py`; tests: `tests/test_screen.py`, `tests/test_screen_notebook.py`.

It answers *can this molecule be placed here without overlapping the protein*. It does **not** rank affinity: on SidA, Vina scores mis-rank
(see the notebook's "How to read the results" and `reports/step4_*`).

---

## Install on a Mac

Everything runs in one conda environment (`sida-screen`) that holds AutoDock Vina, Meeko, RDKit, Open Babel, biotite, PDB2PQR/PROPKA and JupyterLab.
You do not need a GPU, Docker, or any other program.

**Tested:** an Apple-silicon Mac, conda from Homebrew's `miniconda` cask, 2026-11-05: the environment installed, all four external programs were found,
and the notebook ran end to end. **Intel Macs:** the same packages resolve (checked with a dry-run solve only); nobody has run it there.

### 1. Prerequisites: git and conda

```bash
git --version        # if this opens an installer prompt, accept it (Xcode command line tools)
conda --version      # any conda or mamba will do
```

No conda? Install one, then **open a new terminal**:

```bash
brew install --cask miniforge        # or: brew install --cask miniconda
```

(Homebrew itself: https://brew.sh. The environment file asks only for the free `conda-forge` channel.)

### 2. Get the code

```bash
git clone git@github.com:leonardomarino/sida-docking-course.git ~/GitHub/sida-docking-course
cd ~/GitHub/sida-docking-course
git checkout feat/sida-docking-prototype
```

The repository is **private**: you need a GitHub account that has been added to it as a collaborator, and either an SSH key registered on that account
(as above) or the HTTPS address `https://github.com/leonardomarino/sida-docking-course.git` with `gh auth login` / a personal access token.

(Already cloned? `git checkout feat/sida-docking-prototype && git pull`.)

### 3. Create the environment (once, a few minutes, about 250 MB)

```bash
conda env create -f notebooks/environment.yml
```

### 4. Activate it

```bash
conda activate sida-screen
```

If this prints `CondaError: Run 'conda init' before 'conda activate'`, your shell has not been told about conda yet. Either load it for this terminal only:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sida-screen
```

or fix it permanently (this edits your shell startup file; open a new terminal afterwards): `conda init zsh` (or `conda init bash` if your shell is bash; `echo $SHELL` tells you).
The prompt then starts with `(sida-screen)`.

### 5. Install the project code into the environment (once)

Run this from the repository root (the folder that contains `pyproject.toml`):

```bash
pip install -e .
```

### 6. Check the install

```bash
python -c "import sidacore.screen as S; print(S.check_tools())"
```

You should see four paths ending in `.../envs/sida-screen/bin/...` (`vina`, `obabel`, `mk_prepare_receptor.py`, `pdb2pqr30`). A `None` means that program is missing:
see Troubleshooting. Optional quick test (a few seconds; some tests are skipped because the structure cache is not in the repository, which is expected):

```bash
pytest -m "not slow" tests/test_screen.py -q
```

### 7. Run the notebook

```bash
cd notebooks
jupyter lab screen.ipynb
```

In JupyterLab: **Run > Run All Cells**. The first run downloads 4B63 and 6X0J from RCSB (internet needed) and takes several minutes on a laptop
(75 s on the 24-core Linux machine it was built on); later runs reuse the downloads and the prepared receptors.

The last cell compares your table with `reference_results.csv`, produced on Linux with CPU Vina 1.2.7. Vina is deterministic per seed on one machine; small score
differences across machines are possible. A different `reading` or a different number of clashes is worth a look.

### Later

```bash
cd ~/GitHub/sida-docking-course && git pull                       # new code
conda env update -f notebooks/environment.yml --prune             # only if environment.yml changed
conda env remove -n sida-screen                                   # to uninstall everything
```

### Troubleshooting

| Symptom | Fix |
|---|---|
| `CondaError: Run 'conda init' before 'conda activate'` | step 4 above |
| `ModuleNotFoundError: No module named 'sidacore'` | `pip install -e .` from the repository root, with `(sida-screen)` active; in JupyterLab check the kernel (top right) is "Python 3" of `sida-screen`, and start `jupyter lab` from the activated environment |
| `check_tools()` shows `None` | the environment is incomplete: `conda env remove -n sida-screen`, then repeat step 3. For `pdb2pqr30` alone: `pip install pdb2pqr==3.7.1 propka==3.5.1` |
| The notebook stops at the first cell with "missing programs" | same as above; the cell names the missing program |
| A structure cannot be downloaded | RCSB needs internet (and a working proxy if you are behind one). Put the `.cif` file into `notebooks/structures_cache/` yourself, named like `4B63.cif` |
| It is slow | Vina uses every core. Lower `exhaustiveness` or the number of `seeds` in the settings cell, or dock fewer ligands |

---

## Use it with your own targets and ligands

Edit the "parameters" cell (section 2 of the notebook): `TARGETS` (one `S.TargetSpec` per structure), `LIGANDS` (a dict, or the path of a `.csv` with columns `name,smiles`, or a `.smi` file),
`SETTINGS`. Section 3 shows what will be kept and stripped from each structure before anything is docked. Results go to `notebooks/screen_runs/<RUN_NAME>/`:
`results.csv`, `controls.csv`, `run.json` (settings and tool versions), `poses/` (best pose per pair, PDB), `targets/` (prepared receptors).

* `box_size` is the setting that matters most: **12 Å = site-sized** (a ligand that does not fit is forced to clash); **22 Å = roomy** (an oversized ligand docks clash-free at the rim, which is not a fit). Compare both before concluding.
* New targets: `site_ligand` is the residue name of the ligand deposited in the site you care about. For a structure without one, give `center=(x, y, z)` and there is no control.
* Only FAD, FDA and NADP can be kept as cofactors (the only ones with a shipped Meeko template); NADPH (`NDP`), heme, FMN, etc. are stripped, and section 3 says so.

## Tests (for developers)

```bash
pytest -m "not slow" tests/test_screen.py     # 39 tests, ~2 s: clash arithmetic, ligand/target handling, no external programs
pytest -m slow tests/test_screen.py tests/test_screen_notebook.py   # run PDB2PQR, Meeko, Vina and the whole notebook (~2 min on 24 cores)
```

The other test files in `tests/` belong to the SidA docking course and need the `sida-docking` environment (they use `fpocket`, which is not in `sida-screen`); they are not needed to use the notebook.

## What was measured while building it (SidA, ligand-centred box, Linux 2026-11-05)

* The generic preparation reproduces the validated Step 2 receptor for 4B63 atom for atom (test `test_generic_preparation_of_4b63_reproduces_the_validated_receptor`).
* Redock control, 5 seeds: 4B63 **fails** at every box size 12-22 Å (median top-1 RMSD 4.7-4.9 Å); 6X0J chain A passes at every size (about 1 Å);
  6X0J chain C **fails at 22 Å (5.96 Å, as in Step 4) and passes at 12 Å (1.06 Å, 4 of 5 seeds)**: the control gets easier as the box tightens around the deposited ligand.
* Oversized ligand (tetraphenylmethane): at 12 Å it is forced into the pocket (best-pose score +0.4 to +9 kcal/mol, 8-15 clashes, no clash-free pose); at 16-22 Å it docks at the rim with 0 clashes and a score better than the substrate's.
* Clash thresholds (0.4 Å overlap, 2.4 Å N/O) come from five deposited ornithine poses: a small calibration set.
