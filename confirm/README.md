# confirm/ - optional Rosetta check of the docking poses

PyRosetta is free for academic use (students may use it freely) and needs a licence for commercial use. Nothing in `sidacore/`
imports this directory, so the shared layer stays free of non-commercial components; delete `confirm/` and nothing else changes.

Three stages. The middle one runs in an environment that has PyRosetta (here `proteinmpnn_binder_design`); the others run in
`sida-docking`.

```bash
# 1. inputs: receptor (protein + cofactors, heavy atoms), the deposited ligand and every mode of one seed, atom-named, chain L
python -m confirm.make_inputs --config config/default.yaml --prepared prepared --redock redock --out rosetta_in

# 2. Rosetta ref2015: minimise ligand + shell side chains, then interface analysis (PyRosetta env)
for J in $(ls rosetta_in | grep -v manifest); do
  python -m confirm.rosetta_rescore --receptor rosetta_in/$J/receptor.pdb --ligands rosetta_in/$J/poses/*.pdb \
         --out-dir rosetta_out/$J --csv rosetta_out/$J.csv --jobs 12
done

# 3. summary
python -m confirm.summarise --config config/default.yaml --prepared prepared --inputs rosetta_in --scores rosetta_out \
       --json rosetta_out/summary.json --report reports/step4_rosetta.md --notes reports/step4_rosetta_notes.md
```

Tests: `python -m pytest tests/test_confirm_*.py` (sida-docking) and `python -m pytest confirm/tests` (PyRosetta env).

## Step 4b: rescoring the Vina poses (`rescore.py`, no PyRosetta)

Rule: `reports/step4_rescoring_rule.md` (registered in `92d2d0b`, Amendment 1 in `24ad8df`). Scores the existing Vina poses as docked, interaction energy only, with Vina (baseline), AutoDock4 (primary) and
Vinardo (secondary). `vina` comes from `sida-docking`; AutoGrid 4.2.9 (GPL-2.0-only) from `sida-docking-ad4`, run as a separate program on files. Nothing in `sidacore/` imports this.

```bash
python -m confirm.rescore --course . --out rescoring --jobs 12 --report reports/step4_rescoring.md --notes reports/step4_rescoring_notes.md
python -m confirm.rescore --g0-only --out rescoring --course .        # gate G0 only; stops there
python -m confirm.rescore --out rescoring --render-only --report reports/step4_rescoring.md --notes reports/step4_rescoring_notes.md
```

Run from `sida-docking` (`PYTHONPATH=.`); about 90 s on the CPU. Tests: `python -m pytest tests/test_confirm_rescore.py`. `rescoring/` is untracked output.
