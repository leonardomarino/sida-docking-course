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
