#!/bin/bash
# Relax each apo receptor, then dock (10 seeds x 10 decoys) and native-start refine. PyRosetta env; academic licence.
set -u
cd /opt/sida-docking-course
PY=/opt/miniforge3/envs/proteinmpnn_binder_design/bin/python
for J in $(ls rosetta_in | grep -v manifest); do
  echo "=== $J: relax"; env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. $PY -W ignore -m confirm.rosetta_redock relax --receptor rosetta_in/$J/receptor.pdb --out rosetta_in/$J/receptor_relaxed.pdb
  mkdir -p rosetta_dock/$J
  echo "=== $J: dock"; env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. $PY -W ignore -m confirm.rosetta_redock dock --receptor rosetta_in/$J/receptor_relaxed.pdb --starts rosetta_in/$J/starts/start_s*.pdb --mode dock --decoys-per-seed 10 --out-dir rosetta_dock/$J/dock --csv rosetta_dock/$J/dock.csv --jobs 12
  echo "=== $J: native"; env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. $PY -W ignore -m confirm.rosetta_redock dock --receptor rosetta_in/$J/receptor_relaxed.pdb --starts rosetta_in/$J/poses/crystal.pdb --mode native --decoys-per-seed 10 --out-dir rosetta_dock/$J/native --csv rosetta_dock/$J/native.csv --jobs 12
done
echo "=== all done"
