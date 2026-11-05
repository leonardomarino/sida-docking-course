# Licence survey: independent rescoring of the existing Vina poses (hand-written; 2026-11-05)

Purpose: Step 4 left a ranking problem (a near-native 4B63 pose sits 0.4-0.5 kcal/mol behind a stable decoy; see `step4_decoy_characterization.md`). The only experiment aimed at that is to rescore the poses already on
disk with a function other than Vina's. Before anything touches `sidacore`, this records what each candidate costs in licence terms. **Nothing was rescored for a result; no code was added.**

**Policy this is read against** (`provisioning_record.md`): no non-commercial licence and no restricted-training-data weights in any course environment; copyleft (OpenBabel GPL-2.0-only, Meeko LGPL-2.1-only,
OpenMM LGPL-3.0-or-later) is accepted because it matters only if redistributed; `sidacore` never imports PyRosetta. The same rule is applied here: a new tool is run as a separate executable on files, never imported into `sidacore`.

How each licence was established: installed packages from `conda-meta/*.json`; upstream repositories from the GitHub licence API (`gh api repos/<owner>/<repo>/license`) and, for gnina, its README; AutoGrid from the
conda-forge package metadata. Anything not read from one of those is marked *not verified*.

| Candidate | Differs from Vina how? | Licence (source) | Installed? | Verdict for this use |
|---|---|---|---|---|
| **Vina scoring, `--score_only`** (Vina 1.2.7, Uni-Dock 1.2.0) | none: sanity baseline | Apache-2.0 (conda-meta) | yes, both envs | **clean**. Reproduces the docked score (smoke test on the 4B63 decoy: -4.924 vs -4.922 in the docking log) |
| **Vinardo** (`--scoring vinardo` in the same two programs) | re-fitted weights and terms, same family; no electrostatics | Apache-2.0 (same binaries; smina upstream Apache-2.0, GitHub API) | yes | **clean, runs today**; weakly independent of Vina |
| **AutoDock4 scoring** (`vina --scoring ad4` + maps from AutoGrid 4.2.9) | **Coulomb electrostatics, desolvation, directional H-bonds**: the property the external review asked for | Vina Apache-2.0; AutoGrid **GPL-2.0-only** (conda-forge metadata; AutoDock4 upstream GPL-2.0) | Vina yes; AutoGrid not installed (conda-forge `autogrid` 4.2.9 exists) | **clean under the stated policy** if AutoGrid is run as a separate program and only its map files are used; the same GPL-2.0 status as OpenBabel already in the env. Vina refuses `--receptor` with ad4 (smoke test: "No receptor allowed"), so maps are required |
| **gnina** (CNN rescoring) | learned, not physics | Code **dual GPL / Apache-2.0** (README: GPL because of OpenBabel); `gnina/models` repo Apache-2.0 (GitHub API). Training-data terms (CrossDocked2020, and PDBbind if used for the affinity heads) **not verified** | no (needs its own CUDA build; sm_120 support not checked) | **do not adopt until the weights' training data are checked.** Also a leakage question: PDB-trained networks may have seen SidA homologs |
| **ODDT** + RF-Score / NNScore | learned | ODDT BSD-3-Clause (GitHub API); pretrained models trained on PDBbind, whose terms are **not verified**; NNScore licence **not verified** | no | hold: same data-terms question as gnina |
| **OpenMM MM-GBSA** | physics, implicit-solvent electrostatics | OpenMM installed as LGPL-3.0-or-later (conda-meta); OpenFF toolkit MIT (GitHub API); FAD / NADP force-field parameters **not checked** | OpenMM yes in both envs; no parameters for FAD, NADP or the ornithine zwitterion | clean in principle, but days of setup, and MM-GBSA on unrelaxed docking poses is noisy. Not first |
| RosettaLigand / ref2015 | physics + knowledge terms | academic-free, commercial licence required (user statement, this project) | PyRosetta in other envs, not here | **excluded** from `sidacore` by the standing rule; already used as a Step 4 side check |
| GOLD, Glide, MOE, ICM, PLANTS, DOCK6 | various | commercial or academic-only (general knowledge, **not re-verified**) | no | excluded |
| PoseBusters | not a score: physical-plausibility checks | MIT (GitHub API) | no | optional aside, can flag a decoy as unphysical; does not rank |

## Practical findings from the smoke tests

- `vina --scoring {vina,vinardo} --score_only --autobox` works on the existing receptor PDBQT and a single-model ligand PDBQT; a multi-MODEL pose file must be split first ("Unexpected multi-MODEL tag"). Uni-Dock 1.2.0 exposes the same three functions.
- The receptor PDBQT from Meeko already carries AD4 atom types (C, A, OA, N, NA, SA, HD, P) and partial charges, so it can feed AutoGrid without re-preparation. How well AD4 charges treat the buried formal charge of the ornithine
  zwitterion is a modelling caveat, not a licence one, and is untested.
- The smoke runs scored one pose and are not results; nothing from them was recorded as a finding.

## Recommendation (for the user to confirm; nothing registered yet)

1. **Primary: AD4 rescoring**, because it is the only licence-clean candidate with electrostatics and desolvation; **secondary: Vinardo**, which costs nothing. Fix this choice now rather than trying functions until one works.
2. **Register the rule before running**, as for Steps 3-4: rescore the top-20 poses of every seed for 4B63 and for the 6X0J chains; report the rescored top-1 median RMSD and the rank of the best-RMSD pose. The three Vina-passing 6X0J chains (A, B, D) are the controls: a function that
   breaks them is disqualified, whatever it does to 4B63. A pass for 4B63 would be a finding about the scoring function, never a fix, and 4B63 stays primary.
3. **Ceiling**: rescoring can only reorder the 20 poses Vina kept. Seed 1 has two near-native poses (ranks 5 and 8), so the ceiling is not the obstacle there.
4. Implementation if confirmed: a small script under `confirm/` (the PyRosetta precedent), calling `vina` and `autogrid4` as subprocesses; AutoGrid installed into `sida-docking` or a separate env, which is a provisioning change to be approved separately.
5. Defer gnina / ODDT until the training-data terms are read; defer MM-GBSA until a force-field source for FAD / NADP is chosen.
