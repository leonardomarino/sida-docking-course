# Step 4b: rescoring rule (DRAFT, to be registered before any pose is rescored)

**Status: draft for the user's approval.** It becomes the registered rule when the user approves it and it is committed; nothing below has been run. If anything is changed after the first rescoring run, the run
and the change are both reported, as in Steps 3-4.

## Question

Step 4 left a ranking problem: for 4B63 the pose lists contain a near-native pose a few tenths of a kcal/mol behind a stable decoy (`step4_decoy_characterization.md`). **Does a scoring function other than Vina's put a
near-native pose first among the poses Vina already produced, without breaking the receptors that Vina passes?** This is a question about the scoring function. It is not a docking run and it cannot change the Step 4 result.

## Disclosure of what was seen before writing this rule

- Poses and scores from the Step 4 Vina runs, including the decoy characterization (seed-1 poses 5 and 8 on 4B63 are near native, 0.43 and 0.47 kcal/mol behind the top pose).
- One smoke test: the 4B63 seed-1 top pose scored once with `--scoring vina` (-4.924, equal to its docking-log -4.922) and once with `--scoring vinardo` (-3.645). No other pose, no other function, and no AD4 score has been computed.
  Vinardo's absolute value is not comparable to Vina's and says nothing about ranking.

## Inputs (fixed)

| Job | Role | Receptor | Poses |
|---|---|---|---|
| 4B63 chain A, `4B63_A_nadp-keep` | primary benchmark; Vina FAILS (4.80 Å) | `prepared/<job>/receptor.pdbqt` | `redock/<job>/vina/vina_seed{1..10}.pdbqt`, 20 poses each (200) |
| 6X0J chains A, B, D, `6X0J_{A,B,D}_nadp-keep` | **controls**; Vina passes (~1.0 Å) | same | same (200 each) |
| 6X0J chain C, `6X0J_C_nadp-keep` | Vina fails (5.96 Å, 0.04 kcal/mol near-tie); reported, not gating | same | same (200) |

6X0H and 6X0I are excluded (no ornithine). Chimeric receptors are not rescored. Reference ligand for every RMSD: the deposited ornithine of that job (`prepared/<job>/ligand_crystal.pdb`).
Ligand atom types in all pose files: C, HD, N, OA. Receptor atom types: A, C, HD, N, NA, OA, P, SA. All are standard AutoDock4 types.

## Functions (exactly these three; no others may be added after the first run)

- **F0, Vina** (`vina --scoring vina --score_only`): reproducibility gate only.
- **F1, AutoDock4** (`vina --scoring ad4 --score_only` with AutoGrid 4.2.9 maps): **PRIMARY.** Electrostatics, desolvation, directional hydrogen bonds.
- **F2, Vinardo** (`vina --scoring vinardo --score_only`): **secondary**, same functional family as Vina.

## Procedure (fixed)

1. Split each seed file into single-model PDBQT files. Poses are scored **as docked**: no minimization, no re-docking, no flexibility.
2. F0 and F2: `vina --receptor prepared/<job>/receptor.pdbqt --ligand <pose> --scoring <f> --score_only --autobox`.
3. F1: one AutoGrid run per job from a GPF written by the tool with these values and no others: receptor types `A C HD N NA OA P SA`; ligand types `C HD N OA`; `spacing 0.375`; `npts 60 60 60` (22.5 Å);
   `gridcenter` = the box centre in `prepared/<job>/receptor.box.txt`; `smooth 0.5`; `dielectric -0.1465`; receptor charges exactly as in `receptor.pdbqt` (not edited). Scoring: `vina --maps <prefix> --ligand <pose> --scoring ad4 --score_only`.
   **Out-of-map rule, fixed now:** if any atom of any scored pose lies outside the 22.5 Å map of its job, all five jobs are rerun with `npts 80 80 80` (30 Å) and both outcomes are reported. No other grid change is allowed.
4. Score = the program's "Estimated Free Energy of Binding". Within one seed's 20 poses the lowest score ranks first. Two poses within 0.05 kcal/mol at the top are a **tie**: the headline uses the pose that Vina ranked higher, and the number of tied
   seeds is reported.
5. RMSD = `sidacore.redock.pose_rmsd` (symmetry-corrected, no re-alignment) against the deposited ornithine, exactly as in Step 4. **Per-job result = median over the 10 seeds of the rescored top-1 RMSD; pass = median <= 2 Å**, the Step 4 rule.

## Gates, in this order

- **G0, pipeline reproduces Vina.** F0 must reproduce the docking-log score of every one of the 1,000 poses to within 0.05 kcal/mol (the smoke test differed by 0.002). If not, stop: the splitting, receptor or ligand handling is wrong and no F1 / F2 number is read.
- **G1, admissibility (per function).** A function is **admissible** only if, after rescoring, the median top-1 RMSD stays <= 2 Å on **all three controls** (6X0J A, B, D). A function that breaks any control is disqualified, whatever it does to 4B63; its 4B63 numbers are still printed and
  labelled "not interpretable".

## Outcomes, fixed in advance (primary = F1 on 4B63)

| F1 admissible? | F1 4B63 median top-1 | Statement that will be written |
|---|---|---|
| yes | <= 2 Å | AutoDock4 scoring ranks a near-native pose first on 4B63 among Vina's poses and keeps the controls. A finding about scoring functions, not a fix: 4B63 stays primary and its Vina failure stays in the report. |
| yes | > 2 Å | AutoDock4 scoring does not rescue 4B63 either. |
| no | any | AutoDock4 scoring cannot arbitrate this site: it breaks control(s) <names>. |

F2 is filled in the same way and labelled secondary. If F1 is inadmissible or fails and F2 passes, the F2 result is reported with the note that a second function was tried, so a pass is less surprising by chance.
Chain C is reported for every function and gates nothing.

## Also reported for every job and function (descriptive, fixed now)

- the rank, under the rescoring, of the best-RMSD pose in each seed, and the **ceiling**: in how many of the 10 seeds any of the 20 poses is <= 2 Å;
- for 4B63, whether the rescored top-1 pose is the Step 4 decoy (<= 2 Å from the Vina decoy) or a different pose;
- the number of out-of-map poses and tied seeds.

## Not allowed after the first run

Adding a function, changing grid size, spacing, charges or atom types (other than the single fixed `npts 80` fallback), dropping poses or seeds, excluding a control, changing the 2 Å criterion, or swapping primary and secondary receptors.

## What no outcome would show

- The 20 poses were chosen by Vina. A pass shows a better *ranking*, not that the function could *find* the pose; the ceiling above bounds what is possible.
- The 10 seeds of one receptor return essentially the same pose set (4B63 top poses agree to 0.12 Å), and the four 6X0J chains are copies of one protein, so the effective number of independent receptors is small. Seed agreement is not independent evidence.
- AutoDock4 uses the partial charges in the PDBQT (Meeko), rigid receptor, no explicit waters; how it treats the buried formal charge of the ornithine is not tested here.
- Nothing here bears on SidA biology, on the protonation state of the δ-amine, or on whether the 4B63 or 6X0J conformation is the productive one.

## Implementation (not started; proposed)

`confirm/rescore.py` in the course repo, the precedent being the PyRosetta checks: `vina` from `sida-docking` and `autogrid4` from `sida-docking-ad4` are called as subprocesses on files; `sidacore` imports neither and gains no dependency.
Tests first (pose splitting, GPF writer, score parser, tie rule, the gate logic on synthetic inputs, G0 on one real job), then the run, then a generated `reports/step4_rescoring.md`. Cost: 2,000 `--score_only` calls and five AutoGrid runs, minutes on the CPU.
