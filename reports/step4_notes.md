Hand-written. The table above is generated; this is interpretation, and it is an observation, not a conclusion.

**Outcome under the registered rule: the primary benchmark FAILS.** 4B63 (primary): top-1 median 4.80 A (Vina) and 4.82 A
(Uni-Dock), 0 of 10 seeds within 2 A with either engine. 6X0J (secondary) passes in chains A, B and D with both engines
(0.90-1.05 A, 10/10 seeds) and fails in chain C (5.96 A, both engines). Nothing has been changed to improve these numbers,
and the primary / secondary roles are not swapped.

**It is not search noise.** The failing top-1 pose is the same in every seed (4B63: 4.79-4.83 A across 10 seeds in Vina,
4.79-4.85 A in Uni-Dock) and the two engines agree to within 0.05 A. The search finds one stable answer; the question is
why it is not the deposited pose.

**A near-native pose exists in the list but is not ranked first.** Best-in-top-5 median: 4B63 Vina 1.75 A; chain C 1.11 A
(Vina) and 1.09 A (Uni-Dock). For 4B63 Uni-Dock has nothing near-native in its top 5 (4.79 A), but see the pose-count note below. In chain C Uni-Dock's top
pose is near-native in 3 of 10 seeds (3, 4, 10), so there the choice between two poses is close.

**The top-1 scores differ little between a pass and a fail.** Median top-1 score: -4.91 kcal/mol (4B63) against -5.17 to
-5.24 (6X0J passes) and -5.18 to -5.22 (6X0J chain C, a fail).

**Rosetta does not settle it** (`step4_rosetta.md`): it neither endorses the passing redocks nor prefers the near-native pose
in the failing ones.

**Diagnosis done** (`step4_diagnosis.md`, notes in `step4_diagnosis_notes.md`). In the prescribed order (ligand protonation, box,
FAD treatment, receptor protonation) no single-factor change rescues 4B63, so the primary benchmark still fails. 6X0J chain C is
a near-tie (0.04 kcal/mol) that flips with a larger box or with either cofactor removed. Nothing has been adopted as a fix.

**Pose counts differ between engines.** Vina returned exactly 20 poses every seed; Uni-Dock returned 5-12 although `num_modes` was 20
(its own energy-range cut-off applies). Uni-Dock's 'top 5' is therefore often most of its list, and best-in-top-5 comparisons between the
engines are not like for like. The top-1 results, which the pass rule uses, are unaffected.
