Hand-written. The table above is generated; this is interpretation, and it is an observation, not a conclusion.

**Outcome under the registered rule: the primary benchmark FAILS.** 4B63 (primary): top-1 median 4.80 A (Vina) and 4.82 A
(Uni-Dock), 0 of 10 seeds within 2 A with either engine. 6X0J (secondary) passes in chains A, B and D with both engines
(0.90-1.05 A, 10/10 seeds) and fails in chain C (5.96 A, both engines). Nothing has been changed to improve these numbers,
and the primary / secondary roles are not swapped.

**It is not search noise.** The failing top-1 pose is the same in every seed (4B63: 4.79-4.83 A across 10 seeds in Vina,
4.79-4.85 A in Uni-Dock) and the two engines agree to within 0.05 A. The search finds one stable answer; the question is
why it is not the deposited pose.

**A near-native pose exists in the list but is not ranked first.** Best-in-top-5 median: 4B63 Vina 1.75 A; chain C 1.11 A
(Vina) and 1.09 A (Uni-Dock). For 4B63 Uni-Dock has nothing near-native in its top 5 (4.79 A). In chain C Uni-Dock's top
pose is near-native in 3 of 10 seeds (3, 4, 10), so there the choice between two poses is close.

**The top-1 scores differ little between a pass and a fail.** Median top-1 score: -4.91 kcal/mol (4B63) against -5.17 to
-5.24 (6X0J passes) and -5.18 to -5.22 (6X0J chain C, a fail).

**Rosetta does not settle it** (`step4_rosetta.md`): it neither endorses the passing redocks nor prefers the near-native pose
in the failing ones.

**Not yet done, and required by the step.** The prescribed diagnosis, in this order: ligand protonation, box centre and size,
FAD treatment, receptor protonation, each change reported with its effect. Candidates I have not tested and do not claim:
a head-to-tail orientation difference of the ornithine (its two amines are chemically similar), the receptor protonation
(PROPKA's histidine pKa is unreliable; Lys100 is modelled, not observed, in 4B63 and 6X0J), and 4B63's much larger fpocket
cavity (3,325 A^3 against 2,100-2,400 A^3 in 6X0J).
