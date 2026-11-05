## Notes on the RosettaLigand redock (hand-written; numbers are in the generated table or `rosetta_dock/summary.json`)

Registered rule applied as written: per seed the lowest-total-score decoy is the top-1; PASS if the median over 10 seeds is within 2 Å. 100 decoys per receptor, none failed.

- **The protocol fails all five receptors, including the three that Vina passes.** Top-1 medians by total score are 13.35 (4B63), 5.69 (6X0J A), 13.57 (B), 9.62 (C)
  and 13.46 Å (D), with no seed within 2 Å and 0% near-native decoys in any receptor. Ranking by interface delta does not change any verdict.
- **So this protocol cannot arbitrate the 4B63 failure.** I wanted the passing receptors as positive controls: a method that fails them says nothing about whether
  Vina's 4B63 result is a scoring problem. The 4B63 failure here is not evidence about 4B63; every receptor fails the same way.
- **It is not far from the line in some receptors, but never over it.** The best single decoy is 2.09 to 2.83 Å from the deposited pose (4B63 2.16, A 2.53, B 2.83,
  C 2.09, D 2.45), and its rank by total score is 20, 2, 52, 3 and 2. Seed-wise, the smallest top-1 RMSD is 2.09 (C), 2.45 (D) and 2.53 Å (A); in 4B63 it is 4.88 Å.
- **A 2 Å rule is tight for this protocol.** Refining the deposited pose itself (native start, same HighResDocker and FinalMinimizer steps) moves it by a median
  0.65 to 1.50 Å from the crystal (4B63 1.01, A 1.24, B 1.50, C 0.65, D 1.35). So much of the budget is used before any search error.
- **Native-start versus docked scores (search or scoring?)** The refined deposited pose scores lower than the best docked decoy in 6X0J A, C and D (by 3.8, 6.7 and
  0.7 REU), which points to a search miss, and not lower in 4B63 and 6X0J B (the best docked decoy is ahead by 3.0 and 6.1 REU), which points to scoring. These
  margins are small against the spread of decoy scores (interquartile range 62.9 REU in 4B63 and 44.9 in 6X0J A, the two I checked), so I do not read either as a
  finding.
- **What I have not diagnosed.** Why the protocol misses the pocket in the median seed is untested. Candidates, none tested here: sampling (100 decoys per
  receptor, one start conformer per seed), the ref2015 treatment of a buried charged ligand, and the ligand-independent clashes seen in the earlier rescoring.
  Fixing any of these after seeing these numbers would be a new, labelled experiment, not a rescue of this one.
