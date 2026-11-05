Hand-written. The tables above are generated; everything below is interpretation, and it is an observation, not a conclusion.
**Nothing here has been adopted as a fix, and no variant changes the registered verdict for the primary benchmark: 4B63 still
fails after every rung.**

**4B63 (primary) resists every single-factor change.** Top-1 median stays at 4.7-5.1 A through all ligand states, all boxes,
all cofactor edits and both receptor pH values. The only movement is in best-in-top-5, which drops to 1.08 A when FAD and
NADP are both removed (top-1 still 5.08 A). So the failure is not explained by ligand protonation, box, cofactor presence or
receptor pH. The baseline near-native pose (1.71 A) is 0.46 kcal/mol worse than the wrong top-1 and sits at mode 8 (seed 1).

**Ligand protonation is not the cause.** Neutralising either amine, or moving to the -1 state, leaves 4B63 (4.80-4.82 A) and
chain C (5.93-6.06 A) failing and the passing receptors passing (1.0-1.2 A). The +2 state makes everything worse, including the
three receptors that pass at baseline (4.5-6.1 A), so the baseline +1 zwitterion is the better-behaved state of the five tried.
Only one pH-consistent state per charge was tried; other tautomers were not.

**6X0J chain C is a near-tie, not a clear failure.** In seed 1 the near-native pose (1.11 A, mode 2) is only 0.04 kcal/mol
worse than the wrong top-1 (5.97 A), far inside Vina's error. Consistent with that, the verdict flips with small, physically
unremarkable changes: a box 26 or 30 A wide passes (1.06 and 1.05 A, 7/10 and 9/10 seeds), removing either FAD or NADP passes
(1.10 A, 10/10), and a box centre moved by 0.2-0.3 A (the carried box) moves it to 3.52 A (5/10 seeds). Its wrong pose keeps the
delta-amine end within 1.7 A of the deposited one while the carboxylate end is 8-10 A away, so it is not a simple head-to-tail
flip. That it is rescued by removing either cofactor suggests the wrong pose depends on contacts with both; this was not
examined at atom level.

**The 4B63 top-1 pose is displaced along the whole ligand.** Per-atom distances to the deposited pose (both engines): N 5.5,
CA 4.2, CB 2.6, CG 1.5, CD 3.3, NE 4.4, C 5.8, O 5.8, OXT 7.2 A. The passing receptors show 0.2-1.9 A per atom.

**Box findings, with their limits.** The 18 A box contains only 56% of the deposited ligand in 4B63, so that cell says
nothing about 4B63. The oracle box (centred on the deposited ligand) does not rescue 4B63 (4.73 A) or chain C (5.99 A), so
centring is not the problem. The 26 and 30 A boxes pass chain C but not 4B63. A larger box was **not** adopted: the 22 A
size was a starting value, and choosing a size because it makes one receptor pass is the tuning the step forbids. A box rule
fixed in advance (for example the ligand's extent plus a stated margin) would have to be defined first and then all
receptors re-run.

**Two null results that are expected, not informative.** Zeroing FAD's partial charges changes nothing (identical to the
baseline in every receptor), as it should: Vina's scoring function has no electrostatic term, so charges are not used. Receptor
pH 6.0 and 8.5 changed none of the three residues in the Step 2 site table (His91, Lys100, Arg144 in 4B63), but that table does not list Lys107 (and Asp288 in 6X0J), which contact the deposited ornithine and lie 9.4-11.7 A from C4X, outside its 8 A radius, so their states were not examined;
PROPKA's histidine pKa is unreliable (see the Step 2 note), and other His91 tautomers were not tried.

**Not tested.** All crystal waters are stripped (262 in 4B63, the nearest 3.2 A from C4X) in a site that buries a +1 charge, and waters were never tried; Lys107 and Asp288 protonation (see above); Lys100 (modelled, not observed, in 4B63 and 6X0J) in alternative
positions; flexible side chains; other protonation states of the receptor beyond the pH values above; a different scoring
function that includes electrostatics (the Rosetta check of this step was inconclusive); and the large 4B63 cavity (3,161 A^3 against 2,059-2,414 A^3 in 6X0J, fpocket volumes of the run on disk) as a cause.

**What the ladder establishes, stated once.** The failure of 6X0J chain C is resolved by removing either cofactor or by a box of 26-30 A (and
flips with a 0.2 A shift of the box), consistent with a near-tie between two poses. The failure of 4B63, the primary benchmark, survives
every one of the four factors the spec prescribes: no variant brings its top-1 below 4.73 A. The 4B63 and 6X0J ligand positions are the
same (superposed, the deposited ornithines lie 0.34 A apart), so the difference lies in the prepared receptor or box, not in the biology.

**Limits of the ladder.** It runs Vina only; on chain C the engines did not agree at baseline (Uni-Dock found the near-native pose in
3/10 seeds, Vina in 0/10). The 18 A box cell for 4B63 is excluded in the table: that box holds 56% of the deposited ligand, so the
near-native pose cannot be reached. The 'FAD charges zeroed' row is a control, not evidence: Vina ignores partial charges, so it is
identical to the baseline by construction. Rung 1 used hand-written protonation states, not what the pipeline's pH model would give at
other pH values.
