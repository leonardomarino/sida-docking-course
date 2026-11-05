
## Notes (hand-written; numbers are in the tables above or in `rescoring/result.json` / `poses.csv`, untracked and regenerable with `python -m confirm.rescore`)

**Outcome, as registered.** Primary (AutoDock4): admissible, and **does not rescue 4B63** (median top-1 4.90 Å, 2/10 seeds within 2 Å). Secondary (Vinardo): admissible, does not rescue 4B63 (4.90 Å, 0/10). Baseline (Vina, interaction
energy only): 4.89 Å, 0/10, so the qualification sentence of Amendment 1 does not apply. 4B63 stays primary and fails; nothing about the Step 4 result changes.

**Events during the run.** (1) G0 as registered failed once and the rule was amended before any F1 / F2 score existed (Amendment 1, commit `24ad8df`); G0 as amended passed on all 1,000 poses (maximum difference 0.007 kcal/mol).
(2) The out-of-map rule triggered: 24 poses lay outside the 60-point map (4B63: 11, 6X0J_A: 1, 6X0J_B: 3, 6X0J_D: 9, 6X0J_C: 0), so all five jobs were mapped with 80 points per axis (30 Å), as registered; no pose is outside the final maps.
Nothing else departed from the rule.

**4B63, by seed.**
- Dropping the internal-strain term does nothing for 4B63: the interaction-only Vina baseline fails exactly as the docked ranking does (4.89 Å). The 0.4-0.5 kcal/mol gap in the decoy note is therefore not a strain effect at the level of the seed medians.
- AutoDock4 changes which wrong pose wins, not whether one does. Its top-1 is within 2 Å of Vina's own top pose in 4 of 10 seeds (4.86-4.93 Å), is a different pose at 5.06-5.17 Å in 4 seeds, and is near native in 2 seeds (0.99 and 1.78 Å, seeds 9 and 10). Vinardo reproduces Vina's decoy in all 10 seeds.
- The near-native pose is in the 20-pose list in 10 of 10 seeds for every function (the ceiling). Its rank under the rescoring (median, range): AutoDock4 3 (1-11), Vinardo 3 (1-8), Vina interaction-only 5 (3-15).
  In seed 1 the AutoDock4 energy of the best near-native pose (-2.15) is 0.54 kcal/mol behind the top decoy (-2.69), the same order as under Vina.

**Controls.** Vina and Vinardo keep all three at 10/10 seeds. AutoDock4 is admissible by the registered median rule but is the weakest of the three there: 6X0J_B passes with median 1.72 Å and only 5/10 seeds within 2 Å (the others 2.38-4.95 Å), 6X0J_D 7/10, 6X0J_A 10/10.
A 4B63 pass under AutoDock4 would therefore have been weak evidence; its failure is not weak, because a function that is noisy on the controls did not find the pose either.

**Chain C (reported, gates nothing).** Vina and Vina-interaction-only fail (5.96 Å, 0/10). AutoDock4 (1.11 Å) and Vinardo (1.12 Å) pass in 10/10 seeds. So chain C's Vina failure, a 0.04 kcal/mol near-tie in Step 4, is function-specific and goes away
under either alternative, while 4B63's does not.

**What this changes.** The ranking reframe is only partly supported. The near-native pose exists and ranks 1-11 under three functions of different form, but none puts it first on 4B63, while the same functions do fix chain C. The sampling is not the obstacle and the
scoring function alone (among these three) does not explain 4B63. Something in the 4B63 receptor, or in what the model leaves out (explicit waters, the oxidised flavin against the reduced one, the ornithine protonation, the cofactor state), makes the
deposited pose score worse than several alternatives under all of them. This note does not say which.

**Limits.** Interaction-only scores cannot see ligand strain. AutoDock4 uses the Meeko partial charges and a rigid receptor, and its treatment of the buried formal charge is untested. The 10 seeds of a receptor reuse nearly the same pose set, and the four 6X0J chains are
copies of one protein, so the effective number of independent receptors is two (4B63, 6X0J). Only functions registered in advance were run; no further function is permitted by the rule.
