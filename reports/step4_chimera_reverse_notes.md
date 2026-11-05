## Notes on the reverse chimeras (hand-written; numbers are in the generated table or `chimera_reverse/step4_chimera.json`)

Design and reading were registered before the run (see above): each passing 6X0J chain (A, B, D) receives 4B63's side-chain conformation of the named residues as a rotamer transplant onto its own backbone;
a swap that makes a passing chain fail, while the control does not, means those residues matter for the pass. Same pass rule (median top-1 within 2 A over 10 seeds, Vina).

- **The method reproduces the original passes exactly.** The round-trip baselines give 1.00 (A), 1.04 (B) and 0.90 A (D), 10/10 seeds each, the same medians as the original redock.
- **Neither the control nor Asp288 changes anything.** The control (four contact residues taken from 4B63) gives 1.03, 1.08 and 0.94 A, all 10/10; Asp288 alone gives 1.01, 1.04 and 0.81 A, all 10/10.
- **4B63's Gln102 conformation breaks two of the three passing chains.** Chain A fails at 5.07 A (0/10 seeds; per-seed range 5.02-5.10) and chain B at 4.83 A (0/10; 4.79-4.85). Chain D passes with Gln102 alone
  (1.00 A, 8/10 seeds) but fails with Gln102 + Asp288 (5.11 A, 3/10, per-seed range 1.00-5.43). Adding Asp288 to Gln102 changes chain A and B by at most 0.02 A and moves only chain D.
- **It is a ranking change, not a loss of the right pose.** With Gln102 swapped, best-in-top-5 stays at 0.98-1.05 A, and the best pose in the list sits at median rank 3-5.5 (chains A and B; baseline rank 1). The wrong
  pose scores -5.11 (A) and -5.08 kcal/mol (B), against -5.17 and -5.23 for the correct top pose of the baselines: a difference of 0.06-0.15 kcal/mol.
- **The failing distance matches 4B63's.** The wrong top pose is 4.8-5.1 A from the deposited ornithine, as in the 4B63 redock (4.80 A). I compared the distances only, not the poses, so I do not claim it is the same pose.
- **Reading with the forward run.** In the forward run, 6X0J's Gln102 did not fix 4B63 (4.79 A), and it made the near-native pose rank worse. Together: 4B63's Gln102 conformation is enough to break a passing 6X0J chain
  in two of three, but 6X0J's is not enough to fix 4B63, so at least one more difference in 4B63 is involved.
- **Important confound: the Gln102 chimeras are not clean receptors.** 4B63's Gln102 conformation puts its amide oxygen 1.65 A (A), 1.91 A (B) and 1.83 A (D) from the NADP 2'-hydroxyl oxygen of the 6X0J chain, an
  overlap that cannot exist in a real structure. Vina does not score receptor-internal overlaps, but the swapped side chain is in a position this NADP does not allow, so I cannot separate "this conformation changes
  the ranking" from "an impossible receptor changes the ranking". Untested.
- **The all-four-differing variant could not be built in any chain.** 4B63's Arg144 conformation overlaps an NADP nitrogen by 0.92-1.30 A, and Meeko refused all three (excess inter-residue bonds at residue 144 or 179 and at NADP 601).
- **Tests that would settle the confound (not run):** the Gln102 swap with NADP removed from the receptor (if it still fails, the conformation matters on its own); the Gln102 swap together with 4B63's NADP position;
  Arg144 and Arg328 alone. Removing NADP is the treatment that rescued 6X0J chain C in the diagnosis ladder, so it is a natural pairing.
