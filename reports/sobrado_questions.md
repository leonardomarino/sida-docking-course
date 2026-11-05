# SidA ornithine docking: what we see, and what we would like your view on

*Draft from Leonardo Mariño-Ramírez to Pablo Sobrado. For review before sending. First written 4 November 2026; updated 5 November 2026 with the follow-up checks.*

## Why I am writing

I am building a docking benchmark and teaching course around SidA, and I would like the protocol to be right in a way you would accept. You confirmed that the ornithine site is the one to target,
so everything below concerns it. Before using any structure to compare ligands, I redocked the deposited L-ornithine into the structures that contain it. I fixed the pass rule beforehand: the median
top-ranked pose over 10 seeds within 2 Å (symmetry-corrected heavy-atom RMSD) of the deposited ligand. Receptors were prepared with FAD and NADP(H) kept, waters removed,
protein hydrogens from PDB2PQR/PROPKA at pH 7.4. Two engines were used, AutoDock Vina 1.2.7 and Uni-Dock, and their top-ranked poses agree.

## What we see

| Structure | Result | Note |
|---|---|---|
| 6X0J chains A, B, D | pass (0.9-1.1 Å, all seeds) | The deposited pose is recovered at rank 1 |
| 6X0J chain C | **fail with Vina** (top pose 6.0 Å) | A pose at 1.1 Å is in the list at rank 2-4, and in seed 1 the scores differ by only 0.04 kcal/mol. AutoDock4 and Vinardo scoring put the correct pose first (1.1 Å, all 10 seeds). It also passes in a 12 Å box centred on the ornithine (1.06 Å, 4 of 5 seeds) |
| 4B63 (1.9 Å, oxidized FAD, NADP+) | **fail** (top pose 4.8 Å) | It fails with Vina, AutoDock4 and Vinardo scoring (median 4.9 Å under each), and in every box size I tried (12 to 22 Å, 3-5 seeds). RosettaLigand failed on all five receptors, including the three that Vina passes, so it cannot arbitrate |

The deposited ornithine makes the same four protein contacts in every structure: Lys107 NZ to the carboxylate (2.8-3.0 Å), Ser469, Asn293 and Asn323. The wrong top
poses lose Lys107, Asn293 and Ser469 together (Lys107 NZ 8-11 Å from the carboxylate). I then varied one factor at a time: ornithine protonation (including the
neutral delta-amine you argue for), box size and position, removing FAD and/or NADP, and receptor pH. Removing either cofactor, or enlarging the box, rescues
chain C. **Nothing rescues 4B63.** I have not tuned anything to make it pass, and I do not know why 4B63 differs from 6X0J.

**What the wrong 4B63 pose is.** It is one stable pose: all 10 seeds give the same top pose (scores -4.90 to -4.93 kcal/mol). Its carboxylate is 5.8 Å from the deposited one, against Gln102, Asn323 and Leu467,
and 8.9 Å from Lys107 NZ. A near-native pose is in the list (carboxylate 1.4-1.5 Å from the deposited one, against Ser469, Phe296 and Asn293 like the crystal pose), 0.43-0.47 kcal/mol behind the
wrong one. That gap is small next to the usual error of a docking score, so I read this as a ranking failure, not a sampling failure.

**What I checked since, and what did not explain it.**

- *Crystallization additives.* 4B63 has one sulfate and four glycerols that 6X0J lacks. I stripped them. They are far from the ornithine (the sulfate 14.5 Å from its carboxylate, the glycerols 26-29 Å), and
  the wrong pose is not in an anion site. Not the cause.
- *Crystal packing.* 4B63 is I 2 2 2 (one chain); 6X0J is P 1 21 1 (four chains). Gln102, Arg144 and Arg328 have no lattice contact within 4.5 Å in either, and Asp288 touches only a water. Not the cause.
- *Lys107.* I had suspected a Lys107 rearrangement. Its NZ is within 0.1-0.6 Å between 4B63 and the 6X0J chains (chi1 -55 to -60° in all), and it is 2.76-2.98 Å from the carboxylate in all of them. I no longer suspect it.
- *Alternative scoring.* Re-scoring the existing Vina poses with AutoDock4 (with electrostatics) and Vinardo did not rescue 4B63 (4.90 Å for both; Vina's interaction energy alone gives 4.89 Å). AutoDock4 was noisier on the
  controls (6X0J chains B and D pass in only 5 and 7 of 10 seeds).

**What is different between 4B63 and 6X0J.** Four of the 84 residues of the site (the docking box plus everything within 4.5 Å of the ornithine) have side chains that differ by more than the
chains of 6X0J differ among themselves (side-chain heavy-atom RMSD to 4B63, similar for all four 6X0J chains):

| Residue | 4B63 vs 6X0J (Å) | largest difference between two 6X0J chains (Å) | touches the ornithine |
|---|---|---|---|
| Gln102 | 1.9-2.1 | 0.35 | yes |
| Arg144 | 2.7 | 0.24 | no |
| Asp288 | 1.7-1.8 | 0.21 | no |
| Arg328 | 3.5-3.7 | 0.78 | no |

Putting 4B63's Gln102 conformation into 6X0J chains A and B makes both fail (top pose 5.07 and 4.83 Å; chain D fails only when Asp288 is swapped as well; Asp288 alone does nothing). I cannot read this
cleanly, because that conformation puts Gln102 OE1 1.65-1.91 Å from NADP O2D, which a real structure would not. Putting 6X0J's Gln102 and Asp288 into 4B63 does not rescue it. The pose that the Gln102 transplant produces is
not the 4B63 decoy (3.7 Å apart) and resembles the decoy of unmodified chain C instead, so I cannot say that Gln102 explains 4B63. Arg328 CZ is 4.4 Å from the 4B63 sulfate, so its difference may be an
additive effect (see question 3).

## Questions (★ = matters most)

**About the structures**

1. ★ **4B63 versus 6X0J.** Is there anything about 4B63 that I should know: crystallization pH and conditions, why one chain, how confident you are in the ornithine and Lys107 positions?
   Chain C of 6X0J now looks like a near-tie between two poses, not a structural outlier. Do you know of anything that would make it behave differently from A, B and D?
2. ★ **Gln102, Arg144, Asp288 and Arg328.** These side chains sit 1.7-3.7 Å apart in 4B63 and 6X0J, and none is at a crystal contact. Your 2019 Protein Science paper reports that Gln102 interacts with NADPH
   and that Q102A destabilises the C4a-hydroperoxyflavin. The Gln102 side chain differs by about 2 Å between 4B63 and 6X0J, and if I put the 4B63 conformation into a 6X0J chain it collides with the NADP ribose
   (Gln102 OE1 1.65-1.91 Å from NADP O2D). Is the 4B63 conformation of Gln102 (and of the other three) well supported by the density, or is it a modelling choice? Do you see Gln102 in more than one conformation,
   in structures, simulations or from the mutant data?
3. ★ **The 4B63 sulfate.** SO4 1495 has occupancy 1.0 and B of 55-60 Å^2, and sits 4.4 Å from Arg328. Is it real, or a buffer or cryo artifact?
4. **Lys107.** In 4B63 its NZ is 2.9 Å from the carboxylate and it does not move between structures, so I no longer expect it to explain the failure. Do K107 mutant data (Km for ornithine) exist that I can cite?
   The closest-to-crystal 4B63 poses (1.7-2.1 Å) place Lys107 NZ 3.7-4.3 Å from the carboxylate.
5. ★ **Protonation to use.** Your group argues the ornithine delta-amine binds neutral (Protein Sci, "Flavin oxidation in flavin-dependent N-monooxygenases"). Which
   computational studies does that paper cite as refs 7 and 8? For docking, which state would you defend for the delta-amine, the alpha-amine and Lys107? My result did
   not depend on it, but I would like to use the state you consider correct.
6. ★ **Which cofactors belong in the receptor.** Campbell et al. removed NADP(H); the MD of Pierdominici-Sottile et al. (2024) suggests NADPH blocks ornithine in the FAD-out state. For docking into the
   ornithine site, would you keep NADP(H), remove it, or run both? Is the NADP contact to the delta-N in 4B63 (3.3 Å) catalytically meaningful?
7. **The decoy position.** The wrong 4B63 pose and its relatives sit against Gln102, Thr322, Ser286, Asn323 and Leu467, 6-7 Å from the deposited carboxylate. Does that region mean anything
   functionally (substrate entry, a second binding mode, something seen in your simulations), or is it just a pocket that scoring likes?

**About details I handled by default**

8. **Waters.** I remove all of them (262 in 4B63; one is 3.2 Å from FAD C4X). Are any waters in the ornithine site structural?
9. **FAD state.** 4B63 has oxidized FAD; 6X0J has reduced flavin. Do you expect ornithine to bind the same way in both, and which is the relevant state?
10. **FAD-out.** I plan to compare 6X0H (FAD-out, no NADP) with 6X0I (FAD-in, NADP). Because no deposited structure has FAD-out with NADP, I will report that the
    comparison is confounded and not attribute differences to FAD. Do you agree, and is the M101A structure (7JVL) or another structure a better FAD-out reference?

**Data and collaboration**

11. **Ranking test.** On rigid receptors Vina scores arginine (-6.2 to -6.8 kcal/mol) better than lysine (-5.4 to -5.6), and lysine better than ornithine (-4.9 to -5.2), in 4B63 and in 6X0J chains A and C
    (3 seeds, 12 Å box). I treat these scores as unreliable, since a large hydrophobic molecule that does not belong in the site also out-scores ornithine when the box lets it sit at the rim. Do binding data support
    any ordering? I understand the 2012 Franceschini paper has lysine and arginine complexes. Can you point me to the Km/Kd values, and to any inhibitors with a known binding site (for example celastrol) that I may use as a test set?
12. **Mutant structures.** For 5CKU (N323A) and 4NZH (R279A): is ornithine bound in the same pose, and are kinetics available? I would use them as a specificity exercise.
13. **Your own docking.** Have you or your group docked ornithine into 4B63 or 6X0J with any program? If so, would you share the protocol and any redocking numbers?
14. **Review.** Would you be willing to read the protocol section, and tell me where an enzymologist would disagree?

*All numbers above come from runs in the course repository: `reports/step4_redock.md`, `step4_contacts.md`, `step4_site_comparison.md`, `step4_chimera.md`, `step4_chimera_reverse.md`, `step4_decoy_characterization.md`,
`step4_rescoring.md`, `step4_rosetta_redock.md`; the 12 Å box results and the Vina scores of the ligands come from `notebooks/README.md` and `notebooks/reference_results.csv`. Items I have not verified are not stated as fact.*

<!-- EDITOR'S NOTE FOR LEONARDO (delete before sending)
- Sections changed since 4 Nov: "What we see" (rows, wrong-pose, checks since, residue table), questions 2, 3, 4, 7, 11 are new or rewritten; the old 5-11 are renumbered 8-14 and are unchanged except 11.
- You spoke with Pablo after the first draft. If he already answered any of these (4 and 6 are the likeliest), strike them rather than ask twice.
- Question 7 is open-ended on purpose; it asks whether the decoy region means anything, and I do not claim that it does.
-->
