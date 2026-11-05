# SidA ornithine docking: what we see, and what we would like your view on

*Draft from Leonardo Mariño-Ramírez to Pablo Sobrado. For review before sending.*

## Why I am writing

I am building a docking benchmark and teaching course around SidA, and I would like the protocol to be right in a way you would accept. Before using any
structure to compare ligands, I redocked the deposited L-ornithine into the structures that contain it. I fixed the pass rule beforehand: the median top-ranked
pose over 10 seeds within 2 Å (symmetry-corrected heavy-atom RMSD) of the deposited ligand. Receptors were prepared with FAD and NADP(H) kept, waters removed,
protein hydrogens from PDB2PQR/PROPKA at pH 7.4. Two engines were used, AutoDock Vina 1.2.7 and Uni-Dock, and their top-ranked poses agree.

## What we see

| Structure | Result | Note |
|---|---|---|
| 6X0J chains A, B, D | pass (0.9-1.1 Å, all seeds) | The deposited pose is recovered at rank 1 |
| 6X0J chain C | **fail** (top pose 6.0 Å) | A pose at 1.1 Å is in the list at rank 2-4: the scoring ranks a wrong pose above it |
| 4B63 (1.9 Å, oxidized FAD, NADP+) | **fail** (top pose 4.8 Å) | No pose in either list keeps all four protein contacts |

The deposited ornithine makes the same four protein contacts in every structure: Lys107 NZ to the carboxylate (2.8-3.0 Å), Ser469, Asn293 and Asn323. The wrong top
poses lose Lys107, Asn293 and Ser469 together (Lys107 NZ 8-11 Å from the carboxylate). I then varied one factor at a time: ornithine protonation (including the
neutral delta-amine you argue for), box size and position, removing FAD and/or NADP, and receptor pH. Removing either cofactor, or enlarging the box, rescues
chain C. **Nothing rescues 4B63.** I have not tuned anything to make it pass, and I do not know why 4B63 differs from 6X0J.

From the literature I found no published redocking of ornithine into SidA. Campbell et al. (2020) docked with NADP(H), ornithine and waters removed and FAD kept.

## Questions (1-4 matter most)

**About the structures**

1. **4B63 versus 6X0J.** Is there anything about 4B63 that I should know: crystallization pH and conditions, lattice contacts near the ornithine site, why one
   chain, how confident you are in the ornithine and Lys107 positions? Is there any reason chain C of 6X0J would behave differently from A, B and D?
2. **Lys107.** In 4B63 its NZ is 2.9 Å from the carboxylate. Is this residue mobile in your maps or simulations, and do K107 mutant data (Km for ornithine) exist that
   I can cite? The closest-to-crystal 4B63 poses (1.7-2.1 Å) place Lys107 NZ 3.7-4.3 Å from the carboxylate.
3. **Protonation to use.** Your group argues the ornithine delta-amine binds neutral (Protein Sci, "Flavin oxidation in flavin-dependent N-monooxygenases"). Which
   computational studies does that paper cite as refs 7 and 8? For docking, which state would you defend for the delta-amine, the alpha-amine and Lys107? My result did
   not depend on it, but I would like to use the state you consider correct.
4. **Which cofactors belong in the receptor.** Campbell et al. removed NADP(H); the MD of Pierdominici-Sottile et al. (2024) suggests NADPH blocks ornithine in the FAD-out state. For docking into the
   ornithine site, would you keep NADP(H), remove it, or run both? Is the NADP contact to the delta-N in 4B63 (3.3 Å) catalytically meaningful?

**About details I handled by default**

5. **Waters.** I remove all of them (262 in 4B63; one is 3.2 Å from FAD C4X). Are any waters in the ornithine site structural?
6. **FAD state.** 4B63 has oxidized FAD; 6X0J has reduced flavin. Do you expect ornithine to bind the same way in both, and which is the relevant state?
7. **FAD-out.** I plan to compare 6X0H (FAD-out, no NADP) with 6X0I (FAD-in, NADP). Because no deposited structure has FAD-out with NADP, I will report that the
   comparison is confounded and not attribute differences to FAD. Do you agree, and is the M101A structure (7JVL) or another structure a better FAD-out reference?

**Data and collaboration**

8. **Ranking test.** Can docking rank ornithine above lysine and arginine? I understand the 2012 Franceschini paper has lysine and arginine complexes. Can you point me
   to the Km/Kd values, and to any inhibitors with a known binding site (for example celastrol) that I may use as a test set?
9. **Mutant structures.** For 5CKU (N323A) and 4NZH (R279A): is ornithine bound in the same pose, and are kinetics available? I would use them as a specificity exercise.
10. **Your own docking.** Have you or your group docked ornithine into 4B63 or 6X0J with any program? If so, would you share the protocol and any redocking numbers?
11. **Review.** Would you be willing to read the protocol section, and tell me where an enzymologist would disagree?

*All numbers above come from runs in the course repository; the polar-contact table is `reports/step4_contacts.md`. Items I have not verified are not stated as fact.*
