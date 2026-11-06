# Step 4: what the failing 4B63 pose is, and four checks on why (hand-written; probes, not committed code)

Descriptive. Every number below comes from poses and structures already on disk; **no docking was run**. The scripts and their raw output are kept in `/var/tmp/sida-stage/decoy_probes/`
(`decoy_where.py`, `xtal_contacts.py`, `decoy_fingerprint.py`, `decoy_pairwise.py`, `q102.py`, each with its `.out`); they are exploratory probes and are not part of `sidacore`.
Pass rule, benchmark roles and the primary/secondary assignment are unchanged: **4B63 stays primary and fails**, 6X0J chains A, B and D pass, chain C fails.

## 1. Primary result: the failing pose is one stable decoy, and a near-native pose is present just behind it

- **One decoy.** The top-scoring Vina pose of all 10 seeds is the same pose (pairwise RMSD 0.00-0.12 Å; scores -4.90 to -4.93 kcal/mol). Its ligand RMSD to the deposited ornithine is 4.80 Å.
- **Where its carboxylate carbon is** (seed 1; the other seeds differ by under 0.2 Å): 5.8 Å from the deposited carboxylate carbon; 8.9 Å from Lys107 NZ (3.3 Å in the deposited pose);
  9.8 Å from Arg144 CZ; 11.3 Å from Arg328 CZ; 14.2 Å from the crystal sulfate S. Nearest protein atoms: Asn323 OD1 3.9 Å, Gln102 CB 4.0, Leu467 CD1 4.0. It is not at an arginine.
- **The near-native pose is in the list.** Seed 1: poses 5 and 8 (-4.49 and -4.46 kcal/mol) have the carboxylate carbon 1.5 and 1.4 Å from the deposited one, 4.7 Å from Lys107 NZ, against Ser469 / Phe296 / Asn293 like the
  crystal pose. They trail the decoy by 0.43 and 0.47 kcal/mol. Poses 12-15, 18 and 19 also keep the carboxylate within 2.2-4.0 Å of the deposited position.
  So the sampling reaches the right neighbourhood; the failure here is the *ranking* of two sets of poses 0.4 to 0.5 kcal/mol apart.
- **What this does not show.** The 0.4-0.5 kcal/mol gaps are small next to the usual error of a docking score (order of 1-2 kcal/mol; general knowledge, not measured here), so the ranking is not trustworthy either way.
  Earlier in this step the best-in-top-5 median moved from 1.75 to 1.02 Å under 0.1-0.5 Å side-chain changes, so "right pose in the top 5" is itself fragile. Only one seed was read pose by pose (seed 1).

## 2. Crystallization additives: not the cause (negative)

Hypothesis (external review): stripping the extra heteroatoms of 4B63 left a vacated anion site on the Arg144 / Arg328 pair, and the ornithine carboxylate went there.

- 4B63 holds one sulfate (SO4 1495) and four glycerols (GOL 1496-1499) besides FAD, NADP and ornithine; 6X0J has no SO4 or GOL. `prepared/*/removal_log.csv` already lists all of them as
  "crystallization additive / non-cofactor heterogen", with their distance to the site.
- The sulfate S is 14.5 Å from the deposited carboxylate carbon and 14.7 Å from the box centre; offset (x, y, z) from the box centre (-11.5, -5.0, -7.6) Å against a half-edge of 11 Å, so it lies at the edge of the
  box. The four glycerols are 26-29 Å from the site. None is in the substrate cavity.
- The decoy carboxylate is 9.8 and 11.3 Å from Arg144 CZ and Arg328 CZ and 13.3-14.2 Å from the sulfate (nearest O 13.3 Å). It is not in an anion site.
- **One thing survives.** Arg328 CZ is 4.4 Å from the sulfate S in 4B63 (Arg144 CZ: 13.9 Å); 6X0J has no sulfate. This is a candidate explanation of why Arg328's guanidinium differs between the structures (3.5-3.7 Å, Step 4 site
  comparison), i.e. one of the four flagged site residues may reflect an additive rather than biology. It is a distance, not a demonstrated cause, and it is 8.8 Å from the ornithine, so it does not bear on the decoy.
  Whether the sulfate is real (B 55-60 Å^2, occupancy 1.0) is a question for Sobrado.
- The review also attributed the larger fpocket cavity of 4B63 to the stripped heteroatoms. Not tested; given the distances above it is unlikely for the glycerols and unproven for the sulfate.

## 3. Crystal packing: not the cause (negative)

- gemmi contact search, side-chain atoms only, cutoff 4.5 Å, adjacent residues ignored, symmetry mates and other chains of the asymmetric unit included.
- 4B63 (I 2 2 2, one chain) and 6X0J chains A-D (P 1 21 1): **Gln102, Arg144 and Arg328 have no contact within 4.5 Å in any chain; Asp288 touches only a water** (4B63: 4.3 Å to a symmetry-related water; 6X0J chain A: 4.1 Å
  to a water of chain D; none in B-D).
- Lys107 contacts Ile132 of a neighbouring chain (3.7-3.9 Å) in every chain of both structures, so it does not separate them.
- So the Gln102 side-chain difference between 4B63 and 6X0J is not a packing artefact by this test. (4.5 Å is a contact criterion, not a statement about long-range lattice effects.)

## 4. Is the Gln102-chimera decoy the 4B63 decoy? No (negative), and it shows a second, shared decoy

All top poses carried into the 4B63 frame by the Step 3 trimmed CA superposition (fit 0.30-0.35 Å over 394-413 CA); compared with `sidacore.redock.pose_rmsd` (symmetry-corrected, no re-alignment).

- **Chimera vs 4B63 decoy.** 6X0J:A and 6X0J:B with 4B63's Gln102 (top-1 5.07 and 4.83 Å): 3.66 and 3.70 Å from the 4B63 decoy; **0 of 10 seeds within 2 Å**. For scale, the correct poses are 4.7 Å from the 4B63 decoy.
  Carboxylate carbon 6.8-7.0 Å from the deposited position (4B63 decoy: 5.8) against Gln102 / Thr322 / Ser286 (4B63 decoy: Asn323 / Gln102 / Leu467); 9.3-9.5 Å from Arg328 CZ (4B63 decoy: 11.3).
- **Chimera vs chimera vs native chain C** (seed-1 top poses, Å): A+Gln102 vs B+Gln102 **0.27** (the same pose); either vs native 6X0J:C **1.74 and 1.86**; 6X0J:D with Gln102 and Asp288 swapped is 2.39 from the 4B63 decoy,
  3.15-3.19 from A/B, 3.27 from C. Round-trip baselines and the Asp288-only swap (passing controls) sit at 4.68-4.73 from the 4B63 decoy, i.e. on the correct pose.
- **Gln102 in native chain C is not 4B63-like**: side-chain RMSD to 4B63 2.11 Å, the same as chains A, B and D (1.98-2.01 Å). A decoy of the same family as the chimera decoys therefore occurs in an unmodified 6X0J
  receptor whose Gln102 is the 6X0J one.
- **Reading.** There is a decoy region around Gln102 / Thr322 / Ser286 / Asn323 with several nearby poses; the 4B63 decoy, the A/B/C decoy and the D-chimera decoy are different members. Every one of them contacts Gln102, but
  the transplant of 4B63's Gln102 does not reproduce the 4B63 decoy, and chain C has such a decoy without it. This **retires "4B63's Gln102 creates the decoy subsite"** as a mechanism. It does not rule Gln102 out as one of several
  factors in the Gln102 chimera result, which stays confounded by the 1.65-1.91 Å overlap with NADP that the transplant creates.
- Limits: the pairwise table uses seed-1 top poses only. In the D chimera some seeds come closer to the 4B63 decoy (min 2.15-2.54 Å over seeds) and some reach the correct pose (3/10 pass); not analysed seed by seed.

## What this changes

- The Step 4 result is restated as: *the 4B63 top pose is a stable decoy 4.8 Å from the deposited ligand, not explained by crystallization additives, crystal packing or a single transplantable residue; a near-native pose is present at 0.4-0.5 kcal/mol behind it.*
- Site comparison: Arg328 now has a candidate additive explanation, so three of the four flagged residues (Gln102, Arg144, Asp288) remain unexplained.
- The ladder (protonation, box, cofactor, pH) found nothing because it varies the receptor and sampling, and this is a ranking problem; the one experiment aimed at that is an independent rescoring of the existing Vina
  pose lists, after a licence check (see the licence survey).
- Not done: the NADP-removed Gln102 chimera (baseline there is 5.07-5.09 Å, not 4.80); the Asp288 loop backbone; the questions to Sobrado (drafted, not sent).
