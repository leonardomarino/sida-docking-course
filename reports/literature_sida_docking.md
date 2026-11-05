# Literature: docking and modelling studies of SidA (PubMed Central)

Compiled from a PMC/PubMed search run with NCBI E-utilities (key not recorded here). **Verified by me against the retrieved full text:**
Campbell 2020 (docking methods), Robinson 2015 (ornithine contacts), Pierdominici-Sottile 2024 (model building; the 'Orn cannot enter' claim),
Rosy 2022 (receptor preparation, no redocking) and Robinson/Klancher/Rodriguez/Sobrado, Protein Sci (PMC6295899; Orn protonation). Everything else
below is the search agent's reading and has **not** been re-checked line by line; treat it as leads. Two papers (Campbell 2020, Pierdominici-Sottile
2024) have no body in PMC's efetch XML and were read from the article web pages. The reference list of PMC6295899 could not be parsed from the saved
XML, so the 'previous computational studies' it cites for the neutral amine (its refs 7 and 8) are **not identified**.

## What the search found

Only **one** paper docked into SidA structures, and none redocked a crystallographic ligand (ornithine, lysine or NADP) or reported a pose RMSD. No
docking-program benchmark on flavin N-monooxygenases was found. So the redocking benchmark in this repo (4B63 and 6X0J, symmetry-corrected RMSD, 10 seeds)
has no published counterpart in what was retrieved; a null result in a search is not proof that none exists.

| PMCID | Year | Paper | Software and receptor | Relevance |
|---|---|---|---|---|
| PMC7504930 | 2020 | Campbell et al., JBC: new SidA structures (6X0H/I/J/K) | AutoDock Vina 1.1.2, AutoDockTools 1.5.6, Open Babel 2.3.2; 6X0I (FAD-in) and 6X0H (FAD-out) | The only SidA docking: 16 small probes, then 3,035 ChEMBL compounds into the out state. **Verified:** "Solvent molecules, NADP(H), and l-Orn were removed ... FAD was retained"; a large box covering the FAD, NADP(H) and ornithine sites and the channels between them; maximum exhaustiveness. No redocking; the probe control only showed clusters on the NADP site |
| PMC6467063 | 2015 | Robinson et al., Arch Biochem Biophys: ornithine-binding residues | mutagenesis; 4B63 as reference | **Verified:** the ornithine site "includes ionic interaction with Lys107, and hydrogen bonding interactions with Asn293, Ser469, and Asn323" |
| PMC10962481 | 2024 | Pierdominici-Sottile et al., Protein Sci: MD of NADPH and L-Orn binding | AMBER18; 6X0H, 6X0I, 6X0J; Orn models built from the crystal structures (6X0H acetate replaced by Orn; the 6X0J model has NADPH and Orn 'docked', method not stated in the text I read) | **Verified:** umbrella sampling could not remove Orn from FAD-out + NADPH ('the bound NADPH molecule blocked the way'); Orn gets in through the FAD gate in the FAD-out state, or after NADPH has pushed FAD in. The productive complex is FAD-in + NADPH + Orn, which is what 4B63 and 6X0J contain |
| PMC6295899 | 2018/19 | Robinson, Klancher, Rodriguez, Sobrado, Protein Sci: flavin oxidation in N-monooxygenases | kinetics; no docking | **Verified:** argues SidA lowers the pKa of the Orn delta amine so Orn binds with the side chain **neutral**; the alpha-amino group hydrogen-bonds Asn293 and need not be deprotonated. Disagrees with Fredrick et al. (protonated on binding). Cites earlier computational studies (refs not retrieved) |
| PMC12432408 | 2025 | R144 in the flavin oxidation mechanism of a flavin-dependent N-monooxygenase | kinetics | R144A slows flavin dehydration ~60-fold with Orn; no docking or pose |
| PMC3820879 | 2013 | Shirey et al., JBC: Ser257 and NADP(H) sliding | GROMACS; 4B63 | MD only |
| PMC9313252 | 2022 | Rosy et al., Biomolecules: PvdA ornithine-analogue repurposing | Vina, H++ pH 8.3, Desmond; PvdA 3S5W | **Verified:** one chain kept; 'Water molecules and the heteroatoms, including the FAD and NADPH ... were removed'; no redocking or pose validation (RMSD in the text refers to MD). The two-energies claim is not checked |
| PMC12225907 | 2025 | Liu et al., J Adv Res: galangin as a PvdA inhibitor | Vina | Docking is a small part; no validation stated |
| PMC9005580 | 2022 | Shyam Mohan et al., JGEB: putrescine N-monooxygenase | I-TASSER model, AutoDock 4.2.6 | Homology model only; two hits tested in vitro |

Abstract-only, not extractable: Martin del Campo 2016 (celastrol, PMID 27588426; abstract re-read: 'Docking experiments suggest that celastrol binds across the NADPH
and l-Orn pocket'; a noncompetitive inhibitor; worth finding the full text), Franceschini 2012 (SidA with Orn, Lys, Arg), Badieyan 2015 (DFT), Bufkin 2026 (ebselen).

## What it means for this project

- **The closest precedent removes the NADP.** Campbell et al. docked with FAD kept and NADP(H) removed. In our diagnosis ladder, removing NADP (or FAD) is
  exactly what rescues 6X0J chain C. That does not make it the right treatment for the physical receptor, but it shows the published practice and our
  chain C failure are related.
- **No published redocking on SidA** means there is no reference number to compare our 4.8 Å (4B63) and 1.0 Å (6X0J A, B, D) against.
- **The crystallographic binding mode** is a Lys107 salt bridge to the carboxylate with Asn293, Ser469 and Asn323 hydrogen bonds. Our chain C wrong pose puts
  the carboxylate 8-10 Å from the deposited one; whether it loses the Lys107 contact has not been measured.
- **Ornithine protonation is a live question in the literature, and our ladder already covers it.** The Sobrado lab argues for a neutral delta amine (above).
  That state ran in rung 1: 4B63 stayed at 4.80-4.82 Å and chain C at 5.93-6.06 Å. The protonation debate does not explain our failures.
- **Productive-complex logic.** The MD work supports FAD-in + NADPH + Orn as the productive complex and shows NADPH blocking Orn in the FAD-out state. That
  supports the caution already written for Step 5: no deposited structure has FAD-out with NADP, so an FAD-out + NADP arm has no physical reference.
- **Every docking precedent strips the cofactors** (Campbell: NADP(H) and waters removed, FAD kept; Rosy on PvdA: FAD, NADPH and waters removed). None validated
  a pose. The MD papers build Orn in from the crystal structure rather than searching for it.
