# Literature: docking and modelling studies of SidA (PubMed Central)

Compiled from a PMC/PubMed search run with NCBI E-utilities (key not recorded here). **What was verified by me against the retrieved full text:**
Campbell 2020 (docking methods) and Robinson 2015 (ornithine contacts). Everything else below is the search agent's reading of the papers and
has **not** been re-checked line by line; treat it as leads. Two papers (Campbell 2020, Pierdominici-Sottile 2024) have no body in PMC's efetch XML
and were read from the article web pages.

## What the search found

Only **one** paper docked into SidA structures, and none redocked a crystallographic ligand (ornithine, lysine or NADP) or reported a pose RMSD. No
docking-program benchmark on flavin N-monooxygenases was found. So the redocking benchmark in this repo (4B63 and 6X0J, symmetry-corrected RMSD, 10 seeds)
has no published counterpart in what was retrieved; a null result in a search is not proof that none exists.

| PMCID | Year | Paper | Software and receptor | Relevance |
|---|---|---|---|---|
| PMC7504930 | 2020 | Campbell et al., JBC: new SidA structures (6X0H/I/J/K) | AutoDock Vina 1.1.2, AutoDockTools 1.5.6, Open Babel 2.3.2; 6X0I (FAD-in) and 6X0H (FAD-out) | The only SidA docking: 16 small probes, then 3,035 ChEMBL compounds into the out state. **Verified:** "Solvent molecules, NADP(H), and l-Orn were removed ... FAD was retained"; a large box covering the FAD, NADP(H) and ornithine sites and the channels between them; maximum exhaustiveness. No redocking; the probe control only showed clusters on the NADP site |
| PMC6467063 | 2015 | Robinson et al., Arch Biochem Biophys: ornithine-binding residues | mutagenesis; 4B63 as reference | **Verified:** the ornithine site "includes ionic interaction with Lys107, and hydrogen bonding interactions with Asn293, Ser469, and Asn323" |
| PMC10962481 | 2024 | Pierdominici-Sottile et al., Protein Sci: MD of NADPH and L-Orn binding | AMBER18; 6X0H, 6X0I, 6X0J; ornithine placed by superposition, not docked | Agent-reported: ornithine cannot enter when NADPH is bound and FAD is out; it can when FAD is in |
| PMC3820879 | 2013 | Shirey et al., JBC: Ser257 and NADP(H) sliding | GROMACS; 4B63 | MD only |
| PMC9313252 | 2022 | Rosy et al., Biomolecules: PvdA ornithine-analogue repurposing | Vina, H++ pH 8.3, Desmond; PvdA 3S5W | Agent-reported: FAD, NADPH and waters removed; no redocking; the same compound is given two different energies in one paragraph |
| PMC12225907 | 2025 | Liu et al., J Adv Res: galangin as a PvdA inhibitor | Vina | Docking is a small part; no validation stated |
| PMC9005580 | 2022 | Shyam Mohan et al., JGEB: putrescine N-monooxygenase | I-TASSER model, AutoDock 4.2.6 | Homology model only; two hits tested in vitro |

Abstract-only, not extractable: Martin del Campo 2016 (celastrol, PMID 27588426; the abstract says docking suggests it spans the NADPH and ornithine
pockets; worth finding the full text), Franceschini 2012 (SidA with Orn, Lys, Arg), Badieyan 2015 (DFT), Bufkin 2026 (ebselen).

## What it means for this project

- **The closest precedent removes the NADP.** Campbell et al. docked with FAD kept and NADP(H) removed. In our diagnosis ladder, removing NADP (or FAD) is
  exactly what rescues 6X0J chain C. That does not make it the right treatment for the physical receptor, but it shows the published practice and our
  chain C failure are related.
- **No published redocking on SidA** means there is no reference number to compare our 4.8 Å (4B63) and 1.0 Å (6X0J A, B, D) against.
- **The crystallographic binding mode** is a Lys107 salt bridge to the carboxylate with Asn293, Ser469 and Asn323 hydrogen bonds. Our chain C wrong pose puts
  the carboxylate 8-10 Å from the deposited one; whether it loses the Lys107 contact has not been measured.
