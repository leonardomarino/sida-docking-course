## Notes on the Lys107 comparison (hand-written; numbers are in the generated table, except where marked)

- **Lys107 is not different in 4B63.** Its end atom (NZ) sits 0.09-0.59 Å from the 4B63 position in the four 6X0J chains after fitting all CA atoms, and 0.17-0.55 Å
  after fitting only the CA atoms within 12 Å of the deposited ornithine. The fits themselves have RMSDs of 0.55-0.59 Å (all CA) and 0.27-0.28 Å (site CA), so these
  displacements are about the size of the fit's own error: no detectable difference.
- **The side chain has the same shape.** chi1 is -55 degrees in 4B63 and -59 to -60 in 6X0J; chi2 is -161 against -172 to -174; chi3 is -175 against 174 to 178
  (the same trans, either side of the 180 wrap); chi4 is -174 against -152 to -154, the largest difference at about 20 degrees.
- **It is fully ordered and held the same way.** Occupancy 1.0 in every chain, NZ B-factor 18.9-32.2 Å^2 (4B63: 32.2). NZ is 2.76-2.98 Å from the deposited
  carboxylate in all five (4B63: 2.87), and the only other polar neighbours within 3.5 Å are one or two waters.
- **The docking receptors kept these coordinates.** Checked with a one-off comparison, not a committed script: the heavy atoms of Lys107 in each receptor's
  assembled `work/combined.pdb` match the deposited coordinates to 0.000 Å in all five.
- **What this rules out, and what it does not.** A different Lys107 position or order in 4B63 is not the reason the 4B63 redock fails; the receptor side of that
  salt bridge is the same as in the 6X0J chains that pass. It does not say what is different. The failing poses lose the contact because the ligand does not go
  there, and that is still unexplained. The next measurement is the same comparison over all residues of the site, 4B63 against 6X0J, not just Lys107.
