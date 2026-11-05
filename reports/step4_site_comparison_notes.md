## Notes on the all-residue site comparison (hand-written; numbers are in the generated table or `boxes/step4_site_comparison.json` (untracked; regenerate with `python -m sidacore.sitecompare --json boxes/step4_site_comparison.json`), except where marked as a probe)

Descriptive. The flag rule was registered before the run and was not changed afterwards. A flagged residue is a lead, not a cause.

- **Almost the whole site is the same.** 80 of 84 site residues are not flagged (8 of them are glycines, which have no side chain to compare). That includes all four protein contacts of the deposited ornithine
  (Lys107 0.21-0.33 Å, Asn293 0.21-0.26, Asn323 0.19-0.39, Ser469 0.13-0.30) and the other contacting residues (Ile103, Phe296, Thr322, Leu467, each 0.43 Å or less).
- **Four residues differ.** Gln102 (1.93-2.08 Å from the 6X0J chains, against a largest 6X0J-6X0J difference of 0.35), Arg144 (2.67-2.75; 0.24), Asp288 (1.67-1.83; 0.21)
  and Arg328 (3.53-3.68; 0.78). The nearest unflagged residue is Arg409 (0.75 Å, 6X0J spread 0.43).
- **Two of the four are near the substrate** (probe, per-atom shifts after the same fit, 6X0J chains A and C; not a committed script): Gln102 is a contact residue of the
  deposited ornithine (closest side-chain atom 3.9 Å in 4B63) and it is a real reorientation, not an amide flip: CB moves 0.5-0.7 Å, CG 1.0-1.3 and the amide atoms 2.5-3.0.
  Asp288 moves as a whole, 1.4-2.3 Å across its atoms, and is 5.1 Å from the ornithine in 4B63. This agrees with the Step 2 site table, where Asp288 touches the ornithine in all four
  6X0J chains and not in 4B63. Arg144 and Arg328 differ only at the guanidinium end (3-6 Å) and lie 9.3 and 8.8 Å from the ornithine.
- **Cofactors.** Flavin atoms common to both structures differ by 0.55-0.63 Å (the largest 6X0J-6X0J value is 0.28; oxidized FAD in 4B63, reduced in 6X0J), the isoalloxazine ring by
  0.24-0.27 (spread 0.16) and NADP by 0.75-0.85 (spread 0.67). None is flagged: all are under the 1.0 Å floor.
- **The 6X0J chains are close to each other, and this does not explain chain C.** The largest within-6X0J difference is Met101 (1.74 Å), driven by chain B (1.37-1.74 Å from A, C and D, which
  are within 0.64 of each other); chain B passes the redock, so Met101 does not separate pass from fail. For Arg328, chain C differs from A, B and D by 0.73-0.78 Å against at most 0.33 among
  those three; it is the only chain-C-specific difference I looked at, a weak lead 8.8 Å from the ornithine. I checked four of the five residues with a 6X0J spread above 0.5 Å in this way
  (Met101, Leu145, Glu330, Arg328), not Arg152.
- **What this does not show.** Whether Gln102 or Asp288 changes the docking result is untested. A direct test would be a labelled diagnostic: dock into 4B63 with those two residues
  taken from 6X0J. That would show whether the receptor differences matter, not fix anything. Also, 4B63 is one chain from one crystal form, so a difference could reflect crystal packing or
  refinement as much as biology.
