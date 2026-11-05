## Notes on the chimeric-receptor diagnostic (hand-written; numbers are in the generated table or `chimera/step4_chimera.json`)

Design and pass rule were registered before the run. One method change was made before any result was seen: the first attempt copied donor coordinates after a global fit, which stretched the
CA-CB bond to 2.05 A and made Meeko refuse the receptor; it now transplants the donor rotamer onto the 4B63 backbone (CA-CB 1.52-1.54 A in every variant). The first attempts produced no table.
Seeds are fixed, and two independent runs of the five variants that completed in both gave identical medians.

- **The method does not perturb the result.** The round trip with nothing swapped gives 4.80 A (0/10 seeds), best-in-top-5 1.75 A, the same as the original 4B63 redock.
- **No variant rescues 4B63.** None of the six preparable variants has a seed within 2 A. Top-1 is 4.79-4.83 A in five of them and the top-ranked score barely moves (median -4.90 to -4.93 kcal/mol
  against -4.91), so the wrong pose stays the lowest-energy one whatever is swapped.
- **The Gln102 side-chain conformation from 6X0J does not help, and makes the near-native pose rank worse.** With it, the best pose in the list sits at a median rank of 11.5-12 (baseline 7),
  and best-in-top-5 goes from 1.75 to 4.79-4.80 A. This is the opposite of what the hypothesis needed.
- **Asp288's side chain is not the difference.** After placing it on the 4B63 backbone it moves only 0.6 A (0.9 A from chain C), so the Asp288 difference between the structures is mostly backbone, which
  this experiment does not change. The Asp288-only variant is unchanged (4.81 A, best-in-top-5 1.72 A). Whether the Asp288 loop backbone matters is untested.
- **The 6.28 A for Gln102 + Asp288 is not a stable effect.** Its seeds run from 4.79 to 7.81 A, with some seeds finding a second wrong pose near 7.8 A; the median falls between the two groups. The Gln102-only and
  donor-chain-C variants show the same split but medians of 4.79 and 4.81 A.
- **The control shows how fragile the best-in-top-5 number is.** Swapping four contact residues that differ by only 0.1-0.5 A leaves top-1 at 4.83 A but moves the near-native pose from median rank 7 to 4.5
  and best-in-top-5 from 1.75 to 1.02 A (seed range 0.75-4.80). The near-native pose sits at the edge of the top 5, so sub-angstrom changes flip it across that line. The Gln102 effect (rank 7 to about 12) is
  larger than the control's (7 to 4.5) but of the same order, so I do not read it as more than "not helpful".
- **No other trace of a cause in the site.** No swapped side chain came within 2.2 A of another residue in the preparable variants, and the PROPKA protonation states of the site residues did not change in any.
- **The all-four-flagged variant could not be built.** Gln102 and Arg328 in their 6X0J conformations overlap in the 4B63 context (OE1 to NE 1.62 A), so Meeko inferred a spurious bond and refused; it is reported
  as not prepared. That overlap is consistent with the two residues being coupled (the Gln102 rotamer needs the 6X0J position of Arg328), which would also fit both being flagged; this is an inference, not tested.
- **What this does and does not show.** It shows that giving 4B63 the 6X0J side-chain conformation of Gln102 or Asp288 does not make the redock pass, so those two side chains are not sufficient to explain it.
  It does not show they are irrelevant (a reverse chimera, 6X0J with 4B63's residues, tests necessity), does not test the Asp288 backbone, and cannot test Arg144 or Arg328 in combination with Gln102.
