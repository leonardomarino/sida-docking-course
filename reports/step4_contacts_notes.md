## Notes on the contact check (hand-written; every number below is in the generated table above)

Seed 1 only, 5 receptors, both engines. Descriptive: no pass rule was registered for this check and nothing was tuned.

- **The deposited ornithine makes the same four protein contacts in every receptor:** Lys107 NZ to the carboxylate (2.76-2.98 Å), Ser469 OG to the carboxylate,
  Asn293 to the alpha-N and carboxylate, and Asn323 OD1 to the delta-N. This is the binding mode Robinson 2015 describes. 4B63 and 6X0J chain D also have one
  NADP contact (to the delta-N).
- **The registered PASSES keep all four protein contacts** (6X0J chains A and B, both engines; chain D keeps all four protein contacts and loses only the NADP one).
- **The registered FAILS are not near-misses at the top.** The top-ranked pose of 4B63 (both engines) and of chain C (both engines) loses Lys107, Asn293 and Ser469
  together; its Lys107 NZ is 8.0-8.1 Å from the carboxylate in 4B63 and 11.0 Å in chain C, against 2.9 Å deposited. These poses do not keep the salt bridge.
- **Chain C is a ranking problem.** Both engines have a pose at 1.06-1.11 Å (Uni-Dock rank 4, Vina rank 2) that keeps all four contacts (Lys107 NZ 2.98-3.07 Å).
  The scoring function ranks a wrong pose above it.
- **4B63 is not explained by ranking alone.** No pose in either list keeps all four protein contacts (0 of 12 Uni-Dock, 0 of 20 Vina). The closest-to-crystal poses
  (Uni-Dock 2.07 Å, Vina 1.71 Å) lose Lys107 only (Vina's also loses Asn323 and the NADP contact), and their Lys107 NZ is 4.30 and 3.71 Å from the carboxylate, so
  3.71 Å is 0.2 Å past the 3.5 Å convention. Whether that is a docking failure or a cutoff effect is not settled by this table.
- **What it does not show:** why 4B63 differs from 6X0J. Lys107's side-chain position in the prepared 4B63 receptor was not compared with the 6X0J receptors here,
  so I make no claim that a Lys107 rearrangement is the cause. That is the next measurement, and one of the questions for the Sobrado lab.
