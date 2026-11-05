Hand-written. The tables above are generated; everything below is interpretation, and is an observation, not a conclusion.

**The superposition method checks out against itself.** Before measuring, the expectation was registered (in
`tests/test_boxes.py`) that carrying 6X0I's C4X into another FAD-in structure lands within 2 A of that structure's own
C4X. All five other FAD-in receptors (4B63, 6X0J chains A-D) land within 0.2-0.3 A, with 98-99% box overlap.

**FAD-out moves the C4X box by 6.5 A.** In 6X0H the carried box (the 6X0I site expressed in 6X0H's frame) and the box
centred on 6X0H's own C4X overlap by 63%. The protein around the site barely moves between the two (site CA RMSD
0.29 A), but that is over only 5 CA atoms within 8 A of C4X, so it is a coarse measure. Both 6X0H boxes enclose the
same fpocket pocket (86 alpha spheres, all inside either box), so docking in either samples the same cavity. Which box
Step 5 uses for 6X0H is not decided here: the carried box is the "equivalent box" the step asks for, the C4X box is the
sensitivity check, and both will be reported.

**The deposited substrate is inside the box but not at its centre.** In 4B63 and all four 6X0J chains, all 9 ORN heavy
atoms lie inside the 22 A box, with the ORN centroid 7.7-8.0 A from C4X. The nearest alpha sphere to every ORN atom
(no threshold) belongs to fpocket's rank-1 pocket, at a median 0.6-0.9 A and at most 1.7 A. So fpocket independently
finds the substrate site.

**fpocket and the box disagree under the rule fixed in advance, and the rule is left as it was.** The pre-set rule
(site pocket must be rank 1 AND its alpha-sphere centroid within 5 A of the box centre) flags every row. The reason is
not that fpocket looks elsewhere: in every FAD-in receptor the pocket holding the box is rank 1, FAD lines it, and it
holds the substrate. It is that this is a large cleft (about 1,800-3,300 A^3) whose centroid is 9-10.5 A from C4X. The
5 A figure was a guess made before seeing the data. A different criterion can be chosen, but it would be post hoc and
must be labelled so; it has not been tuned to make the table pass.

**FAD-out presents a much smaller pocket in the box, but the structures cannot attribute that to FAD.** fpocket's pocket
in the 6X0H box is 453 A^3 (druggability 0.66, 2 FAD atoms in the lining) against about 1,800 A^3 (0.996, 17 FAD atoms)
in 6X0I and 2,100-2,400 A^3 in 6X0J. No deposited SidA has FAD-out with NADP bound, so FAD conformation and NADP
occupancy are confounded. Stripping NADP from 6X0I does not remove the confound (below). Docking scores in Step 5 will
carry the same limitation.

**NADP-stripped 6X0I (unexplained, not investigated).** The site pocket is geometrically unchanged (227 alpha spheres,
146 inside the box, volume within fpocket's own run-to-run noise), but its fpocket score falls from 0.804 to -0.348 and
its rank from 1 to 34. The stripped arm is not a controlled comparison; this is one more sign of that.

**Two fpocket behaviours worth knowing.**
- fpocket silently ignores hetero residues whose names it does not know and logs nothing. It keeps `FAD` and `NAP` and
  drops the reduced flavin `FDA` (6X0J). Before the relabel, 6X0J's rank-1 site pocket had 158 alpha spheres and no FAD
  in its lining; relabelled `FAD` it has 257 and FAD lines it. The relabel is applied to fpocket's input only and is
  recorded per receptor. It is covered by a regression test.
- The reported volume is Monte Carlo and varies by about 2-3% between identical runs. Score and sphere count are
  deterministic.
