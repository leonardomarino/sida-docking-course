Hand-written. The tables above are generated; everything below is interpretation, and is an observation, not a conclusion.

**The superposition method checks out against itself.** Before measuring, the expectation was registered (in
`tests/test_boxes.py`) that carrying 6X0I's C4X into another FAD-in structure lands within 2 A of that structure's own
C4X. All five other FAD-in receptors (4B63, 6X0J chains A-D) land within 0.2-0.3 A, with 98-99% box overlap.

**FAD-out moves the C4X box by 6.5 A.** In 6X0H the carried box (the 6X0I site expressed in 6X0H's frame) and the box
centred on 6X0H's own C4X overlap by 63%. The site itself moves more than in any FAD-in structure: over the 80 residues that lie in the docking box, the CA RMSD to 6X0I is 0.61 A with one residue displaced 3.98 A, against 0.22-0.26 A (largest displacement 0.4-0.6 A) for 4B63 and the four 6X0J chains. Both 6X0H boxes enclose the
same fpocket pocket (86 alpha spheres, all inside either box), so docking in either samples the same cavity. That is a statement about the cavity, not about the residues: the box centred on 6X0H's own C4X
does not contain Lys107, Asp108 or Glu146, which are inside the 6X0I box and inside the carried 6X0H box (same composition as 6X0I's). The carried box is
the "equivalent box" the step asks for and is the comparable one; the C4X box would be a sensitivity check. Both will be reported in Step 5.

**The deposited substrate is inside the box but not at its centre.** In 4B63 and all four 6X0J chains, all 9 ORN heavy
atoms lie inside the 22 A box, with the ORN centroid 7.7-8.0 A from C4X. The nearest alpha sphere to every ORN atom
(no threshold) belongs to fpocket's rank-1 pocket, at a median 0.6-0.8 A and at most 1.7 A. So fpocket independently
finds the substrate site.

**fpocket and the box disagree under the rule fixed in advance, and the rule is left as it was.** The pre-set rule
(site pocket must be rank 1 AND its alpha-sphere centroid within 5 A of the box centre) flags every row. The reason is
not that fpocket looks elsewhere: in every FAD-in receptor the pocket holding the box is rank 1, FAD lines it, and it
holds the substrate. It is that this is a large cleft (about 1,800-3,300 A^3 across repeat runs) whose centroid is 9-10.5 A from C4X. The
5 A figure was a guess made before seeing the data. A different criterion can be chosen, but it would be post hoc and
must be labelled so; it has not been tuned to make the table pass.

**FAD-out presents a much smaller pocket in the box, but the structures cannot attribute that to FAD.** fpocket's pocket
in the 6X0H box is about 440 A^3 (439-453 over three runs; druggability 0.66, 2 FAD atoms in the lining) against about 1,800 A^3 (1,774-1,892 over three runs; 0.996, 17 FAD atoms) in 6X0I and 2,059-2,414 A^3 in 6X0J. No deposited SidA has FAD-out with NADP bound, so FAD conformation and NADP
occupancy are confounded. Stripping NADP from 6X0I does not remove the confound (below). Docking scores in Step 5 will
carry the same limitation.

**NADP-stripped 6X0I (unexplained, not investigated).** The site pocket is geometrically unchanged (227 alpha spheres,
146 inside the box; its volume differs between the keep and strip receptors by less than fpocket's own run-to-run spread: keep 1,774-1,892 over three runs, strip 1,760-1,835 over two, and the sign of the difference flipped between runs), but its fpocket score falls from 0.804 to -0.348 and
its rank from 1 to 34. The stripped arm is not a controlled comparison; this is one more sign of that.

**Two fpocket behaviours worth knowing.**
- fpocket silently ignores hetero residues whose names it does not know and logs nothing. It keeps `FAD` and `NAP` and
  drops the reduced flavin `FDA` (6X0J). Before the relabel, 6X0J's rank-1 site pocket had 158 alpha spheres and no FAD
  in its lining; relabelled `FAD` it has 257 and FAD lines it (a one-off check; the pre-relabel output was not kept). The relabel is applied to fpocket's input only and is
  recorded per receptor. It is covered by a regression test.
- The reported volume is Monte Carlo and varies by up to about 6% between identical runs (three runs of one receptor: 1,774-1,892 A^3). Score and sphere count are
  deterministic.
**Limits of these measures.** Volumes are quoted as ranges over repeat runs because fpocket's volume is Monte Carlo; the numbers in the
generated tables are from the run on disk (`boxes/step3_report.json`). The site used for the site RMSD is every protein residue with an atom
in the docking box (80 residues). An earlier version used the five CA atoms within 8 A of C4X, a sphere centred on the cofactor that missed
the substrate site (the deposited ornithine's centroid is 7.7-8.0 A from C4X).
