# Construct check: 4B63 vs the 6X0* batch (Amendment 1, mandatory)

**Verdict: 4B63 remains the primary benchmark; 6X0J the secondary.** The sequences of the resolved residues are
identical and the author numbering is shared, so no residue mapping or renumbering is needed. Nothing was renumbered.

Method: biotite optimal alignment of the resolved polymer residues of chain A (`sidacore.mapping`), pairing each
residue of the reference with one of the query. Full per-residue tables: `reports/construct_check/*.csv`
(columns `ref_num, ref_aa, qry_num, qry_aa, identical, same_number`). Regression-tested in `tests/test_mapping.py`.

| reference | query | aligned | identical | mismatches | numbers differ | resolved only in reference | resolved only in query |
|---|---|---|---|---|---|---|---|
| 4B63 | 6X0H | 448 | 448 | 0 | 0 | 29, 490-491 | 71-74, 384-385, 391 |
| 4B63 | 6X0I | 443 | 443 | 0 | 0 | 29, 68-70, 75, 392, 490-491 | none |
| 4B63 | 6X0J | 447 | 447 | 0 | 0 | 29, 392, 490-491 | 71-74 |
| 6X0H | 6X0I | 443 | 443 | 0 | 0 | 68-75, 384-385, 391-392 | none |

- **No point mutations** anywhere in the aligned residues (0 mismatches of 1,338 aligned pairs across the three
  comparisons to 4B63). The differences are only in *which terminal and loop residues each crystal resolved*.
- **SEQRES differs, resolved sequence does not:** 6X0H and 6X0J were crystallised with an N-terminal His-tag
  (SEQRES 494), 4B63 and 6X0I without (501). The tag is unresolved in every structure, so it never reaches a receptor.
- **Missing stretches are far from the site.** Distances from C4X to the residues flanking each gap or terminus
  (biotite): loops 70-75 and 383-394 are 29-41 A away; the resolved termini are 26-38 A away. A 22 A box centred on
  C4X reaches 11 A. So not rebuilding missing residues has no effect on the pocket.
- 6X0I chain A lacks residues 68-75, which the others (except 6X0H's 68-75 gap in the pair comparison) resolve:
  this is the loop near residue 71-74, also ~30 A from the site.
