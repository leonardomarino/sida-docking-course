# Step 2 — receptor preparation

Everything below was produced by one command, on CPU, in 32 s for 8 receptors (about 4 s each):

    python -m sidacore.prep --config config/default.yaml --out prepared

`sidacore.prep` is thin glue over installed tools; it implements no chemistry or geometry of its own.

| Stage | Tool (version) |
|---|---|
| read mmCIF, alternate conformers, atom selection, every distance and gap | biotite 1.6.0 |
| protein hydrogens and pKa at `protonation_ph` | PDB2PQR 3.7.1 + PROPKA 3.5.1 (the `af2` env's binaries) |
| receptor -> PDBQT, cofactor templates, polar H, charges, Vina box | Meeko 0.8.0 |

Each receptor directory (`prepared/<ID>_<chain>_nadp-<keep|strip>/`) holds `receptor.pdbqt`, `receptor.box.txt`,
`receptor_prepared.pdb`, `ligand_crystal.pdb` (the substrate pose, for the redocking RMSD), `removal_log.csv`
(every removed atom with a reason), `protonation_states.csv`, `prep_report.json` (versions, exact commands,
checksums, counts) and `resolved_config.yaml`. A directory that already exists is never overwritten.

## Policies (all set in `config/default.yaml`)

- **FAD stays in the receptor.** It is a cofactor, not the ligand, and is never co-docked. Resolved from the config
  atom `FAD:C4X`; `FAD` also matches the reduced flavin `FDA` (6X0J).
- **The substrate (ORN) is removed from the receptor** and written to `ligand_crystal.pdb`.
- **NADP:** kept or stripped per arm (`nadp`). Stripping it leaves an empty cavity: it does **not** remove the
  FAD-conformation/NADP confound between 6X0H and 6X0I (Amendment 1).
- **Waters:** stripped by default. `keep-within-pocket` (within `pocket_radius_angstrom` = 8 A of C4X) was run end to
  end and the waters reach the PDBQT (5 in 4B63, 1 in 6X0I).
- **Ions:** stripped by default (the Ca2+ sits 15-30 A from the site in 6X0H/6X0I and was crystallisation additive);
  `ions=keep` was run end to end and the ion reaches the PDBQT.
- **Other heterogens** (acetate, PEG, glycerol, sulfate): removed as crystallisation additives, each logged.
- **Alternate conformers:** highest occupancy, chosen by biotite per atom. Checked: no residue in any of the four
  structures ends up with mixed conformer labels. 6X0I has 20 residues with alternates (131 atoms removed over its four
  chains, including Met 101); 4B63 has 4.
- **Missing residues are not rebuilt.** All gaps are 29-41 A from the site (see `construct_check.md`).
- **Chain:** `receptor_chain` = A. 6X0I chain C has a split FAD (0.47/0.53 occupancy), so it is avoided. All four
  6X0J chains are prepared (one redock each).

## Results (default config, pH 7.4)

| receptor | site (C4X, A) | FAD / NADP charge | Meeko-deleted | PDB2PQR rebuilt (atoms / residues) | nearest rebuilt residue to site |
|---|---|---|---|---|---|
| 4B63 A, NADP kept | -25.4, -31.4, 14.5 | -2 / -3 | 0 | 103 / 29 | Lys100 6.3 A |
| 6X0J A..D, NADP kept | per chain | FDA -2 / -3 | A: Asp177, Ser179 (43, 47 A); B-D: 0 | 114-134 / 29-32 | Lys100 6.6-6.7 A |
| 6X0H A (as deposited) | 112.9, 2.5, 45.5 | -2 / none | 0 | 97 / 25 | Asn323 6.8 A |
| 6X0I A, NADP kept | 112.2, -0.1, 47.2 | -2 / -3 | 0 | 57 / 15 | 13.7 A |
| 6X0I A, NADP stripped | 112.2, -0.1, 47.2 | -2 / none | 0 | 57 / 15 | 13.7 A |

Site coordinates are in each crystal's own frame; Step 3 superposes the conformers before comparing them.

### Ionizable residues of the substrate site (state; PROPKA pKa)

The site is every protein residue with an atom inside the docking box (the 22 A cube about C4X) plus any residue within
`site_contact_angstrom` = 4.5 A of the deposited ligand (benchmark structures only). That gives 19-24 ionizable residues per receptor;
the full tables are `prepared/<receptor>/protonation_states.csv`, with `in_box` and `contacts_ligand` columns. The docking inputs
(receptor PDBQT, box) are byte-identical to those made with the earlier 8 A-sphere table: only this report changed.

- **Residues that touch the deposited ornithine:** Lys107 in 4B63 and all four 6X0J chains (pKa 9.1-9.3, charged); Asp288 in the four 6X0J
  chains but not in 4B63 (pKa 3.2-3.9, ionized). Both are in ordinary charged states in every receptor.
- **His91** is HID everywhere (pKa 4.1-4.5, so neutral; PROPKA's histidine pKa is no better than assuming 6.5, treat as uncertain). **His228** differs
  by structure (HID in 4B63, 6X0H and 6X0J chain B; HIE elsewhere). Glu260 has pKa 4.2 in 4B63 and 6.0-6.3 in the 6X0 structures; it is ionized
  at pH 7.4 in all and neutral at pH 6.0 in 6X0I and 6X0J.

**Retraction.** An earlier version of this section listed three to four ionizable residues from an 8 A sphere about C4X and concluded that the
site residues differ between the conformers (Tyr324, Tyr407, Arg409 in 6X0H against Lys100, Arg144, Arg328 in 6X0I). That was an artefact of which
residues sat at the edge of the sphere: with the corrected site, Lys100, Tyr324, Tyr407 and Arg409 are in every receptor. The sphere was centred on
the cofactor, but the substrate sits 7.7-8.0 A from C4X, so it also missed Lys107 and Asp288.

**What does differ between conformers (observation, not a conclusion):** the box centred on 6X0H's own C4X does not contain Lys107, Asp108 or
Glu146, which are in the 6X0I box; the box carried from 6X0I into 6X0H does contain them (and has the same composition as 6X0I's own box). The
carried box is therefore the comparable one for 6X0H. The two structures also differ in NADP occupancy, so the comparison cannot isolate FAD
position.

## Caveats that Step 4's diagnosis must keep in view

1. **Lys100's side chain is modelled, not observed, in 4B63 and all four 6X0J chains** (PDB2PQR rebuilt CE/NZ, or
   CG/CD/CE/NZ, 6.3-6.7 A from the site); in 6X0H it is Asn323. A modelled, charged group lining the pocket is directly
   relevant to docking a cationic amino-acid substrate. 6X0I has nothing rebuilt within 13.7 A. If a redock fails,
   this is the first place to look under "receptor protonation / side chains".
2. **His91 is HID (neutral) in every receptor**, with PROPKA pKa 4.1-4.5 (model 6.5). PROPKA's absolute histidine
   pKa is no better than assuming 6.5 (Hogues et al. 2025), and PDB2PQR's HID-versus-HIE choice is a tautomer
   assignment, so treat His91's state as uncertain. Every histidine row carries this note in `protonation_states.csv`.
3. **Cofactor protonation does not follow `protonation_ph`.** Meeko builds the FAD/FDA/NADP templates from the PDB
   component dictionary and gives FAD -2 (8 polar H), FDA -2 (10 polar H), NADP -3 (7 polar H); only the protein
   changes with pH (tested: pH 4.0 vs 7.4 changes site residue states). The reduced flavin FDA keeps its N1-H and N5-H as
   in the dictionary; whether N1 is deprotonated in the enzyme is not resolved.
4. **PDB2PQR rebuilds truncated side chains** (57-134 atoms per receptor, all listed in `prep_report.json`).
   It is a heuristic placement with no energy minimisation; in 6X0J chain A it put Ser179's hydroxyl H 1.14 A from the
   Asp177 carboxylate (a clash).

## Tool behaviours found while building this (each cost real time; the code now guards or documents them)

- **`FAD:C4A` silently selects the wrong atom.** FAD has both `C4A` (the adenine carbon) and `C4X` (isoalloxazine C4a,
  the hydroxylation site), 18.0 A apart in 4B63. The prompt's default would have centred the box on the adenine without
  an error. The default is now `FAD:C4X` and `resolve_site_center` refuses a flavin atom not bonded to N5.
- **biotite's `filter_amino_acids` counts L-ornithine as an amino acid**, so the substrate would have stayed in the
  "protein". The protein filter is `filter_canonical_amino_acids`; a test asserts ORN never reaches the receptor.
- **Meeko refuses residues it thinks have extra bonds** (here the Ser179 clash). `--delete_bad_res_from_box_radius 8`
  deletes such residues only if their centre is more than 8 A outside the box and raises for any nearer one, so
  failures near the site stay loud and far ones are logged (6X0J-A: Asp177, Ser179).
- **Meeko downloads cofactor templates from rcsb.org** at run time. `sidacore/data/meeko_cofactor_templates.json`
  (12.7 KB: FAD, FDA, NAP only) makes it offline. Verified: with network access blocked, output was byte-identical to a
  clean online run for 4B63 and 6X0J. **Do not use Meeko's full `--cache_templates` file in its place:** loading it
  replaces Meeko's own defaults and swapped proline CB/CD atom names in 25 residues (names only; charges identical).
- **PDB2PQR prints "Missing atom" warnings** that are not errors: they announce atoms it is about to model.

## Not tested / not done

- Vina has not yet consumed these receptors (Step 4). Meeko's PDBQT and box file are the hand-off.
- Waters and ions were run through the pipeline but not through a docking run.
- The student-laptop route needs PDB2PQR/PROPKA in the pinned env; today they come from the `af2` env on this host.
