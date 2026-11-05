Hand-written. **This check does not corroborate the docking ranking, and it is too unreliable in this form to arbitrate.** It is
reported as that, not as support for either outcome.

**What it shows.** Rosetta's own top-ranked pose is often far from the deposited one (5-18 A in most receptors) and its rank correlation with the engine score ranges from -0.69 to +0.43 (dG view) and from -0.74 to +0.29 (total-energy view). It does not endorse the passing redocks either (6X0J chain A: engine top-1 at 1.0 A, Rosetta's pick at 5.0 A for Vina and 18.4 A for Uni-Dock, as docked; 4.5 and 19.9 A after Rosetta's own minimisation). Among docked poses it therefore does not
prefer the near-native one in the failing receptors (4B63: near-native pose ranked 10th-19th by dG; chain C: 9th-20th).

**Why it cannot be trusted here.**
- The baseline energies carry clashes unrelated to the ligand: the minimised crystal complex scores -337, +64, +278, -232 and -8 REU in 4B63 and 6X0J chains A, B, C and D (chain D is the unremarkable one). For the crystal pose of 6X0J chains A, B and D the shell-only minimisation changed the total by 0.0-0.1 REU (it moved 4B63 by -37 and chain C by -20 REU), so in A and B the clashes sit where the allowed degrees of freedom cannot reach them. Their cause was not investigated (candidates:
  unrelaxed backbone, Rosetta-built hydrogens, the cofactor residue types built from the component dictionary).
- dG_separated is positive for most in-site poses of this +1 ligand (buried charge), so a pose at the edge of the box with
  little interface can beat a deeply buried one.
- "The crystal pose scores better than the docked one" holds in 8 of 10 cases by both measures, including poses within 1.0 A
  of the deposited one (41 REU for 6X0J chain A). The crystal ligand is pre-organised with its own side chains, so this is
  expected for any scoring function and is not evidence about Vina.

**Provenance.** A first version of this table was wrong for 6X0J chains B-D: the interface was hard-coded against chain A, so
every dG was 0.0. That version was discarded; the chain letters are now read from the file and a batch with no contact fails.
The total-energy view was added after the dG results were seen and is labelled as such in the table.

**Licence.** PyRosetta is free for academic use and needs a licence for commercial use. This check lives in `confirm/`, which
`sidacore` never imports and which can be removed without touching the shared layer.

**Pose counts.** Uni-Dock returned 5-12 poses per seed although 20 were requested, Vina always 20, so "best RMSD in the list" and its
Rosetta rank compare lists of different length across engines and are not like for like. The RMSD columns in the table show the pose as
docked and after Rosetta's minimisation; the second can differ materially (6X0J chain D, Vina's pick: 1.00 -> 2.01 A).
