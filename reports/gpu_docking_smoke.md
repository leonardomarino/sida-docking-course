# GPU docking arm — tooling smoke test (NOT the Step 4 redock)

Run 2026-11-04 on mj-ai (RTX 5090, compute capability 12.0, driver 620.70.06). Purpose: show that a GPU docking engine
is available, commercially licensed and fast enough to matter, before Step 4 builds a GPU arm. These runs use one ligand
and no acceptance criteria; they say nothing about whether docking reproduces the crystal pose.

## Tool and licence

**Uni-Dock 1.2.0** (conda-forge, CUDA 12.9 build, env `sida-docking-gpu`, defined in protein-binders `environments/`).
It is Vina's scoring function and search on the GPU. Licence read from the package metadata: **1.2.0 is Apache-2.0**;
versions up to 1.1.3 are LGPL-3.0-only. Both permit commercial use; 1.2.0 is the cleanest and matches Vina's Apache-2.0.
It pulls NVIDIA-EULA runtime libraries (`cuda-nvrtc`, `libcufft`, `libcurand`, `cuda-version`): commercial use permitted,
redistribution restricted. Not available on a laptop route (linux-64 + NVIDIA GPU).

Other GPU docking tools searched on conda-forge/bioconda and **not found**: AutoDock-GPU, GNINA, Vina-GPU, DiffDock. They
were not evaluated (no package to inspect), so nothing is claimed about their licences.

## Setup

Receptor: `prepared/4B63_A_nadp-keep/receptor.pdbqt` (Step 2). Ligand: L-ornithine from the PDB component dictionary
(biotite), protonated at pH 7.4 by Open Babel (net charge +1: both amines protonated, carboxylate negative), PDBQT by Meeko
(15 atoms, 6 torsions). Box 22 x 22 x 22 A centred on FAD C4X (-25.372, -31.449, 14.541). Identical settings for both
engines: Vina scoring, exhaustiveness 32, 20 modes, seed 1. The GPU run went through `gpu-pipeline-wrapper` as marino.

## Results

**One ligand.** Vina 1.2.7 on the CPU (24 threads): **1.6 s**, top score **-4.922** kcal/mol. Uni-Dock on the GPU: top score
**-4.914** (0.01 kcal/mol apart, as expected for the same scoring function); its own batch timer reports 8.3 s, almost all
GPU start-up, and 16.7 s wall including the wrapper stopping and restoring `llama-swap`. For a single ligand the CPU wins.
The two engines' mode lists differ slightly because the stochastic searches differ.

**200 copies of the same ligand** (throughput; every copy and seed identical):

| engine | total | per ligand | notes |
|---|---|---|---|
| Vina 1.2.7, 24 parallel single-thread runs | 185.7 s | 0.93 s | 200 outputs |
| Uni-Dock 1.2.0, GPU | kernel 32.0 s, batch 32.6 s, wall 40.9 s incl. wrapper | **0.16 s** | 200 outputs; peak VRAM 1,470 MiB; GPU utilisation 100% |

So **about 5.7x faster than a well-parallelised 24-core CPU for this batch** (185.7 s / 32.6 s), at about 1.5 GB of VRAM, which
fits alongside anything the wrapper leaves running.

## Limits of these numbers

- One small ligand repeated, not a library: a larger or more flexible ligand costs more on both engines, and I have not
  measured how the ratio moves. A real library would also amortise the ~8 s start-up.
- No pose-quality comparison against the crystal ligand was made; that is Step 4's redock, with its own protocol.
- The wrapper's stop/restore overhead (about 8 s here) applies to every GPU job on this host.
- The GPU arm is an addition for this host. Steps 1-5 of the teaching pipeline stay CPU-only so students can run them.
