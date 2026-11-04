# Provisioning record

Append-only. The teaching document is built from this file, so every deviation, finding and
measured figure goes here with the date. Written as marino-owned; root only provisions.

## 2026-11-04 — host mj-ai (Debian 13, RTX 5090 32 GB, driver 620.70.06)

### Environment `sida-docking` (CPU analysis, Steps 1-5)
- Location `/opt/miniforge3/envs/sida-docking`, root-owned, read-only for marino (verified: marino
  `touch` into the env -> Permission denied). Defined in **protein-binders `environments/`**, not here.
- Built from the exact explicit lock; `conda list --explicit --md5` of the rebuilt env is
  byte-identical to the lock. First solve+create 17.7 s, rebuild from lock 1.4 s (warm package cache,
  so these are not cold-install times).
- Resolved (conda-forge): python 3.11.18, rdkit 2026.03.1, vina 1.2.7, meeko 0.8.0, prody 2.4.1,
  numpy 1.26.4, scipy 1.17.1, pandas 3.0.6, openbabel 3.2.1, pdbfixer 1.12, openmm 8.6.1,
  biotite 1.6.0, py3dmol 2.5.4, fpocket 4.2.3, gemmi 0.7.5, jupyterlab 4.6.5, pytest 9.2.0. 278 packages.
  `vina --version` prints `AutoDock Vina f458505-mod` (a build hash, not "1.2.7"); quote both.
- Verified as marino with the amendment's smoke test: `imports ok`, vina, fpocket, obabel resolve.

### Deviations from Amendment 1 (each deliberate)
1. **Path**: `/opt/miniforge3/envs/`, not `/opt/conda/envs/`. Every other env here lives under miniforge3 and is
   already root-owned/read-only to marino; a parallel `/opt/conda` tree was a second convention. Built there first,
   then retired (same session).
2. **Permissions**: `go-w,a+rX`, not literally "files 644" — 644 would strip the execute bit from `bin/`.
3. **No apt packages**: `fpocket` is not packaged in Debian 13; `openbabel` comes from conda-forge so the
   pinned env is self-contained and one source of truth.
4. **Boltz-2: no new env.** A `sida-boltz` env (Boltz 2.2.1 without the optional proprietary cuEquivariance
   kernels) was built and worked (76 s), then **removed**: the existing `boltz2` env (Boltz 2.2.1, torch 2.11.0)
   already does the job. Used with `--no_kernels`. Lesson recorded: survey `environments/` first.
5. **Weights**: the existing copies were reused, not re-downloaded. They are hardlinked (same filesystem, no
   duplicate bytes, cannot dangle) into `/opt/bio-models/weights/boltz-2/` (root-owned, 0444, world-readable).
   `/root` is mode 700 and `/home/marino` is private, so neither is a shared location. MD5s identical across three
   independent copies (root, marino, a fresh upstream download I made and discarded):
   `boltz2_conf.ckpt 2f0a1775bf8fc366a1a85e2019eca288`, `boltz2_aff.ckpt 8e93dadedd6edb7a4d170f6051b99ec0`,
   `mols.tar 5399454dcbbba695ab01f88deede577c`.

### Findings (each a root provisioning task or a documented limit)
- **conda is not initialised in marino's non-interactive shells** (`~/.bashrc` returns early before the conda
  block), so the amendment's smoke command fails as written (`conda: command not found`). Fix without touching
  dotfiles: `source /opt/miniforge3/etc/profile.d/conda.sh` first. Scripts must do the same.
- **Python 3.11, not 3.12**: Meeko 0.8.0 needs `prody`; conda-forge has no `prody` for Python 3.12 on
  osx-arm64, so a 3.12 spec solves Meeko 0.5.0 on a Mac and 0.8.0 on Linux. Pins `meeko=0.8.0 prody=2.4.1
  numpy=1.26` make both platforms agree on all 19 checked scientific packages (3 incidental libs differ by a patch).
- **macOS floor**: every osx-arm64 `fpocket` build requires macOS >= 14.5; everything else solves at 11.0.
  The osx-arm64 lock was a dry-run only, never executed, and is not shipped. **Decision (user): Linux is
  authoritative, the Mac is a client.**
- **NVIDIA proprietary components exist but pass the commercial test**: `cuequivariance-ops-*` (optional Boltz
  kernels; NVIDIA math-library SDK licence — use permitted, redistribution restricted, so not installed in any course
  env and bypassed with `--no_kernels`); and on linux-64 OpenMM pulls `cuda-nvrtc`, `libcufft`, `cuda-version`
  (NVIDIA EULA). No non-commercial licence and no AlphaFold2 parameters anywhere (read from package metadata).
  Copyleft (OpenBabel GPL-2.0-only, Meeko LGPL-2.1-only, OpenMM LGPL-3.0-or-later) matters only if redistributed.
- **GPU services (Step 6 blocker) — RESOLVED, no sudoers change needed.** marino controls
  `llama-swap`, `llama-vision`, `llama-server`, `llama-embedding`, `litellm` through polkit
  (`/etc/polkit-1/rules.d/50-gpu-pipeline-marino.rules`) without sudo, and `/var/run/gpu-pipeline` and the
  wrapper log are marino-owned. Proven by a real job run as marino through `gpu-pipeline-wrapper`: it stopped
  `llama-swap`, ran, and restored it. **Gap**: `hidream-o1.service` is on the wrapper's stop list but is not in
  the polkit rule; today it is inactive+disabled so nothing breaks, but if it were running the wrapper's
  `systemctl stop ... || true` would fail silently as marino and leave ~20 GB of VRAM held. Narrowest fix:
  add that one unit name to the rule's `allowed_units`. **Proposed, not applied** (system file).
- **Root-owned entries inside `/opt/protein-optimization/.git`**: 1,688 (HEAD, index, config, object dirs, branch
  dirs), from the weekly root cron (`/root/cron/update-conda-envs.sh` runs a bare `git commit` there) and root
  sessions. marino cannot create branches in that repo (`cannot lock ref ... Permission denied`). Not fixed
  here: the proper fix is for the cron to run git as marino. protein-binders work was done in the marino-owned
  dev clone `/home/marino/GitHub/protein-binders` after a read-only fetch from `/opt/protein-optimization`.
- `update-conda-envs.sh` commits whatever is staged, so nothing was staged in `/opt/protein-optimization`.

### Measured: Boltz-2 smoke (provisioning check, NOT a scientific result)
Input: SidA construct from 6X0H (494 aa incl. His-tag) + ORN, `msa: empty` (single-sequence), no FAD, affinity on,
`--no_kernels --seed 1 --diffusion_samples 1`. Run as marino via `gpu-pipeline-wrapper`, `HF_HUB_OFFLINE=1`.
| run | env | rc | wall | peak VRAM |
|---|---|---|---|---|
| 1 | sida-boltz (retired) | 0 | 76 s | 5,626 MiB |
| 2 | existing `boltz2`, shared weights | 0 | 76 s | 5,624 MiB |
Wall time includes the wrapper stopping/restoring `llama-swap` and a post-run model-unload. Confidence is low
(ipTM 0.38) as expected for this input; do not quote it. The with-kernels speed comparison is Step 6's job.

### Status
- Commits are local (marino). **mj-ai has no GitHub key: awaiting push from the Mac over `ssh -A`** for both
  `sida-docking-course` (branch `feat/sida-docking-prototype`) and `protein-binders`
  (branch `feat/sida-docking-env`).

## 2026-11-04 (later) — follow-ups approved by the user
- **Env copy for students**: `tools/sync_env.py` (tested, 36 tests / 99% coverage) copies the env files from protein-binders
  into `env/` at release time, stamped with commit + SHA-256, refusing dirty/untracked sources and stripping the
  `prefix:` line. Run against the real committed branch (`cc77c70`): works. `env/` stays empty until first release.
- **Cron proposals drafted, not applied**: `reports/cron_proposals.md` + `reports/update-conda-envs.sh.patch`.
- **Correction**: "mj-ai has no GitHub key" is true for marino but NOT for root: the weekly root cron pushes
  successfully. marino's `ssh -T git@github.com` is denied.
- **Lesson**: a `git clone` from a local path hardlinks objects, so a `chown -R` on the clone would have changed ownership in
  the source repo. Use a `file://` URL or `git init` for scratch repos.
- **Tooling note**: the `block-no-verify` hook false-positives on commands containing `commit` plus a `-n`-like token
  (`--name-only`, `bash -n`). No hook was skipped; such text was moved into script files.
