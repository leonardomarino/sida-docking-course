# Environment

**The environment is not defined in this repo.** All conda environments are defined in
`leonardomarino/protein-binders`, `environments/` (single source of truth, one convention):

- `environments/sida-docking.yml` — `conda env export` snapshot (exact versions + builds)
- `environments/sida-docking.linux-64.explicit.txt` — exact reproduction (URLs + md5)
- `environments/patches/sida_docking_caches.sh` — cache hook (not captured by `conda env export`)
- `environments/README.md`, note "sida-docking" — why each pin exists

On mj-ai it is `/opt/miniforge3/envs/sida-docking` (root-owned, read-only for marino).
Boltz-2 (Step 6) uses the existing `boltz2` env from the same place, not a new one.

Linux is the authoritative platform; macOS is a client only (see that README note).
