# Environment

**The environment is defined once, in `leonardomarino/protein-binders`, `environments/`** (single
source of truth, one convention for every env on mj-ai):

- `environments/sida-docking.yml` — `conda env export` snapshot (exact versions + builds)
- `environments/sida-docking.linux-64.explicit.txt` — exact reproduction (URLs + md5)
- `environments/patches/sida_docking_caches.sh` — cache hook (not captured by `conda env export`)
- `environments/README.md`, note "sida-docking" — why each pin exists

On mj-ai it is `/opt/miniforge3/envs/sida-docking` (root-owned, read-only for marino).
Boltz-2 (Step 6) uses the existing `boltz2` env from the same place, not a new one.
Linux is the authoritative platform; macOS is a client only (see that README note).

## For students: a generated copy

Students without access to protein-binders get a **generated copy** in this directory at release
time. Do not edit it by hand:

```bash
python -m tools.sync_env --source /path/to/protein-binders     # writes env/ + env/SOURCE.txt
```

It refuses (writes nothing) unless every source file is tracked and committed, stamps the copy
with the source commit and per-file SHA-256, and removes the `prefix:` line from the `.yml` so
`conda env create -f` builds at the user's own location instead of `/opt/miniforge3/...`.
Until the first release this directory holds only this file.
