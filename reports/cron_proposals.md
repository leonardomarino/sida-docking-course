# Proposal for `/root/cron/update-conda-envs.sh` (NOT applied)

Drafted 2026-11-04 at the user's request. Nothing under `/root` was modified. `/root/cron` is not under
version control, so today the script has no history; if preferred this belongs in `mj-ai-config`.

Patch: `reports/update-conda-envs.sh.patch` (2 hunks). Verified: `bash -n` syntax ok; `patch --dry-run` applies
cleanly to the live script; the live script was byte-identical before and after.

    cd /root/cron && patch update-conda-envs.sh < /opt/sida-docking-course/reports/update-conda-envs.sh.patch

## Smoke-test `sida-docking` from the weekly cron

Adds the env to the "ML env sanity (no updates)" section using the script's own `smoke_py` / `smoke_cli`
helpers: rdkit, vina, meeko, openmm, biotite, pdbfixer imports plus `vina --version` and `obabel -V`.
The env is deliberately NOT added to any update list: it is frozen and exact-pinned.

**Tested** by running the script's real helper functions plus the proposed lines as root: 8/8 PASS.

**The `PYTHONDONTWRITEBYTECODE=1` guard is load-bearing.** The cron runs as root, and root can write into the
read-only env. A negative control without the guard wrote 5 `.pyc` files into `/opt/miniforge3/envs/sida-docking`
(removed; the env was re-verified byte-identical to the committed lock). With the guard: 0 files changed.

**Not tested**: the live script end to end (it runs real `conda update`s). First scheduled run after applying is
Sunday 2026-11-08 06:00; check `/var/log/update-conda-envs.log` for the `sida-docking` PASS lines.

## Withdrawn: running the cron's git as marino

An earlier draft changed the cron's git handling and required `chown -R marino:marino /opt/protein-optimization`.
**Decision (user, 2026-11-04): `/opt/protein-optimization` stays root-owned** -- it is the authoritative
repo and holds production pieces -- so root running git there is by design, not a defect. Changes reach it the
normal way: committed in a marino-owned dev clone, pushed from the Mac, pulled into the production checkout.
Observation worth keeping, not a request: the script never checks the exit status of its push, so a broken key
would only show up as a repo that quietly stops updating.
