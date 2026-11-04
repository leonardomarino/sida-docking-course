# Proposals for `/root/cron/update-conda-envs.sh` (NOT applied)

Drafted 2026-11-04 at the user's request. Nothing under `/root` was modified. These concern server
infrastructure rather than the course; if preferred they belong in `mj-ai-config` (`/root/cron` is not
under version control at all, so today the script has no history).

Patch: `reports/update-conda-envs.sh.patch` (3 hunks). Verified: `bash -n` syntax ok; `patch --dry-run`
applies cleanly to the live script; the live script was byte-identical before and after.

    cd /root/cron && patch update-conda-envs.sh < /opt/sida-docking-course/reports/update-conda-envs.sh.patch

## A. Smoke-test `sida-docking` (hunks 1-2)

Adds the env to the "ML env sanity (no updates)" section, using the script's own `smoke_py` / `smoke_cli`
helpers: rdkit, vina, meeko, openmm, biotite, pdbfixer imports plus `vina --version` and `obabel -V`.
The env is deliberately NOT added to any update list: it is frozen and exact-pinned.

**Tested** by running the script's real helper functions plus the proposed lines as root: 8/8 PASS.

**The `PYTHONDONTWRITEBYTECODE=1` guard is load-bearing.** The cron runs as root, and root can write into the
read-only env. A negative control without the guard wrote 5 `.pyc` files into `/opt/miniforge3/envs/sida-docking`
(removed; the env was re-verified byte-identical to the committed lock). With the guard: 0 files changed.

## B. Run git as marino, not root (hunk 3)

Today the script runs `git add/commit/push` as root inside a marino-owned repo. Consequences measured on
2026-11-04: `/opt/protein-optimization/.git` has 1,688 root-owned entries and the working tree has 687
root-owned files, so marino cannot commit or branch there (`insufficient permission for adding an object`,
`cannot lock ref`). Also, a bare `git commit` takes whatever is staged, so anything staged in that checkout
before Sunday 06:00 would be swept into the automated commit.

The patch: commits **as marino** by **pathspec** (only the two yml files), and pushes **as root by URL**.
Root must still push because **marino has no GitHub key** (`ssh -T git@github.com` as marino: Permission denied;
as root: authenticates, which is how the weekly push works). A push by URL, not `origin`, writes no `refs/remotes/*`.

**One-time root step, required before the patch can work** (marino cannot write the root-owned index):

    find /opt/protein-optimization ! -user marino | wc -l      # look first: expect ~2,400
    chown -R marino:marino /opt/protein-optimization

**Tested** in a scratch repo made by root (root-owned `.git`, marino-owned directories, like the real state):
- reproduced the real failure: marino's commit -> `insufficient permission for adding an object to repository database`
- after the one-time chown, the proposed block run as root: commit authored by marino containing only
  `feasibility.yml`; an unrelated staged file was left staged, not swept in; the push reached the (local, bare)
  remote; **0 root-owned entries** afterwards; a second run with no changes does nothing.
- on the REAL repo, read-only: `git push --dry-run <url> master` as root authenticates ("Everything up-to-date")
  and left the root-owned count unchanged at 1,688.

**Not tested**: the live script end to end (it runs real `conda update`s). First scheduled run after applying is
Sunday 2026-11-08 06:00; check `/var/log/update-conda-envs.log` for `Pushed:` and that
`find /opt/protein-optimization ! -user marino` is still 0.

**Caveats**
- Any root session that runs git in that checkout re-pollutes it. Policy, not code.
- Root still needs its `safe.directory` entry for `/opt/protein-optimization` (it already has one, among ~20).
- **Option B2, your call:** register marino's `~/.ssh/id_ed25519.pub` with GitHub. Then the whole job (and every
  marino commit on mj-ai) could run as marino, the root `safe.directory` entries for marino's repos could be
  removed, and "awaiting a push from the Mac" disappears. The cost is a private key on an internet-facing host,
  which is presumably why it is not registered; note root's key already lives there and already pushes.
- Observation, not in the patch: a failed push is silent today (the script never checks the push's exit status),
  so a broken key would only show up as a repo that quietly stops updating.
