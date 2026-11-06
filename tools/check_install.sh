#!/usr/bin/env bash
# Check, from scratch and on this machine, that the documented installation works.
#
#   tools/check_install.sh            # follow the README steps in a new throwaway environment, then run every test
#   tools/check_install.sh --quick    # same, without the slow tests (docking and the whole notebook)
#   tools/check_install.sh --keep     # keep the throwaway folder even when everything passes
#
# It uses a copy of this working tree, a new conda environment in a temporary folder and a fresh package cache, so it never touches your
# existing environments or this repository, and nothing left over from an earlier installation can hide a problem. It needs conda and
# git (README, Step 1) and an internet connection. The same steps run in GitHub Actions (.github/workflows/ci.yml).
# Exit status: 0 everything passed, 1 a step failed, 2 conda or git missing / bad option.

# shellcheck disable=SC2329   # the functions below are run indirectly, through step()
set -o pipefail

QUICK=0
KEEP=0
for arg in "$@"; do
  case "$arg" in
    --quick) QUICK=1 ;;
    --keep) KEEP=1 ;;
    -h|--help) sed -n '2,11p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $arg (try --help)" >&2; exit 2 ;;
  esac
done

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
command -v git >/dev/null 2>&1 || { echo "git not found: see README, Step 1." >&2; exit 2; }
command -v conda >/dev/null 2>&1 || { echo "conda not found: see README, Step 1 (then open a new terminal)." >&2; exit 2; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/sida-check-install.XXXXXX")" || exit 2
LOG="$WORK/log.txt"
START=$SECONDS
FAILED=""

export CONDA_PKGS_DIRS="$WORK/pkgs"       # a fresh package cache
export SIDA_CI=1                          # a skipped test counts as a failure (tests/conftest.py)
export PYTHONDONTWRITEBYTECODE=1

step() {  # step "title" command...   runs the command, shows and logs its output, stops the run at the first failure
  local title="$1"; shift
  [ -n "$FAILED" ] && return 1
  echo; echo "=== $title"
  local t0=$SECONDS
  "$@" 2>&1 | tee -a "$LOG"
  local rc=${PIPESTATUS[0]}
  if [ "$rc" = 0 ]; then
    echo "--- ok ($((SECONDS - t0)) s)"
  else
    echo "--- FAILED ($((SECONDS - t0)) s)"; FAILED="$title"
  fi
  return "$rc"
}

copy_tree() {      # the working tree as it is now (tracked and new files, not ignored ones) into $WORK/repo
  mkdir -p "$WORK/repo" || return 1
  git -C "$REPO" ls-files -co --exclude-standard | while IFS= read -r f; do
    [ -f "$REPO/$f" ] || continue
    mkdir -p "$WORK/repo/$(dirname "$f")" && cp -p "$REPO/$f" "$WORK/repo/$f" || exit 1
  done
}

create_env() { ( cd "$WORK/repo" && conda env create -f notebooks/environment.yml -p "$WORK/sida-screen" ); }

install_code() { ( cd "$WORK/repo" && pip install -e . ); }

check_tools() { ( cd "$WORK/repo" && python -c "import sidacore.screen as S; t = S.check_tools(); print(t); assert all(t.values()), t" ); }

fetch_structures() {
  ( cd "$WORK/repo" && for attempt in 1 2 3; do
      python -c "from sidacore.prep import fetch_structure; [fetch_structure(i, 'structures_cache') for i in ('4B63', '6X0H', '6X0J')]" && break
      echo "download failed (attempt $attempt)"; sleep 20
    done && ls -l structures_cache )
}

quick_tests() { ( cd "$WORK/repo" && python -m pytest -p no:cacheprovider --disable-warnings -m "not slow" tests/test_screen.py tests/test_install.py tests/test_docs.py -q -rs ); }

slow_tests() { ( cd "$WORK/repo" && python -m pytest -p no:cacheprovider --disable-warnings -m slow tests/test_screen.py tests/test_install.py tests/test_screen_notebook.py -q -rs ); }

echo "Checking the installation of $REPO"
echo "machine: $(uname -sm), $(conda --version), working folder: $WORK"

step "Copy the working tree" copy_tree &&
step "Step 3: conda env create -f notebooks/environment.yml" create_env

# Activation must happen in this shell (step() runs its command in a pipe), so it is not a step() call.
if [ -z "$FAILED" ]; then
  echo; echo "=== Step 4: conda activate"
  if eval "$(conda shell.bash hook)" && conda activate "$WORK/sida-screen"; then
    echo "python: $(command -v python)"; echo "--- ok"
  else
    echo "--- FAILED"; FAILED="Step 4: conda activate"
  fi
fi

step "Step 5: pip install -e ." install_code &&
step "Step 6: check the installation" check_tools &&
step "Download the test structures from RCSB" fetch_structures &&
step "Quick tests" quick_tests
if [ "$QUICK" = 0 ]; then
  step "Slow tests (PDB2PQR, Meeko, Vina, the whole notebook)" slow_tests
fi

echo
if [ -z "$FAILED" ]; then
  echo "PASS: the installation works on this machine ($((SECONDS - START)) s)."
  if [ "$KEEP" = 1 ]; then
    echo "kept: $WORK"
  else
    case "$WORK" in */sida-check-install.*) rm -rf "$WORK" ;; esac
  fi
  exit 0
fi
echo "FAIL at: $FAILED ($((SECONDS - START)) s)."
echo "The folder was kept for inspection: $WORK"
echo "Full output: $LOG"
exit 1
