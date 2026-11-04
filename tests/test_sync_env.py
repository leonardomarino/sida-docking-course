"""Release-time copy of the env definition from protein-binders into this repo's env/."""
import hashlib
import subprocess
from pathlib import Path

import pytest

from tools.sync_env import FILES, SyncError, sync_env

YML = "name: sida-docking\nchannels:\n  - conda-forge\ndependencies:\n  - vina=1.2.7=h1\nprefix: /opt/miniforge3/envs/sida-docking\n"
LOCK = "@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/vina-1.2.7-h1.conda#abc\n"
HOOK = "export MPLCONFIGDIR=x\n"


def _git(repo, *args):
    return subprocess.run(
        ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-C", str(repo), *args],
        check=True, capture_output=True, text=True,
    ).stdout.strip()


@pytest.fixture
def source(tmp_path):
    repo = tmp_path / "protein-binders"
    (repo / "environments" / "patches").mkdir(parents=True)
    (repo / "environments" / "sida-docking.yml").write_text(YML)
    (repo / "environments" / "sida-docking.linux-64.explicit.txt").write_text(LOCK)
    (repo / "environments" / "patches" / "sida_docking_caches.sh").write_text(HOOK)
    _git(repo, "init", "-q")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "init")
    return repo


def test_copies_files_and_strips_the_prefix_line(source, tmp_path):
    dest = tmp_path / "env"
    sync_env(source, dest)
    out = (dest / "sida-docking.yml").read_text()
    assert "prefix:" not in out                       # would make conda write to /opt on a student machine
    assert out == YML.replace("prefix: /opt/miniforge3/envs/sida-docking\n", "")
    assert (dest / "sida-docking.linux-64.explicit.txt").read_text() == LOCK
    assert (dest / "sida_docking_caches.sh").read_text() == HOOK


def test_records_source_commit_and_hashes(source, tmp_path):
    dest = tmp_path / "env"
    sync_env(source, dest)
    rec = (dest / "SOURCE.txt").read_text()
    assert _git(source, "rev-parse", "HEAD") in rec
    assert hashlib.sha256(LOCK.encode()).hexdigest() in rec
    assert "prefix line removed" in rec


def test_refuses_modified_source_file(source, tmp_path):
    (source / "environments" / "sida-docking.yml").write_text(YML + "# edited\n")
    with pytest.raises(SyncError, match="uncommitted"):
        sync_env(source, tmp_path / "env")
    assert not (tmp_path / "env").exists()            # nothing half-written


def test_refuses_untracked_source_file(source, tmp_path):
    _git(source, "rm", "-q", "--cached", "environments/sida-docking.linux-64.explicit.txt")
    _git(source, "commit", "-q", "-m", "untrack")
    with pytest.raises(SyncError, match="not tracked"):
        sync_env(source, tmp_path / "env")


def test_refuses_missing_file_and_non_repo(source, tmp_path):
    (source / "environments" / "patches" / "sida_docking_caches.sh").unlink()
    with pytest.raises(SyncError):
        sync_env(source, tmp_path / "env")
    plain = tmp_path / "plain"; plain.mkdir()
    with pytest.raises(SyncError, match="git"):
        sync_env(plain, tmp_path / "env2")


def test_resync_overwrites_generated_files(source, tmp_path):
    dest = tmp_path / "env"
    sync_env(source, dest)
    (dest / "sida-docking.yml").write_text("stale")
    sync_env(source, dest)
    assert "stale" not in (dest / "sida-docking.yml").read_text()


def test_file_list_is_exactly_the_three_env_files():
    assert [Path(f).name for f in FILES] == [
        "sida-docking.yml", "sida-docking.linux-64.explicit.txt", "sida_docking_caches.sh"]


def test_cli_success_and_refusal_exit_codes(source, tmp_path, capsys):
    from tools.sync_env import main
    assert main(["--source", str(source), "--dest", str(tmp_path / "env")]) == 0
    (source / "environments" / "sida-docking.yml").write_text(YML + "# edited\n")
    assert main(["--source", str(source), "--dest", str(tmp_path / "env2")]) == 1
    assert "refused" in capsys.readouterr().err
