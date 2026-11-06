"""Does the documented installation produce a working environment?

These are the checks behind the README's "Check the installation" step, made automatic: the pinned versions in
notebooks/environment.yml are what is installed, every external program starts, the package is importable from outside the
repository (so `pip install -e .` really happened), and a built wheel carries the package data. CI runs them on a fresh runner
after following the README, so a README that stops working fails the build.
"""
import importlib.metadata as md
import importlib.resources as resources
import json
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest
from packaging.version import Version

from sidacore import screen as S

REPO = Path(__file__).resolve().parents[1]
ENV_YML = REPO / "notebooks" / "environment.yml"
PYTHON_PACKAGES = ("vina", "meeko", "biotite", "rdkit", "numpy", "pdb2pqr", "propka")


def pins():
    found = {}
    for line in ENV_YML.read_text(encoding="utf-8").splitlines():
        m = re.match(r"\s+- ([A-Za-z0-9_.-]+)(?:==|=)([0-9][^=\s]*)\s*$", line)
        if m:
            found[m.group(1).lower()] = m.group(2)
    return found


def test_the_environment_file_is_parsed():
    assert set(PYTHON_PACKAGES) | {"python", "openbabel"} <= set(pins())      # guards the parser itself


@pytest.mark.parametrize("package", PYTHON_PACKAGES)
def test_installed_version_matches_the_environment_file(package):
    assert Version(md.version(package)) == Version(pins()[package])


def test_python_version_matches_the_environment_file():
    assert "%d.%d" % sys.version_info[:2] == pins()["python"]


def test_every_external_program_is_found_and_starts():
    tools = S.check_tools()
    assert all(tools.values()), "missing programs: %s" % [k for k, v in tools.items() if not v]
    for name, args in (("vina", ["--version"]), ("obabel", ["-V"]), ("pdb2pqr30", ["--version"]), ("mk_prepare_receptor.py", ["--help"])):
        run = subprocess.run([tools[name], *args], capture_output=True, text=True, timeout=120)
        assert run.returncode == 0, "%s failed to start: %s" % (name, run.stderr[-300:])
        if name == "obabel":
            assert re.search(r"Open Babel %s\b" % re.escape(pins()["openbabel"]), run.stdout + run.stderr)
    if Path(sys.prefix).name == "sida-screen":                                  # all four come from this environment, not a stray system copy
        assert {Path(p).resolve().parent for p in tools.values()} == {Path(sys.prefix).resolve() / "bin"}


def test_the_package_is_importable_from_outside_the_repository(tmp_path):
    """-I ignores the current directory and PYTHONPATH, so this only passes after `pip install -e .`."""
    run = subprocess.run([sys.executable, "-I", "-c", "import sidacore.screen as S; print(S.__file__)"],
                         cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, "run `pip install -e .` from the repository root: %s" % run.stderr[-300:]


def test_package_data_is_present_and_valid_json():
    root = resources.files("sidacore")
    assert json.loads((root / "config_schema.json").read_text(encoding="utf-8"))
    templates = [p for p in (root / "data").iterdir() if p.name.endswith(".json")]
    assert templates
    for p in templates:
        assert json.loads(p.read_text(encoding="utf-8"))


@pytest.mark.slow
def test_a_built_wheel_ships_the_package_data(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    shutil.copy(REPO / "pyproject.toml", src)
    shutil.copytree(REPO / "sidacore", src / "sidacore", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    out = tmp_path / "wheel"
    run = subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation", "-w", str(out), str(src)],
                         capture_output=True, text=True, timeout=600)
    assert run.returncode == 0, run.stderr[-500:]
    names = zipfile.ZipFile(next(out.glob("sidacore-*.whl"))).namelist()
    assert "sidacore/config_schema.json" in names
    assert any(n.startswith("sidacore/data/") and n.endswith(".json") for n in names)
