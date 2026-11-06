"""The README is the installation contract: keep it consistent with the files it describes."""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
README = (REPO / "README.md").read_text(encoding="utf-8")
NB_README = (REPO / "notebooks" / "README.md").read_text(encoding="utf-8")
ENV = (REPO / "notebooks" / "environment.yml").read_text(encoding="utf-8")

INSTALL_COMMANDS = (
    "conda env create -f notebooks/environment.yml",
    "conda activate sida-screen",
    "pip install -e .",
    "jupyter lab screen.ipynb",
    "conda env remove -n sida-screen",
    'pytest -m "not slow" tests/test_screen.py',
)


def slug(heading):
    return re.sub(r"\s", "-", re.sub(r"[^\w\s-]", "", heading.strip().lower()))


def test_both_readmes_give_the_same_install_commands():
    for cmd in INSTALL_COMMANDS:
        assert cmd in README and cmd in NB_README, cmd


def test_readme_versions_match_the_environment_file():
    assert "python=%s" % re.search(r"python=([0-9.]+)", ENV).group(1) in README.replace("python-", "python=")
    assert "Vina%%20%s" % re.search(r"vina=([0-9.]+)", ENV).group(1) in README
    for pin in re.findall(r"(?:pdb2pqr|propka)==[0-9.]+", ENV):
        assert pin in README, pin


def test_in_page_links_resolve_and_files_exist():
    heads = {slug(m.group(2)) for m in re.finditer(r"^(#{1,6})\s+(.*)$", README, re.M)}
    assert [a for a in re.findall(r"\]\(#([^)]+)\)", README) if a not in heads] == []
    for target in re.findall(r"\]\((?!http|#)([^)]+)\)", README):
        assert (REPO / target).exists(), target


def test_layout_tree_lists_only_existing_entries():
    tree = re.search(r"```text\n(.*?)```", README, re.S).group(1)
    for name in re.findall(r"[├└]── (\S+)", tree):
        n = name.rstrip("/")
        assert any((REPO / base / n).exists() for base in ("", "notebooks", "sidacore")), name


def test_code_fences_are_balanced():
    for text in (README, NB_README):
        assert len(re.findall(r"^\s*```", text, re.M)) % 2 == 0


def test_licence_file_matches_the_declared_licence():
    assert (REPO / "LICENSE").read_text(encoding="utf-8").startswith("MIT License")
    assert 'license = { text = "MIT" }' in (REPO / "pyproject.toml").read_text(encoding="utf-8")


def test_no_stale_branch_or_private_repository_wording():
    for text in (README, NB_README):
        assert not re.search(r"feat/sida-docking-prototype|collaborator|\*\*private\*\*", text, re.I)
