"""Execute notebooks/screen.ipynb end to end and compare its table with the committed reference run."""
import re
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
NB = REPO / "notebooks" / "screen.ipynb"
CACHE = REPO / "structures_cache"


def _have(name):
    return bool(shutil.which(name)) or (Path(sys.executable).parent / name).exists()


@pytest.mark.slow
@pytest.mark.skipif(not (NB.exists() and (CACHE / "4B63.cif").exists() and (CACHE / "6X0J.cif").exists()), reason="notebook or cached structures missing")
@pytest.mark.skipif(not (_have("vina") and _have("pdb2pqr30")), reason="vina / pdb2pqr30 not available")
def test_the_notebook_runs_and_reproduces_the_reference_table(tmp_path):
    import nbformat
    import pandas as pd
    from nbclient import NotebookClient

    nb = nbformat.read(NB, as_version=4)
    work = tmp_path / "work"
    for c in nb.cells:
        if "parameters" in c.metadata.get("tags", []):
            src = re.sub(r"(?m)^WORK_DIR = .*$", "WORK_DIR = Path(r'%s')" % work, c.source)
            src = re.sub(r"(?m)^CACHE_DIR = .*$", "CACHE_DIR = Path(r'%s')" % CACHE, src)
            assert src != c.source
            c.source = src
    NotebookClient(nb, timeout=900, resources={"metadata": {"path": str(NB.parent)}}).execute()

    got = pd.read_csv(work / "results.csv").set_index(["target", "ligand"])
    ref = pd.read_csv(NB.parent / "reference_results.csv").set_index(["target", "ligand"])
    assert set(got.index) == set(ref.index)
    got, ref = got.loc[ref.index], ref
    # Vina is deterministic per seed on one machine, not bit-identical across machines: an Apple-silicon Mac differed from this reference by at most
    # 0.02 kcal/mol in score and by one clash, with identical readings. The tolerances leave room for that and no more.
    assert (got["best_score"] - ref["best_score"]).abs().max() < 0.1
    assert ((got["n_clashes"] - ref["n_clashes"]).abs() <= 2).all()
    assert (got["reading"] == ref["reading"]).all()
    stderr = [o for c in nb.cells if c.cell_type == "code" for o in c.get("outputs", []) if o.get("output_type") == "stream" and o.get("name") == "stderr"]
    assert stderr == []                                                       # no tracebacks or warnings leaking into the notebook
