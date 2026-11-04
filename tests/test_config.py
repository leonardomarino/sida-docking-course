"""Step 1b: one versioned config, strictly validated, resolved copy written next to every run."""
import copy
import hashlib
from pathlib import Path

import pytest
import yaml

from sidacore.config import ConfigError, load_config, write_resolved

DEFAULT = Path(__file__).resolve().parent.parent / "config" / "default.yaml"


def _raw():
    return yaml.safe_load(DEFAULT.read_text())


def _write(tmp_path, data, name="c.yaml"):
    p = tmp_path / name
    p.write_text(yaml.safe_dump(data))
    return p


def test_default_config_matches_documented_defaults():
    c = load_config(DEFAULT)
    assert c.protonation_ph == 7.4
    assert c.box_center_atom == "FAD:C4A"
    assert c.box_size_angstrom == (22.0, 22.0, 22.0)
    assert c.exhaustiveness == 32
    assert c.num_modes == 20
    assert c.seeds == tuple(range(1, 11))
    assert c.waters == "strip"
    assert c.alt_conformers == "highest-occupancy"
    assert c.redock_rmsd_acceptance == 2.0


def test_default_config_encodes_amendment_1_structures_and_arms():
    c = load_config(DEFAULT)
    assert (c.structures.benchmark_primary, c.structures.benchmark_secondary) == ("4B63", "6X0J")
    assert (c.structures.fad_out, c.structures.fad_in) == ("6X0H", "6X0I")
    assert c.redock_nadp == "keep"
    arms = {a.name: (a.structure, a.nadp) for a in c.cross_docking_arms}
    assert arms == {
        "6X0H_as_deposited": ("6X0H", "keep"),
        "6X0I_with_nadp": ("6X0I", "keep"),
        "6X0I_nadp_stripped": ("6X0I", "strip"),
    }


def test_config_is_immutable():
    c = load_config(DEFAULT)
    with pytest.raises(Exception):
        c.protonation_ph = 5.0


def test_unknown_key_is_rejected_not_ignored(tmp_path):
    d = _raw(); d["exhaustivness"] = 8          # typo must not silently fall back to a default
    with pytest.raises(ConfigError, match="exhaustivness"):
        load_config(_write(tmp_path, d))


def test_missing_key_is_rejected_no_silent_default(tmp_path):
    d = _raw(); del d["protonation_ph"]
    with pytest.raises(ConfigError, match="protonation_ph"):
        load_config(_write(tmp_path, d))


@pytest.mark.parametrize("key,bad", [
    ("protonation_ph", 15.0), ("protonation_ph", -1), ("protonation_ph", "7.4"),
    ("exhaustiveness", 0), ("exhaustiveness", 8.5), ("num_modes", 0),
    ("box_size_angstrom", [22, 22]), ("box_size_angstrom", [22, 22, 0]),
    ("box_center_atom", "C4A"), ("box_center_atom", "FAD:"),
    ("seeds", []), ("seeds", [1, 1, 2]), ("seeds", [1, 2.5]),
    ("waters", "maybe"), ("alt_conformers", "random"),
    ("redock_rmsd_acceptance", 0), ("redock_nadp", "sometimes"),
])
def test_invalid_values_are_rejected(tmp_path, key, bad):
    d = _raw(); d[key] = bad
    with pytest.raises(ConfigError):
        load_config(_write(tmp_path, d))


def test_arm_must_reference_a_known_structure_and_names_must_be_unique(tmp_path):
    d = _raw(); d["cross_docking_arms"][0]["structure"] = "9ZZZ"
    with pytest.raises(ConfigError, match="9ZZZ"):
        load_config(_write(tmp_path, d))
    d = _raw(); d["cross_docking_arms"][1]["name"] = d["cross_docking_arms"][0]["name"]
    with pytest.raises(ConfigError, match="unique"):
        load_config(_write(tmp_path, d))


def test_override_is_validated_and_recorded(tmp_path):
    c = load_config(DEFAULT, overrides={"protonation_ph": 6.5})
    assert c.protonation_ph == 6.5
    assert c.overrides == {"protonation_ph": 6.5}
    with pytest.raises(ConfigError):
        load_config(DEFAULT, overrides={"protonation_ph": 99})
    with pytest.raises(ConfigError, match="no_such_key"):
        load_config(DEFAULT, overrides={"no_such_key": 1})


def test_resolved_config_records_values_overrides_and_source_hash(tmp_path):
    c = load_config(DEFAULT, overrides={"protonation_ph": 6.5})
    out = write_resolved(c, tmp_path / "run1")
    assert out.name == "resolved_config.yaml"
    r = yaml.safe_load(out.read_text())
    assert r["config"]["protonation_ph"] == 6.5
    assert r["overrides"] == {"protonation_ph": 6.5}
    assert r["source_sha256"] == hashlib.sha256(DEFAULT.read_bytes()).hexdigest()
    assert r["source_path"] == str(DEFAULT)
    # the resolved values reload to an identical config (traceable, round-trips)
    again = load_config(out)
    assert again.protonation_ph == 6.5 and again.seeds == c.seeds


def test_resolved_config_is_never_silently_overwritten(tmp_path):
    c = load_config(DEFAULT)
    write_resolved(c, tmp_path / "run1")
    with pytest.raises(ConfigError, match="exists"):
        write_resolved(c, tmp_path / "run1")


def test_malformed_yaml_and_missing_file_raise_config_error(tmp_path):
    p = tmp_path / "bad.yaml"; p.write_text("a: [unclosed")
    with pytest.raises(ConfigError):
        load_config(p)
    with pytest.raises(ConfigError):
        load_config(tmp_path / "nope.yaml")


def test_non_mapping_yaml_is_rejected(tmp_path):
    p = tmp_path / "list.yaml"; p.write_text("- just\n- a list\n")
    with pytest.raises(ConfigError, match="mapping"):
        load_config(p)
