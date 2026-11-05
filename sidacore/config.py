"""Load, validate and record the single versioned config (Step 1b).

Strict by design: unknown keys, missing keys and out-of-range values raise ConfigError.
There are no silent defaults -- the YAML file is the only source of values.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, fields
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping

import yaml
from jsonschema import Draft202012Validator

SCHEMA_PATH = Path(__file__).with_name("config_schema.json")
RESOLVED_NAME = "resolved_config.yaml"


class ConfigError(ValueError):
    """The config is invalid or cannot be read. Never downgraded to a default."""


@dataclass(frozen=True)
class Structures:
    benchmark_primary: str
    benchmark_secondary: str
    fad_out: str
    fad_in: str


@dataclass(frozen=True)
class Arm:
    name: str
    structure: str
    nadp: str


@dataclass(frozen=True)
class Config:
    config_version: int
    protonation_ph: float
    box_center_atom: str
    box_size_angstrom: tuple
    exhaustiveness: int
    num_modes: int
    seeds: tuple
    waters: str
    alt_conformers: str
    redock_rmsd_acceptance: float
    ligand_resname: str
    structures: Structures
    redock_nadp: str
    cross_docking_arms: tuple
    receptor_chain: str
    ions: str
    pocket_radius_angstrom: float
    superposition_outlier_threshold: float
    pocket_agreement_offset_angstrom: float
    site_contact_angstrom: float
    overrides: Mapping[str, Any]
    source_path: str
    source_sha256: str

    def to_plain(self) -> dict:
        """Parameter values only (no provenance), as plain YAML-safe types."""
        skip = {"overrides", "source_path", "source_sha256"}
        out = {}
        for f in fields(self):
            if f.name in skip:
                continue
            v = getattr(self, f.name)
            if f.name == "structures":
                v = {g.name: getattr(v, g.name) for g in fields(v)}
            elif f.name == "cross_docking_arms":
                v = [{g.name: getattr(a, g.name) for g in fields(a)} for a in v]
            elif isinstance(v, tuple):
                v = list(v)
            out[f.name] = v
        return out


def _schema_errors(data: dict) -> list:
    schema = json.loads(SCHEMA_PATH.read_text())
    errs = sorted(Draft202012Validator(schema).iter_errors(data), key=lambda e: list(e.path))
    return ["%s: %s" % (".".join(str(p) for p in e.path) or "<root>", e.message) for e in errs]


def _cross_checks(data: dict) -> list:
    errs = []
    known = set(data["structures"].values())
    names = [a["name"] for a in data["cross_docking_arms"]]
    if len(set(names)) != len(names):
        errs.append("cross_docking_arms: names must be unique, got %s" % names)
    for a in data["cross_docking_arms"]:
        if a["structure"] not in known:
            errs.append("cross_docking_arms[%s]: structure %s is not listed under structures"
                        % (a["name"], a["structure"]))
    return errs


def load_config(path, overrides: Mapping[str, Any] | None = None) -> Config:
    path = Path(path)
    try:
        raw_bytes = path.read_bytes()
        data = yaml.safe_load(raw_bytes)
    except (OSError, yaml.YAMLError) as e:
        raise ConfigError("cannot read config %s: %s" % (path, e)) from e
    if not isinstance(data, dict):
        raise ConfigError("config %s must be a YAML mapping" % path)

    recorded = {}
    if "source_sha256" in data and "config" in data:      # a resolved_config.yaml: reload it
        recorded = dict(data.get("overrides") or {})
        data = data["config"]
    merged = {**data, **dict(overrides or {})}
    errs = _schema_errors(merged)
    if not errs:
        errs = _cross_checks(merged)
    if errs:
        raise ConfigError("invalid config %s:\n  " % path + "\n  ".join(errs))

    return Config(
        **{k: merged[k] for k in ("config_version", "protonation_ph", "box_center_atom",
                                  "exhaustiveness", "num_modes", "waters", "alt_conformers",
                                  "redock_rmsd_acceptance", "ligand_resname", "redock_nadp",
                                  "receptor_chain", "ions", "pocket_radius_angstrom",
                                  "superposition_outlier_threshold", "pocket_agreement_offset_angstrom",
                                  "site_contact_angstrom")},
        box_size_angstrom=tuple(float(x) for x in merged["box_size_angstrom"]),
        seeds=tuple(merged["seeds"]),
        structures=Structures(**merged["structures"]),
        cross_docking_arms=tuple(Arm(**a) for a in merged["cross_docking_arms"]),
        overrides=MappingProxyType({**recorded, **dict(overrides or {})}),
        source_path=str(path),
        source_sha256=hashlib.sha256(raw_bytes).hexdigest(),
    )


def write_resolved(config: Config, run_dir) -> Path:
    """Write the resolved config next to a run's outputs. Refuses to overwrite."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    out = run_dir / RESOLVED_NAME
    doc = {
        "source_path": config.source_path,
        "source_sha256": config.source_sha256,
        "overrides": dict(config.overrides),
        "config": config.to_plain(),
    }
    try:
        with open(out, "x") as fh:                       # 'x': fail if it already exists
            yaml.safe_dump(doc, fh, sort_keys=False)
    except FileExistsError as e:
        raise ConfigError("resolved config already exists, refusing to overwrite: %s" % out) from e
    return out
