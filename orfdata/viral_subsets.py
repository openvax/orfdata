"""Read declarative viral-reference subsets from package YAML."""

from functools import lru_cache
from importlib.resources import files

import yaml

_RESOURCE = "data/viral/subsets.yaml"
_SCHEMA = "orfdata.viral_subsets.v1"


def _strings(values):
    return isinstance(values, list) and bool(values) and all(
        isinstance(value, str) and value for value in values)


@lru_cache(maxsize=1)
def _registry():
    data = yaml.safe_load(files("orfdata").joinpath(_RESOURCE).read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema") != _SCHEMA:
        raise ValueError("Unsupported viral subset metadata")
    if not isinstance(data.get("subsets"), dict) or not isinstance(data.get("aliases"), dict):
        raise ValueError("Invalid viral subset metadata")
    for name, record in data["subsets"].items():
        required = {"title", "description", "abbreviations", "source_urls", "limitations"}
        if not isinstance(name, str) or not isinstance(record, dict) or required - set(record):
            raise ValueError("Invalid viral subset record")
        if not all(isinstance(record[field], str) and record[field]
                   for field in ("title", "description", "limitations")):
            raise ValueError("Invalid viral subset text")
        if not _strings(record["abbreviations"]):
            raise ValueError("Invalid viral subset abbreviations")
        if not _strings(record["source_urls"]) or not all(
                value.startswith("https://") for value in record["source_urls"]):
            raise ValueError("Invalid viral subset sources")
    if not all(isinstance(alias, str) and isinstance(name, str) and name in data["subsets"]
               for alias, name in data["aliases"].items()):
        raise ValueError("Invalid viral subset aliases")
    return data


def viral_subsets():
    """List dictionaries describing named, reproducible viral subsets."""
    registry = _registry()
    return tuple(dict(registry["subsets"][name], name=name,
                      references=len(registry["subsets"][name]["abbreviations"]))
                 for name in sorted(registry["subsets"]))


def normalize_viral_subsets(subsets):
    """Normalize named viral subsets and their concise aliases."""
    if subsets is None:
        return ()
    subsets = (subsets,) if isinstance(subsets, str) else tuple(subsets)
    if not subsets or any(not isinstance(value, str) or not value.strip() for value in subsets):
        raise ValueError("subsets must contain one or more nonempty names")
    registry = _registry()
    normalized = tuple(dict.fromkeys(
        registry["aliases"].get(value.strip().lower(), value.strip().lower()) for value in subsets))
    unknown = sorted(set(normalized) - set(registry["subsets"]))
    if unknown:
        raise ValueError(
            f"Unknown viral subset(s) {', '.join(unknown)}; available: "
            f"{', '.join(sorted(registry['subsets']))}"
        )
    return normalized


def viral_subset_abbreviations(subsets):
    """Return GDC abbreviations in one or more named subsets."""
    names = normalize_viral_subsets(subsets)
    if not names:
        return ()
    records = _registry()["subsets"]
    return tuple(sorted({abbreviation for name in names for abbreviation in records[name]["abbreviations"]}))
