"""Read declarative virus--tumour mechanism selectors from package YAML."""

from functools import lru_cache
from importlib.resources import files

import yaml

_RESOURCE = "data/viral/mechanisms.yaml"
_SCHEMA = "orfdata.viral_mechanisms.v1"


def _strings(values):
    return isinstance(values, list) and bool(values) and all(
        isinstance(value, str) and value for value in values)


@lru_cache(maxsize=1)
def _registry():
    data = yaml.safe_load(files("orfdata").joinpath(_RESOURCE).read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema") != _SCHEMA:
        raise ValueError("Unsupported viral mechanism metadata")
    for section in ("mechanisms", "subsets", "aliases"):
        if not isinstance(data.get(section), dict):
            raise ValueError("Invalid viral mechanism metadata")
    for name, record in data["mechanisms"].items():
        required = {"title", "description", "abbreviations", "source_urls", "limitations"}
        if not isinstance(name, str) or not isinstance(record, dict) or required - set(record):
            raise ValueError("Invalid viral mechanism record")
        if not all(isinstance(record[field], str) and record[field] for field in
                   ("title", "description", "limitations")):
            raise ValueError("Invalid viral mechanism text")
        if not _strings(record["abbreviations"]):
            raise ValueError("Invalid viral mechanism abbreviations")
        if not _strings(record["source_urls"]) or not all(
                value.startswith("https://") for value in record["source_urls"]):
            raise ValueError("Invalid viral mechanism sources")
        for setting in record.get("tumour_settings", []):
            required = {"abbreviations", "tumour_types", "viral_feature", "evidence", "source_urls"}
            if not isinstance(setting, dict) or required - set(setting):
                raise ValueError("Invalid viral tumour setting")
            if not _strings(setting["abbreviations"]) or not _strings(setting["tumour_types"]):
                raise ValueError("Invalid viral tumour setting references")
            if not all(isinstance(setting[field], str) and setting[field]
                       for field in ("viral_feature", "evidence")):
                raise ValueError("Invalid viral tumour setting text")
            if (not _strings(setting["source_urls"])
                    or not all(value.startswith("https://") for value in setting["source_urls"])):
                raise ValueError("Invalid viral tumour setting sources")
    for name, record in data["subsets"].items():
        required = {"title", "description"}
        if not isinstance(name, str) or not isinstance(record, dict) or required - set(record):
            raise ValueError("Invalid viral subset record")
        if not all(isinstance(record[field], str) and record[field]
                   for field in ("title", "description")):
            raise ValueError("Invalid viral subset text")
        if "mechanisms" not in record and "abbreviations" not in record:
            raise ValueError("Invalid viral subset mechanisms")
        if "mechanisms" in record and not _strings(record["mechanisms"]):
            raise ValueError("Invalid viral subset mechanisms")
        if "abbreviations" in record and not _strings(record["abbreviations"]):
            raise ValueError("Invalid viral subset abbreviations")
        if "source_urls" in record and (not _strings(record["source_urls"])
                                        or not all(value.startswith("https://")
                                                   for value in record["source_urls"])):
            raise ValueError("Invalid viral subset sources")
        if "mechanisms" in record and record.get("mechanism_match", "all") not in ("any", "all"):
            raise ValueError("Invalid viral subset mechanism match")
        if "mechanisms" in record and set(record["mechanisms"]) - set(data["mechanisms"]):
            raise ValueError("Unknown viral subset mechanism")
    return data


def _public_record(name, record):
    result = dict(record)
    result["name"] = name
    return result


def viral_mechanisms():
    """List dictionaries describing tumour-scoped viral mechanism selectors."""
    registry = _registry()
    return tuple(dict(_public_record(name, registry["mechanisms"][name]),
                      references=len(registry["mechanisms"][name]["abbreviations"]))
                 for name in sorted(registry["mechanisms"]))


def normalize_viral_mechanisms(mechanisms):
    """Normalize mechanism names and their concise aliases."""
    if mechanisms is None:
        return ()
    mechanisms = (mechanisms,) if isinstance(mechanisms, str) else tuple(mechanisms)
    if not mechanisms or any(not isinstance(value, str) or not value.strip() for value in mechanisms):
        raise ValueError("mechanisms must contain one or more nonempty names")
    registry = _registry()
    aliases = registry["aliases"].get("mechanisms", {})
    normalized = tuple(dict.fromkeys(
        aliases.get(value.strip().lower(), value.strip().lower()) for value in mechanisms))
    unknown = sorted(set(normalized) - set(registry["mechanisms"]))
    if unknown:
        raise ValueError(
            f"Unknown viral mechanism(s) {', '.join(unknown)}; available: "
            f"{', '.join(sorted(registry['mechanisms']))}"
        )
    return normalized


def viral_mechanism_abbreviations(mechanisms, *, match="any"):
    """Return GDC abbreviations selected by mechanism names."""
    names = normalize_viral_mechanisms(mechanisms)
    if match not in ("any", "all"):
        raise ValueError("match must be 'any' or 'all'")
    if not names:
        return ()
    records = _registry()["mechanisms"]
    groups = [set(records[name]["abbreviations"]) for name in names]
    selected = set.union(*groups) if match == "any" else set.intersection(*groups)
    return tuple(sorted(selected))


def viral_subsets():
    """List dictionaries describing named, reproducible viral subsets."""
    registry = _registry()
    return tuple(dict(_public_record(name, registry["subsets"][name]),
                      references=len(viral_subset_abbreviations((name,))))
                 for name in sorted(registry["subsets"]))


def normalize_viral_subsets(subsets):
    """Normalize named viral subsets and their concise aliases."""
    if subsets is None:
        return ()
    subsets = (subsets,) if isinstance(subsets, str) else tuple(subsets)
    if not subsets or any(not isinstance(value, str) or not value.strip() for value in subsets):
        raise ValueError("subsets must contain one or more nonempty names")
    registry = _registry()
    aliases = registry["aliases"].get("subsets", {})
    normalized = tuple(dict.fromkeys(
        aliases.get(value.strip().lower(), value.strip().lower()) for value in subsets))
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
    selected = set()
    for name in names:
        record = records[name]
        selected.update(record.get("abbreviations", ()))
        if "mechanisms" in record:
            selected.update(viral_mechanism_abbreviations(
                record["mechanisms"], match=record.get("mechanism_match", "all")))
    return tuple(sorted(selected))
