"""Discover, download, inspect and query quantitative datasets."""

import sqlite3
import sys
from contextlib import redirect_stdout

from .cli_output import manifest_summary, print_fields, print_record_stream, print_result
from .cli_parser import UsageError, parse_source_paths
from .command_parsers import datasets_parser
from .datasets import DatasetRegistry
from .downloads import normalize_status
from .reference import ScopedTarget
from .registry import load_registry


def _download_options(args):
    """What the download, download-reference and build commands pass on to DataCache."""
    return dict(show_progress=not args.quiet, timeout=args.timeout)


def _reference(references, key, option):
    """The dataset reference ``option`` names, or a usage error listing the dataset's references."""
    if key not in references:
        raise UsageError(f"{option} {key!r} is not a reference of this dataset; choose from "
                         f"{', '.join(references)}")
    return references[key]


def validate_search_reference_options(registry, name, version=None, *, reference=None, asset_id=None):
    """Validate ``--reference`` and ``--asset`` individually and as a pair against a dataset's search references.

    Raises :class:`~orfdata.cli_parser.UsageError` (CLI exit status 2) when
    either value is unknown or when both are supplied and do not match the same
    declared search reference.

    Parameters
    ----------
    registry : DatasetRegistry
        Dataset registry containing ``name``.
    name : str
        Dataset identifier.
    version : str, optional
        Dataset release version (defaults to the registered default).
    reference : str, optional
        Exact native reference context key supplied via ``--reference``.
    asset_id : str, optional
        Declared source asset ID supplied via ``--asset``.
    """
    specs = registry.search_references(name, version)
    for option, value, choices in (("--reference", reference, [s.reference.key for s in specs]),
                                   ("--asset", asset_id, [s.source.asset_id for s in specs])):
        if value is not None and value not in choices:
            raise UsageError(f"{option} {value!r} is not a search reference of this dataset; choose from "
                             f"{', '.join(dict.fromkeys(choices)) or 'none'}")
    if reference is not None and asset_id is not None:
        if not any(s.reference.key == reference and s.source.asset_id == asset_id for s in specs):
            pairs = ", ".join(f"{s.reference.key} ({s.source.asset_id})" for s in specs) or "none"
            raise UsageError(
                f"--reference {reference!r} and --asset {asset_id!r} do not match the same search reference "
                f"of this dataset; available pairs: {pairs}"
            )


def _check_search_reference(registry, args):
    """A usage error for an unknown or mismatched --reference / --asset pair."""
    validate_search_reference_options(
        registry, args.name, args.version, reference=args.reference, asset_id=args.asset,
    )


def _check_query(registry, args, parser):
    """The dataset's references, after checking the query or pairs options its manifest settles.

    This needs only the registered manifest, so a command-line mistake exits 2
    whether or not the dataset is installed.
    """
    references = {r.key: r for r in registry.info(args.name, args.version).references}
    if args.reference is not None:
        _reference(references, args.reference, "--reference")
    if args.command == "pairs":
        for side in ("left", "right"):
            values = [getattr(args, f"{side}_{field}") for field in ("target", "kind", "namespace", "reference")]
            if any(values) and not all(values):
                parser.error(f"{side.title()} target queries require --{side}-target, --{side}-kind, "
                             f"--{side}-namespace and --{side}-reference")
            if values[3] is not None:
                _reference(references, values[3], f"--{side}-reference")
        return references
    if any((args.target, args.kind, args.namespace)) and not all((args.target, args.kind, args.namespace,
                                                                   args.reference)):
        parser.error("Target queries require --target, --kind, --namespace and --reference")
    if args.mapped and not args.target:
        parser.error("--mapped requires an exact target and reference")
    if not args.mapped and (args.mapping_status or args.require_sequence):
        parser.error("Mapping filters require --mapped")
    return references


def run(argv=None, prog=None):
    parser = datasets_parser(prog)
    args = parser.parse_args(argv)
    try:
        registry = load_registry(DatasetRegistry, args.registry, cache_root=args.cache_root)
        if args.command in ("inspect-reference", "download-reference", "query-reference"):
            _check_search_reference(registry, args)
        if args.command == "list":
            result = [m.to_dict() for m in registry.list(
                species=args.species, taxon_id=args.taxon_id, assay=args.assay)]
        elif args.command == "info":
            result = registry.info(args.name, args.version).to_dict()
            result["download"] = registry.download_metadata(args.name, args.version)
        elif args.command == "references":
            result = [spec.to_dict() for spec in registry.search_references(args.name, args.version)]
        elif args.command == "inspect-reference":
            inspected = registry.inspect_reference(args.name, args.version, reference=args.reference, asset_id=args.asset)
            result = dict(
                status=normalize_status(inspected.status),
                datacache_status=inspected.status,
                path=str(inspected.path),
                error=str(inspected.error) if inspected.error else None,
            )
        elif args.command == "download-reference":
            with redirect_stdout(sys.stderr):
                result = registry.download_reference(args.name, args.version, reference=args.reference, asset_id=args.asset,
                                                     force=args.force, source_path=args.source_path, **_download_options(args))
        elif args.command == "query-reference":
            reference = registry.open_reference(args.name, args.version, reference=args.reference, asset_id=args.asset)
            target = ScopedTarget(reference.reference, args.kind, args.target, args.namespace or reference.spec.namespace)
            result = dict(source=reference.spec.to_dict(), membership=reference.membership(target).to_dict())
            if not args.json:
                membership = result["membership"]
                result = dict(target=args.target, namespace=target.namespace, status=membership["status"],
                              reason=membership["reason"], reference=target.reference.key,
                              source_asset=reference.spec.source.asset_id)
        elif args.command == "inspect":
            inspect = registry.inspect_build if args.built else registry.inspect
            inspected = inspect(args.name, args.version)
            result = dict(
                status=normalize_status(inspected.status),
                datacache_status=inspected.status,
                path=str(inspected.path),
                error=str(inspected.error) if inspected.error else None,
            )
            if args.built or not registry.published(args.name, args.version):
                result.update(origin="source_build", verified=inspected.verified)
        elif args.command == "download":
            sources = {"dataset.sqlite": args.source_path} if args.source_path else None
            with redirect_stdout(sys.stderr):
                result = registry.download(args.name, args.version, force=args.force,
                                           source_paths=sources, **_download_options(args))
        elif args.command == "build":
            sources = parse_source_paths(args.source_path, metavar="ASSET=PATH")
            with redirect_stdout(sys.stderr):
                dataset = registry.build(args.name, args.version, force=args.force,
                                         source_paths=sources, download_options=_download_options(args))
            result = dict(path=str(dataset.path), counts=dataset.counts,
                          dataset_id=dataset.manifest.dataset_id, version=dataset.manifest.version)
        elif args.command == "pairs":
            references = _check_query(registry, args, parser)
            open_dataset = registry.open_build if args.built else registry.open
            dataset = open_dataset(args.name, args.version)
            selections = {}
            for side in ("left", "right"):
                identifier, kind, namespace, reference = (
                    getattr(args, f"{side}_{field}") for field in ("target", "kind", "namespace", "reference"))
                if identifier:
                    selections[f"{side}_target"] = ScopedTarget(references[reference], kind, identifier, namespace)
                selections[f"{side}_filters"] = {
                    field: getattr(args, f"{side}_{field}")
                    for field in ("quantity", "status", "unit", "scale", "log_base")
                    if getattr(args, f"{side}_{field}") is not None}
            matches = dataset.matched_measurements(
                allowed_bases=args.basis, allowed_temporal=args.temporal or ["concurrent"],
                reference=args.reference, pair_ids=args.pair_id, relation_ids=args.relation_id,
                **selections,
            )
            print_record_stream((match.to_table_row() for match in matches), format=args.format, limit=args.limit,
                                columns=(("PAIR", "pairing.pair_id"), ("BASIS", "pairing.basis"),
                                         ("TEMPORAL", "pairing.temporal_relation"),
                                         ("RELATION", "target_relation.relation"),
                                         ("LEFT REF", "target_relation.left.reference"),
                                         ("LEFT", "target_relation.left.identifier"),
                                         ("RIGHT REF", "target_relation.right.reference"),
                                         ("RIGHT", "target_relation.right.identifier")))
            return
        else:
            references = _check_query(registry, args, parser)
            open_dataset = registry.open_build if args.built else registry.open
            dataset = open_dataset(args.name, args.version)
            target = None
            if args.target:
                target = ScopedTarget(references[args.reference], args.kind, args.target, args.namespace)
            filters = dict(sample_id=args.sample, assay_id=args.assay_id,
                           quantity=args.quantity, status=args.status,
                           unit=args.unit, scale=args.scale, log_base=args.log_base)
            if args.mapped:
                measurements = dataset.mapped_measurements(
                    target, mapping_statuses=args.mapping_status or ["unique"],
                    require_sequence=args.require_sequence, **filters)
            else:
                measurements = dataset.measurements(target=target, reference=args.reference, **filters)
            columns = (("REF", "target.reference"), ("TARGET", "target.identifier"), ("SAMPLE", "sample_id"),
                       ("QUANTITY", "quantity"), ("VALUE", "value"),
                       ("UNIT", "unit"), ("SCALE", "scale"), ("STATUS", "status"))
            if args.mapped:
                columns = tuple((label, "measurement." + path) for label, path in columns) + (
                    ("MAPPING", "attribution"),)
                if args.format == "table":
                    print_fields(dict(requested_target=target.identifier, namespace=target.namespace,
                                      requested_reference=target.reference.key))
                    print("Values retain their original targets; group membership supplies context only.\n")
            print_record_stream((measurement.to_dict() for measurement in measurements),
                                format=args.format, limit=args.limit, columns=columns)
            return
        columns = (("dataset_id", "version", "assay_kinds") if args.command == "list"
                   else ("source.asset_id", "namespace", "reference.species", "reference.genome_build"))
        if args.command == "info" and not args.json:
            result = dict(manifest_summary(result), download=result["download"])
        print_result(result, json_output=args.json, columns=columns, limit=getattr(args, "limit", 25))
        if args.command in ("info", "query-reference") and not args.json:
            print("\nUse --json for full source, reference, and provenance metadata.")
        if not args.json and args.command == "list":
            print("\nDetails: orfdata datasets info NAME\nInstall: orfdata data fetch NAME")
    except (ValueError, OSError, ImportError, RuntimeError, sqlite3.Error) as exc:
        parser.fail(exc)
