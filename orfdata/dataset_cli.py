"""Discover, download, inspect and query quantitative datasets."""

import json
import sqlite3
import sys
from contextlib import redirect_stdout

from .cli_output import manifest_summary, print_fields, print_record_stream, print_result
from .cli_parser import UsageError
from .command_parsers import datasets_parser
from .datasets import DatasetRegistry
from .reference import ScopedTarget


def _download_options(args):
    """What the download, download-reference and build commands pass on to DataCache."""
    return dict(show_progress=not args.quiet, timeout=args.timeout)


def _reference(references, key, option):
    """The dataset reference ``option`` names, or a usage error listing the dataset's references."""
    if key not in references:
        raise UsageError(f"{option} {key!r} is not a reference of this dataset; choose from "
                         f"{', '.join(references)}")
    return references[key]


def _check_search_reference(registry, args):
    """A usage error for a --reference or --asset that none of the dataset's search references has."""
    specs = registry.search_references(args.name, args.version)
    for option, value, choices in (("--reference", args.reference, [s.reference.key for s in specs]),
                                   ("--asset", args.asset, [s.source.asset_id for s in specs])):
        if value is not None and value not in choices:
            raise UsageError(f"{option} {value!r} is not a search reference of this dataset; choose from "
                             f"{', '.join(dict.fromkeys(choices)) or 'none'}")


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
        definitions = None
        if args.registry:
            with open(args.registry) as handle:
                try:
                    definitions = json.load(handle)
                except ValueError as error:
                    raise ValueError(f"Invalid dataset registry {args.registry}: {error}") from error
        try:
            registry = DatasetRegistry(definitions, cache_root=args.cache_root)
        except (KeyError, TypeError, AttributeError) as error:
            if not args.registry:  # The packaged registry: a bug, which keeps its traceback.
                raise
            raise ValueError(f"Invalid dataset registry {args.registry}: {type(error).__name__}: {error}") from error
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
            result = dict(status=inspected.status, path=str(inspected.path),
                          error=str(inspected.error) if inspected.error else None)
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
            result = dict(status=inspected.status, path=str(inspected.path),
                          error=str(inspected.error) if inspected.error else None)
            if args.built or not registry.published(args.name, args.version):
                result.update(origin="source_build", verified=inspected.verified)
        elif args.command == "download":
            sources = {"dataset.sqlite": args.source_path} if args.source_path else None
            with redirect_stdout(sys.stderr):
                result = registry.download(args.name, args.version, force=args.force,
                                           source_paths=sources, **_download_options(args))
        elif args.command == "build":
            sources = {}
            for value in args.source_path:
                asset, separator, path = value.partition("=")
                if not separator or not asset or not path or asset in sources:
                    parser.error("Source overrides require unique ASSET=PATH entries")
                sources[asset] = path
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
