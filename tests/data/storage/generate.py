"""Run from released 0.3.0/0.4.0/0.5.0 to generate genuine legacy databases."""

import hashlib
import json
import sys
from pathlib import Path

import orfdata
from orfdata import (
    Assay,
    DatasetManifest,
    Measurement,
    QuantitativeDataset,
    ReferenceContext,
    Sample,
    ScopedTarget,
    SourceAsset,
)

output = Path(sys.argv[1])
version = int(sys.argv[2])
references = tuple(ReferenceContext(taxon, name, assembly, "ensembl", "110")
                   for taxon, name, assembly in ((9606, "human", "GRCh38"),
                                                 (10090, "mouse", "GRCm39"),
                                                 (9615, "dog", "ROS_Cfam_1.0"),
                                                 (4932, "yeast", "R64-1-1")))
source = SourceAsset("input.tsv", "https://example.invalid/storage-fixture.tsv",
                     hashlib.sha256(b"synthetic fixture\n").hexdigest(), 18,
                     "Synthetic test data; Apache-2.0", redistributable=True)
manifest = DatasetManifest("legacy-fixture", "1", "Synthetic legacy storage validation",
                           references, (source,), ("Synthetic, not experimental evidence",),
                           ("rna_seq", "digest_ms"), {"scope": "four deliberately colliding species"},
                           {"writer": "orfdata", "writer_version": orfdata.__version__})
samples, assays, observations, pairs, relations, mappings = [], [], [], [], [], []
for index, reference in enumerate(references):
    samples.append(Sample("sample", reference, biological_unit_id="culture"))
    assays.extend((Assay("rna", "sample", "rna_seq", reference), Assay("protein", "sample", "digest_ms", reference)))
    gene = ScopedTarget(reference, "gene", "same-id", "fixture")
    protein = ScopedTarget(reference, "protein", "same-id", "fixture")
    observations.extend((
        Measurement("rna-value", gene, "sample", "rna", "rna_abundance", index, "TPM",
                    manifest.dataset_id, manifest.version, provenance={"source_asset": "input.tsv", "row": index}),
        Measurement("protein-value", protein, "sample", "protein", "protein_abundance",
                    None if index == 1 else 10.25 + index, "iBAQ", manifest.dataset_id, manifest.version,
                    status="below_detection" if index == 1 else "measured", detection_limit=20 if index == 1 else None,
                    uncertainty=0.5, uncertainty_type="standard_deviation",
                    provenance={"source_asset": "input.tsv", "row": index, "native": "<20" if index == 1 else "observed"}),
        Measurement("derived", gene, "sample", "rna", "rna_abundance", 1, "fixture-units",
                    manifest.dataset_id, manifest.version, origin="derived", input_measurements=("rna-value",),
                    derivation={"method": "synthetic lineage fixture"}),
    ))
    if version >= 2:
        from orfdata import AssayPair, TargetRelation

        pairs.append(AssayPair("pair", reference, manifest.dataset_id, manifest.version,
                               "rna", "protein", "same_material", {"source_asset": "input.tsv"}, "concurrent"))
        relations.append(TargetRelation("relation", manifest.dataset_id, manifest.version, gene, protein,
                                        "canonical_translation", "one_to_one", {"source_asset": "input.tsv"}))
    if version >= 3:
        from orfdata import ReferenceMapping

        mappings.append(ReferenceMapping("map", manifest.dataset_id, manifest.version,
                                         protein, reference, (ScopedTarget(reference, "protein", "mapped-id", "fixture"),),
                                         "unique", "synthetic identifier match", {"sequence": "not validated"},
                                         provenance={"source_asset": "input.tsv"}))
kwargs = dict(samples=samples, assays=assays, measurements=observations)
if version >= 2:
    kwargs.update(pairings=pairs, relations=relations)
if version >= 3:
    kwargs.update(reference_mappings=mappings)
dataset = QuantitativeDataset.create(output / f"v{version}.sqlite", manifest, **kwargs)
expected = dict(writer_version=orfdata.__version__, manifest=manifest.to_dict(),
                measurements=[m.to_dict() for m in observations], counts=dataset.counts)
(output / f"v{version}.json").write_text(json.dumps(expected, indent=2) + "\n")
