# Source references and PyEnsembl mappings

A proteomics experiment may search a known protein FASTA without naming a
genomic assembly. An older RNA study may report its annotation provider but
omit the release. ORFdata keeps these source identities and their unknown fields.

```python
from orfdata import ReferenceContext

source = ReferenceContext(
    taxon_id=9615, species="dog", genome_build=None,
    annotation_name="uniprot", annotation_version=None,
    sequence_sha256=searched_fasta_sha256,
)
```

Unknown assembly/provider fields require a pinned `sequence_sha256` or
`source_sha256`. The former identifies a sequence set; the latter can pin an
input table when its reference release is unreported. These source contexts
cannot bind a genomic transcript. A `ReferenceORF` still requires a complete
genomic model validated against its actual PyEnsembl transcript. Existing
complete-reference keys are unchanged.

## Query through an explicit mapping

`ReferenceMapping` connects a native `ScopedTarget` to all candidate targets in
one destination reference. It records the dataset/version, mapping method,
validation, source provenance, status and unresolved/unmapped reason. Importers
must perform the stated checks against the pinned inputs. Unique identifier
resolution is separate from sequence validation.

```python
from orfdata import ScopedTarget

target = ScopedTarget.from_pyensembl(gene, kind="gene", namespace="ensembl:110")
for view in dataset.mapped_measurements(target):
    print(view.measurement.value, view.measurement.unit)
    print(view.measurement.target, view.requested_target)
    print(view.mapping.method, view.mapping.validation)
```

The view contains the **original observation**, including its source reference,
target, sample, assay, unit and provenance. Nothing is relabeled or copied onto
the requested target. Direct `dataset.measurements(target=...)` still returns
only observations native to that exact target/reference.

By default, mapped queries return uniquely resolved assignments. Use
`mapping_statuses=["ambiguous", "unresolved"]` to inspect other candidates, and
`require_sequence=True` to require recorded exact sequence validation. A group
value returned for a candidate member has `attribution="group_membership_context"`;
it does not independently quantify that protein. Parent gene/transcript RNA
cannot map directly onto an ORF. Use a separate parent-context relationship.

The CLI exposes the same query:

```sh
orfdata datasets query DATASET --mapped --reference ref:SHA256 \
    --kind protein --namespace ensembl:110 --target PROTEIN_ID --require-sequence
```

## Importer contract

Pass `reference_mappings=` to `QuantitativeDataset.create`. Declare source and
destination contexts in the manifest. Each source/destination-reference pair
has one mapping containing all candidates, preventing competing unique claims.
Unmapped records remain visible through `dataset.reference_mappings()`.

The status is `unique`, `ambiguous`, `unmapped` or `unresolved`. A single reported
protein representative with incomplete group membership is unresolved; its
measurement cannot become a uniquely measured member. Identifier-based mappings
must say what was checked and leave sequence validation absent when unavailable.
For an exact sequence search, record:

```python
validation = {"sequence": {
    "source_sha256": source_protein_sequence_sha256,
    "target_sha256": matching_protein_sequence_sha256,
    "all_candidates_checked": True,
}}
```

These are hashes of the compared sequences, distinct from the whole searched
FASTA's hash in the reference context. Both hashes must agree to count as an
exact match. ORF assignments additionally require their own `frame` and
`position` evidence. Known reference-strain differences need explicit
`strain_relationship` evidence. Cross-species mappings use the separate
[orthology API](pairing.md#across-species).

Storage version 3 preserves these mappings. Earlier versions remain readable
and report no declared reference mappings; no mapping is inferred on load.
