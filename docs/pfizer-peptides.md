# Pfizer dog peptide and protein-group evidence

`pfizer-atlas-dog-peptide-evidence` adds the underlying peptide, PSM, protein-group
and reporter evidence from the [Oziolor et al. atlas](https://doi.org/10.1038/s41597-026-07512-9),
MassIVE **MSV000099744 / PXD070290**. Its companion
[`pfizer-atlas-dog-protein`](pfizer-protein.md) retains the published iBAQ/rIBAQ
summary. They describe the same experiment.

Release **2026-10-27.1**, available with ORFdata **0.40.0** or later, downloads
6.90 GB in 4 compressed parts and installs a 24.71 GB SQLite
artifact. The full build and native-source acceptance receipts are included in
[the data release](https://github.com/openvax/orfdata/releases/tag/data-pfizer-dog-peptide-evidence-2026-10-27.1).

## Coverage and interpretation

The selected source contains 70,280,798 records across 42 tables, including
49,192,521 reported membership links. It retains target and decoy partitions,
confidence definitions, modifications, acquisition identifiers and quantitative
vectors. Spectral peak arrays, chromatograms and workflow/UI payloads are excluded.
The source lineage covers 96 runs and 44 channels: 38 tissue pools, four bridges
and two empty channels. A positive pool does not identify individual positive dogs.

The deposited FASTA has 20,649 searched sequences. The independent sequence audit
found 953,980 distinct target peptide sequences and 1,171,840 compatible placements:
1,132,603 exact and 39,237 additional isoleucine/leucine-equivalent placements.
Reported source memberships and sequence-compatible candidates are separate
queries. A candidate placement does not establish unique protein attribution.
The FASTA excludes spliced isoforms; it is not a complete canine isoform background.

No confidence threshold or new normalization is applied. Native confidence,
exclusion, modification and ambiguity fields remain available for an explicit
downstream policy. Decoy records remain controls even when their native confidence
label is high. Bulk tissue MS supplies no normal-tissue DLA-presentation evidence
and does not add noncanonical ORFs to the dog proteome.

## Explicit download and offline queries

```sh
orfdata datasets download pfizer-atlas-dog-peptide-evidence
```

The download installs a checksummed SQLite artifact from ordered compressed parts.
Normal imports and queries do not acquire the 86 GB vendor database or a separate
source FASTA. After installation:

```python
from itertools import islice
from orfdata import DatasetRegistry

dataset = DatasetRegistry().open("pfizer-atlas-dog-peptide-evidence")
target = next(dataset.measurements(status="measured")).target

for value in islice(dataset.measurements(target=target), 10):
    print(value.sample_id, value.assay_id, value.unit, value.value, value.status)

# The same TargetMapping model used by other peptide datasets.
assignments = list(dataset.mappings(reference=target.reference,
                                   source_identifier=target.identifier))
for assignment in assignments:
    print(assignment.evidence["basis"], assignment.status, assignment.targets)

# The same TargetRelation model used for other group memberships.
protein = next(t for assignment in assignments for t in assignment.targets)
groups = list(dataset.relations(right_target=protein))
```

These are the same high-level queries and return types used for other sources
and species. A reported assignment is tied to its original peptide-group or PSM
record; a `sequence_compatibility` mapping retains all compatible proteins and
positions separately. Selecting one destination with `mappings(target=protein)`
still returns every candidate in each selected mapping. An incomplete reported
assignment has status `unresolved`, a reason and any known candidates. A unique
mapping does not independently establish protein detection.

Measurement units name the original column: `PD2.4 TargetPeptideGroups.Abundances`,
`PD2.4 TargetPeptideGroups.AbundancesNormalized`,
`PD2.4 TargetPeptideGroups.AbundancesScaled` or `PD2.4 TargetPsms.Abundances`.
These are distinct series, not calibrated copies per cell or the summary's iBAQ.
Numeric zeros remain measured zeros. Null and nonfinite values remain explicitly
missing, with their native state preserved in provenance. The low-level adapter’s
`missing_acquisitions()` reports requested acquisition keys with no source row, which differs from a
present row with missing quantities. Queries stream observations from compact
vectors rather than loading a dense matrix.

Group relations retain the pinned same-study protein-summary context. They
neither transfer a group's abundance to an individual protein nor supply another
independent measurement. Publication and deposit instrument reports remain
explicit where unresolved.

For low-level inspection of vendor table IDs, decoy controls, confidence
definitions or individual acquisition records, use the source adapter explicitly:

```python
from orfdata.sources.pfizer_attribution import CanineEvidence

native = CanineEvidence(dataset.path)
context = native.summary_group("A0A5F4CET6")
```

This adapter exposes native records without changing the source-independent
measurement, mapping and relationship APIs.

## Reproducibility and source integrity

The published source component retains every selected typed value. Full-source
acceptance compares the consumer artifact against the native export, checks
membership endpoints and quantitative axes, and verifies the independently
established candidate digest. Source and artifact archives retain separate hashes
for compressed parts and their installed logical file.

Rebuilding from the pinned normalized sources is explicit:

```python
dataset = DatasetRegistry().build("pfizer-atlas-dog-peptide-evidence")
```

The returned dataset uses the same query API. Use `open_build()` to reopen a
local rebuild, whose own output checksum is recorded in its build receipt.

The builder creates one private database and publishes it only after complete
integrity checks. With the 32 GiB output cap, an already expanded source requires
34.4 GB of additional free space; this is a safety reservation, not the measured
artifact size. This release measured 24.71 GB of additional build files at peak.
Compressed source acquisition also needs space for its expansion.
Interrupted builds remove their own private output and preserve existing datasets.
See the [extraction and build guide](https://github.com/openvax/orfdata/blob/main/tools/README-dog-results.md)
for the selected tables, resource receipts and provenance.

Selective FTPS acquisition checked remote identity, offsets, lengths and acquired
chunk hashes. There is no authoritative checksum for the entire vendor database;
the recorded integrity claim covers the selected acquired facts. Source records
are CC0 under the MassIVE deposit, and searched FASTA records retain UniProt
attribution and CC-BY-4.0 terms. Compiled lineage retains factual source locations;
the publication's original supplementary workbook is not redistributed.
