# Viral reference data

The metadata-only GDC `GRCh38.d1.vd1` inventory table is built from the
[published viral inventory](https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt).
The source file's SHA256 is
`f221e758f47e4784714b4a89198b89536ddacdc905da9c778a6802d3e1e51c56`.
The unusual `GRCh83` spelling belongs to GDC's file name; the reference name
is `GRCh38.d1.vd1`. Rebuild the inventory with:

```sh
python tools/build_gdc_viral_coverage.py GRCh83.d1.vd1_virus_decoy.txt \
  orfdata/data/viral/gdc-grch38-d1-vd1-coverage.tsv
```

This inventory supplies names and source accession leads, not qualified coding
annotations. Only the released HPV16 GDC placement is marked as bundled.

## Source ledgers

Each source the collection draws on has one packaged ledger, all in the
`orfdata.reference_ledger.v2` format: one row per inventory entry, recording
the identifier that source resolved it to, the molecule's length, MD5 and
SHA-256, how many coding features it annotates, its NCBI taxon and its own
description. The ledgers hold checksums and counts, never sequence.

`gdc-grch38-d1-vd1-genbank.tsv` records all 200 GenBank accessions in the
inventory, checked 2026-10-24: every one resolved, and its versioned record is
what `download_viral_references()` fetches and checks what it received
against. NCBI E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK25501/ .

`gdc-grch38-d1-vd1-pave.tsv` records the PaVE reference the inventory names
for each of the 189 papillomavirus rows. 187 are `published`; the other two
name types PaVE does not publish (`HPV-mRTRX7` and `HPV-mSD2`) and are
recorded `unavailable`. A row whose published coordinates did not reproduce
its published protein would be recorded `inconsistent` and never used.

Rebuild either with:

```sh
python tools/build_reference_ledger.py --source genbank \
  orfdata/data/viral/gdc-grch38-d1-vd1-genbank.tsv
python tools/build_reference_ledger.py --source pave \
  orfdata/data/viral/gdc-grch38-d1-vd1-pave.tsv
```

A ledger pins what ORFdata's own readers counted when it was built, so rebuild
it when a source is re-checked and when those readers change. Each ledger's
SHA-256 is pinned in `orfdata/viral_registry.py`, and the GenBank ledger's is
also the collection's cache identity.

## On-demand GenBank source-record bundle

`download_viral_references()` and
`orfdata download viral --species human` use the
version-pinned accessions in the GenBank ledger to retrieve all 200 GenBank
records from NCBI EFetch, and check each one against the length, checksums and
coding-feature count that ledger pins. The bundle is installed in the ORFdata
cache and keeps three layers: viral genomic DNA, one coding transcript per source CDS
feature with a GTF placing it on the genome, and the source translation of each
transcript. It also keeps the tables describing them, each source GenBank
record exactly as NCBI sent it, compressed, and a SHA256 manifest. These are
source records for alignment and annotation workflows; downloading them does
not make the entries bundled, curated coding references in ORFdata.

`human` is the target analysis species. It is not a claim that every source
record declared a human host. `gdc-viruses` is the public collection identity
for the 200-reference GDC viral screen and `GRCh38.d1.vd1` is the source
inventory name. The
`cancer-transcriptome-screening` context is its current analysis context. The named
`cancer-causing` and `cancer-associated` subsets select IARC Group 1 and IARC
Groups 1, 2A and 2B reference sets, respectively. GDC membership is not a
claim that each virus causes cancer, and contexts do not determine cache
identity.

## Which molecule a collection carries, and from which source

The inventory names both a GenBank accession and a PaVE reference for each
papillomavirus row, and they are not always the same molecule. GDC's decoy
FASTA carries the PaVE one: PaVE's `HPV16REF` reproduces the GDC HPV16 DNA MD5
`f2463721818c29998ef5ab917e6df74a` pinned in `tools/build_hpv16_cds.py`
exactly, while GenBank's `K02718.1` is a different, two-base-shorter sequence
with a sparser 1985 annotation. 24 of the 189 papillomavirus rows differ this
way, among them HPV16, HPV18, HPV56 and HPV59; the inventory's own `PaVE Notes`
column marks 18 of them `Revised`.

The source ledgers above make that choice offline, and a fetched molecule is
verified against the pinned length and checksums before anything is written.
Every genome and transcript records which source it came from.

PaVE: https://pave.niaid.nih.gov/ (NIAID), API at https://pave.niaid.nih.gov/api .

## Annotation for records whose own source leaves them unannotated

Three inventory rows — HPV2 (`X55964.1`), HPV8 (`M12737.1`) and HPV51
(`M62877.1`) — point at GenBank submissions that carry no CDS feature at all.
NCBI publishes no RefSeq for HPV8 or HPV51, and no other GenBank record holds
the same sequence with annotation, so the coding features come from the
alternate reference the inventory itself names for every papillomavirus row.

The fill is used only for a genome whose own record annotates nothing, only
when the alternate publishes the identical sequence, and only when each
feature's coordinates reproduce the protein that source publishes for it. Each
transcript and each genome records `annotation_source`, so a model taken from
elsewhere is never mistaken for the record's own. `--no-annotation-fill`
(``annotate_from_alternate=False``) keeps strictly to each entry's own
record.

## HPV16 coding reference

Source: NCBI RefSeq **NC_001526.4**, complete GenBank record:
https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_001526.4&db=nuccore&report=genbank&retmode=text

Source record SHA256:
`282be442f226cc809fefbda616ac136b8d20dfde23fad09558f3f653a3de96a4`.
Source authors and annotation provenance remain in the versioned record,
preserved as a compressed test fixture. NCBI database usage policies apply:
https://www.ncbi.nlm.nih.gov/home/about/policies/ . No endorsement is implied.

The catalogue contains nine source-annotated CDSs, including two joined paths.
It does not claim all reported HPV16 ORFs or complete transcript models.
RNA support does not establish translation or protein expression.

The released catalogue also contains the full forward reference bases for both
qualified placements. HPV16's documented human host is sourced from PaVE's
HPV16REF record (`/host="Human (Homo sapiens)"`):
https://pave.niaid.nih.gov/locus_viewer?feature_id=HPV16_E7&seq_id=HPV16REF .
This host association applies to HPV16 as a type, not to all papillomaviruses,
and does not claim an exhaustive host range. ICTV lists HPV16 in the species
*Alphapapillomavirus 9*:
https://ictv.global/report/chapter/papillomaviridae/taxonomy/papillomaviridae .
These are distinct from the RefSeq DNA/CDS source and the GDC coordinate
placement. Consulted 2026-10-16.

Native DNA MD5: `08d2b3849496d12cc45fd4df713d5c2f`.
GDC GRCh38.d1.vd1 HPV16 DNA MD5: `f2463721818c29998ef5ab917e6df74a`.
GDC DNA is the source circle rotated by 864 bases. Both placements retain
their own SHA256 and ordered, 1-based inclusive coding blocks.
GDC reference files: https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files .

Rebuild with `orfdata` installed:

```sh
python tools/build_hpv16_cds.py tests/data/viral/NC_001526.4.gb.gz \
  orfdata/data/viral/hpv16-refseq-CDS.json.gz
```

The builder checks source DNA identity, the exact GDC rotation and all nine
source translations. The resource checksum is pinned in `orfdata/viral.py`.
Broader coverage: https://github.com/openvax/orfdata/issues/91 through 98.
