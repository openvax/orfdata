# Viral reference data

The metadata-only GDC `GRCh38.d1.vd1` coverage ledger is built from the
[published viral inventory](https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt).
The source file's SHA256 is
`f221e758f47e4784714b4a89198b89536ddacdc905da9c778a6802d3e1e51c56`.
The unusual `GRCh83` spelling belongs to GDC's file name; the reference name
is `GRCh38.d1.vd1`. Rebuild the ledger with:

```sh
python tools/build_gdc_viral_coverage.py GRCh83.d1.vd1_virus_decoy.txt \
  orfdata/data/viral/gdc-grch38-d1-vd1-coverage.tsv
```

This inventory supplies names and source accession leads, not qualified coding
annotations. Only the released HPV16 GDC placement is marked as bundled.

## GenBank record summaries

`gdc-grch38-d1-vd1-genbank-summaries.tsv` records NCBI ESummary metadata for
all 200 GenBank accessions in the inventory, checked 2026-10-22. It contains
the resolved accession.version, title, sequence length, taxon ID and record
dates. The lookup matched each inventory accession; it did not download or
store nucleotide sequences or CDS annotations. NCBI E-utilities:
https://www.ncbi.nlm.nih.gov/books/NBK25501/ .

Summary manifest SHA256:
`5d6c7cbc3f77ba5236acc2e94bc49cf78e3ff486c12535d2b51e80f3ab875498`.

## On-demand GenBank source-record bundle

`download_viral_references()` and
`orfdata download viral --species human` use the
version-pinned accessions in the summary manifest to retrieve all 200 GenBank
records from NCBI EFetch. The bundle is installed in the ORFdata cache and
includes genome, source-CDS nucleotide and protein FASTAs, manifest and feature tables,
individual compressed GenBank records, and a SHA256 manifest. These are source
records for alignment and annotation workflows; downloading them does not make
the entries bundled, curated coding references in ORFdata.

`human` is the target analysis species. It is not a claim that every source
record declared a human host. `gdc-grch38-d1-vd1` is the collection identity;
`cancer-transcriptome-screening` and `cancer-causing` are nonexclusive
contexts recorded for it. The latter selects viral agents that IARC classifies
as Group 1 carcinogenic to humans; GDC membership is not a claim that each virus causes cancer,
and contexts do not determine cache identity. `orfdata viral download` remains
a compatibility alias.

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

Rebuild with `orfdata[build-data]` installed:

```sh
python tools/build_hpv16_cds.py tests/data/viral/NC_001526.4.gb.gz \
  orfdata/data/viral/hpv16-refseq-CDS.json.gz
```

The builder checks source DNA identity, the exact GDC rotation and all nine
source translations. The resource checksum is pinned in `orfdata/viral.py`.
Broader coverage: https://github.com/openvax/orfdata/issues/91 through 98.
