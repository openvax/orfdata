# HPV16 coding reference

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
