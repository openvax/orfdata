# Dog peptides: MHC ligands and the Canine PeptideAtlas

Two dog peptide datasets record which peptides of known dog proteins mass
spectrometry has observed. Neither searched for noncanonical ORFs, so they
don't add ORFs to the dog proteome
([why dog has no ncORF catalogue](catalogs.md#dog)).

| Dataset | What it holds | How to get it |
|---|---|---|
| `kaabinejadian-dog-dla-ligands` | DLA class I ligands from four dog tumours and three DLA-88 monoallelic cell lines | `orfdata datasets download` (CC0, 25 MB) |
| `peptideatlas-dog-peptides` | 149,546 tryptic peptides from 18 Labrador tissue, plasma and urine samples, with genome locations | `orfdata datasets build` (built locally from 71 MB of sources) |

## DLA class I ligands

[Kaabinejadian et al. (2026)](https://doi.org/10.1016/j.isci.2026.115975)
eluted peptides from DLA class I molecules with the BB7.6 or H58A antibody. PRIDE
[PXD074485](https://www.ebi.ac.uk/pride/archive/projects/PXD074485) deposits
their PEAKS Studio 11.5 peptide tables (1% FDR), which this dataset imports in
full: 76,550 rows, each with its PEAKS area, intensity and spectrum count.

| Sample | Tumour | Breed | DLA-88 typing | Antibodies |
|---|---|---|---|---|
| Lola | osteosarcoma | Black Russian Terrier | \*501:01 | BB7.6, H58A |
| 163828A | anal sac adenocarcinoma | English Cocker Spaniel | \*004:02 | BB7.6 |
| Lily | anal sac adenocarcinoma | English Cocker Spaniel | \*501:01, \*501:02, \*004:02 | H58A |
| Bogey | oral melanoma | American Cocker Spaniel | \*003:02, \*006:01, \*017:01 | BB7.6 |
| HCT116-DLA88… | human HCT 116 cells, one DLA-88 allele each | — | \*003:02, \*012:01 or \*501:01 | BB7.6 |

Each tumour's spectra were searched four times: against UniProt's dog
proteome, Basenji and Great Dane (Zoey) Ensembl 104 proteins, and a merged
database of four dog assemblies. Each search is its own assay
(`Lola-BB76:UniProt9615`, `Lola-BB76:MergedCanFams`, ...). A peptide found in
two searches of one sample is one observation, not two. The HCT 116 lines were
searched against human SwissProt, and their peptides are human proteins, so
they sit in the dataset's separate human reference context.

Database search finds 11,363 distinct dog peptides in the UniProt searches and
12,948 across all four databases, leaving out iRT standards (below). For the
cell lines, the counts of unique 8–14-mers from database search (591 for
DLA-88\*003:02, 1,465 for \*501:01) match the article's. The article's tumour counts come from its curated supplementary
lists, which aren't imported (they are CC BY-NC-ND), and match no single table.

Every row keeps its PEAKS identification type in provenance as `found_by`:
`DB Search`, `DeepNovo` (de novo sequencing) or `Homolog` (SPIDER homology
search). Only database-search peptides are guaranteed to occur in a listed
protein, so filter on it for a conservative ligand list.

### iRT standards

Every run was spiked with Biognosys iRT retention-time peptides, and their
spectra appear in every search, often in disguise. 10,629 rows are marked
`irt_standard`, and `irt_rule` names the first rule that marked them:

| `irt_rule` | Rows | Marks a row that |
|---|---|---|
| `irt_sequence` | 1,179 | is a kit peptide or a fragment of one, with I and L equal |
| `kit_only` | 1,097 | PEAKS attributes only to the kit (homology reads of kit peptides) |
| `irt_precursor` | 1,276 | sits at the precursor of a row marked by the first two rules in the same immunoprecipitation: same charge, m/z within 20 ppm, retention time within a minute |
| `de_novo_irt` | 7,077 | is a de novo sequence that PEAKS matches to the kit, or that shares 6 or more consecutive residues with a kit peptide |

The precursor rule applies to de novo and homology rows. It applies to
database-search rows only in the merged database, which lacks the kit, and only
when no search with the kit found the same peptide. Elsewhere PEAKS could
choose the kit peptide and didn't.

- **De novo rows:** half of them (8,448 of 16,830) are iRT reads, often a
  kit peptide with residues added or changed (`KATLFLQFGAQGSPFLK`). Don't take
  an unmarked de novo sequence for a new ligand without checking its spectrum.
- **Homology rows:** 8 are marked by precursor, such as `TGFLLDPVLNR`, which
  has exactly the mass of `GTFIIDPGGVIR` and co-elutes with it.
- **Database-search rows:** apart from 1,075 kit peptides, only one is marked:
  `IFLAFQQKAYVAP`, matched to dog proteins in the merged database at the
  precursor of `LFLQFGAAGSPFLK`.
- **Chance matches:** few. Shuffled kit sequences share 6 residues with about
  3 de novo rows, and shifting the kit's retention times by 10–20 minutes puts
  about 14 de novo rows at a kit precursor.

A mapping is `unmapped` when any of its rows is marked. Release 2026-10-16.1
(importer `kaabinejadian-dla-v1`) marked only exact kit sequences and kit-only
attributions (2,242 rows), so it missed the disguised reads.

```python
from orfdata import DatasetRegistry, ScopedTarget

registry = DatasetRegistry()
ligands = registry.open("kaabinejadian-dog-dla-ligands")  # After downloading it.
dog = next(r for r in ligands.manifest.references if r.taxon_id == 9615)

# Where a peptide was eluted, and how it was identified.
peptide = ScopedTarget(dog, "peptide", "AIKKGLEEL", "amino-acid-sequence")
for m in ligands.measurements(target=peptide, quantity="spectrum_count"):
    print(m.sample_id, m.assay_id, m.value, m.provenance["found_by"])

# Lola's database-search ligands (BB7.6, UniProt search).
lola = {m.target.identifier for m in ligands.measurements(
            reference=dog, assay_id="Lola-BB76:UniProt9615", quantity="spectrum_count")
        if m.provenance["found_by"] == "DB Search" and not m.provenance["irt_standard"]}

# Why a row is marked as an iRT standard.
for m in ligands.measurements(target=peptide, quantity="spectrum_count"):
    print(m.assay_id, m.provenance["irt_standard"], m.provenance["irt_rule"])
```

Targets are unmodified sequences; each row's modified form (for example
`M(+15.99)`) is in its provenance as `peptide`. `mappings()` lists, for each
search and identification type, the proteins PEAKS reported for a peptide:
`unique` or `ambiguous` for database-search and homology matches, and
`unmapped` with a reason for iRT standards and de novo sequences, whose listed
proteins stay in the mapping's evidence. The search FASTAs weren't deposited,
so no accession can be checked against a sequence.

### Runs missing from the deposit

PEAKS reports, for each peptide, the raw file and scan of its best spectrum,
and sums its spectra and area over all of the sample's runs. The authors' PEAKS
projects for two immunoprecipitations included runs that the deposit lacks, so
each row records `raw_file_deposited`. For 1,695 rows the best spectrum comes
from a missing run:

- **Lily (H58A):** a second injection of the remaining eluate, 22-1107 "1uL
  Residue" (fractions 12–30), 299 rows. The first injection, 22-1104, is
  deposited.
- **DLA-88\*012:01:** the 20240929 acquisition (fractions 14–25), 1,396 of its
  1,635 rows. The line's archive holds only its 20241013 runs (fractions 26–28
  in triplicate, and blanks), though the SDRF names the 20240929 run as the
  assay.

Those identifications come from the authors' tables, but their spectra can't
be checked against the deposit. The raw files of 163828A's H58A run have no
deposited results and aren't included.

## Canine PeptideAtlas

The [Canine PeptideAtlas](https://peptideatlas.org/builds/canine/)
([Kusebauch et al. 2026](https://doi.org/10.1038/s41597-026-06647-z))
reprocessed 138 runs from Labrador Retrievers (13 tissues, plasma and urine;
PRIDE PXD052932) at a PSM FDR of 0.001. This dataset imports its 2024-01 build
(573):

- PSM counts for 149,546 peptides in each of the 18 samples where they were
  observed (528,272 values);
- every protein each peptide maps to, with its position and, for Ensembl
  proteins, its ROS_Cfam_1.0 coordinates (exon-spanning peptides have one
  location per exon);
- the build's 17,795 protein calls with presence levels (6,423 canonical). Of
  its 17,307 UniProt and Ensembl calls, 16,959 map to one Ensembl 110 protein
  with an identical sequence, 151 to several and 197 to none; contaminant and
  decoy calls aren't mapped.

PeptideAtlas states no licence for its downloads, so orfdata never
redistributes this dataset. `build` downloads the build's 61 MB export and
Ensembl 110's dog proteins, checks their SHA-256, and builds the dataset
locally (about a minute; 540 MB on disk). After that, `open` works offline.

```sh
orfdata datasets build peptideatlas-dog-peptides
```

```python
from pyensembl import EnsemblRelease
from orfdata import DatasetRegistry, ScopedTarget

atlas = DatasetRegistry().open("peptideatlas-dog-peptides")  # After building it.
source, ensembl = atlas.manifest.references

# PSMs of one peptide in each sample.
peptide = ScopedTarget(source, "peptide", "AAHEEICTTNEGVMYR", "amino-acid-sequence")
for m in atlas.measurements(target=peptide):
    print(m.sample_id, m.value)

# The PeptideAtlas call for a PyEnsembl dog protein (albumin).
transcript = EnsemblRelease(110, species="dog").transcript_by_id("ENSCAFT00845010439")
protein = ScopedTarget(ensembl, "protein", transcript.protein_id, "ensembl")
for mapped in atlas.mapped_measurements(protein, require_sequence=True):
    call = mapped.measurement
    print(call.target.identifier, call.value, call.provenance["presence_level"])
```

The `ensembl` context is exactly PyEnsembl's dog release 110
(`ReferenceContext.from_genome(EnsemblRelease(110, species="dog"))`), so
PyEnsembl protein IDs work directly. The [dog RNA dataset](pfizer-rna.md) uses
Ensembl 110 IDs too, in a context that also pins its GTF. Ensembl 105–111 have
identical sequences for all 25,710 Ensembl proteins in PeptideAtlas's search
set.

Each peptide's mapping (`mappings()`, keyed by PeptideAtlas accession) keeps its
locations as `[protein, start, end, chromosome, strand, chromosome_start,
chromosome_end, exact]`. `exact` is false for the 4,767 locations where the
protein has isoleucine for leucine or the reverse, which PeptideAtlas maps as
equal. 204 peptides map only to decoy sequences; they pass the FDR threshold but
are probably false identifications, and are `unmapped` with that reason.

A protein call's PSM count includes peptides it shares with other proteins, so
calls don't add up to the build's 13,781,115 PSMs. PeptideAtlas's web page
reports 245,884 peptides and 13,716,253 PSMs for the build; its export, which
this dataset reproduces, has 149,546 peptides (197,285 distinct modified forms)
and 13,781,115 PSMs.
