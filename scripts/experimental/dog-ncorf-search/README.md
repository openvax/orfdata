# Dog ncORF search (experimental)

> **EXPERIMENTAL.** These scripts aren't part of the orfdata package and aren't installed
> with it. Their outputs aren't an orfdata dataset, and nobody has reviewed their results.
> Treat every ORF they report as a candidate, not as evidence of translation.

No public Ribo-seq data or noncanonical ORF (ncORF) mass-spectrometry search exists for
dog. So these scripts search the Canine PeptideAtlas's spectra again (PRIDE
[PXD052932](https://www.ebi.ac.uk/pride/archive/projects/PXD052932), Kusebauch et al.,
Sci Data 2026) against a database of Ensembl 110 dog proteins plus candidate ncORF
translations, and report the ncORFs whose peptides survive a class-specific FDR and a set
of filters.

## Steps

| Step | Script | Time | Network |
|---|---|---|---|
| 1. Get reference files | (manual) | minutes | about 200 MB |
| 2. Build the search database | `build_db.py` | a few minutes | none |
| 3. List the runs | `list_runs.py` | seconds | a few MB of ZIP directories |
| 4. Search every run | `driver.sh` | about 11–13 h | about 204 GB, streamed |
| 5. Aggregate | `aggregate.py` | minutes; can run any time | none |

Run every Python script with `python -I` (isolated mode). Each script finds its helper
modules (`standards.py`, `remotezip.py`) in its own directory.

### Requirements

- **Sage** v0.14.7, from <https://github.com/lazear/sage/releases/tag/v0.14.7>. The
  macOS arm64 build used here is `sage-v0.14.7-aarch64-apple-darwin.tar.gz` (SHA-256
  `6dc2270f311c6f3405bd7baff7ee1d478422ffc19574504af9c5dc5562fa2ac4`). Its binary reports
  version 0.14.6.
- **Python 3** with pandas, pyarrow, numpy, biopython, pyensembl and inflate64. The PXD052932
  archives compress with Deflate64, which Python's `zipfile` can't decompress, so
  `fetch_member.py` uses inflate64. Tested with pandas 3.0.6, pyarrow 26.0.0, numpy 2.5.4,
  biopython 1.84, pyensembl 2.27.2 and inflate64 1.0.4.
- About 10 GB of free memory while Sage runs, and at least `MIN_FREE_GB` + 4 GB of free disk
  (see [Resources](#resources)).

### 1. Get the reference files

Dog, Ensembl 110 (ROS_Cfam_1.0). `pyensembl install --release 110 --species dog` puts them
in PyEnsembl's cache, or download them from
`https://ftp.ensembl.org/pub/release-110/`:

| File | Path under release-110/ | SHA-256 used in the pilot |
|---|---|---|
| GTF | `gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.110.gtf.gz` | `f0f226cc42b8579d1f6d1668344f60804b6df6f68cce74de73604e074b160651` |
| cDNA | `fasta/canis_lupus_familiaris/cdna/Canis_lupus_familiaris.ROS_Cfam_1.0.cdna.all.fa.gz` | `a76e86c8e36bc7bf8ed20bc39aa90b69fac2815883dd62713cd9ee3c2b8929d2` |
| ncRNA | `fasta/canis_lupus_familiaris/ncrna/Canis_lupus_familiaris.ROS_Cfam_1.0.ncrna.fa.gz` | `2dd613c83b743483ce804b45ade281657ec996e3f274e877c141ad9a3ae2a417` |
| Proteins | `fasta/canis_lupus_familiaris/pep/Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa.gz` | `9daa065155157ab337fbdb882914acfd87aefaf8a7425bfcb63a3e1243b0a070` |

PeptideAtlas's search FASTA (62,456,242 bytes, SHA-256
`e13ab1dbd80853ea430b44c967f4a9f9cc08fb439e098a1be2207e2a1daafac3`): UniProt UP000805418,
Ensembl ROS_Cfam_1.0 proteins and 499 contaminants, with decoys:
<https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/12/PXD052932/UP000805418_2023-09-25_ENS-ROS_Cfam_1.0_CONTAM_targetdecoy.fasta>

Human and mouse proteins, used only to judge conservation (the pilot used Ensembl 115):
`https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz`
and `https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz`.

### 2. Build the search database

```sh
python -I build_db.py --gtf GTF --cdna CDNA --ncrna NCRNA --pep PEP \
    --peptideatlas-fasta UP000805418_2023-09-25_ENS-ROS_Cfam_1.0_CONTAM_targetdecoy.fasta \
    --pyensembl-cache WORK/pyensembl-cache --out-dir WORK/db
```

| Option | Meaning | Default |
|---|---|---|
| `--gtf`, `--cdna`, `--ncrna`, `--pep` | Ensembl 110 dog files (step 1) | required |
| `--peptideatlas-fasta` | PeptideAtlas's search FASTA (step 1) | required |
| `--pyensembl-cache` | Directory where PyEnsembl writes its index of the GTF (about 500 MB) | required |
| `--out-dir` | Output directory | required |
| `--min-length` | Shortest candidate ncORF, in amino acids, stop codon excluded | 10 |

The database has two classes.

- **Canonical** (`CANON000000`, ...): Ensembl 110 proteins, then PeptideAtlas's targets
  (UniProt, Ensembl, contaminants; its `DECOY_` entries are dropped), then one entry per
  retention-time standard peptide (`std|iRT|LGGNEQVTR`, ...; see
  [Retention-time standards](#retention-time-standards)). Identical sequences become one
  entry, whose header lists up to five source IDs (`ens|`, `pa|`, `contam|`, `std|`).
- **ncORF** (`NCORF000000`, ..., ordered by contig and start): ATG-initiated ORFs that end at
  a stop codon inside an Ensembl 110 transcript.

  | Kind | Where it starts |
  |---|---|
  | `uORF` | 5′ UTR, ending before the CDS starts |
  | `uoORF` | 5′ UTR, out of frame, ending inside or after the CDS |
  | `dORF` | 3′ UTR |
  | `lncORF` | any frame of a transcript of a lncRNA gene |
  | `ncTxORF` | any frame of a `retained_intron` or `processed_transcript` transcript |

  For each stop codon only the most upstream ATG counts, since in-frame downstream ATGs give
  suffixes of the same ORF. In-frame upstream ATGs with no stop before the CDS (N-terminal
  extensions) are skipped. Candidates are deduplicated by sequence, and dropped when they are
  identical to or contained in a canonical sequence, or when they contain a standard peptide
  that no canonical protein contains (I and L count as equal). An ORF ID is only meaningful
  within one build of the database.

Outputs:

- `dog_db.fasta`: the search database. The pilot's had 44,120 canonical and 310,527 ncORF
  sequences; this version adds 26 standard peptides to the canonical class and is otherwise
  identical (checked for the canonical half; no pilot candidate contains a standard).
- `ncorf_candidates.tsv`: one row per ncORF. Columns: `orf_id`, `length` (aa), `kinds`
  (comma-separated), `genes`, `gene_names`, `transcripts` (Ensembl IDs, comma-separated),
  `locus` (`contig:start-stop_end:strand` per placement, semicolon-separated; 1-based
  genomic positions of the A of the ATG and of the last base of the stop codon), `sequence`.
- `build_db_summary.json`: inputs with their SHA-256, counts per step (canonical inputs,
  transcript classes, candidates removed and why, ncORF kinds), the standards that canonical
  proteins already contain, and the SHA-256 of `dog_db.fasta`.

### 3. List the runs

```sh
python -I list_runs.py WORK/runs.tsv
```

reads the ZIP central directories of PXD052932's four mzML archives over HTTP range requests
and writes one row per `.mzML` member. Other members, the six Heart S-Trap runs' `.pin` files
and a directory entry, are skipped and reported. The archives are, under
`https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/12/PXD052932/`:

| Archive (`group`) | mzML runs | Instrument and samples |
|---|---|---|
| `mzML_QExactive_plasma_urine.zip` | 12 | Q Exactive; plasma and urine |
| `mzML_timsTOF_tissue_PAC_DDA.zip` | 78 | timsTOF, DDA; tissues, protein aggregation capture |
| `mzML_timsTOF_tissue_STrap_DDA.zip` | 24 (+6 `.pin` only) | timsTOF, DDA; tissues, S-Trap |
| `mzML_timsTOF_plasma_DDA.zip` | 18 | timsTOF, DDA; plasma |

That's 132 of the atlas's 138 runs: 203.9 GB compressed, 354.7 GB of mzML.

| Option | Meaning | Default |
|---|---|---|
| `OUT_TSV` | runs.tsv to write | required |
| `--base-url` | Directory URL holding the archives | the PXD052932 URL above |
| `--group NAME` | Archive name without `.zip`; repeatable | all four, in the order above |
| `--listing GROUP=PATH` | Use a saved listing (`remotezip.py list` output) instead of the network | none |
| `--save-listings DIR` | Also save each archive's listing as `DIR/GROUP.tsv` | off |
| `--match REGEX` | Keep only runs whose name matches, for example `Brain` for a trial | all |

`runs.tsv` has no header. Columns: `run` (member file name without `.mzML`), `group`,
`zip_url`, `member` (path inside the archive), `mzml_bytes` (uncompressed),
`compressed_bytes`. The driver searches runs in file order; edit the file to reorder or
drop runs.

### 4. Search every run

```sh
export DB_FASTA=WORK/db/dog_db.fasta SAGE=/path/to/sage PYTHON=python3
bash driver.sh --dry-run WORK                      # check the plan; touches nothing
nohup bash driver.sh WORK > WORK/driver.out 2>&1 &
```

A background fetch loop streams each run's mzML out of its archive with parallel HTTP range
requests (`fetch_member.py`), checks its CRC-32 and size, records its SHA-256 and MS2 count,
and stores it gzip-compressed in `WORK/work/`. It stays at most `MAX_RUNS_ON_DISK` runs ahead
of the search loop and waits while free disk is below `MIN_FREE_GB`. The search loop runs Sage
on each run, keeps its rank-1 PSMs (`post_run.py`), appends a receipt to the run manifest and
deletes the mzML. Downloads, Sage and post-processing all run under `nice`.

`--dry-run` checks the tools and inputs, then prints the settings, how many runs are done,
failed and pending, how much is left to download, free disk, and the exact commands for the
first pending run. It downloads, writes and deletes nothing.

| Setting | Meaning | Default |
|---|---|---|
| `WORK_DIR` | Work directory; the first argument, or this variable | required |
| `RUNS_TSV` | Run list from `list_runs.py` | `$WORK_DIR/runs.tsv` |
| `DB_FASTA` | Search database from `build_db.py` | `$WORK_DIR/db/dog_db.fasta` |
| `SAGE` | Sage binary | `sage` on `PATH` |
| `PYTHON` | Python with pandas, pyarrow and inflate64 | `python3` |
| `SAGE_PARAMS` | Sage parameter template (the driver fills in `database.fasta`) | `sage_params.json` here |
| `NCONN` | Parallel range requests per download (16 MiB each) | 16 |
| `SAGE_THREADS` | Sage threads (`RAYON_NUM_THREADS`) | 6 |
| `MIN_FREE_GB` | Don't start a download while the work directory's disk has less free space (GB) | 5 |
| `MAX_RUNS_ON_DISK` | Most runs' mzML on disk at once, including one being downloaded | 2 |
| `NICE` | `nice` increment for downloads, Sage and post-processing | 19 |
| `RETRY_FAILED` | `1` clears failure markers and retries those runs | 0 |

A work directory holds one search. On its first start the driver writes
`search_config.json` (the resolved Sage parameters and the database's SHA-256) and
`sage_params.resolved.json`. Once any run has finished, it refuses to restart with different
parameters or a different database (compared by checksum, so moving the file is fine); use a
new work directory instead. Only one driver can use a work directory at a time
(`.driver.lock`; remove it by hand if a driver was killed with `kill -9`).

#### Sage parameters (`sage_params.json`)

Fully tryptic peptides (cleave after K or R, not before P) with up to 2 missed cleavages, 7–30
residues and 500–5,000 Da. Carbamidomethyl cysteine is fixed; methionine oxidation and protein
N-terminal acetylation (`[`) are variable, at most 2 per peptide. Sage generates reversed
decoys (`rev_` prefix). Precursor and fragment tolerances are ±20 ppm, with precursor isotope
errors −1 to +3 and charges 2–4 where a spectrum has none. Singly charged b and y ions from
b2 and y2, between m/z 150 and 2,000; at most 150 peaks per spectrum, at least
15 peaks and 4 matched fragments. No deisotoping, chimeric or wide-window (DIA) search.
Retention time prediction is on, as a scoring feature. Sage reports two PSMs per spectrum;
only rank 1 is kept. `bucket_size` 32768 sets Sage's fragment index bucket size.

#### Resuming and failures

Stop the driver at any time with Ctrl-C or `kill`; it stops its download, removes its lock and
exits (after `kill`, once the current Sage search ends; that run is searched again later).
Start it again with the same command: runs with a row in `run_manifest.tsv` are never
searched again, a run that finished downloading is searched without downloading it again,
and partial downloads are deleted and restarted.

| Failure | What the driver does | Marker in `WORK/work/` | Log |
|---|---|---|---|
| Download fails 3 times (a minute apart; each try retries every range up to 10 times) | Skips the run | `RUN.fetch_failed` | `logs/fetch_RUN.log` |
| Sage exits non-zero | Deletes the mzML, keeps Sage's output directory `RUN_sage/` | `RUN.sage_failed` | `logs/sage_RUN.log` |
| `post_run.py` fails | Deletes the mzML, keeps `RUN_sage/` | `RUN.post_failed` | `logs/post_RUN.log` |

A restart skips runs with a marker. To retry them, fix the cause and restart with
`RETRY_FAILED=1`, which deletes the markers and leftover `RUN_sage/` directories; or delete
one run's marker by hand. `logs/driver.log` records every start, download (with its
receipt), finished run (with its manifest row), retry, low-disk wait and error.

### 5. Aggregate

```sh
python -I aggregate.py config.json
```

pools the PSM files of every finished run, so it can run during the search.
`aggregate_config.example.json` shows every key; replace `WORK` and `ENSEMBL` with real paths.

| Key | Meaning | Default |
|---|---|---|
| `psm_dir` | Directory of per-run PSM parquet files (`WORK/results/psms`) | required |
| `runs_tsv` | The run list, for the "runs aggregated: n of N" count | required |
| `out_dir` | Output directory | required |
| `db_fasta` | The database that was searched | required |
| `candidates_tsv` | Its `ncorf_candidates.tsv` | required |
| `human_pep`, `mouse_pep` | Ensembl human and mouse protein FASTAs, for conservation | required |
| `gtf`, `cdna`, `ncrna`, `pep` | The Ensembl 110 dog files used to build the database, for ORF placement | required |
| `pyensembl_cache` | PyEnsembl index directory (reuse `build_db.py`'s) | required |
| `manifest_tsv` | `run_manifest.tsv`; reports runs searched but missing from `psm_dir` | optional |
| `extra_dog_fastas` | More dog protein FASTAs (for example newer UniProt or RefSeq), used only to filter peptides | `[]` |
| `fdr` | FDR threshold for PSMs and peptides | 0.01 |
| `homology_min_score` | Smallest local alignment score for an ORF to count as homologous | 100 |
| `homology_min_identity` | Smallest identity over aligned positions for the same | 0.4 |
| `isobaric_tolerance_da` | Mass tolerance (Da) for near-isobaric 2-residue differences | 0.02 |

### Self-test

```sh
python -I check_standards.py
```

checks the retention-time standard handling of `standards.py` and `build_db.py` on tiny
made-up inputs, with no network or Ensembl files.

## Resources

- **Download:** about 204 GB. The pilot streamed at 4.3–5.1 MB/s with 16 connections, so
  about 11–13 h for all runs. A timsTOF run is 1.2–2.1 GB compressed (2.3–3.6 GB of mzML), a
  Q Exactive run about 0.2 GB (0.7 GB). The driver stores each as gzip level 1.
- **Disk:** at most `MAX_RUNS_ON_DISK` runs' mzML (about 4 GB with the default 2), plus Sage's
  output for the run being searched (about 100 MB). Kept results: about 20 MB of PSMs per run,
  about 3 GB in all.
- **Memory:** Sage held 10.3 GB resident (14.5 GB peak footprint) with the full database,
  whose index has 14.2 million peptides and 316 million fragments. It rebuilds the index for
  each run (about 10 s).
- **CPU:** the pilot searched 3 timsTOF runs in 154 s with all cores, so searching takes far
  less time than downloading.

## Outputs

`WORK/` after a search:

| Path | Contents |
|---|---|
| `runs.tsv` | The run list (step 3) |
| `search_config.json`, `sage_params.resolved.json` | The Sage parameters and database checksum this search uses |
| `results/run_manifest.tsv` | One row per finished run |
| `results/psms/RUN.parquet` | The run's rank-1 PSMs (zstd parquet) |
| `results/runs/RUN.receipt.json` | Its download receipt |
| `results/runs/RUN.sage_results.json` | Sage's `results.json`: version, the parameters it used, inputs and outputs |
| `logs/` | `driver.log`, and per run `fetch_RUN.log`, `sage_RUN.log`, `post_RUN.log` |
| `work/` | Transient: downloads in progress, the run being searched, failure markers |

**`run_manifest.tsv`** (tab-separated, with header):

| Column | Meaning |
|---|---|
| `run`, `member`, `zip_url` | The run, its path in the archive, the archive's URL |
| `compressed_bytes`, `mzml_bytes` | Size in the archive and uncompressed |
| `mzml_crc32`, `crc_ok` | The archive's CRC-32 of the member, and whether the download matched it and its size |
| `mzml_sha256` | SHA-256 of the uncompressed mzML |
| `n_ms2` | MS2 spectra (count of `ms level` 2 tags) |
| `download_s`, `download_MBps` | Download time and speed |
| `sage_s` | Sage's wall time, including building its index |
| `sage_rows`, `rank1_psms` | PSMs Sage reported (ranks 1 and 2), and rank-1 PSMs kept |
| `target_psms_q01`, `target_peptides_q01`, `target_proteins_q01` | Sage's own per-run counts at q ≤ 0.01 (spectrum, peptide and protein q), over the whole database |
| `canonical_rank1_targets`, `canonical_rank1_decoys` | Rank-1 PSMs to a canonical-class target or decoy (any protein `CANON…`) |
| `ncorf_rank1_targets`, `ncorf_rank1_decoys` | Rank-1 PSMs only to ncORF-class targets or decoys |
| `psm_parquet` | The PSM file's name |
| `completed_utc` | When the run finished |
| `irt_psms_q01` | Target PSMs to iRT standard peptides at spectrum q ≤ 0.01 |

**`psms/RUN.parquet`**: `run`, then every column of Sage's `results.sage.tsv` for rank-1 PSMs
(`psm_id`, `peptide` with modifications such as `M[+15.9949]`, `proteins` (semicolon-separated),
`num_proteins`, `filename`, `scannr`, `rank`, `label` (1 target, −1 decoy), `expmass`,
`calcmass`, `charge`, `peptide_len`, `missed_cleavages`, `semi_enzymatic`, `isotope_error`,
`precursor_ppm`, `fragment_ppm`, `hyperscore`, `delta_next`, `delta_best`, `rt`, `aligned_rt`,
`predicted_rt`, `delta_rt_model`, `ion_mobility`, `predicted_mobility`, `delta_mobility`,
`matched_peaks`, `longest_b`, `longest_y`, `longest_y_pct`, `matched_intensity_pct`,
`scored_candidates`, `poisson`, `sage_discriminant_score`, `posterior_error`, `spectrum_q`,
`peptide_q`, `protein_q`, `ms2_intensity`). See Sage's documentation for each.

**`RUN.receipt.json`**: `member`, `zip_url`, `compressed_bytes`, `mzml_bytes`, `mzml_crc32`,
`crc_ok`, `mzml_sha256`, `gz_bytes` (stored size), `n_ms2`, `download_s`, `MBps`.

`aggregate.py` writes to `out_dir`:

- **`aggregate_summary.txt`**: one line per step. Runs aggregated; per class, rank-1 targets
  and decoys, PSMs and peptides passing; retention-time standard PSMs; ncORF peptides passing;
  the size of the filtering reference; filter outcomes; supported ORFs by classification and
  kind; and candidate novel ncORFs with 2 or more peptides, in 2 or more runs, and without
  ORF-level homology.
- **`aggregate_summary.json`**: `summary` (those lines), `canonical` (canonical-class PSMs
  and peptides passing), `per_run_canonical_psms_q01`, `runtime_s`.
- **`run_qc.tsv`**: per run, target rank-1 PSMs at class q ≤ `fdr`: `tissue` (from the run
  name), `canonical_psms`, `ncorf_psms`, `standard_psms` and `standard_peptides` (spike-in
  standards in the canonical class), `standard_psms_ncorf_class` (should be 0 with a
  database from `build_db.py`), `natural_standard_psms` (standard sequences that dog proteins
  contain, which can't be told from endogenous peptides).
- **`ncorf_peptides_all_candidates.tsv`**: one row per ncORF-class target peptide at
  peptide-level q ≤ `fdr`, with the PSM columns: `peptide` (unmodified), `q_peptide`,
  `best_q_psm`, `n_psms_q01` (PSMs at class q ≤ `fdr`), `n_psms_rank1` (all its rank-1 PSMs),
  `n_runs`, `tissues`, `per_run_psms_q01` (JSON), `spectrum_ids` (`file:scan`,
  semicolon-separated), `best_hyperscore`, `best_score` (Sage discriminant score), and the
  best PSM's `best_spectrum`, `best_modified`, `best_charge`, `best_precursor_ppm`,
  `best_isotope_error`, `best_delta_next`, `best_matched_peaks`. Mapping: `orf_ids`,
  `orf_unique`, `kinds`. Filters: `IL_in_dog_reference`, `dog_1aa` and `n_dog_1aa` (up to 5
  reference sequences one residue away), `dog_near_isobaric_2aa` and `n_dog_near_isobaric`
  (two residues away with a mass difference explained by a modification), `human_exact`,
  `human_1aa`, `mouse_exact`, `mouse_1aa`, `fail_reason`, `passes`, `conserved_peptide`,
  `classification`, `orf_homologous_to_annotated_protein` (empty for filtered peptides).
- **`ncorf_peptides_passing.tsv`**: the same, only peptides that pass every filter.
- **`supported_orfs.tsv`**: one row per ncORF with a passing peptide. Identity: `orf_id`,
  `kind`, `orf_length_aa`, `protein_sequence`, `start_codon` (ATG). Placement on one
  transcript: `transcript_id`, `gene_id`, `gene_name`, `transcript_biotype`, `contig`,
  `strand`, `coding_blocks` (JSON list of 1-based inclusive genomic intervals, ascending,
  stop codon included), `start_codon_position`, `n_genomic_placements`, `all_placements`
  (JSON, every transcript the ORF sits on). Constants: `genome_build` ROS_Cfam_1.0,
  `annotation` Ensembl 110, `taxon_id` 9615, `species`, `translation_table` 1. Evidence:
  `n_peptides`, `n_psms_q01`, `n_runs`, `tissues`, `peptides`, `peptide_evidence` (JSON per
  peptide: `method` proteomics, `sequence`, `protein_intervals` (1-based), `source`
  sage_reanalysis, `accession` PXD052932, `mapping_count`, `spectrum_count`, `datasets` (runs),
  `sequence_match`, `attribution`, `quality`, `spectrum_ids`, `q_peptide`, `best_q_psm`,
  `best_hyperscore`, `classification`), `classification`. Homology:
  `orf_homologous_to_annotated_protein`, and per species `human_best_hit` (protein ID and gene),
  `human_score`, `human_identity`, `human_aln_len`, and the `mouse_` equivalents. Sorted by
  classification, then PSMs.
- **`supported_orfs.jsonl`**: the same rows as JSON, with `coding_blocks`, `all_placements`
  and `peptide_evidence` as JSON values.

The peptide and ORF files are written only when there is at least one row.

## Statistics and filters

1. **Classes.** A PSM is canonical-class if any of its proteins is a canonical entry (`CANON…`,
   or its `rev_` decoy), otherwise ncORF-class. The ncORF class has about seven times as many
   sequences as the canonical one but far fewer true matches, so a shared FDR would let
   ncORF matches ride on canonical ones.
2. **Class-specific FDR.** Within each class, rank-1 PSMs from every run are pooled and ranked
   by Sage's discriminant score; q-values come from target-decoy counts (decoys / targets,
   made monotone). Peptide-level q-values use each unmodified peptide's best score, separately
   for targets and decoys.
3. **ncORF peptides.** ncORF-class target peptides at peptide-level q ≤ `fdr`.
4. **Filters**, in this order (`fail_reason` names the first that applies). The dog
   reference is every canonical-class sequence in the database (Ensembl, UniProt, contaminants,
   standards), `extra_dog_fastas`, and the standard peptides; I and L count as equal.
   - `I/L-identical to dog reference`: the peptide occurs in a reference sequence.
   - `1 aa from dog reference`: a reference window differs at one position.
   - `near-isobaric (2 aa) to dog reference`: a window differs at two positions and their mass
     difference is within `isobaric_tolerance_da` of 0, oxidation (15.995), deamidation
     (0.984), acetylation (42.011) or one isotope (1.003).
   - `not ORF-unique`: the ncORFs containing the peptide don't share a stop codon, so it can't
     be assigned to one ORF.
5. **Classification** of passing peptides: `conserved_coding_probable_annotation_gap` if the
   peptide occurs, or nearly occurs (one residue away), in a human or mouse Ensembl protein,
   which suggests a protein dog Ensembl misses or mis-models; otherwise
   `candidate_novel_ncORF`. An ORF takes the first class if any of its peptides has it.
6. **ORF-level homology** (secondary evidence): each supported ORF is aligned locally
   (BLOSUM62, gap open −11, extend −1) against the 25 human and the 25 mouse proteins sharing
   the most 6-mers with it. It is homologous to an annotated protein if its best alignment in
   either species scores at least `homology_min_score` with identity at least
   `homology_min_identity`.

## Retention-time standards

Labs often spike synthetic retention-time peptides into runs. A database without them can
match their intense spectra to the nearest sequence it has, including an ncORF. So
`build_db.py` adds the 11 Biognosys iRT peptides (the same list as
`orfdata.sources.kaabinejadian.IRT_PEPTIDES`) and the 15 Pierce PRTC peptides
(`standards.py`) to the canonical class, and drops ncORF candidates that contain a standard
no canonical protein contains (none in the pilot database did). `aggregate.py` also puts
them in the filtering reference and counts their PSMs per run (`run_qc.tsv`);
`post_run.py` counts iRT PSMs per run.

PXD052932 most likely wasn't spiked: the Canine PeptideAtlas build has no iRT peptides and
only two PRTC peptides (GILFVGSGVSGGEEGAR, NGFILDGFPR), both ordinary dog sequences. The
standards stay in the database as a safeguard. Four PRTC sequences (DIPVPKPK, IGDYAGIK,
GILFVGSGVSGGEEGAR, NGFILDGFPR) occur in dog proteins, so their PSMs can't be told from
endogenous peptides and aren't counted as spike-ins. The PRTC kit's peptides also carry heavy
C-terminal lysine or arginine, which these Sage parameters don't search for.

## Pilot

Three timsTOF PAC runs (brain, liver and spleen), searched together with the same parameters
and the pilot database (without standards):

- canonical class: 249,156 PSMs and 37,744 peptides at q ≤ 0.01;
- ncORF class: 194 PSMs and 43 peptides at q ≤ 0.01. Filters removed 9 (4 I/L-identical, 3
  one residue away, 2 near-isobaric), leaving 34 peptides on 17 ORFs (7 dORF, 5 lncORF,
  3 uoORF, 2 uORF);
- 11 ORFs are `conserved_coding_probable_annotation_gap` and 6 `candidate_novel_ncORF`. Five
  of the six are homologous to a human or mouse protein at ORF level, and none has 2 peptides.
  The one without homology is NCORF049411, a 131-aa lncORF on chromosome 13
  (ENSCAFT00845026041) resembling Ly6-family proteins (best mouse hit Ly6g: score 129,
  identity 0.27), supported by one peptide (TCVQSCETLNPMPSR, 9 PSMs in liver and spleen).

`aggregate.py` reproduces those tables byte for byte from the pilot's PSMs and database. No
run has been searched with `driver.sh` yet.

## Caveats

- **Pooled per-run scores.** The driver runs Sage once per run, so each run's discriminant
  score comes from its own model, and `aggregate.py` pools scores from different models as if
  comparable. The pilot searched its three runs together. Rescoring all PSMs together (for
  example with mokapot) would be sounder.
- **One search engine, no validation.** Nothing confirms a match with a second engine,
  predicted spectra or synthetic peptides. A single-peptide ORF is weak evidence.
- **Candidates come from Ensembl 110 models only**: ATG starts, at least 10 residues, on
  annotated transcripts. Non-ATG starts, unannotated transcripts and other assemblies aren't
  searched.
- **Annotation gaps dominate.** Most supported ORFs in the pilot look like conserved proteins
  that dog Ensembl lacks or mis-models, not new ncORFs.
- **Contaminants** are only the 499 in PeptideAtlas's FASTA, plus the standards above.
- **Untested at scale.** `driver.sh` was tested with a stand-in Sage on a local archive, and
  `aggregate.py` on three runs. Memory and run time for all 132 runs are estimates.
