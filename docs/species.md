# Species and references

Resolve names and IDs offline with `orfdata.reference.species_taxon(name)` and
`species_name_for_taxon(taxon_id)`. The latter returns a canonical species name,
or `None` for an unknown positive ID, without requiring an Ensembl annotation.

ORFdata uses the same objects for human, mouse, dog, zebrafish, C. elegans, yeast and custom species.
Every ORF states its NCBI species taxon, assembly and genetic code. A
`ReferenceContext` also pins the annotation provider/version and can retain a
versioned assembly accession, strain taxon/label and sequence-set SHA-256.

```python
from orfdata import ReferenceContext, ScopedTarget

yeast = ReferenceContext(
    taxon_id=4932, species="saccharomyces_cerevisiae",
    genome_build="R64-1-1", assembly_accession="GCA_000146045.2",
    annotation_name="ensembl", annotation_version="115",
    strain="S288C", strain_taxon_id=559292,
)
target = ScopedTarget(yeast, "protein", "YFL039C", "ensembl:115")
print(target.key)
```

Target keys include the reference, source namespace, kind and source identifier.
Identical local IDs in different species, reference versions or datasets remain
separate. Known aliases such as `canis_familiaris` and
`canis_lupus_familiaris` resolve to dog taxon 9615 through a versioned mapping;
the original name stays in exports. Assembly names are never automatically
equated, and strain identity is never discarded. Unknown species require their
own explicit taxonomy context; they do not require a schema change.

`ReferenceContext.from_genome(genome)` resolves available PyEnsembl metadata
offline. A custom Genome without species/taxonomy needs an explicit context:

```python
bound = record.bind(transcript, context=yeast)
# Or: catalog.for_transcript(transcript, context=yeast)
```

The caller declares metadata absent from a custom genome; any metadata already
present must agree. The binder still checks transcript version, strand, complete
spliced path, coding sequence and protein. An ORF model with a
`reference_context` requires that exact context. Older catalogues without a
pinned annotation context can bind an identical transcript in another release;
the actual resolved context and transcript/coding/protein hashes are exported.
No binding downloads a reference or performs liftover.

`context.validate_genome(genome)` checks available metadata against an explicit
context. Genomic alleles can use `include_annotation=False`: their coordinates
are independent of a gene annotation release, while taxonomy, strain and
assembly checks still apply. This does not relax ORF or RNA-product binding.

`ScopedTarget` supports genes, transcripts, ORFs, proteins, ambiguous protein
groups and peptides. A parent's gene or RNA measurement remains at its original
resolution; this identity layer does not copy it onto every ORF. Orthology and
equal peptide sequences do not transfer observations between species.

## Genetic codes

Each ORF declares its NCBI `translation_table`. The binder supports complete
stopped ORFs under table 1 (standard), table 2 (vertebrate mitochondrial) and
table 3 (yeast mitochondrial), each with that table's NCBI initiation codons.
An initiation codon reads as methionine only at the start. Other codes and
initiation codons fail explicitly; there is no fallback to the standard code.
A contig named `MT` never selects a code implicitly: a human mitochondrial ORF
must declare table 2, and a yeast one table 3. See
[NCBI's code definitions](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

`bound.genetic_code` returns the declared code. Tools that build modified
products, such as Varcode and Isovar, must translate them with that code.
For example, OpenProt's IP_306399 on human MT-ND4 begins ATA and ends AGG. It
binds under table 2; table 1 would read those codons as isoleucine and arginine.
