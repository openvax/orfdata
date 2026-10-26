# Transcript-bound ORFs and incomplete coding boundaries

Scope: ORFdata #61/#64 and Varcode #638. A TranscriptORF is a coding placement on one RNA, including an annotated CDS. A catalogue record is optional provenance, not a prerequisite for consequence calculation.

## Representation

Keep the original RNA, exact coding interval and coding blocks. Add an explicit first-complete-codon offset for bounded fragments; retain leading/trailing partial-codon bases. Start and stop completeness are independent true/false/unknown values. Reading frame and boundary claims are separate: a fragment can have a known local frame without a known initiation site. Complete catalogue ORFs keep their strict validation.

Support bounded source intervals through TranscriptORF.at without replacing their parent RNA with a trimmed transcript. Source-supported non-AUG starts retain their declaration. Reject an incompatible phase, false stop declaration, reference mismatch or contradictory supplied protein. Missing sequence is never padded. Preserve the definition and evidence in a checked serialization that can be restored against an archived RNA.

## Consumers

Varcode uses common ORF fields rather than assuming record/model are present. Complete, ordinary coding placements use its existing transcript-model engine. Partial boundaries retain explicit unavailable start-loss, stop-loss or full-product results; a validated local coding change may still be calculated in the declared frame. Unsupported recoding, circular or unmapped models remain accessible with a structured reason.

Native JSON and Proteome exports retain the source ORF, parent RNA, source/host IDs, phase, boundaries, evidence and unsupported outcomes. No fabricated catalogue evidence or annotated transcript identity is introduced.

## Validation

Reuse pinned SNURF/SNRPN, MIEF1/MIURF, non-AUG and alternative-start examples. Add native partial annotations and synthetic both-strand/splice-junction controls for phase and exact boundaries. Exercise the ordinary public effect API, source-ID/reference isolation, JSON reconstruction and joint producer export. Check human and yeast declared-code fixtures independently through Varcode and Isovar.

Primary conventions: INSDC feature-table v11.4 (partial feature locations and codon_start), https://www.insdc.org/submitting-standards/feature-table/; NCBI genetic codes, https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.
