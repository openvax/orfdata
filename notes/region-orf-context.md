# Region context for candidate and partial ORFs

The #61/#64 closure audit found that ScopedTarget.from_orf still dereferences
the optional catalogue model. Use the common ORF ID so annotated and source
placements can participate in quantitative source relationships.

Add region_orf_context to the existing TargetRelation contract. It links a
source region to a candidate/partial ORF while preserving both identities and
explicit mapping evidence. It does not transfer region measurements, establish
translation, or authorize ORF-specific ratios. Multiple candidate relations and
unknown cardinality remain valid. Require left region/right ORF, as with the
other typed context relations; keep reference scope and provenance intact.

Test partial and annotated placements, native IDs, serialization, dataset
storage/query, multiple candidate relationships and refusal of inferred ratios.
Document the evidence semantics and release a compatible patch.
