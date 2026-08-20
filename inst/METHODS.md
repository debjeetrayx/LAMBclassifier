# LAMBclassifier: Locked Method and Validation Policy

This document describes classifier specification **1.0.0**, as shipped in
package **1.0.0**. The external evaluation distributed with this release was run
against specification 1.0.0.

LAMB is a research classifier. Its expression-defined classes are not
clinical diagnoses and must not replace pathology, ER/PR/HER2 testing, or a
validated clinical assay.

## Development data and independence

The training set is an in-house PDX cohort. It is the sole training set, and it
cannot be published yet.

Because that cohort was used for signature-gene selection, any performance
measured on it is internal and is never reported as independent validation. All
endpoints reported for this release come from external cohorts only.

TCGA was not used to choose signature genes or train decision rules. It was one
of six cohorts used only to define gene availability in the 13,468-gene scoring
universe, together with METABRIC, ICGC, EORTC and TRANSTART. Consequently, TCGA
results are characterization rather than training performance. Because TCGA
contributed to the universe, it is not fully independent evaluation of the
complete frozen pipeline.

## Locked scoring specification

For every sample, genes in the fixed 13,468-gene universe are ranked and each
signature is scored with `singscore::simpleScore()`. Because ranking is within
sample, a sample receives the same score whether processed alone or with other
samples, provided gene identifiers and preprocessing are unchanged.

The decision rule is:

1. If `max(chondroid, squamous)` exceeds
   `max(basal_like, luminal, apocrine)`, call the higher metaplastic signature.
2. Otherwise, if `basal_like` exceeds both `luminal` and `apocrine`, call
   `Basal-like`.
3. Otherwise, call `LHS` when `luminal > apocrine`, and `Apocrine` otherwise.
4. In step 3 only, call `MA-Lo (HER2)` when `malo_her2` exceeds both luminal
   lineage scores.

The MA-Lo rule detects ERBB2-amplicon **expression**, not DNA amplification and
not clinical HER2 status. Chondroid and Squamous are expression labels, not
histopathologic diagnoses. An Apocrine call is not evidence of ER negativity or
AR protein positivity without the corresponding assays.

With fewer than `min_genes_per_sig` present genes (default 3), that signature is
not scored.

A call names a class and implies every exclusion made to reach it. If any of
`chondroid`, `squamous`, `basal_like`, `luminal` or `apocrine` is unscored, the
corresponding gate cannot be evaluated, and the sample is returned as `NC` with
the `incomplete_gates` flag. The `Unscored_Signatures` column records which
signatures were missing. Missing MA-Lo genes skip only the MA-Lo check, because
that check refines a call rather than excluding a class.

Scoring is performed strictly within the locked universe, which is what the
universe hash on every result refers to. At least 100 universe genes must be
present; below that, scoring is refused.

`close_call` (margin below 0.02) and `low_coverage` (any decision-relevant
signature below 60% coverage) are quality flags, not post-hoc exclusion
criteria.

### What the identifiers cover

The two SHA-256 values identify the **ordered signature genes** and the
**reference universe**. They do not cover the decision rule, the gate order or
the thresholds, which are identified by the **specification version**. Results
are comparable only when the version and both hashes match.

The default ordered gene lists and fixed universe have SHA-256 identifiers:

- signature SHA-256:
  `8626709901e967f04edfaaa6034e8d67c088da4ab6aa8e4e30e5a0d548989ffa`
- universe SHA-256:
  `81454b9a99665ed0ff17892eb66db028d67123d7b1aeb5779e39752694c3c7c0`

Call `lamb_specification()` to retrieve the lists and identifiers. Every
classification row records them, and `export_results()` writes companion
provenance and ordered-signature files.

## Independent validation policy

The classifier, signatures, order of gates, thresholds, preprocessing, class
mappings, and endpoints must be frozen before viewing an external cohort's
outcomes. A cohort used to alter any of these is thereafter development data
and needs a new untouched cohort for independent validation.

Validation is performed per cohort; estimates must not be naively pooled across
platforms or overlapping studies. Patient/sample identifiers are checked for
overlap before synthesis, and biological replicates—not expression columns—are
the unit used for clinical performance. Technical replicates are reserved for
repeatability analysis.

Primary analyses include every evaluable sample. Missing reference labels are
reported and excluded endpoint-by-endpoint, with denominators shown. `NC`,
`close_call`, and `low_coverage` calls are retained in the primary analysis;
flag-restricted analyses may be reported only as sensitivity analyses.

### PAM50 comparison

PAM50 and LAMB are different taxonomies. Their comparison is a
cross-classification, not multiclass accuracy with PAM50 treated as truth.
Report the complete table, row/column proportions, Cramer's V, and pre-specified
compatibility summaries such as LHS versus Luminal A/B, Basal-like versus
Basal, and MA-Lo versus HER2-enriched. These summaries test convergent validity,
not interchangeability.

### Clinicopathologic concordance

Clinical annotations are assessed with explicit binary mappings, cohort-level
2x2 counts, sensitivity, specificity, PPV, NPV, balanced accuracy, and 95%
Wilson intervals. Pre-specified mappings are:

- LHS prediction versus ER-positive reference;
- MA-Lo prediction versus clinical HER2-positive reference;
- Basal-like prediction versus ER-/PR-/HER2- triple-negative reference; and
- Apocrine prediction versus ER-negative reference (supportive only).

Because these endpoints are correlated and class prevalences differ across
cohorts, raw agreement alone is insufficient and confirmatory claims require a
pre-specified multiplicity strategy and sample-size justification.

### Histology

A metaplastic endpoint requires an actual pathology annotation. GSE58812 is a
triple-negative breast cancer cohort; its `meta` field is metastatic outcome,
not metaplastic histology, so it must not be used as a metaplastic validation
set. The current external workflow contains no cohort with component-level
Chondroid or Squamous pathology and therefore makes no external validation
claim for those two LAMB classes.

## Required reporting

For every evaluation, retain accession, platform, eligibility rules, sample
flow, missingness, biological and technical replicate handling, preprocessing,
classifier hashes, full prediction file, all 2x2 cells, uncertainty intervals,
and all planned sensitivity analyses. Report unfavorable results and failed
endpoints. Do not select cohorts or endpoints after inspecting performance.
