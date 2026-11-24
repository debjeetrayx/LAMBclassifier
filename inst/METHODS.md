# LAMBclassifier: Locked Method and Validation Policy

This document describes classifier specification **2.0.0**, as shipped in
package **1.0.1**.

Limits: see the package README.

## Development data and independence

The training set is not public. Because it was used for signature-gene
selection, any performance measured on it is internal and is never reported as
independent validation. All endpoints reported for this release come from
external cohorts only.

The 13,468-gene scoring universe ships with the package (`lamb_universe()`).

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
4. In step 3 only, call `MA-Lo` when `malo` exceeds both luminal lineage
   scores.

**MA-Lo** is a luminal-branch tumour in which neither the ER-driven (LHS) nor the
AR-driven (Apocrine) programme dominates, and the MA-Lo programme does (step 4):
ER activity is well below LHS and AR activity well below Apocrine. HER2
amplification is reported separately (`HER2_Amplicon`). The MA-Lo programme is
`TBX10, BMP7, TDO2, SCGB2A1, PHGDH, ACY3, RASD1, CEACAM5, MAPK8IP1, NDUFA4L2,
ACOX2, DEFB1, RNF145, SCGB2A2, DUSP6`. Each MA-Lo call is reported as
LHS-leaning or apocrine-leaning (`MA_Lo_Lean`: the higher of the LHS and
Apocrine scores), or balanced when they differ by less than 0.02, the default
class margin. The name follows the molecular apocrine MAlo group of
LABClassifier.

Chondroid and Squamous are expression labels, not histopathologic diagnoses. An Apocrine call
is not evidence of ER negativity or AR protein positivity without the
corresponding assays.

### Read-out beside the class (specification 2.0.0)

One read-out is reported next to the class. It is not part of the decision
tree and never changes a class.

**HER2 amplicon modifier (`HER2_Amplicon`), every class.** Module:
`ERBB2, GRB7, MIEN1, PGAP3, STARD3`. The sample is restricted to universe
genes; each measured gene gets a within-sample percentile from the top,
`100 * (1 - (rank(v, ties = "average") - 1) / (n_measured - 1))`, where 0 is the
highest and missing values are excluded from `n_measured`. The sample is
`high` when at least 3 of the 5 module genes have percentile `<= q`, with
`q = 9` (fixed, stored in the calibration as `her2_amplicon$q`); `not
assessable` when fewer than 3 module genes are measured; `not high` otherwise.
The counts behind the call are in `HER2_Amplicon_N_High` and
`HER2_Amplicon_N_Measured`. This is an expression read-out, not DNA
amplification and not clinical HER2 status.

**Chondroid coverage warning.** On a
targeted panel (fewer than 1,000 universe genes measured), when fewer than all 9
of its genes are measured, reports warn:
"Chondroid rests on X of 9 genes on this platform; read Chondroid calls (and
their absence) with care."

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

`close_call` (margin below 0.02) and `low_coverage` (any of the six class
signatures below 60% coverage in that sample) are quality flags, not post-hoc
exclusion criteria.

Both thresholds were set before the evidence below was computed and were not
changed by it (`inst/validation/threshold_evidence.R` and `.tsv`, public SCAN-B
data).

- **Margin (0.02).** In 136 tumours sequenced twice, the margin changed by
  0.006 at the median, 0.017 at the 90th percentile, 0.035 at the 95th and
  0.160 at the 99th. The LHS minus Apocrine difference behind the MA-Lo lean
  changed by 0.022 at the 90th percentile, so about 1 lean in 10 can cross the
  0.02 band on re-sequencing. Separate biopsies vary more than this.
- **Coverage (60%).** The same random genes were removed from 600 SCAN-B
  profiles (21 runs, 12,600 calls) and each call compared with the call on
  full data. With every class signature at 30% coverage or more, no call was
  refused and 1.3% (70% or more), 2.7% (60-69%), 3.3% (50-59%) and 3.8%
  (30-49%) changed class. Below 30%, 56% were refused (NC) and 2.9% changed
  class. There is no step at 60%: among scored calls, flagged calls kept their
  class 94.9% of the time and unflagged calls 98.2%.

### What the identifiers cover

The two SHA-256 values identify the **ordered signature genes** and the
**reference universe**. The signature hash covers the
six scored signatures (`chondroid, squamous, basal_like, luminal, apocrine,
malo`) followed by the ordered HER2 amplicon module. They do not cover the
decision rule, the gate order or the thresholds (including `q`), which
are identified by the **specification version**. Results are comparable only
when the version and both hashes match.

Specification 2.0.0 has SHA-256 identifiers:

- signature SHA-256:
  `fe7c5dc714d7c6e72b1b70042076b6d067164d18f8b95648a299fe4c9be0336d`
- universe SHA-256 (unchanged from 1.0.0):
  `81454b9a99665ed0ff17892eb66db028d67123d7b1aeb5779e39752694c3c7c0`

The files in `inst/validation/` were produced by specification 2.0.0. For
reference, the earlier specification 1.0.0 had signature SHA-256
`8626709901e967f04edfaaa6034e8d67c088da4ab6aa8e4e30e5a0d548989ffa` and the same
universe.

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
compatibility summaries: LHS versus Luminal A/B and Basal-like versus Basal.
These summaries test convergent validity,
not interchangeability.

### Clinicopathologic concordance

Clinical annotations are assessed with explicit binary mappings, cohort-level
2x2 counts, sensitivity, specificity, PPV, NPV, balanced accuracy, and 95%
Wilson intervals. Pre-specified mappings are:

- LHS prediction versus ER-positive reference;
- HER2-amplicon high versus clinical HER2-positive
  reference, per cohort, in GSE20194, GSE20271, GSE25066, GSE41998 and
  GSE96058; full 2x2, sensitivity, specificity, PPV, NPV and Wilson
  95% intervals. GSE20194 and GSE25066 overlap and are never pooled.
- Basal-like prediction versus ER-/PR-/HER2- triple-negative reference; and
- Apocrine prediction versus ER-negative reference (supportive only).

Because these endpoints are correlated and class prevalences differ across
cohorts, raw agreement alone is insufficient and confirmatory claims require a
pre-specified multiplicity strategy and sample-size justification.

### Histology

A metaplastic endpoint requires an actual pathology annotation. GSE58812 is a
triple-negative breast cancer cohort; its `meta` field is metastatic outcome,
not metaplastic histology, so it must not be used as a metaplastic validation
set.

## Required reporting

For every evaluation, retain accession, platform, eligibility rules, sample
flow, missingness, biological and technical replicate handling, preprocessing,
classifier hashes, full prediction file, all 2x2 cells, uncertainty intervals,
and all planned sensitivity analyses. Report unfavorable results and failed
endpoints. Do not select cohorts or endpoints after inspecting performance.
