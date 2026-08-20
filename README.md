# LAMBclassifier

Luminal-Apocrine-Metaplastic-Basal classifier. It assigns a breast cancer
expression profile to one of six classes: `LHS`, `Apocrine`, `MA-Lo (HER2)`,
`Basal-like`, `Chondroid` or `Squamous`.

Each sample is scored on its own, so the same specification runs on bulk
RNA-seq, microarray, single-cell cluster pseudobulk and spatial panels. External
evaluation covers bulk cohorts only.

---

## How it works

Signatures are built from transcription factors expressed in the normal mammary
epithelium. Each signature is scored with `singscore` rank enrichment over a
fixed 13,468-gene universe. Genes are ranked within one sample, so a sample gets
the same score whether it is run alone or with 3,000 others. There is no cohort
fitting and no batch centring.

The classes are assigned by a four-step decision tree:

1. If the higher of `chondroid` / `squamous` is above all of `basal_like`,
   `luminal` and `apocrine`, call that metaplastic class.
2. Otherwise, if `basal_like` is above both `luminal` and `apocrine`, call
   `Basal-like`.
3. Otherwise call `LHS` if `luminal > apocrine`, and `Apocrine` if not.
4. Inside step 3 only, call `MA-Lo (HER2)` if `malo_her2` is above both
   `luminal` and `apocrine`.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("debjeetrayx/LAMBclassifier")
```

Requires R >= 4.1. `singscore` and `genefu` come from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("singscore", "genefu"))
```

---

## Quick start

```r
library(LAMBclassifier)

# `expr` : genes x samples matrix, HGNC symbols as rownames, log2-transformed
results <- LAMBclassifier(expr)

head(results[, c("Sample", "Pred_Class", "Margin", "Flag")])

export_results(results, "lamb_calls.csv")
```

`export_results()` writes three files: the calls, a `_provenance.tsv` with the
specification version and hashes, and a `_signatures.tsv` with the ordered gene
lists used.

Input must be genes x samples, with HGNC symbols as rownames, on a log2 scale.
Values above 100 are treated as untransformed and `log2(x + 1)` is applied with
a warning. Use `prepare_lamb_data()` to convert raw counts to TPM first.

### Single-cell and spatial data

There is no separate single-cell version. Aggregate cells to pseudobulk, then
classify. Because genes are ranked within a sample, a cluster pseudobulk is
scored on the same scale as a bulk tumour.

**Not externally validated.** No single-cell or spatial cohort with independent
class labels has been evaluated. Calls on such data are exploratory. Check
per-signature coverage and the `Flag` column before interpreting them.

```r
# counts: genes x cells,  clusters: factor of cluster labels
pseudobulk <- vapply(
  levels(clusters),
  function(cl) Matrix::rowSums(counts[, clusters == cl, drop = FALSE]),
  numeric(nrow(counts))
)
log_cpm <- log2(t(t(pseudobulk) / colSums(pseudobulk)) * 1e6 + 1)

cluster_calls <- LAMBclassifier(log_cpm, run_pam50 = FALSE)
```

Classifying each cluster and weighting by cell count gives a class composition
per sample. Aggregating a whole heterogeneous section into one pseudobulk mixes
populations and is not recommended.

### Plots

Each function takes the result of `LAMBclassifier()` and writes files to
`output_prefix`:

```r
plot_comprehensive_heatmap(expr, results, output_prefix = "figures/heatmap")
plot_lineage_tree(results, output_prefix = "figures/tree")
plot_decision_boundaries(results, output_prefix = "figures/boundaries")
plot_expression_scatter(expr, results, "GATA3", "AR", output_prefix = "figures/g_a")
plot_confusion_matrix(results, truth_col = "PAM50", output_prefix = "figures/cm")
generate_standard_heatmaps(expr, results, output_prefix = "figures/standard")
```

Worked examples are in the vignette: `vignette("LAMBclassifier")`.

---

## The six classes

| Class | Description |
| :--- | :--- |
| `Chondroid` | Metaplastic, cartilage-like differentiation |
| `Squamous` | Metaplastic, squamous differentiation |
| `Basal-like` | Basal-like tumours |
| `LHS` | Luminal hormone-sensing |
| `Apocrine` | Molecular apocrine (AR/FOXA1, ER-low) |
| `MA-Lo (HER2)` | Apocrine lineage with ERBB2-amplicon expression |

---

## Reproducibility

Every classification carries a SHA-256 hash of the ordered signature genes and
of the reference universe. The hashes are computed from the gene lists and
universe themselves, not from a version string.

```r
spec <- lamb_specification()
spec$version           # "1.0.0"
spec$signature_sha256  # 8626709901e967f04edfaaa6034e8d67c088da4ab6aa8e4e30e5a0d548989ffa
spec$universe_sha256   # 81454b9a99665ed0ff17892eb66db028d67123d7b1aeb5779e39752694c3c7c0
spec$universe_size     # 13468
```

Two version numbers are tracked separately. The **package version** in
`DESCRIPTION` is the software release. The **specification version** from
`lamb_specification()` covers the signatures, universe, gate order and
thresholds. Both are `1.0.0` here. See [`NEWS.md`](NEWS.md).

The two hashes identify the **signature genes** and the **universe**. They do
not cover the decision rule, which the specification version identifies. Results
are comparable only when the version and both hashes match.

### Quality flags

| Flag | Meaning |
| :--- | :--- |
| `close_call` | Winning margin below 0.02 |
| `low_coverage` | A scored signature had under 60% of its genes present |
| `NC` | A gate needed to reach a call could not be evaluated |
| `incomplete_gates` | Names why a sample is `NC`; see `Unscored_Signatures` |

`low_coverage` considers every signature the decision depended on. Flags do not
remove a sample from the results.

A call names a class and implies every exclusion made to reach it. Assigning a
class therefore requires every gate on the path to be evaluable. If any of
`chondroid`, `squamous`, `basal_like`, `luminal` or `apocrine` is unscored, the
sample returns `NC`. Missing MA-Lo genes skip only the MA-Lo check, which
refines a call rather than excluding a class.

### Targeted panels

Panels carry fewer genes than the full universe. Coverage is recorded per
signature. A signature with fewer than `min_genes_per_sig` genes present is left
unscored, and any sample whose decision path needed it returns `NC`.

Scoring is performed strictly within the locked universe, which is what the
universe hash on every result refers to. At least 100 reference-universe genes
must be present; below that, scoring is refused.

A panel must therefore carry enough of each signature to evaluate every gate.
Coverage thresholds have not been set by systematic dropout experiments; the
current requirement is `min_genes_per_sig` per signature plus 100 universe genes
overall.

### Input checks

Non-numeric matrices, empty dimensions, missing rownames and all-missing data
raise an error. Duplicated gene symbols raise a warning and only the first
occurrence of each is scored.

### Testing

`R CMD check` passes on Linux, macOS and Windows via GitHub Actions. The tests
cover each gate of the decision tree, coverage fallback, quality flags, input
checks, single-sample invariance, and a check that the computed hashes match the
ones in `inst/METHODS.md`.

---

## Validation

Specification 1.0.0 was evaluated retrospectively on six public cohorts, 4,759
expression profiles. The specification and endpoints were set in code before
these runs, but there is no time-stamped preregistration, so this is reported as
**retrospective external evaluation**, not prospectively locked validation.
Confirmatory validation on an untouched cohort has not been done.

The training set is an in-house PDX cohort. It is the sole training set, and it
cannot be published yet. The six cohorts below were used to define which genes
are available in the scoring universe. They were not used to set any decision
rule.

| Cohort | Profiles |
| :--- | ---: |
| GSE20194 | 278 |
| GSE20271 | 178 |
| GSE25066 | 508 |
| GSE41998 | 279 |
| GSE58812 | 107 (TNBC; described only) |
| GSE96058 (SCAN-B) | 3,409 (3,273 primary + 136 technical replicates) |

**Results**

| Endpoint | Cohort | Sensitivity | Specificity | Balanced accuracy |
| :--- | :--- | ---: | ---: | ---: |
| LHS vs ER-positive | SCAN-B (n=2,783) | 0.926 | 0.921 | 0.923 |
| LHS vs ER-positive | GSE41998 (n=279) | 0.852 | 0.918 | 0.885 |
| Basal-like vs triple-negative | GSE41998 (n=278) | 0.814 | 0.848 | 0.831 |
| Basal-like vs triple-negative | SCAN-B (n=2,564) | 0.730 | 0.970 | 0.850 |
| Basal-like vs PAM50 Basal | SCAN-B (n=2,969) | 0.838 | 0.996 | 0.917 |

PPV for LHS vs ER-positive in SCAN-B is 0.993.

**Technical repeatability.** 135 of 136 SCAN-B technical replicate pairs got the
same call (99.3%, Cohen's kappa 0.985).

**PAM50 association.** Cramér's V 0.523 (GSE25066) and 0.587 (SCAN-B). LAMB and
PAM50 are different taxonomies, so this is reported as association, not
accuracy.

**Endpoints that did not perform**

- MA-Lo vs clinical HER2-positive: balanced accuracy 0.49–0.83 across cohorts.
  In GSE25066 it identified 0 of 6 clinically HER2-positive cases.
- Apocrine vs ER-negative: sensitivity 0.03–0.10 in every cohort tested.
- Chondroid and Squamous: no external validation, because no available public
  cohort carries component-level metaplastic pathology.

**Cohort overlap.** A 300-gene correlation screen found 187 near-identical
profiles between GSE20194 and GSE25066 (Pearson r > 0.99), and 47 shared
normalised identifiers between GSE20271 and GSE25066. These cohorts are not
pooled as independent observations. The `meta` field in GSE58812 is metastatic
outcome, not metaplastic histology, and is not used as a metaplastic endpoint.

Per-cohort calls, 2x2 tables, Wilson intervals, cross-classifications, overlap
screens and hashes ship in [`inst/validation/`](inst/validation/). The protocol
is in [`inst/METHODS.md`](inst/METHODS.md).

---

## Limits

- **Training data.** The training set is an in-house PDX cohort. It is the sole
  training set, it cannot be published yet, and it is excluded from every
  endpoint reported above.
- **Metaplastic classes.** No external validation is claimed for `Chondroid` or
  `Squamous`.
- **MA-Lo.** It scores ERBB2-amplicon expression. It does not measure DNA
  amplification or clinical HER2 status.
- **Apocrine.** Sensitivity against an ER-negative reference is low in every
  cohort tested.
- **PAM50.** LAMB and PAM50 are different taxonomies. Comparisons are reported
  as cross-classification and Cramér's V.
- **Research use.** The six classes are expression-defined labels. They are not
  clinical diagnoses and do not replace pathology or ER/PR/HER2 testing.
- **Single-cell and spatial.** Supported mechanically, not externally validated.
  No single-cell or spatial cohort with independent labels has been evaluated.
- **Evidence status.** This is retrospective external evaluation. There is no
  preregistration and no untouched confirmatory cohort. The complete six-class
  system is not a validated classifier; LHS and Basal-like have the strongest
  external support.

---

## Predecessor

LAMBclassifier extends `LABClassifier` (Luminal-Apocrine-Basal), a four-class
classifier (Luminal / Basal / MAlo / MAhi) by Richard Iggo and Élodie Darbo. The
changes are: chondroid and squamous classes added at the first gate,
single-sample rank enrichment in place of cohort-dependent mixture modelling,
and the addition of hashes and quality flags.

---

## Citation

See `citation("LAMBclassifier")`.

## License

GPL (>= 3). See [LICENSE.md](LICENSE.md).

Includes code adapted from `LABClassifier` by Élodie Darbo and Richard Iggo,
licensed under GPL-3. See [NOTICE](NOTICE).
