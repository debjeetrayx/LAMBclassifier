# Locked external-validation report

Classifier specification: 1.0.0; signature SHA-256 `8626709901e967f04edfaaa6034e8d67c088da4ab6aa8e4e30e5a0d548989ffa`.

The training set is an in-house PDX cohort. It is the sole training set, and it cannot be published yet. TCGA was not training data and contributed only to gene availability in the fixed scoring universe.

PAM50 results are reported as cross-classification/association between different taxonomies. Clinical results are pre-specified binary concordance endpoints with all 2x2 cells and Wilson intervals.

No external histology cohort with component-level Chondroid or Squamous pathology is included; the current workflow therefore makes no external metaplastic-class validation claim.

SCAN-B technical repeatability was 135/136 (99.3%; Cohen's kappa 0.985). Technical replicates were not counted in clinical concordance.

Overlap screens found 187 near-identical GSE20194/GSE25066 profiles (maximum Pearson r > 0.99 on a deterministic 300-gene screen) and 47 shared normalized GSE20271/GSE25066 identifiers; these cohorts must not be pooled as independent observations.

GSE58812 is used only to describe calls within TNBC. Its `meta` field denotes metastatic outcome, not metaplastic histology.

See the TSV and per-cohort call files in this directory for complete denominators, cells, intervals, hashes, and sample-level results.
