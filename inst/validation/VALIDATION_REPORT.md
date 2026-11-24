# External evaluation report

Classifier specification: 2.0.0; signature SHA-256 `fe7c5dc714d7c6e72b1b70042076b6d067164d18f8b95648a299fe4c9be0336d`.

The training set is not public and is excluded from every endpoint.

PAM50 results are reported as cross-classification/association between different taxonomies. Clinical results are pre-specified binary concordance endpoints with all 2x2 cells and Wilson intervals.

HER2 is reported as an ERBB2-amplicon read-out on every class. The pre-declared HER2 endpoint is HER2-amplicon high (at least 3 of ERBB2, GRB7, MIEN1, PGAP3, STARD3 in the top 9% of the sample) versus the clinical HER2-positive reference.

SCAN-B technical repeatability was 136/136 (100.0%; Cohen's kappa 1.000). Technical replicates were not counted in clinical concordance.

Overlap screens found 187 near-identical GSE20194/GSE25066 profiles (maximum Pearson r > 0.99 on a deterministic 300-gene screen) and 47 shared normalized GSE20271/GSE25066 identifiers; these cohorts must not be pooled as independent observations.

GSE58812 is used only to describe calls within TNBC. Its `meta` field denotes metastatic outcome, not metaplastic histology.

See the TSV and per-cohort call files in this directory for complete denominators, cells, intervals, hashes, and sample-level results.
