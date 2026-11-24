# Validation outputs — classifier specification 2.0.0

These are the complete outputs behind the numbers quoted in the package README.
They were produced by specification 2.0.0 and carry its hashes, so any file here
can be matched against the classifier that generated it.

Every `*_calls.csv` records, per sample, the predicted class, the winning
margin, any quality flag, the metaplastic score and all six signature scores,
the HER2 amplicon read-out, and the specification version and signature and
universe SHA-256 hashes. Sample identifiers are the GEO sample accessions or,
for GSE96058, the sample titles used by that series.

## Files

| File | Contents |
| :--- | :--- |
| `cohort_inventory.tsv` | Profiles, metadata rows, replicate counts and preprocessing per cohort |
| `binary_concordance_metrics.tsv` | All pre-specified endpoints: full 2x2 cells, sensitivity, specificity, PPV, NPV, balanced accuracy, Wilson 95% intervals |
| `pam50_association.tsv` | Cramér's V per cohort, reported as association between non-equivalent taxonomies |
| `*_pam50_counts.tsv`, `*_pam50_by_lamb.tsv` | Full LAMB x PAM50 cross-classifications and row proportions |
| `GSE96058_technical_repeatability.tsv` | Agreement and Cohen's kappa over 136 technical replicate pairs |
| `GSE96058_technical_repeat_pairs.tsv` | The individual replicate pairs |
| `cohort_overlap_checks.tsv` | Near-duplicate and shared-identifier screens between cohorts |
| `threshold_evidence.tsv`, `threshold_evidence.R` | Score change between technical repeats, and same / different / refused calls under gene loss, behind the 0.02 margin and the 60% coverage flag |
| `GSE58812_TNBC_distribution.tsv` | Class distribution within TNBC; descriptive only |
| `*_calls.csv` | Per-sample predictions, scores, flags and hashes |
| `*_calls_provenance.tsv` | Specification version, hashes, universe size, thresholds, HER2 amplicon genes and q, export timestamp |
| `*_calls_signatures.tsv` | The exact ordered gene lists of the six class signatures used for that run |
| `run_provenance.tsv` | Specification version, hashes, universe size, run time and R version |
| `VALIDATION_REPORT.md` | Narrative summary of the evaluation |

## Reproducing

The protocol is in `inst/METHODS.md`. Evaluation requires the six public GEO
series listed in `cohort_inventory.tsv`, which are not redistributed here.

## Not included

The training set is not public. It is excluded from every endpoint reported in
this directory.
