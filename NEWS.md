# LAMBclassifier 1.0.0

First public release.

## What this release provides

A classifier that assigns breast cancer expression profiles to six classes:
`LHS`, `Apocrine`, `MA-Lo (HER2)`, `Basal-like`, `Chondroid` and `Squamous`.

Scoring uses `singscore` rank enrichment within each sample. A profile gets the
same score whether it is run alone or with thousands of others. The same
specification is used for bulk RNA-seq, microarray, single-cell cluster
pseudobulk and targeted spatial panels.

## Classifier specification 1.0.0

- Signature SHA-256: `8626709901e967f04edfaaa6034e8d67c088da4ab6aa8e4e30e5a0d548989ffa`
- Universe SHA-256: `81454b9a99665ed0ff17892eb66db028d67123d7b1aeb5779e39752694c3c7c0`
- Universe size: 13,468 genes
- `min_margin` 0.02, `min_gene_coverage` 0.60, `min_genes_per_sig` 3

| Signature | Genes |
| :--- | ---: |
| chondroid | 9 |
| squamous | 9 |
| basal_like | 14 |
| luminal | 15 |
| apocrine | 14 |
| malo_her2 | 9 |

Every signature gene is present in the reference universe. Hashes are computed
from the gene lists and universe themselves, not from a version string.
`lamb_specification()` returns them, and `export_results()` writes them with
every set of calls.

The hashes identify the gene lists and the universe. They do not cover the
decision rule, the gate order or the thresholds, which the specification version
identifies. Results are comparable only when the version and both hashes match.

## Calls and gates

A call names a class and implies every exclusion made to reach it. Assigning a
class therefore requires every gate on the path to be evaluable. If any of
`chondroid`, `squamous`, `basal_like`, `luminal` or `apocrine` is unscored, the
sample returns `NC` with the `incomplete_gates` flag, and the
`Unscored_Signatures` column names the signatures that were missing. Missing
MA-Lo genes skip only the MA-Lo check, which refines a call rather than
excluding a class.

`low_coverage` applies when any signature the decision depended on is below the
coverage threshold.

Scoring is performed strictly within the locked universe, which is what the
universe hash on every result refers to. At least 100 universe genes must be
present; below that, scoring is refused.

## Version numbering

Two version numbers are tracked separately:

- **Package version** (`DESCRIPTION`) — the software release. Bug fixes, new
  plotting functions and refactors change this number without changing any call.
- **Specification version** (`lamb_specification()$version`) — the ordered
  signature genes, reference universe, gate order and thresholds. A change here
  changes classifications.

Both are `1.0.0` at this release.

## Validation

Evaluated retrospectively on six public cohorts, 4,759 expression profiles. The
specification and endpoints were set in code before these runs, but there is
no time-stamped preregistration, so this is reported as retrospective external
evaluation rather than prospectively locked validation. No untouched
confirmatory cohort has been evaluated. Per-sample calls, 2x2 cells, Wilson
intervals, cross-classifications, overlap screens and hashes ship in
`inst/validation/`.

Endpoints that did not perform are reported with those that did. Two classes,
`Chondroid` and `Squamous`, have no external validation. Single-cell and spatial
pseudobulk are supported mechanically and are not externally validated. See the
README and `inst/METHODS.md` for the limits.

The training set is an in-house PDX cohort. It is the sole training set, and it
cannot be published yet.

## Quality

`R CMD check` passes on Linux, macOS and Windows via GitHub Actions. The tests
cover each gate of the decision tree, coverage fallback, quality flags, input
checks, single-sample invariance, and a check that the computed hashes match the
ones in `inst/METHODS.md`.

Non-numeric matrices, empty dimensions, missing rownames and all-missing data
raise an error. Duplicated gene symbols raise a warning and only the first
occurrence of each is scored. Values above 100 are treated as untransformed and
`log2(x + 1)` is applied with a warning.

## Scope

LAMBclassifier is a research tool. Its six classes are expression-defined
labels, not clinical diagnoses, and do not replace pathology or ER/PR/HER2
testing. See `inst/METHODS.md` for the locked specification and the stated
limits of the supporting evidence.
