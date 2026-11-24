# LAMBclassifier

Luminal-Apocrine-Metaplastic-Basal classifier. It assigns a breast cancer
expression profile to one of six classes: `LHS`, `Apocrine`, `MA-Lo`,
`Basal-like`, `Chondroid` or `Squamous`, and reports a HER2 amplicon read-out
beside every class.

Each sample is scored on its own, so the same specification runs on bulk
RNA-seq, microarray, single-cell cluster pseudobulk and spatial panels.

---

## How it works

Each class has a fixed signature of 9 to 15 genes. Each signature is scored
with `singscore` rank enrichment over a fixed 13,468-gene universe. Genes are
ranked within one sample, so a sample gets the same score whether it is run
alone or with thousands of others. There is no cohort fitting and no batch
centring.

The classes are assigned by a four-step decision tree:

1. If the higher of `chondroid` / `squamous` is above all of `basal_like`,
   `luminal` and `apocrine`, call that metaplastic class.
2. Otherwise, if `basal_like` is above both `luminal` and `apocrine`, call
   `Basal-like`.
3. Otherwise call `LHS` if `luminal > apocrine`, and `Apocrine` if not.
4. Inside step 3 only, call `MA-Lo` if `malo` is above both `luminal` and
   `apocrine`.

Beside the class, **`HER2_Amplicon`** is reported for every class and never
changes it: "high" when at least 3 of `ERBB2, GRB7, MIEN1, PGAP3, STARD3` rank
in the sample's top 9% of universe genes; "not assessable" when fewer than 3 of
them are measured; "not high" otherwise. Expression only, not a HER2 IHC/ISH
result.

---

## Installation

Requires R >= 4.1. In R:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("singscore", "genefu"))
remotes::install_github("debjeetrayx/LAMBclassifier")
```

Optional, only for some inputs: `readxl` (Excel), `rhdf5` or `hdf5r` (10x `.h5`,
`.h5ad`), `arrow` (Visium HD), `SeuratObject` / `SingleCellExperiment`, and
`png`, `jpeg` or `terra` (H&E image).

A `Dockerfile` with pinned R and Bioconductor versions is in the repository.

---

## Run it without programming

```r
LAMBclassifier::lamb_run()
```

With nothing in the brackets, it asks whether your data is a file or a folder,
opens your computer's own file picker, runs, and opens the report. Nothing else
needs to be set. The data stay on your computer.

For a desktop icon that does the same without opening R, run this once:

```r
LAMBclassifier::lamb_install_launcher()
```

To give the path directly, or from a terminal:

```r
LAMBclassifier::lamb_run("path/to/your/data")
```
```sh
Rscript -e 'LAMBclassifier::lamb_run("path/to/your/data")'
```

Try it on the bundled public example (one patient):

```r
LAMBclassifier::lamb_run(system.file("extdata", "example_patient.csv", package = "LAMBclassifier"),
                         out_dir = tempfile("lamb_"))
```

Without `out_dir`, results go into a new folder next to the input
(`LAMB_results_<name>`).

### What it accepts

Formats are recognised from the file contents, never from column positions. When
an input is ambiguous, it stops and says why in plain words rather than guessing.

| Give it | Recognised as |
| :--- | :--- |
| A table: `.csv`, `.tsv`, `.txt` (gzipped or not; `,` `;` or tab separated; decimal commas), `.xlsx`, or an `.rds` matrix | bulk, one patient or many; genes in rows **or** columns; HGNC symbols or Ensembl IDs; extra annotation columns ignored |
| One quantifier file: Salmon `quant.sf`, kallisto `abundance.tsv`, RSEM `genes.results`, featureCounts, STAR `ReadsPerGene.out.tab`, HTSeq-count | bulk, one patient |
| A folder with one such file per patient | bulk cohort |
| A 10x Cell Ranger folder or `.h5`; an `.h5ad` file with raw counts; a Seurat or SingleCellExperiment `.rds` | single-cell |
| A Xenium, Visium, Visium HD, CosMx (flat files) or MERSCOPE output folder | spatial |


Every decision LAMB makes about a file (which column holds the genes, whether the
table is turned round, whether values are counts or log, what happens to duplicates
and missing values) is listed with its threshold in `vignette("input-rules")`, with
the run-log message it writes when it applies.
Raw counts, TPM / FPKM and log2 values are detected and converted. PAM50 is added
only for 20 or more samples (it centres each gene across the samples run
together); a single patient is classified by LAMB alone.

Single-cell and spatial data are scored per group of cells: the clusters the
instrument wrote (10x, Xenium, Visium), the active identities (Seurat),
`colLabels()` (SingleCellExperiment) or a `leiden` / `louvain` / `cluster`
column (`.h5ad`). `clusters =` can instead name a cell annotation column
(for example `clusters = "cell_type"`) or a CSV file. CosMx and MERSCOPE carry
no clusters; give your own with `lamb_run("folder", clusters = "my_clusters.csv")`
(two columns: cell id, group).

**Xenium with the H&E (optional).** If the post-run H&E was aligned in Xenium
Explorer, give the image and the alignment file Explorer exported:

```r
lamb_run("xenium_output_folder",
         he_image     = "HE.ome.tif",
         he_alignment = "HE_alignment_files/matrix.csv")
```

PNG / JPEG and standard TIFF / OME-TIFF images are read. For a JPEG2000-compressed
OME-TIFF, export a smaller PNG (QuPath: File > Export images, downsample 16) and
add `he_downsample = 16`. Without an H&E the classes are
drawn on a plain background. The run log reports how many cells land on the
image; below 90 % it warns that the alignment may belong to another section.

### What you get

| File | Bulk | Single-cell | Spatial |
| :--- | :---: | :---: | :---: |
| `LAMB_report.pdf`: one page per sample — the class and what it means, how sure the call is, gene coverage, the decision path that explains it, and a marker strip | ✓ | | |
| `LAMB_genes.pdf`: one page per sample or classified cluster — every gene LAMB scores, grouped by signature and the HER2 amplicon | ✓ | ✓ | ✓ |
| `LAMB_calls.csv`: the class per sample or per cluster, the HER2 amplicon read-out, how sure the call is, and a plain-language note (for cells also the cluster size, `Cell_type` and `Epithelial_share`) | ✓ | ✓ | ✓ |
| `LAMB_summary.pdf`: the slide-level class and the percentage of epithelial cells in each class | | ✓ | ✓ + tissue map when cell positions are present (on the H&E for Xenium, if given) |
| `LAMB_clusters.pdf`: why each epithelial cluster got its class | | ✓ | ✓ |
| `LAMB_slide.csv`: the slide-level class and the percentages | | ✓ | ✓ |
| `LAMB_xenium_explorer.csv`: load in Xenium Explorer (Cells → Cell groups) | | | Xenium |
| `technical/`: scores, signature coverage, provenance hashes, run log | ✓ | ✓ | ✓ |
| `technical/LAMB_gene_percentiles.csv`: each scored gene's rank within the sample or cluster | ✓ | ✓ | ✓ |

**Marker strip.** Each bulk page shows ESR1, PGR, AR, ERBB2, FOXA1, GATA3,
KRT5, TP63 and MKI67 as their rank within the sample (100 = highest of the
measured genes): RNA levels comparable across platforms, not IHC results.

**How sure is the call?** *Clear* when the winning class led by at least 0.02 at
the deciding step; *close call* when two classes scored almost the same;
*not classifiable* when a required signature could not be scored on the
platform.

**Single-cell and spatial.** LAMB scores groups of cells (pseudobulk: summed
counts, log2(CPM + 1)), never single cells. The classes describe breast
epithelium, so each cluster first passes a **cell-type gate**: in its
pseudobulk, gene modules for epithelium, immune cells, vessels, adipocytes and
fibroblasts are compared (mean within-cluster percentile rank of each module's 3
highest-ranked genes; a module needs 3 measured genes). A cluster is classified
unless another module ranks above the epithelial one; such clusters are named by
what dominates (immune, vascular, adipocyte, fibroblast-like), marked "Not
epithelial" and left out of the summary. Metaplastic tumour cells can lose
epithelial genes: name those clusters with `tumour_clusters =` to classify them.
`Epithelial_share` (the share of cells detecting an epithelial gene) is reported
but not used to decide. Clusters under 50 cells are listed but not classified.
The **slide-level class** is the class holding most classified epithelial
cells; when no class reaches 50 % it reads "mixed", naming every class at 10 %
or more (for example "mixed: LHS 33% / MA-Lo 33% / Basal-like 33%"). A Visium spot holds several cells, often of
different types, so a Visium class describes a cluster of spots. On targeted
panels some signatures are partly missing, so calls can carry a coverage note.

For scripted use, `LAMBclassifier(expr)` classifies a genes x samples log2
matrix directly, `lamb_report(results, "report.pdf")` draws the per-sample
report, and `export_results()` writes the full output.


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
| `Apocrine` | Molecular apocrine (AR / FOXA1 programme) |
| `MA-Lo` | A luminal-branch tumour in which neither the ER-driven (LHS) nor the AR-driven (Apocrine) programme dominates, and the MA-Lo programme does: ER activity well below LHS, AR activity well below Apocrine. Reported as LHS-leaning, apocrine-leaning, or balanced when the two scores are within 0.02; named after LABClassifier's MAlo group |

The HER2 amplicon is not a class: `HER2_Amplicon` is reported for every class.

---

## Reproducibility

Every classification carries a SHA-256 hash of the ordered signature genes and
of the reference universe. The hashes are computed from the gene lists and
universe themselves, not from a version string.

```r
spec <- lamb_specification()
spec$version           # "2.0.0"
spec$signature_sha256  # fe7c5dc714d7c6e72b1b70042076b6d067164d18f8b95648a299fe4c9be0336d
spec$universe_sha256   # 81454b9a99665ed0ff17892eb66db028d67123d7b1aeb5779e39752694c3c7c0
spec$universe_size     # 13468
```

Two version numbers are tracked separately. The **package version** in
`DESCRIPTION` is the software release. The **specification version** from
`lamb_specification()` covers the signatures, universe, gate order and
thresholds. The package is `1.0.0`; the specification is `2.0.0`. See
[`NEWS.md`](NEWS.md).

The two hashes identify the **signature genes** and the **universe**. They do
not cover the decision rule, which the specification version identifies. Results
are comparable only when the version and both hashes match.

### Quality flags

| Flag | Meaning |
| :--- | :--- |
| `close_call` | Winning margin below 0.02 |
| `low_coverage` | One of the six class signatures had under 60% of its genes measured in that sample |
| `incomplete_gates` | The sample is `NC` because a gate could not be evaluated; see `Unscored_Signatures` |

`NC` (not classifiable) is returned as the class when a gate needed to reach a
call could not be evaluated. `low_coverage` checks all six class signatures,
whatever path the decision took. Flags do not remove a sample from the results.

A call names a class and implies every exclusion made to reach it. Assigning a
class therefore requires every gate on the path to be evaluable. If any of
`chondroid`, `squamous`, `basal_like`, `luminal` or `apocrine` is unscored, the
sample returns `NC`. Missing MA-Lo genes skip only the MA-Lo check, which
refines a call rather than excluding a class. The HER2 amplicon is a read-out
and never makes a sample `NC`.

### Targeted panels

Panels carry fewer genes than the full universe. Coverage is recorded per
signature. A signature with fewer than `min_genes_per_sig` genes present is left
unscored, and any sample whose decision path needed it returns `NC`.

Scoring is performed strictly within the locked universe, which is what the
universe hash on every result refers to. At least 100 reference-universe genes
must be present; below that, scoring is refused.

A panel must therefore carry at least `min_genes_per_sig` genes of each
signature, and 100 universe genes overall.

### Input checks

In `LAMBclassifier()`, non-numeric matrices, empty dimensions, missing
rownames and all-missing data raise an error. Duplicated gene symbols raise a
warning and only the first occurrence of each is scored. `lamb_run()` collapses
duplicates before that (see `vignette("input-rules")`).

---

## Validation

Specification 2.0.0 was evaluated on six public cohorts, 4,759 expression
profiles. The specification and endpoints were set in code before these runs.

The training set is not public. It is excluded from every endpoint below.

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
| LHS vs ER-positive | SCAN-B (n=2,783) | 0.957 | 0.893 | 0.925 |
| LHS vs ER-positive | GSE41998 (n=279) | 0.843 | 0.883 | 0.863 |
| Basal-like vs triple-negative | GSE41998 (n=278) | 0.814 | 0.848 | 0.831 |
| Basal-like vs triple-negative | SCAN-B (n=2,564) | 0.730 | 0.970 | 0.850 |
| Basal-like vs PAM50 Basal | SCAN-B (n=2,969) | 0.838 | 0.996 | 0.917 |
| HER2 amplicon high vs clinical HER2-positive | SCAN-B (n=2,868) | 0.685 | 0.993 | 0.839 |
| HER2 amplicon high vs clinical HER2-positive | GSE20194 (n=278) | 0.576 | 0.995 | 0.786 |

PPV for LHS vs ER-positive in SCAN-B is 0.991; for HER2 amplicon high vs
clinical HER2-positive, 0.938 (SCAN-B) and 0.971 (GSE20194).

**Technical repeatability.** All 136 SCAN-B technical replicate pairs got the
same call (100%, Cohen's kappa 1.000).

**PAM50 association.** Cramér's V 0.516 (GSE25066) and 0.541 (SCAN-B). LAMB and
PAM50 are different taxonomies, so this is reported as association, not
accuracy.

**Cohort overlap.** GSE20194 and GSE25066 share 187 near-identical profiles,
and GSE20271 and GSE25066 share 47 identifiers. Each cohort is reported on its
own; they are not pooled.

Per-cohort calls, 2x2 tables, Wilson intervals, cross-classifications, overlap
screens and hashes ship in [`inst/validation/`](inst/validation/). The protocol
is in [`inst/METHODS.md`](inst/METHODS.md).

---

## Limits

- **HER2 amplicon on microarray data.** Sensitivity against clinical HER2
  status was 0.36–0.58 on microarray cohorts (0 of 6 in GSE25066), against
  0.685 on RNA-seq (SCAN-B). A "not high" read-out on microarray data is weak
  evidence. Specificity was 0.974–0.995 in all five cohorts with HER2 status.
- **LHS on ER-negative tumours.** In one cohort (GSE20271), 26 of 80
  clinically ER-negative tumours were called LHS (8–12% in the four other
  cohorts with ER status). These tumours expressed FOXA1, GATA3 and ER target
  genes such as TFF1, with low ESR1: LHS identifies a luminal hormone-sensing
  programme, not ER protein. An LHS call does not replace ER
  immunohistochemistry.
- **Apocrine versus ER-negative.** Sensitivity was low, 0.04–0.15.
- **MA-Lo, Chondroid and Squamous.** No external validation.
- **Single-cell and spatial.** Calls are not externally validated.
- **Evaluation.** The six cohorts had been seen under the earlier
  specification (1.0.0), so the evaluation is retrospective and not on
  untouched data.
- **Research use only.** A class is not a clinical diagnosis and does not
  replace pathology or ER/PR/HER2 testing.
- **Two inputs that change results.** Zeros used to mean "missing" are scored
  as silent genes; store missing values as empty or `NA`. Whole numbers that
  are not counts (rounded TPM or intensities) are read as counts; say what they
  are with `values = "tpm"` or `values = "log2"`.

---

## Predecessor and acknowledgement

LAMBclassifier is built on the concept of
[LABClassifier](https://github.com/ElodieDarbo/LABClassifier)
(Luminal-Apocrine-Basal), a four-class classifier (Luminal / Basal / MAlo /
MAhi). Credit for that concept goes to Richard Iggo. Code in LAMBclassifier is
adapted from LABClassifier under the GNU General Public License, version 3 or
later.

The changes are: chondroid and squamous classes added at the first gate;
singscore rank enrichment with fixed score comparisons in place of ssGSEA
scores and a cohort-fitted mixture-model cut-off (which LABClassifier uses to
separate MAhi from MAlo); the HER2 amplicon read-out; and hashes and quality
flags.

---

## Citation

See `citation("LAMBclassifier")`.

## License

GPL (>= 3). See [LICENSE.md](LICENSE.md).

Includes code adapted from
[LABClassifier](https://github.com/ElodieDarbo/LABClassifier), copyright
Richard Iggo, licensed under GPL-3. See [NOTICE](NOTICE).
