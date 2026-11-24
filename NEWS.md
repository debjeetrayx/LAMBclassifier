# LAMBclassifier 1.0.1

Classes, scores and the classifier specification (2.0.0) are unchanged. Every
reader was run end to end on real public data.

## Inputs

- **Seurat v5:** a merged object with split count layers is read in full (layers
  joined); before, only the first layer was classified. An SCT or integrated
  default assay is passed over for RNA (or Spatial); log-normalised values in
  the counts slot are refused; the default `orig.ident` identities are no
  longer taken as clusters.
- **AnnData (`.h5ad`):** gene symbols are found when `raw/var` holds only IDs;
  nullable integer and boolean columns are read; repeated cell names keep their
  own clusters.
- **SingleCellExperiment:** cluster columns are found in `colData`; Ensembl row
  names are mapped through `rowData` symbols.
- **Visium:** the folder exactly as downloaded (prefixed `.h5` beside
  `spatial/`) is recognised, and pointing at the matrix itself keeps the spot
  positions. **Visium HD:** the finest available bin other than 2 um is used and
  named. **CosMx and MERSCOPE:** prefixed file names are recognised, files are
  paired by run, and a folder holding several runs asks for one.
- Spots and bins are called spots and bins in the log; `Epithelial_share` is
  left empty when the panel has no epithelial marker.
- Installing needs no separate Bioconductor step (`biocViews` in DESCRIPTION);
  a missing optional package is named with its install command.

## Reports

- The dot legend reads "wins the step / does not win; further right = higher
  score". "How sure" is now "How clear-cut", and its line states the 0.02
  margin rule without implying a probability.
- Marker strip: "percentile among the genes measured in this sample";
  values below 1 print as "<1".
- Decision steps name the strongest program generically ("Strongest
  luminal"), never a class that contradicts the page. Pages without a class
  say "Why no class" and name the unscored signature once.

## Documentation

- The README lists every `lamb_run()` argument; its plotting example runs as
  written; `?lamb_run` has a runnable example.

# LAMBclassifier 1.0.0

First public release. Classifier specification 2.0.0.

## The classifier

Six expression-defined classes, `LHS`, `Apocrine`, `MA-Lo`, `Basal-like`,
`Chondroid` and `Squamous`, and a HER2 amplicon read-out beside every class.
Scoring uses `singscore` rank enrichment within each sample, so a profile gets
the same score whether it is run alone or with thousands of others.

- **HER2 amplicon read-out** (`HER2_Amplicon`): "high" when at least 3 of
  ERBB2, GRB7, MIEN1, PGAP3 and STARD3 rank in the sample's top 9%; "not
  assessable" with fewer than 3 of them measured; "not high" otherwise. It does
  not change the class. It is an expression read-out, not a DNA amplification
  test or a clinical HER2 status.
- **MA-Lo** has ER activity well below LHS and AR activity well below Apocrine,
  and is named after LABClassifier's molecular apocrine MAlo group. Each MA-Lo
  call is reported as LHS-leaning, apocrine-leaning, or balanced when the two scores are within the class margin (`MA_Lo_Lean`) (15 genes: TBX10, BMP7,
  TDO2, SCGB2A1, PHGDH, ACY3, RASD1, CEACAM5, MAPK8IP1, NDUFA4L2, ACOX2, DEFB1,
  RNF145, SCGB2A2, DUSP6). Inside the luminal branch, a sample is MA-Lo when its
  MA-Lo score exceeds both the LHS and Apocrine scores.
- Reports warn when the Chondroid signature is only partly measured on a
  platform (targeted panels).

### Specification 2.0.0

- Signature SHA-256: `fe7c5dc714d7c6e72b1b70042076b6d067164d18f8b95648a299fe4c9be0336d`
  (the six class signatures and the ordered amplicon module)
- Universe SHA-256: `81454b9a99665ed0ff17892eb66db028d67123d7b1aeb5779e39752694c3c7c0`
  (13,468 genes)
- `min_margin` 0.02, `min_gene_coverage` 0.60, `min_genes_per_sig` 3,
  HER2 amplicon q = 9

| Signature | Genes |
| :--- | ---: |
| chondroid | 9 |
| squamous | 9 |
| basal_like | 14 |
| luminal (LHS) | 15 |
| apocrine | 14 |
| malo | 15 |
| HER2 amplicon module | 5 |

The hashes identify the gene lists and the universe; the specification version
identifies the decision rule, gate order and thresholds. Results are comparable
only when the version and both hashes match.

## Use it without programming

- `lamb_run()` with nothing in the brackets asks whether your data is a file or
  a folder, opens the computer's own file picker, runs, and opens the report.
- `lamb_install_launcher()` puts a double-click icon on the desktop (macOS and
  Linux `.command`, Windows `.bat`; on Windows a file can be dropped onto it).
  No extra software beyond R; the data never leave the computer.
- Plain-language messages throughout: what was recognised, what was converted,
  and why an input was refused.

## Inputs

One reader per data type, recognising formats by their contents, never by
column position:

- **Bulk, one patient or a cohort:** CSV / TSV / TXT (gzipped or not; `,` `;`
  tab separators; decimal commas), `.xlsx`; genes in rows or in columns; extra
  annotation columns ignored; Salmon `quant.sf`, kallisto `abundance.tsv`, RSEM
  `genes.results`, featureCounts, STAR `ReadsPerGene.out.tab`, HTSeq-count; a
  folder with one such file per patient; an `.rds` matrix. Counts, TPM / FPKM
  and log2 values are detected and converted.
- **Single-cell:** 10x Cell Ranger folder or `.h5`, `.h5ad` (raw counts),
  Seurat or SingleCellExperiment `.rds`.
- **Spatial:** Xenium, Visium, Visium HD, CosMx (flat files), MERSCOPE. A tissue
  map coloured by class for every platform with positions; H&E overlay and a
  Xenium Explorer cell-group file for Xenium.
- PAM50 runs only with 20 or more samples (`run_pam50 = "auto"`), because it
  centres genes across the samples run together; a single patient is classified
  by LAMB alone.

## Reports

- `LAMB_report.pdf` (bulk, via `lamb_report()`): one page per sample with the
  class and what it means, how sure the call is (margin against the 0.02
  threshold), signature gene coverage, and the decision path showing this
  sample's scores at each step; a cohort overview first.
- Marker strip on every `LAMB_report.pdf` sample page: ESR1, PGR, AR, ERBB2,
  FOXA1, GATA3, KRT5, TP63 and MKI67 as within-sample percentile ranks
  (100 = highest).
- `LAMB_genes.pdf` (via `lamb_gene_report()`): one page per sample (bulk) or
  classified cluster, with every signature gene and the HER2 amplicon genes as
  within-sample percentile ranks; `LAMBclassifier()` attaches them as the
  `"gene_percentiles"` attribute and `lamb_run()` writes
  `technical/LAMB_gene_percentiles.csv`.
- `LAMB_summary.pdf` and `LAMB_clusters.pdf` (single-cell and spatial): the
  slide-level class, class composition, tissue map (when cell positions are
  present), and why each epithelial cluster got its class.
- Single-cell and spatial: a cell-type gate leaves out clusters where an
  immune, vascular, adipocyte or fibroblast gene module outranks the epithelial
  module (top 3 genes per module); `tumour_clusters =` classifies named
  clusters anyway; `clusters =` may name a cell annotation column.
  `Epithelial_share` is reported, not used to decide.
- `LAMB_calls.csv` carries a `Confidence` column (clear / close call / not
  classifiable) and a plain-language note.
- Squamous and Chondroid are described as expression programs, not a
  histological diagnosis of metaplastic carcinoma; the HER2 amplicon read-out
  as an expression signal, not a HER2 IHC/ISH result.

## Example data

`inst/extdata/` ships public data only: six GEO GSE20194 profiles (one per
class) and a single-patient file; a 235-cell subset of GEO GSE176078 (Wu et
al. 2021) in Cell Ranger layout; and a Xenium-layout folder *simulated* from
GSE20194 profiles. Scripts that rebuild them are in
`data-raw/` of the source repository.

## Missing values and messy files

- A missing value (`NA`) in one sample is treated like a gene the platform did
  not measure: the sample is scored on its own measured genes under the same
  3-gene and 60% rules, and flagged when it falls below them.
- Data centred per gene is refused (it changes the order of genes inside a
  sample). Thousands separators are read; summary rows (`Total`, `__…`) are
  dropped; negative (log) data is never logged again; `values =` declares
  what the numbers are. GEO series matrices and non-upper-case gene names are
  refused with the reason.
- Whole numbers that add up to about one million in every sample are read as
  rounded TPM, not counts. Row names that are mostly not human gene symbols or
  Ensembl IDs stop the run with examples of the unrecognised names.
- Every rule is listed with its threshold in `vignette("input-rules")`.

## Calls and gates

A call names a class and implies every exclusion made to reach it, so every gate
on the path must be evaluable. If any of `chondroid`, `squamous`, `basal_like`,
`luminal` or `apocrine` is unscored, the sample returns `NC` with the
`incomplete_gates` flag, and `Unscored_Signatures` names the missing
signatures. Missing MA-Lo genes skip only the MA-Lo check. `low_coverage`
applies when any of the six class signatures is below the coverage threshold
in that sample.
Scoring is restricted to the locked universe; at least 100 universe genes must
be present.

## Version numbering

- **Package version** (`DESCRIPTION`), the software release: `1.0.0`.
- **Specification version** (`lamb_specification()$version`), the signature
  genes, universe, gate order and thresholds: `2.0.0`. A change here changes
  classifications.

## Validation

Evaluated on six public cohorts (4,759 profiles); the specification and
endpoints were fixed in code before these runs. Highlights:

- LHS vs ER-positive in SCAN-B: sensitivity 0.957, specificity 0.893.
- Basal-like vs triple-negative: sensitivity 0.53-0.81, specificity 0.85-0.97.
- HER2 amplicon high vs clinical HER2-positive in SCAN-B: sensitivity 0.685,
  specificity 0.993.
- SCAN-B technical repeatability: 136 of 136 pairs (kappa 1.000).

Full tables are in `inst/validation/`; the protocol is in `inst/METHODS.md`.
The training set is not public and is excluded from every endpoint.
