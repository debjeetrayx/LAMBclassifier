# Example data shipped with LAMBclassifier

Every file here comes from a public dataset, or is simulated from one and says
so in its name. No patient-identifiable or unpublished data are
included. The scripts that rebuild each example are in `data-raw/` of the
source package.

Find a file after installation with
`system.file("extdata", "<name>", package = "LAMBclassifier")`.

| Path | What it is | Real or simulated | Size |
|---|---|---|---|
| `example_bulk.csv` | 6 breast tumour microarray profiles, 10,104 genes x 6 samples, log2 | Real (GEO GSE20194) | ~550 KB |
| `example_patient.csv` | One of those profiles alone (GSM505471) | Real (GEO GSE20194) | ~145 KB |
| `example_10x/` | One ER+ tumour, 235 cells, Cell Ranger layout | Real (GEO GSM5354516 / GSE176078) | ~820 KB |
| `xenium_sim/` | Xenium-layout folder, 450 cells on a 2 x 1.5 mm field | **SIMULATED** from GSE20194 profiles | ~140 KB |

These examples show how the software runs.

## example_bulk.csv and example_patient.csv

- **Source:** GEO series [GSE20194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20194)
  (MAQC-II breast cancer, fine-needle aspirates, Affymetrix HG-U133A / GPL96,
  MAS5, log2).
- **Citation:** Popovici V, et al. Effect of training-sample size and
  classification difficulty on the accuracy of genomic predictors. *Breast
  Cancer Res* 2010;12:R5. doi:10.1186/bcr2468. MAQC Consortium, *Nat
  Biotechnol* 2010;28:827-838. doi:10.1038/nbt.1665.
- **Terms:** GEO data are public; NCBI places no restrictions on their use or
  distribution ([GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)).
  Cite the papers above.
- **Subset:** probes collapsed to genes (first GPL96 symbol, mean of probes,
  15 significant digits), then restricted to the genes of `lamb_universe()`.
  Scoring happens only inside that universe, so the calls and every score are
  identical to a run on all 13,237 genes. `data-raw/example_bulk.R` checks
  this. Do not round the values further: rounding creates ties that change ranks.
- **Samples** (LAMB call with `run_pam50 = FALSE`; receptor status from the GEO record):

| Sample | LAMB class | Margin | ER / PR / HER2 (GEO) |
|---|---|---|---|
| GSM505471 | LHS | 0.4718 | + / + / - |
| GSM505546 | MA-Lo | 0.1593 | - / - / + |
| GSM505331 | Apocrine | 0.4948 | - / - / - |
| GSM505423 | Basal-like | 0.5859 | - / - / - |
| GSM505491 | Squamous | 0.2339 | - / - / - |
| GSM505558 | Chondroid | 0.0644 | - / - / - |

  `example_patient.csv` is GSM505471 alone. In bulk data, Squamous and
  Chondroid mean an expression programme, not a histological diagnosis.

## example_10x/

- **Source:** GEO sample [GSM5354516](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5354516)
  (patient CID3941, ER+ primary breast cancer, 10x Chromium), series
  [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078).
- **Citation:** Wu SZ, et al. A single-cell and spatially resolved atlas of
  human breast cancers. *Nat Genet* 2021;53:1334-1347.
  doi:10.1038/s41588-021-00911-1 (open access, CC BY 4.0).
- **Terms:** as for GEO above. Cite Wu et al. 2021.
- **Subset:** raw UMI counts, unchanged. 120 of the 196 cancer epithelial cells,
  60 of the 266 T cells and all 55 B cells (sampled with `set.seed(20211)`).
  Genes: the `lamb_universe()` genes and the `Epithelial_share` genes
  present in the source (13,410 genes).
- **Clusters:** `analysis/clustering/gene_expression_graphclust/clusters.csv`
  holds the **authors' cell-type labels** (`celltype_major`), not a Cell
  Ranger clustering. Only the file layout follows Cell Ranger.
- **Result:** cancer epithelial cluster LHS. The T and B cell clusters are
  left out as immune by the cell-type gate, so their classes are not reported.

## xenium_sim/  (SIMULATED)

- **No tissue was measured.** Cells were Poisson-sampled from three public
  GSE20194 profiles in `example_bulk.csv`: GSM505471 (LHS), GSM505423
  (Basal-like) and GSM505546 (MA-Lo), 110 cells each, about 300
  transcripts per cell. A fourth group of 120 cells has the `Epithelial_share`
  genes and LAMB signature genes set close to zero and the immune gate genes
  (PTPRC, CD3E, ...) among its highest. It stands in for immune cells and is
  **not** a real immune profile.
- **Panel:** 392 genes (the LAMB class-signature, HER2-amplicon,
  `Epithelial_share` and immune gate genes present in `example_bulk.csv`, plus
  300 random universe genes), plus 5 negative-control features.
- **Layout:** each epithelial group occupies its own region of the field, and
  the non-epithelial cells fill the space between them. Clusters 1-4 are the
  simulated groups. `experiment.xenium` has `"pixel_size": 0.2125` and a
  note saying the data are simulated. `SIMULATED_DATA.txt` in the folder says
  the same.
- **Result:** clusters 1/2/3 are LHS / Basal-like / MA-Lo, and cluster 4 is
  left out as immune. HER2 amplicon is "high" on cluster 3 (its source profile
  is clinically HER2-positive) and "not high" on clusters 1 and 2. The slide-level call is "mixed" because no
  class holds half of the epithelial cells.
- **Script:** `data-raw/xenium_sim.R` (`set.seed(20194)`).
