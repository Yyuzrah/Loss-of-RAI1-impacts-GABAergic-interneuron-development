# Visium HD isocortex, SMS versus CTR

R scripts for the 16 µm Visium HD comparison of Rai1 haploinsufficient (SMS) and control (CTR) isocortex. Space Ranger 4.0.1 wrote the inputs. Each capture area holds two animals of the same genotype.

The count matrices and the GSE124952 reference are not in this repository.

## Setup

R 4.4.1, Seurat 5. Other packages used here: arrow, DESeq2, apeglm, ggplot2, ggrepel, dplyr, stringr, clusterProfiler, org.Mm.eg.db, FNN, rstan, CellChat, slingshot, SingleCellExperiment, viridis.

In `R/00_config.R`, set `proj_dir` to the directory that contains:

- `Huang_Minza__Ctr_1_A1/square_016um/`
- `Huang_Minza__SMS_1_D1/square_016um/`
- `Task2/DataFromPaper/GSE124952_expression_matrix.csv.gz`
- `Task2/DataFromPaper/GSE124952_meta_data.csv.gz`

Each `square_016um` folder needs its `spatial/` directory. `SMS_vs_CTR_isocortex.ipynb` is the same code, one step per cell, with no saved output. Set the working directory to this folder before running it. From this directory in R:

```r
source("R/01_import_isocortex.R")
source("R/02_pseudobulk.R")
source("R/03_label_transfer.R")
source("R/04_celltype_de.R")
source("R/05_rai1_dotplot.R")
source("R/06_interneuron_subtypes.R")
source("R/07_cellchat.R")
source("R/08_slingshot.R")
source("R/09_cortical_depth.R")
source("R/10_distect_pixel.R")
source("R/11_distect_array.R")
```

Outputs go to `results/`. `01` stops if the crop is not 27,586 / 28,126 / 26,691 / 25,575 bins. `09` stops if the L1 counts are not 5,397 / 4,540 / 3,404 / 4,308.

## Scripts

`01` keeps the four hand-drawn boxes, joins the samples, and log-normalizes the joined counts. In `GetTissueCoordinates()` for these captures, x is the full-resolution row and y is the full-resolution column. Pseudobulk uses the raw counts, so it does not depend on that normalization.

`02` sums one profile per animal (`spatial_replicate`: one box, one animal) and fits DESeq2, `~ condition`, with apeglm shrinkage. The unshrunken p-values are kept. Bins below 100 UMIs or 50 genes are a second run of the same model. Whole-isocortex GO uses genes with BH adjusted p < 0.05, split by the sign of the shrunken log2 fold change.

`03` transfers labels from Bhattacherjee et al. 2019 (GSE124952). `04` is the bin-level Wilcoxon test within each predicted class (log2 fold change 0.1, min.pct 0.1). The class volcano colors genes with absolute log2 fold change above 0.25 and Bonferroni p below 0.05. Cell-class GO uses that Bonferroni column, then BH inside `enrichGO`. Genes called in four or more classes are written out separately. `05` is the Rai1 dot plot. A group with fewer than 20 detected bins is drawn at the color midpoint.

`06` takes bins with Gad2 > 0, clusters them (15 PCs, resolution 0.4), and names a cluster PV, SST, or VIP when that marker has the highest mean and the mean is above 0.1. Subtype differential expression uses log2 fold change 0.05 and min.pct 0.05. The volcano and the subtype GO use the unadjusted p-value. `07` fits CellChat separately in CTR and SMS on excitatory, PV, SST, and VIP bins (`type = "triMean"`, `min.cells = 3`). If that filter returns no interactions, the fit is repeated without it. `08` fits Slingshot on CTR and SMS together. Inhibitory bins use the subtype labels. Excitatory bins are reclustered at resolution 0.3.

`09` sets depth as distance to the pial contour, divided by the piece-wise 99.5th percentile. The cuts 0.12, 0.42, and 0.72 are display bands, as a fraction of cortical thickness. They were not fit from the marker curves. Upper cortex is depth at or below 0.42. Deep cortex is depth at or above 0.55.

`10` and `11` fit the DiSTect autologistic model (Zhao, Deng, and Zhang, Bioinformatics, 2025) to 200 variable genes. `set.seed(42)` draws 500 bins per piece for all cortex and 400 for upper and deep. The variational fit is mean-field, seed 128, and each model takes about 15 minutes. A gene is written to the selected table when its 95% credible interval excludes 0. `10` passes image pixels as the neighbor coordinates. Neighboring 16 µm bins are about 56 px apart, and the neighbor radius is 1 coordinate unit, so that term is off. `11` passes Space Ranger array row and column on the same draw. Full-rank variational inference did not adapt on these data. Pareto k from the fits that were run is about 4, so the selected genes are preliminary.

Label transfer, the clustering, and Slingshot were not seeded when those plots were first drawn. These scripts set seed 1 there, so a rerun will not match those plots cell for cell. The DiSTect subsample uses the seeds above.

## Not in this repository

The comparison with Huang et al. 2016 and Garay et al. 2020 has no script here. The Wilcoxon test of subtype abundance across depth was not saved as a script. Drafts that measured depth along a straight line, or that treated an oligodendrocyte-rich interior as corpus callosum, are not part of this analysis.
