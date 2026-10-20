# Label transfer from Bhattacherjee et al. 2019, GSE124952.
# The score is label confidence, not a cell count.
source("R/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

if (!file.exists(path_ref_counts) || !file.exists(path_ref_meta)) {
  stop("GSE124952 counts or metadata not found under Task2/DataFromPaper.")
}

ref_counts <- read.csv(path_ref_counts, row.names = 1, check.names = FALSE)
ref_meta <- read.csv(path_ref_meta, row.names = 1)
ref <- CreateSeuratObject(counts = ref_counts, meta.data = ref_meta)
if (!"CellType" %in% colnames(ref[[]])) stop("Reference metadata has no CellType column.")
Idents(ref) <- ref$CellType
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ref <- ScaleData(ref, verbose = FALSE)
ref <- RunPCA(ref, verbose = FALSE)

cortex <- load_cortex()
DefaultAssay(cortex) <- "Spatial"
set.seed(1)
cortex <- FindVariableFeatures(cortex, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
cortex <- ScaleData(cortex, verbose = FALSE)
cortex <- RunPCA(cortex, verbose = FALSE)

anchors <- FindTransferAnchors(
  reference = ref, query = cortex,
  normalization.method = "LogNormalize", dims = 1:30, verbose = FALSE
)
pred <- TransferData(
  anchorset = anchors, refdata = ref$CellType,
  prediction.assay = TRUE, weight.reduction = cortex[["pca"]],
  dims = 1:30, verbose = FALSE
)
cortex[["prediction"]] <- pred
probs <- GetAssayData(cortex, assay = "prediction")
# TransferData appends a feature named max. It is not a class.
if ("max" %in% rownames(probs)) {
  probs <- probs[rownames(probs) != "max", , drop = FALSE]
}
cortex$predicted_celltype <- rownames(probs)[apply(probs, 2, which.max)]
cortex$prediction_score <- apply(probs, 2, max)
Idents(cortex) <- "predicted_celltype"

out <- out_dir("celltypes")
scores <- data.frame(
  bins = ncol(cortex),
  score_at_least_0.5 = sum(cortex$prediction_score >= 0.5)
)
write.csv(scores, file.path(out, "prediction_score_summary.csv"), row.names = FALSE)
write.csv(
  as.data.frame(table(cortex$predicted_celltype, cortex$spatial_replicate)),
  file.path(out, "celltype_counts_by_replicate.csv"),
  row.names = FALSE
)

counts <- as.data.frame(table(cortex$predicted_celltype, cortex$spatial_replicate))
names(counts) <- c("celltype", "replicate", "n")
counts$condition <- sub("_Rep[12]$", "", counts$replicate)
p <- ggplot(counts, aes(celltype, n, fill = condition, group = replicate)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(CTR = "grey80", SMS = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Bins", fill = NULL, title = "Predicted class by animal")
ggsave(file.path(out, "celltype_counts_by_replicate.png"), p, width = 9, height = 5, dpi = 140)

save_cortex(cortex)
message(sprintf(
  "Score >= 0.5: %d / %d",
  scores$score_at_least_0.5, scores$bins
))
