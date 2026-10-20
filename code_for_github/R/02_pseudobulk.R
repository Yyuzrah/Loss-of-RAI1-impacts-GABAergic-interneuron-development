# One pseudobulk profile per animal. CTR is the reference, so a positive
# log2 fold change is higher in SMS. Drd2 and Adora2a stay in the model
# and are removed only from the saved table.
source("R/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

cortex <- load_cortex()
out <- out_dir("pseudobulk")

medians <- data.frame(
  replicate = names(tapply(cortex$nCount_Spatial, cortex$spatial_replicate, median)),
  median_umis = as.numeric(tapply(cortex$nCount_Spatial, cortex$spatial_replicate, median)),
  median_genes = as.numeric(tapply(cortex$nFeature_Spatial, cortex$spatial_replicate, median))
)
write.csv(medians, file.path(out, "per_replicate_medians.csv"), row.names = FALSE)

run_pseudobulk <- function(obj) {
  # AggregateExpression already sums the counts layer. Passing layer=
  # again is an error in Seurat 5.3.
  pseudo <- AggregateExpression(
    obj,
    group.by = "spatial_replicate",
    assays = "Spatial",
    return.seurat = FALSE
  )
  cts <- round(pseudo$Spatial)
  coldata <- data.frame(
    row.names = colnames(cts),
    condition = factor(
      ifelse(grepl("^CTR", colnames(cts)), "CTR", "SMS"),
      levels = c("CTR", "SMS")
    )
  )
  dds <- DESeqDataSetFromMatrix(cts, coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  dds <- DESeq(dds)
  coef_name <- "condition_SMS_vs_CTR"
  if (!coef_name %in% resultsNames(dds)) {
    stop("DESeq2 coefficient ", coef_name, " was not found.")
  }
  res_base <- results(dds, contrast = c("condition", "SMS", "CTR"))
  shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res_base)
  df <- as.data.frame(shrunk)
  df$gene <- rownames(df)
  list(dds = dds, table = df)
}

full <- run_pseudobulk(cortex)
res <- full$table[!full$table$gene %in% c("Drd2", "Adora2a"), ]
write.csv(res, file.path(out, "pseudobulk_deseq2_shrunken.csv"), row.names = FALSE)

ieg <- c("Rai1", "Fos", "Junb", "Fosl2", "Egr1", "Egr2", "Nr4a1", "Nr4a2", "Nr4a3", "Arc", "Dusp5")
res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 0 & res$padj < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange < 0 & res$padj < 0.05] <- "DOWN"
res$label <- ifelse(res$gene %in% ieg, res$gene, NA)

volcano <- ggplot(res, aes(log2FoldChange, -log10(padj), col = diffexpressed, label = label)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c(DOWN = "blue", NO = "grey", UP = "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = subset(res, !is.na(label)),
    max.overlaps = Inf, box.padding = 0.8, color = "black", size = 3
  ) +
  theme_minimal() +
  labs(
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p",
    color = NULL,
    title = "Isocortex pseudobulk, SMS versus CTR"
  ) +
  theme(legend.position = "none")
ggsave(file.path(out, "pseudobulk_volcano.png"), volcano, width = 8, height = 6, dpi = 140)

vsd <- vst(full$dds, blind = TRUE, nsub = min(1000, nrow(full$dds)))
pca <- plotPCA(vsd, intgroup = "condition") +
  geom_text_repel(aes(label = name), size = 3.5) +
  scale_color_manual(values = c(CTR = "grey40", SMS = "red")) +
  theme_minimal() +
  labs(title = "Pseudobulk PCA")
ggsave(file.path(out, "pseudobulk_pca.png"), pca, width = 6, height = 5, dpi = 140)

# Sensitivity: drop bins below 100 UMIs or 50 genes and repeat the same model.
pass <- cortex$nCount_Spatial >= 100 & cortex$nFeature_Spatial >= 50
removed <- !pass
qc <- data.frame(
  group = c("CTR", "SMS", "all"),
  bins_removed = c(
    sum(removed & cortex$condition == "CTR"),
    sum(removed & cortex$condition == "SMS"),
    sum(removed)
  ),
  bins_total = c(
    sum(cortex$condition == "CTR"),
    sum(cortex$condition == "SMS"),
    ncol(cortex)
  ),
  umis_removed = c(
    sum(cortex$nCount_Spatial[removed & cortex$condition == "CTR"]),
    sum(cortex$nCount_Spatial[removed & cortex$condition == "SMS"]),
    sum(cortex$nCount_Spatial[removed])
  ),
  umis_total = c(
    sum(cortex$nCount_Spatial[cortex$condition == "CTR"]),
    sum(cortex$nCount_Spatial[cortex$condition == "SMS"]),
    sum(cortex$nCount_Spatial)
  )
)
write.csv(qc, file.path(out, "qc_sensitivity_counts.csv"), row.names = FALSE)

filtered <- run_pseudobulk(subset(cortex, cells = colnames(cortex)[pass]))
filt <- filtered$table
filt <- filt[filt$gene %in% ieg, c("gene", "log2FoldChange", "lfcSE", "pvalue", "padj")]
write.csv(filt, file.path(out, "qc_sensitivity_ieg.csv"), row.names = FALSE)

# Whole-isocortex GO uses these pseudobulk calls (BH padj < 0.05), split by sign.
for (direction in c("UP", "DOWN")) {
  genes <- res$gene[res$diffexpressed == direction]
  if (length(genes) < 5) next
  ego <- enrichGO(
    genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.05
  )
  if (is.null(ego) || nrow(ego) == 0) next
  write.csv(as.data.frame(ego), file.path(out, paste0("go_", tolower(direction), ".csv")), row.names = FALSE)
  p <- dotplot(ego, showCategory = 15) + ggtitle(paste("Isocortex", direction))
  ggsave(file.path(out, paste0("go_", tolower(direction), ".png")), p, width = 7, height = 6, dpi = 140)
}

p <- DotPlot(cortex, features = ieg, group.by = "condition") +
  ggtitle("Selected genes, all isocortical bins")
ggsave(file.path(out, "selected_genes_dotplot.png"), p, width = 8, height = 4, dpi = 140)
message("Pseudobulk tables are in ", out)
