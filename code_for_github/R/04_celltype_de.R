# Wilcoxon tests within each predicted class. p_val_adj is Seurat's Bonferroni.
# Drd2 and Adora2a are removed after the test.
source("R/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

cortex <- load_cortex()
if (!"predicted_celltype" %in% colnames(cortex[[]])) {
  stop("predicted_celltype is missing. Run R/03_label_transfer.R.")
}
DefaultAssay(cortex) <- "Spatial"
Idents(cortex) <- "predicted_celltype"
out <- out_dir("celltypes")
blacklist <- c("Drd2", "Adora2a")

# Initial whole-isocortex screen. The reported isocortex test is the pseudobulk.
screen <- FindMarkers(
  cortex, ident.1 = "SMS", ident.2 = "CTR", group.by = "condition",
  logfc.threshold = 0.1
)
screen$gene <- rownames(screen)
screen <- screen[!screen$gene %in% blacklist, ]
screen <- screen[order(screen$avg_log2FC, decreasing = TRUE), ]
write.csv(screen, file.path(out, "whole_isocortex_wilcox.csv"), row.names = FALSE)

tables <- list()
for (ctype in unique(cortex$predicted_celltype)) {
  sub <- subset(cortex, idents = ctype)
  if (length(unique(sub$condition)) < 2) next
  markers <- tryCatch(
    FindMarkers(
      sub, ident.1 = "SMS", ident.2 = "CTR", group.by = "condition",
      logfc.threshold = 0.1, min.pct = 0.1
    ),
    error = function(e) NULL
  )
  if (is.null(markers) || nrow(markers) == 0) next
  markers$gene <- rownames(markers)
  markers$cell_type <- ctype
  markers <- markers[!markers$gene %in% blacklist, ]
  tables[[ctype]] <- markers
  write.csv(markers, file.path(out, paste0("markers_", gsub(" ", "_", ctype), ".csv")), row.names = FALSE)

  # |log2FC| > 0.25 and Bonferroni p < 0.05. The test itself used a threshold of 0.1.
  markers$call <- "NS"
  markers$call[markers$avg_log2FC > 0.25 & markers$p_val_adj < 0.05] <- "UP"
  markers$call[markers$avg_log2FC < -0.25 & markers$p_val_adj < 0.05] <- "DOWN"
  p <- ggplot(markers, aes(avg_log2FC, -log10(p_val_adj), color = call)) +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_color_manual(values = c(DOWN = "blue", NS = "grey80", UP = "red")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    labs(
      title = paste(ctype, "SMS versus CTR"),
      x = "log2 fold change",
      y = "-log10 Bonferroni p"
    ) +
    theme(legend.position = "none")
  ggsave(
    file.path(out, paste0("volcano_", gsub(" ", "_", ctype), ".png")),
    p, width = 7, height = 6, dpi = 140
  )

  for (direction in c("up", "down")) {
    if (direction == "up") {
      genes <- markers$gene[markers$avg_log2FC > 0 & markers$p_val_adj < 0.05]
    } else {
      genes <- markers$gene[markers$avg_log2FC < 0 & markers$p_val_adj < 0.05]
    }
    if (length(genes) <= 5) next
    ego <- enrichGO(genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05)
    if (is.null(ego) || nrow(ego) == 0) next
    n_show <- min(15, nrow(ego))
    p <- dotplot(ego, showCategory = n_show) +
      ggtitle(paste(ctype, direction)) +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45))
    ggsave(
      file.path(out, paste0("go_", direction, "_", gsub(" ", "_", ctype), ".png")),
      p, width = 7, height = 4 + 0.25 * n_show, dpi = 140
    )
  }
}

all_markers <- bind_rows(tables)
write.csv(all_markers, file.path(out, "markers_all_celltypes.csv"), row.names = FALSE)

# Bonferroni p < 0.05 in four or more classes.
sig <- all_markers[all_markers$p_val_adj < 0.05, c("gene", "cell_type", "avg_log2FC")]
n_classes <- tapply(sig$cell_type, sig$gene, function(x) length(unique(x)))
shared_genes <- names(n_classes[n_classes >= 4])
if (length(shared_genes)) {
  shared <- sig[sig$gene %in% shared_genes, ]
  write.csv(shared, file.path(out, "genes_significant_in_four_or_more_classes.csv"), row.names = FALSE)
  p <- ggplot(shared, aes(cell_type, gene, fill = avg_log2FC)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = "log2 FC", title = "Genes called in four or more classes")
  ggsave(file.path(out, "genes_shared_across_classes.png"), p, width = 8, height = 0.28 * length(shared_genes) + 2, dpi = 140)
}

message("Cell-type tables: ", nrow(all_markers), " rows")
