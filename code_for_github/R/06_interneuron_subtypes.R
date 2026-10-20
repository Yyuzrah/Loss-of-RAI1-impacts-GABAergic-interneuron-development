# Gad2-detected bins, clustered at resolution 0.4 on 15 PCs.
# A cluster is PV, SST, or VIP when that marker has the highest mean and the mean exceeds 0.1.
source("R/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

cortex <- load_cortex()
DefaultAssay(cortex) <- "Spatial"
out <- out_dir("subtypes")

inhib <- subset(cortex, subset = Gad2 > 0)
inhib <- NormalizeData(inhib, verbose = FALSE)
inhib <- FindVariableFeatures(inhib, verbose = FALSE)
inhib <- ScaleData(inhib, verbose = FALSE)
inhib <- RunPCA(inhib, verbose = FALSE)
set.seed(1)
inhib <- RunUMAP(inhib, dims = 1:15, verbose = FALSE)
inhib <- FindNeighbors(inhib, dims = 1:15, verbose = FALSE)
inhib <- FindClusters(inhib, resolution = 0.4, verbose = FALSE)

avg <- AverageExpression(inhib, features = c("Pvalb", "Sst", "Vip"), assays = "Spatial")$Spatial
if (all(grepl("^g[0-9]+$", colnames(avg)))) {
  colnames(avg) <- sub("^g", "", colnames(avg))
}

cluster_map <- data.frame(cluster = colnames(avg), subtype = "Other", stringsAsFactors = FALSE)
for (i in seq_len(ncol(avg))) {
  top <- rownames(avg)[which.max(avg[, i])]
  if (max(avg[, i]) <= 0.1) next
  cluster_map$subtype[i] <- switch(top, Pvalb = "PV", Sst = "SST", Vip = "VIP", "Other")
}
write.csv(cluster_map, file.path(out, "cluster_to_subtype.csv"), row.names = FALSE)

inhib$subtype <- cluster_map$subtype[match(as.character(Idents(inhib)), cluster_map$cluster)]
if (all(is.na(inhib$subtype))) stop("Cluster ids did not match the marker table.")
Idents(inhib) <- "subtype"
inhib <- subset(inhib, idents = c("PV", "SST", "VIP"))
print(table(inhib$subtype, inhib$condition))

p <- DimPlot(inhib, reduction = "umap", group.by = "subtype", label = TRUE) +
  ggtitle("Inhibitory bins")
ggsave(file.path(out, "umap_subtypes.png"), p, width = 6, height = 5, dpi = 140)
p <- DimPlot(inhib, reduction = "umap", group.by = "subtype", split.by = "condition") +
  ggtitle("Inhibitory bins")
ggsave(file.path(out, "umap_subtypes_by_condition.png"), p, width = 10, height = 5, dpi = 140)

# FindMarkers Bonferroni uses every gene. Subtype GO and volcanoes use the unadjusted p.
tables <- list()
for (st in c("PV", "SST", "VIP")) {
  sub <- subset(inhib, idents = st)
  if (length(unique(sub$condition)) < 2) next
  markers <- FindMarkers(
    sub, ident.1 = "SMS", ident.2 = "CTR", group.by = "condition",
    logfc.threshold = 0.05, min.pct = 0.05, verbose = FALSE
  )
  if (nrow(markers) == 0) next
  markers$gene <- rownames(markers)
  markers$cell_type <- st
  markers$direction <- ifelse(markers$avg_log2FC > 0, "UP", "DOWN")
  tables[[st]] <- markers
  write.csv(markers, file.path(out, paste0("markers_", st, ".csv")), row.names = FALSE)

  for (direction in c("UP", "DOWN")) {
    genes <- markers$gene[markers$direction == direction & markers$p_val < 0.05]
    if (length(genes) < 5) next
    ego <- enrichGO(genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05)
    if (is.null(ego) || nrow(ego) == 0) next
    gp <- dotplot(ego, showCategory = 10) + ggtitle(paste(st, direction))
    ggsave(file.path(out, paste0("go_", tolower(direction), "_", st, ".png")), gp, width = 6, height = 6, dpi = 140)
  }
}
all_sub <- bind_rows(tables)
write.csv(all_sub, file.path(out, "markers_all_subtypes.csv"), row.names = FALSE)

# Rai1 is labeled on every subtype volcano, with the immediate-early genes.
down_core <- c("Fos", "Junb", "Fosl2", "Egr1", "Egr2", "Nr4a1", "Nr4a2", "Nr4a3", "Arc", "Dusp5", "Rai1")
up_list <- list(
  VIP = c("Scg5", "Rbfox1", "Plch2", "Btbd9", "Camk2n1"),
  SST = c("Camk2n1", "Glo1", "mt-Co2", "Nrn1", "Scn1b"),
  PV = c("Glo1", "mt-Co2", "Scn1b", "Camk2n1", "Nrn1")
)
for (st in names(tables)) {
  df <- tables[[st]]
  df$color <- "NS"
  df$color[df$avg_log2FC > 0.05 & df$p_val < 0.05] <- "UP"
  df$color[df$avg_log2FC < -0.05 & df$p_val < 0.05] <- "DOWN"
  label_genes <- intersect(df$gene, c(down_core, up_list[[st]]))
  df$label <- ifelse(df$gene %in% label_genes, df$gene, NA)
  p <- ggplot(df, aes(avg_log2FC, -log10(p_val), color = color, label = label)) +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_color_manual(values = c(DOWN = "blue", NS = "grey80", UP = "red")) +
    geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(max.overlaps = 50, force = 2, size = 3, min.segment.length = 0) +
    theme_minimal() +
    labs(title = paste(st, "SMS versus CTR"), x = "log2 fold change", y = "-log10 p") +
    theme(legend.position = "none")
  ggsave(file.path(out, paste0("volcano_", st, ".png")), p, width = 7, height = 6, dpi = 140)
}

cortex$subtype <- NA_character_
cortex$subtype[match(colnames(inhib), colnames(cortex))] <- as.character(inhib$subtype)
save_cortex(cortex)
saveRDS(inhib, file.path(results_dir, "inhib.rds"))
message("Subtype bins: ", ncol(inhib))
