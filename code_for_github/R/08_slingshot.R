# Shared pseudotime for CTR and SMS. Inhibitory bins use the subtype labels.
# Excitatory bins are reclustered at resolution 0.3 to give Slingshot anchor groups.
source("R/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(ggplot2)
  library(viridis)
})

draw_trajectory <- function(df, curves, title, path) {
  p <- ggplot() +
    geom_point(data = df, aes(UMAP_1, UMAP_2, color = Pseudotime), size = 0.8, alpha = 0.7) +
    scale_color_viridis_c(option = "magma", na.value = "grey90") +
    geom_path(data = curves, aes(UMAP_1, UMAP_2, group = Lineage),
              color = "black", linewidth = 0.8,
              arrow = arrow(length = unit(0.25, "cm"))) +
    facet_wrap(~ Condition) +
    theme_minimal() +
    labs(title = title, color = "Pseudotime")
  ggsave(path, p, width = 10, height = 5, dpi = 140)

  dens <- ggplot(df, aes(Pseudotime, fill = Condition)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c(CTR = "grey60", SMS = "red")) +
    theme_minimal() +
    labs(title = title, x = "Pseudotime", y = "Density")
  ggsave(sub("trajectory", "density", path), dens, width = 7, height = 4, dpi = 140)
}

curve_frame <- function(sce) {
  curves <- slingCurves(sce)
  bind <- lapply(seq_along(curves), function(i) {
    xy <- curves[[i]]$s[curves[[i]]$ord, , drop = FALSE]
    data.frame(UMAP_1 = xy[, 1], UMAP_2 = xy[, 2], Lineage = paste0("Lineage", i))
  })
  do.call(rbind, bind)
}

out <- out_dir("trajectory")

inhib_path <- file.path(results_dir, "inhib.rds")
if (!file.exists(inhib_path)) stop("inhib.rds is missing. Run R/06_interneuron_subtypes.R.")
inhib <- readRDS(inhib_path)
set.seed(1)
sce <- slingshot(as.SingleCellExperiment(inhib), clusterLabels = "subtype", reducedDim = "UMAP")
pt <- slingPseudotime(sce)
umap <- reducedDims(sce)$UMAP
df <- data.frame(
  UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
  Pseudotime = pt[, 1],
  Condition = inhib$condition,
  Subtype = inhib$subtype
)
draw_trajectory(df, curve_frame(sce), "Inhibitory subtypes", file.path(out, "inhibitory_trajectory.png"))

cortex <- load_cortex()
Idents(cortex) <- "predicted_celltype"
excit <- subset(cortex, idents = "Excitatory")
DefaultAssay(excit) <- "Spatial"
excit <- FindVariableFeatures(excit, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
excit <- ScaleData(excit, verbose = FALSE)
excit <- RunPCA(excit, verbose = FALSE)
set.seed(1)
excit <- RunUMAP(excit, dims = 1:30, verbose = FALSE)
excit <- FindNeighbors(excit, dims = 1:30, verbose = FALSE)
excit <- FindClusters(excit, resolution = 0.3, verbose = FALSE)

set.seed(1)
sce <- slingshot(as.SingleCellExperiment(excit), clusterLabels = "seurat_clusters", reducedDim = "UMAP")
pt <- slingPseudotime(sce)
umap <- reducedDims(sce)$UMAP
df <- data.frame(
  UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
  Pseudotime = pt[, 1],
  Condition = excit$condition,
  Cluster = excit$seurat_clusters
)
draw_trajectory(df, curve_frame(sce), "Excitatory bins", file.path(out, "excitatory_trajectory.png"))
message("Trajectory plots are in ", out)
