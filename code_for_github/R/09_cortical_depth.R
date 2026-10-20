# Depth is distance to the pial contour, divided by the piece-wise 99.5th percentile.
# Faces of the crop rectangle are not treated as pia.
source("R/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(FNN)
})

dilate4 <- function(z) {
  nr <- nrow(z)
  nc <- ncol(z)
  z |
    rbind(z[-1, , drop = FALSE], FALSE) |
    rbind(FALSE, z[-nr, , drop = FALSE]) |
    cbind(z[, -1, drop = FALSE], FALSE) |
    cbind(FALSE, z[, -nc, drop = FALSE])
}

flood_exterior <- function(tissue) {
  nr <- nrow(tissue)
  nc <- ncol(tissue)
  ext <- matrix(FALSE, nr, nc)
  ext[1, ] <- !tissue[1, ]
  ext[nr, ] <- !tissue[nr, ]
  ext[, 1] <- !tissue[, 1]
  ext[, nc] <- !tissue[, nc]
  for (i in seq_len(4000)) {
    nxt <- dilate4(ext) & !tissue
    if (identical(nxt, ext)) break
    ext <- nxt
  }
  ext
}

cortex <- load_cortex()
cortex$depth_from_pia <- NA_real_
cortex$layer_band <- NA_character_
cortex$layer_class <- NA_character_
pia_list <- list()
bin_list <- list()

for (rep in sort(unique(cortex$spatial_replicate))) {
  cells <- colnames(cortex)[cortex$spatial_replicate == rep]
  xy <- data.frame(x = cortex$pxl_x[cells], y = cortex$pxl_y[cells])
  step <- median(FNN::get.knnx(as.matrix(xy), as.matrix(xy), k = 2)$nn.dist[, 2])
  if (!is.finite(step) || step <= 0) step <- 56.25
  x0 <- min(xy$x)
  y0 <- min(xy$y)
  jj <- as.integer(round((xy$x - x0) / step)) + 3L
  ii <- as.integer(round((xy$y - y0) / step)) + 3L
  tissue <- matrix(FALSE, max(ii) + 2L, max(jj) + 2L)
  tissue[cbind(ii, jj)] <- TRUE
  outer_b <- tissue & dilate4(flood_exterior(tissue))
  ob <- which(outer_b, arr.ind = TRUE)
  bound <- data.frame(
    x = x0 + (ob[, 2] - 3) * step,
    y = y0 + (ob[, 1] - 3) * step
  )
  rx <- range(xy$x)
  ry <- range(xy$y)
  pad <- 2.2 * step
  on_crop <- bound$x <= rx[1] + pad | bound$x >= rx[2] - pad |
    bound$y <= ry[1] + pad | bound$y >= ry[2] - pad
  pia <- bound[!on_crop, , drop = FALSE]
  if (nrow(pia) < 2) stop("Pial contour is empty for ", rep)
  dist_pia <- FNN::get.knnx(as.matrix(pia), as.matrix(xy), k = 1)$nn.dist[, 1]
  depth <- pmin(pmax(dist_pia / (as.numeric(quantile(dist_pia, 0.995)) + 1e-6), 0), 1)
  band <- ifelse(depth <= cut_L1, "L1",
          ifelse(depth <= cut_L23, "L2/3",
          ifelse(depth <= cut_L5, "L5", "L6")))
  cls <- ifelse(depth <= upper_max, "upper", ifelse(depth >= deep_min, "deep", NA_character_))
  cortex$depth_from_pia[cells] <- depth
  cortex$layer_band[cells] <- band
  cortex$layer_class[cells] <- cls
  pia_list[[rep]] <- data.frame(pia, replicate = rep)
  bin_list[[rep]] <- data.frame(x = xy$x, y = xy$y, depth = depth, band = band, replicate = rep)
  message(rep, "  L1 ", sum(band == "L1"), "  L2/3 ", sum(band == "L2/3"),
          "  L5 ", sum(band == "L5"), "  L6 ", sum(band == "L6"))
}

# These L1 counts are from the executed run. A mismatch means the contour moved.
expected_l1 <- c(CTR_Rep1 = 5397L, CTR_Rep2 = 4540L, SMS_Rep1 = 3404L, SMS_Rep2 = 4308L)
got_l1 <- sapply(names(expected_l1), function(rep) {
  sum(cortex$layer_band == "L1" & cortex$spatial_replicate == rep, na.rm = TRUE)
})
if (!identical(as.integer(got_l1), as.integer(expected_l1))) {
  stop("L1 counts changed: ", paste(got_l1, collapse = ", "))
}

out <- out_dir("depth")
pia_df <- do.call(rbind, pia_list)
bin_df <- do.call(rbind, bin_list)
write.csv(
  as.data.frame(table(cortex$layer_band, cortex$spatial_replicate)),
  file.path(out, "layer_band_counts.csv"),
  row.names = FALSE
)

theme_sp <- theme_minimal() + theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
p_pia <- ggplot() +
  geom_point(data = bin_df, aes(x, y), color = "grey82", size = 0.12) +
  geom_point(data = pia_df, aes(x, y), color = "#00BFC4", size = 0.4) +
  scale_y_reverse() +
  facet_wrap(~replicate, scales = "free") +
  theme_sp +
  labs(title = "Pial contour", x = NULL, y = NULL)
ggsave(file.path(out, "pial_contour.png"), p_pia, width = 12, height = 9, dpi = 140)

p_band <- ggplot(bin_df, aes(x, y, color = band)) +
  geom_point(size = 0.15) +
  scale_color_manual(values = c("L1" = "#F4A261", "L2/3" = "#2C7BB6", "L5" = "#D7191C", "L6" = "#6A3D9A")) +
  scale_y_reverse() +
  facet_wrap(~replicate, scales = "free") +
  theme_sp +
  labs(title = "Display bands", x = NULL, y = NULL, color = NULL)
ggsave(file.path(out, "layer_bands.png"), p_band, width = 12, height = 9, dpi = 140)

markers <- c("Reln", "Lamp5", "Cux2", "Rorb", "Ndnf", "Sncg", "Pvalb", "Nefm", "Bcl11b", "Fezf2", "Foxp2", "Tbr1")
markers <- intersect(markers, rownames(cortex))
expr <- as.matrix(GetAssayData(cortex, layer = "data")[markers, , drop = FALSE])
long <- do.call(rbind, lapply(markers, function(g) {
  data.frame(
    replicate = cortex$spatial_replicate,
    depth = cortex$depth_from_pia,
    gene = g,
    expr = expr[g, ],
    stringsAsFactors = FALSE
  )
}))
p_san <- ggplot(long, aes(depth, expr, color = gene)) +
  geom_smooth(method = "loess", span = 0.35, se = FALSE, linewidth = 0.7, formula = y ~ x) +
  facet_wrap(~replicate, ncol = 2) +
  theme_minimal() +
  labs(x = "Depth from pia", y = "Log-normalized expression", title = "Laminar markers")
ggsave(file.path(out, "markers_vs_depth.png"), p_san, width = 12, height = 8, dpi = 140)

save_cortex(cortex)
message("Depth written. L1 counts match the executed run.")
