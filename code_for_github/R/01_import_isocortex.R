# Four isocortical pieces. Counts are joined, then log-normalized.
source("R/00_config.R")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(arrow))

crop_piece <- function(obj, boxes, sample_name, positions) {
  coords <- GetTissueCoordinates(obj)
  if ("cell" %in% colnames(coords)) rownames(coords) <- coords$cell
  x_col <- grep("col|x", colnames(coords), value = TRUE, ignore.case = TRUE)[1]
  y_col <- grep("row|y", colnames(coords), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(x_col) || is.na(y_col) || x_col == y_col) {
    stop("Could not tell the coordinate columns apart.")
  }

  replicate <- rep(NA_character_, ncol(obj))
  names(replicate) <- colnames(obj)
  keep <- character()
  for (i in seq_len(nrow(boxes))) {
    ids <- rownames(coords[
      coords[[x_col]] >= boxes$xmin[i] & coords[[x_col]] <= boxes$xmax[i] &
        coords[[y_col]] >= boxes$ymin[i] & coords[[y_col]] <= boxes$ymax[i],
    ])
    if (length(ids)) {
      replicate[ids] <- paste0(sample_name, "_Rep", i)
      keep <- c(keep, ids)
    }
  }
  keep <- unique(keep)
  if (!length(keep)) stop("No bins in the ", sample_name, " boxes.")

  obj <- subset(obj, cells = keep)
  obj$spatial_replicate <- replicate[colnames(obj)]
  obj$condition <- sample_name
  obj$pxl_x <- coords[colnames(obj), x_col]
  obj$pxl_y <- coords[colnames(obj), y_col]

  hit <- match(colnames(obj), positions$barcode)
  if (anyNA(hit)) stop("Array coordinates missing for ", sum(is.na(hit)), " ", sample_name, " bins.")
  obj$array_row <- positions$array_row[hit]
  obj$array_col <- positions$array_col[hit]
  obj
}

read_positions <- function(sample_dir) {
  arrow::read_parquet(file.path(sample_dir, "spatial", "tissue_positions.parquet"))
}

ctr <- Load10X_Spatial(path_ctr, filename = "filtered_feature_bc_matrix.h5")
sms <- Load10X_Spatial(path_sms, filename = "filtered_feature_bc_matrix.h5")
ctr <- crop_piece(ctr, ctr_boxes, "CTR", read_positions(path_ctr))
sms <- crop_piece(sms, sms_boxes, "SMS", read_positions(path_sms))

cortex <- merge(ctr, y = sms, add.cell.ids = c("CTR", "SMS"))
rm(ctr, sms)
cortex <- JoinLayers(cortex)
DefaultAssay(cortex) <- "Spatial"
cortex <- NormalizeData(cortex, normalization.method = "LogNormalize", verbose = FALSE)
Idents(cortex) <- "condition"

counts <- table(cortex$spatial_replicate)
print(counts)
if (!identical(as.integer(counts[names(expected_bins)]), as.integer(expected_bins))) {
  stop("Crop bin counts do not match the analysis set.")
}
message("Bins: ", ncol(cortex), "  Genes: ", nrow(cortex))
save_cortex(cortex)
