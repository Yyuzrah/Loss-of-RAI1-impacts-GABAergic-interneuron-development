# Same bins and the same seeds as R/10_distect_pixel.R.
# Neighbor distance is Space Ranger array row and column instead of image pixels.
source("R/00_config.R")
source("R/distect_model.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(rstan)
  library(ggplot2)
})
rstan_options(auto_write = TRUE)

cortex <- load_cortex()
if (!"depth_from_pia" %in% colnames(cortex[[]])) {
  stop("Depth is missing. Run R/09_cortical_depth.R.")
}
if (anyNA(cortex$array_row) || anyNA(cortex$array_col)) {
  stop("Array coordinates are missing.")
}
cortex$coord_x <- cortex$array_row
cortex$coord_y <- cortex$array_col
fit_three(cortex, out_dir("distect_array"))
