# Pixel-coordinate DiSTect. This is the fit behind the current credible-interval figure.
# Each model is mean-field ADVI and takes roughly 15 minutes.
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
cortex$coord_x <- cortex$pxl_x
cortex$coord_y <- cortex$pxl_y
fit_three(cortex, out_dir("distect"))
