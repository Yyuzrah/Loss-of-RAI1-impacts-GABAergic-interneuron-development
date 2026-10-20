# Edit proj_dir, then source the numbered scripts from this directory.
#   source("R/01_import_isocortex.R")

# Analysis-machine project directory. The escapes keep the path valid
# regardless of how R reads this file.
proj_dir <- "F:/H.M.\u6253\u706b\u673a"

repo_dir <- normalizePath(getwd(), winslash = "/")
results_dir <- file.path(repo_dir, "results")

path_ctr <- file.path(proj_dir, "Huang_Minza__Ctr_1_A1", "square_016um")
path_sms <- file.path(proj_dir, "Huang_Minza__SMS_1_D1", "square_016um")
path_ref_counts <- file.path(proj_dir, "Task2", "DataFromPaper", "GSE124952_expression_matrix.csv.gz")
path_ref_meta <- file.path(proj_dir, "Task2", "DataFromPaper", "GSE124952_meta_data.csv.gz")

ctr_boxes <- data.frame(
  xmin = c(4000, 15500), xmax = c(14500, 27500),
  ymin = c(7000, 6000), ymax = c(17500, 16500)
)
sms_boxes <- data.frame(
  xmin = c(15000, 2000), xmax = c(26000, 14500),
  ymin = c(10000, 12500), ymax = c(19500, 21500)
)

# Drawn on the H&E before differential expression. In GetTissueCoordinates()
# for these captures, x is full-resolution row and y is full-resolution column.
expected_bins <- c(CTR_Rep1 = 27586L, CTR_Rep2 = 28126L, SMS_Rep1 = 26691L, SMS_Rep2 = 25575L)

cut_L1 <- 0.12
cut_L23 <- 0.42
cut_L5 <- 0.72
upper_max <- 0.42
deep_min <- 0.55

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
options(warn = 1)
options(future.globals.maxSize = 8000 * 1024^2)

out_dir <- function(...) {
  d <- file.path(results_dir, ...)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

load_cortex <- function() {
  p <- file.path(results_dir, "cortex.rds")
  if (!file.exists(p)) stop("cortex.rds is missing. Run R/01_import_isocortex.R from the repository root.")
  readRDS(p)
}

save_cortex <- function(obj) {
  saveRDS(obj, file.path(results_dir, "cortex.rds"))
  invisible(obj)
}
