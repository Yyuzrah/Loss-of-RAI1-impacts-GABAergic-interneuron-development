# Rai1 by predicted class. Groups with fewer than 20 detected bins are drawn at the midpoint.
source("R/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

cortex <- load_cortex()
if (!"predicted_celltype" %in% colnames(cortex[[]])) {
  stop("predicted_celltype is missing. Run R/03_label_transfer.R.")
}

dat <- FetchData(cortex, vars = c("Rai1", "predicted_celltype", "condition"), layer = "data")
stats <- dat %>%
  group_by(predicted_celltype, condition) %>%
  summarise(
    bins = n(),
    bins_detected = sum(Rai1 > 0),
    percent_detected = 100 * sum(Rai1 > 0) / n(),
    mean_count = mean(expm1(Rai1)),
    .groups = "drop"
  )
stats$scaled <- as.numeric(scale(log1p(stats$mean_count)))
stats$color_value <- stats$scaled
# Fewer than 20 detected bins is drawn as the midpoint, not as a z-score.
stats$color_value[stats$bins_detected < 20] <- 0

out <- out_dir("celltypes")
write.csv(stats, file.path(out, "rai1_by_class.csv"), row.names = FALSE)

p <- ggplot(stats, aes(condition, predicted_celltype)) +
  geom_point(aes(size = percent_detected, color = color_value)) +
  scale_size_continuous(range = c(3, 9), limits = c(0, 100), name = "% bins detected") +
  scale_color_gradient2(
    low = "blue", mid = "grey95", high = "red", midpoint = 0,
    name = "Relative expression"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, title = "Rai1")
ggsave(file.path(out, "rai1_dotplot.png"), p, width = 6, height = 7, dpi = 140)
