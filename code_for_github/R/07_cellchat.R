# Ligand-receptor inference for excitatory, PV, SST, and VIP bins.
# CTR and SMS are fit separately. min.cells = 3 because several subtypes are small.
source("R/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
})

cortex <- load_cortex()
if (!"subtype" %in% colnames(cortex[[]])) stop("subtype is missing. Run R/06_interneuron_subtypes.R.")

cortex$fine_label <- as.character(cortex$predicted_celltype)
has_subtype <- !is.na(cortex$subtype)
cortex$fine_label[has_subtype] <- cortex$subtype[has_subtype]

keep <- cortex$fine_label %in% c("Excitatory", "PV", "SST", "VIP")
mini <- subset(cortex, cells = colnames(cortex)[keep])
mini$fine_label <- factor(mini$fine_label, levels = c("Excitatory", "PV", "SST", "VIP"))
Idents(mini) <- "fine_label"

fit_cellchat <- function(obj) {
  obj$fine_label <- droplevels(factor(
    obj$fine_label, levels = c("Excitatory", "PV", "SST", "VIP")
  ))
  obj <- JoinLayers(obj)
  expr <- GetAssayData(obj, assay = "Spatial", layer = "data")
  cc <- createCellChat(object = expr, meta = obj@meta.data, group.by = "fine_label")
  cc@DB <- CellChatDB.mouse
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 3)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  if (sum(cc@net$count) == 0) {
    message("No interactions at min.cells = 3. Refitting without that filter.")
    cc <- computeCommunProb(cc, type = "triMean")
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
  }
  cc
}

ctr <- fit_cellchat(subset(mini, subset = condition == "CTR"))
sms <- fit_cellchat(subset(mini, subset = condition == "SMS"))
saveRDS(list(CTR = ctr, SMS = sms), file.path(results_dir, "cellchat.rds"))

out <- out_dir("cellchat")
colors <- scPalette(4)
names(colors) <- c("Excitatory", "PV", "SST", "VIP")
for (nm in c("CTR", "SMS")) {
  cc <- if (nm == "CTR") ctr else sms
  message(nm, " interactions: ", sum(cc@net$count))
  if (sum(cc@net$count) == 0) next
  png(file.path(out, paste0("circle_", nm, ".png")), width = 2000, height = 2000, res = 300)
  netVisual_circle(
    cc@net$count,
    vertex.weight = as.numeric(table(cc@idents)),
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = colors[levels(cc@idents)],
    title.name = nm
  )
  dev.off()
}
