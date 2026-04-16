# ==================== CellChat Analysis Script (Run Separately) ====================
# Save as: run_cellchat_analysis.R
# Run: source("run_cellchat_analysis.R")

cat("========== CellChat Analysis ==========\n")

# Set working directory (modify according to your actual path)
work.dir <- "D:/文章/生信分析/小鼠肝脏单细胞RNA测序分析/20260325"
setwd(work.dir)

# Load RData
rdata_file <- file.path(work.dir, "results", "2_analysis", "high_quality_cells_v4.RData")
if(!file.exists(rdata_file)) {
  stop("File not found: ", rdata_file)
}
load(rdata_file)
cat("✓ Data loading completed:", ncol(high_quality_cells), "cells\n")

# Install/Load CellChat
if(!requireNamespace("CellChat", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("sqjin/CellChat")
}
library(CellChat)
library(Seurat)
library(dplyr)

# Create directory
cellchat.dir <- file.path(work.dir, "results", "2_analysis", "cellchat")
dir.create(cellchat.dir, showWarnings = FALSE, recursive = TRUE)

# Prepare data
cat("\nExtracting HSC and Hepatocytes...\n")
hsc_hep_cells <- subset(high_quality_cells, 
                        subset = celltype_reannotated %in% c("HSC_Act", "HSC_Qui", "Hepatocyte") |
                          celltype %in% c("HSC", "Hepatocyte"))

cat("Cell composition:\n")
print(table(hsc_hep_cells$celltype_reannotated, hsc_hep_cells$group))

# Function to create CellChat object
create_cellchat_object <- function(seurat_obj, group_name) {
  cat(sprintf("\nProcessing %s group...\n", group_name))
  
  # Extract data
  data.input <- GetAssayData(seurat_obj, layer = "data")
  meta <- data.frame(
    labels = seurat_obj$celltype_reannotated,
    row.names = colnames(seurat_obj)
  )
  
  # Create object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "labels")
  
  # Set database
  cellchat@DB <- CellChatDB.mouse
  
  # Analysis workflow
  cat("  Preprocessing...\n")
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cat("  Computing communication probabilities...\n")
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Save
  saveRDS(cellchat, file.path(cellchat.dir, sprintf("cellchat_%s.rds", group_name)))
  cat(sprintf("  ✓ Saved: cellchat_%s.rds\n", group_name))
  
  return(cellchat)
}

# Process both groups
hsc_hep_ccl4 <- subset(hsc_hep_cells, subset = group == "CCl4")
hsc_hep_ctrl <- subset(hsc_hep_cells, subset = group == "Control")

cellchat_ccl4 <- create_cellchat_object(hsc_hep_ccl4, "CCl4")
cellchat_ctrl <- create_cellchat_object(hsc_hep_ctrl, "Control")

cat("\n========== CellChat Analysis Completed ==========\n")
cat("Output files:\n")
cat(sprintf("  - %s/cellchat_CCl4.rds\n", cellchat.dir))
cat(sprintf("  - %s/cellchat_Control.rds\n", cellchat.dir))
cat("\nNow you can run the main analysis script to generate Figure 3\n")
