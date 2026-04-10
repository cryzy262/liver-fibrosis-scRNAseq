# ==================== CellChat分析脚本（单独运行）====================
# 保存为: run_cellchat_analysis.R
# 运行: source("run_cellchat_analysis.R")

cat("========== CellChat Analysis ==========\n")

# 设置工作目录（根据您的实际路径修改）
work.dir <- "D:/文章/生信分析/小鼠肝脏单细胞RNA测序分析/20260325"
setwd(work.dir)

# 加载RData
rdata_file <- file.path(work.dir, "results", "2_analysis", "high_quality_cells_v4.RData")
if(!file.exists(rdata_file)) {
  stop("找不到文件: ", rdata_file)
}
load(rdata_file)
cat("✓ 加载数据完成:", ncol(high_quality_cells), "个细胞\n")

# 安装/加载CellChat
if(!requireNamespace("CellChat", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("sqjin/CellChat")
}
library(CellChat)
library(Seurat)
library(dplyr)

# 创建目录
cellchat.dir <- file.path(work.dir, "results", "2_analysis", "cellchat")
dir.create(cellchat.dir, showWarnings = FALSE, recursive = TRUE)

# 准备数据
cat("\n提取HSC和肝细胞...\n")
hsc_hep_cells <- subset(high_quality_cells, 
                        subset = celltype_reannotated %in% c("HSC_Act", "HSC_Qui", "Hepatocyte") |
                          celltype %in% c("HSC", "Hepatocyte"))

cat("细胞组成:\n")
print(table(hsc_hep_cells$celltype_reannotated, hsc_hep_cells$group))

# 创建CellChat对象的函数
create_cellchat_object <- function(seurat_obj, group_name) {
  cat(sprintf("\n处理 %s 组...\n", group_name))
  
  # 提取数据
  data.input <- GetAssayData(seurat_obj, layer = "data")
  meta <- data.frame(
    labels = seurat_obj$celltype_reannotated,
    row.names = colnames(seurat_obj)
  )
  
 
  
  # 设置数据库
  cellchat@DB <- CellChatDB.mouse
  
  # 分析流程
  cat("  预处理...\n")
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cat("  计算通讯概率...\n")
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # 保存
  saveRDS(cellchat, file.path(cellchat.dir, sprintf("cellchat_%s.rds", group_name)))
  cat(sprintf("  ✓ 保存: cellchat_%s.rds\n", group_name))
  
  return(cellchat)
}

# 分别处理两组
hsc_hep_ccl4 <- subset(hsc_hep_cells, subset = group == "CCl4")
hsc_hep_ctrl <- subset(hsc_hep_cells, subset = group == "Control")

cellchat_ccl4 <- create_cellchat_object(hsc_hep_ccl4, "CCl4")
cellchat_ctrl <- create_cellchat_object(hsc_hep_ctrl, "Control")

cat("\n========== CellChat分析完成 ==========\n")
cat("输出文件:\n")
cat(sprintf("  - %s/cellchat_CCl4.rds\n", cellchat.dir))
cat(sprintf("  - %s/cellchat_Control.rds\n", cellchat.dir))
cat("\n现在可以运行主分析脚本生成Figure 3\n")
