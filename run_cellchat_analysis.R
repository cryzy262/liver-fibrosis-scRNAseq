
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
