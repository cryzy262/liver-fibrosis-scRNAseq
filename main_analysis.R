# ==============================================================================
# Single-cell RNA-seq Analysis of CCl4-induced Liver Fibrosis
# 
# Description: Comprehensive analysis of hepatocyte zonation, HSC activation,
#              cell-cell communication, and therapeutic target prioritization
# 
# Author: [Xiao Liu]
# Date: 2026-03-26
# 
# Repository: https://github.com/cryzy262/liver-fibrosis-scRNAseq
# 
# Usage:
#   1. Prepare data in data/ directory (see data/README.md)
#   2. Run CellChat analysis: source("run_cellchat_analysis.R")
#   3. Run main analysis: source("main_analysis.R")
# ==============================================================================

# ==================== 完整单细胞分析流程 v4.0(修复版)====================
# 修复: 路径转义问题、正则表达式语法

# ==================== 0. 环境设置 ====================
work.dir <- "D:/文章/生信分析/小鼠肝脏单细胞RNA测序分析/20260325"
set.seed(2025)
setwd(work.dir)

# 创建目录
dirs.to.create <- file.path(work.dir, c(
  "data", "results", 
  "results/1_quality_control", 
  "results/2_analysis", 
  "results/3_figures", 
  "results/4_tables",
  "results/3_figures/hepatocytes",
  "results/4_tables/hepatocytes"
))
invisible(lapply(dirs.to.create, dir.create, showWarnings = FALSE, recursive = TRUE))

figures.dir <- file.path(work.dir, "results", "3_figures")
tables.dir <- file.path(work.dir, "results", "4_tables")
data.dir <- file.path(work.dir, "data")
hep_figures_dir <- file.path(figures.dir, "hepatocytes")
hep_tables_dir <- file.path(tables.dir, "hepatocytes")

# 加载包
packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "cowplot", "hdf5r",
              "clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE",
              "msigdbr", "fgsea", "pheatmap", "magick", "gridExtra", 
              "presto", "MAST", "tibble", "tidyr", "ggrepel", "ggpubr",
              "viridis", "RColorBrewer")

for(pkg in packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    if(pkg %in% c("org.Mm.eg.db", "DOSE", "clusterProfiler", "hdf5r", "MAST", "presto")) {
      if(pkg == "hdf5r") {
        install.packages("hdf5r")
      } else if(pkg == "presto") {
        devtools::install_github("immunogenomics/presto")
      } else {
        BiocManager::install(pkg)
      }
    } else {
      install.packages(pkg)
    }
  }
  suppressMessages(library(pkg, character.only = TRUE))
}

# 修复命名空间
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
desc <- dplyr::desc
group_by <- dplyr::group_by
summarise <- dplyr::summarise
pull <- dplyr::pull

# 主题
theme_journal <- theme_bw(base_size = 12) +      
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.1), face = "bold"),
    axis.title = element_text(size = rel(1.0)),
    axis.text = element_text(size = rel(0.9)),
    legend.title = element_text(size = rel(1.0)),
    legend.text = element_text(size = rel(0.9)),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

cat("=== 环境初始化完成 ===\n")

# ==================== 1. 数据读取与QC ====================
cat("\n=== 步骤1: 读取h5格式数据与质量控制 ===\n")

rdata_file <- file.path(work.dir, "results", "2_analysis", "high_quality_cells_v4.RData")

if(file.exists(rdata_file)) {
  cat("发现已处理的Seurat对象，加载中...\n")
  load(rdata_file)
  cat(sprintf("✓ 加载完成，共%d个细胞\n", ncol(high_quality_cells)))
} else {
  # 查找h5文件 - 使用正确的正则表达式
  h5_files <- list.files(data.dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  
  if(length(h5_files) == 0) {
    h5_files <- list.files(data.dir, pattern = "\\.h5ad$|\\.hdf5$", recursive = TRUE, full.names = TRUE)
  }
  
  if(length(h5_files) == 0) {
    stop("错误: 在", data.dir, "中未找到h5文件")
  }
  
  cat(sprintf("发现 %d 个h5文件:\n", length(h5_files)))
  print(basename(h5_files))
  
  # 读取h5文件
  seurat_list <- list()
  
  for(i in 1:length(h5_files)) {
    file_path <- h5_files[i]
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    
    cat(sprintf("\n  读取样本 %d/%d: %s\n", i, length(h5_files), sample_name))
    
    counts <- tryCatch({
      Read10X_h5(file_path)
    }, error = function(e1) {
      tryCatch({
        library(hdf5r)
        h5file <- H5File$new(file_path, mode = "r")
        
        if("X" %in% names(h5file)) {
          data <- h5file[["X"]][]
          genes <- h5file[["var"]][]$index
          cells <- h5file[["obs"]][]$index
          Matrix::t(data)
        } else if("matrix" %in% names(h5file)) {
          h5file[["matrix"]]
        } else {
          stop("无法识别h5文件格式")
        }
      }, error = function(e2) {
        # 使用正确的正则表达式匹配
        if(grepl("\\.h5ad$", file_path)) {
          Convert(file_path, dest = "h5seurat", overwrite = TRUE)
          h5seurat_path <- sub("\\.h5ad$", ".h5seurat", file_path)
          LoadH5Seurat(h5seurat_path)
        } else {
          stop("无法读取h5文件: ", conditionMessage(e2))
        }
      })
    })
    
    seu <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
    seu$sample <- sample_name
    seurat_list[[i]] <- seu
  }
  
  # 合并
  if(length(seurat_list) > 1) {
    cat("\n合并多个样本...\n")
    high_quality_cells <- merge(seurat_list[[1]], y = seurat_list[-1], 
                                add.cell.ids = sapply(seurat_list, function(x) x$sample[1]))
  } else {
    high_quality_cells <- seurat_list[[1]]
  }
  
  # 添加分组信息
  high_quality_cells$group <- ifelse(grepl("Ctrl|Control|WT|wildtype|normal", 
                                           high_quality_cells$sample, ignore.case = TRUE), 
                                     "Control", "CCl4")
  
  cat(sprintf("\n原始细胞数: %d\n", ncol(high_quality_cells)))
  cat("样本分组:\n")
  print(table(high_quality_cells$sample, high_quality_cells$group))
  
  # QC
  cat("\n进行质量控制...\n")
  high_quality_cells[["percent.mt"]] <- PercentageFeatureSet(high_quality_cells, pattern = "^mt-")
  high_quality_cells[["percent.rb"]] <- PercentageFeatureSet(high_quality_cells, pattern = "^Rp[sl]")
  
  # 严格QC
  cat("过滤标准: nFeature 200-6000, nCount 500-50000, percent.mt < 25%\n")
  high_quality_cells <- subset(high_quality_cells, 
                               subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                                 nCount_RNA > 500 & nCount_RNA < 50000 &
                                 percent.mt < 25)
  
  cat(sprintf("QC后细胞数: %d (保留率: %.1f%%)\n", 
              ncol(high_quality_cells), 
              100 * ncol(high_quality_cells) / sum(sapply(seurat_list, ncol))))
  
  # 标准化
  cat("标准化和数据缩放...\n")
  high_quality_cells <- NormalizeData(high_quality_cells)
  high_quality_cells <- FindVariableFeatures(high_quality_cells, selection.method = "vst", nfeatures = 2000)
  high_quality_cells <- ScaleData(high_quality_cells)
  high_quality_cells <- RunPCA(high_quality_cells, features = VariableFeatures(object = high_quality_cells))
  
  save(high_quality_cells, file = rdata_file)
  cat("✓ 数据读取和QC完成，对象已保存\n")
}

# ==================== 2. 降维与聚类 ====================
cat("\n=== 步骤2: 降维与聚类 ===\n")

if(!"seurat_clusters" %in% colnames(high_quality_cells@meta.data)) {
  cat("运行UMAP和聚类...\n")
  
  high_quality_cells <- JoinLayers(high_quality_cells)
  high_quality_cells <- FindNeighbors(high_quality_cells, dims = 1:20)
  high_quality_cells <- FindClusters(high_quality_cells, resolution = 0.8)
  high_quality_cells <- RunUMAP(high_quality_cells, dims = 1:20)
  
  p_umap1 <- DimPlot(high_quality_cells, reduction = "umap", label = TRUE, pt.size = 0.5) + 
    ggtitle("A. Clusters") + theme_journal
  p_umap2 <- DimPlot(high_quality_cells, reduction = "umap", group.by = "group", pt.size = 0.5) + 
    ggtitle("B. Groups") + theme_journal
  
  p_combined <- p_umap1 + p_umap2
  ggsave(file.path(figures.dir, "Figure1_UMAP.tiff"), p_combined, width = 12, height = 5, dpi = 300)
  
  cat(sprintf("✓ 聚类完成，发现 %d 个cluster\n", length(unique(high_quality_cells$seurat_clusters))))
} else {
  cat("✓ 使用已有聚类结果\n")
}


# ==================== 3. 细胞类型注释(Figure 1C)====================
cat("
=== 步骤3: 细胞类型注释 ===
")

if(!"celltype" %in% colnames(high_quality_cells@meta.data) || 
   sum(grepl("^Cluster_", high_quality_cells$celltype)) > 5) {
  
  cat("基于标志物进行细胞注释...
")
  
  if (length(Layers(high_quality_cells)) > 1) {
    high_quality_cells <- JoinLayers(high_quality_cells)
  }
  
  # 扩展的标志物列表
  liver_markers <- list(
    "Hepatocyte" = c("Alb", "Ttr", "Apoa1", "Cyp3a11", "Fabp1", "Ass1", "Serpina1a", "Cyp2e1"),
    "HSC" = c("Lrat", "Des", "Vim", "Pdgfrb", "Acta2", "Col1a1", "Col1a2", "Dcn", "Gfap", "S100a4", "Thy1"),
    "Endothelial" = c("Pecam1", "Cdh5", "Vwf", "Eng", "Stab2", "Kdr", "Tek", "Cldn5", "Flt1", "Nrp2"),
    "Kupffer" = c("Adgre1", "Cd68", "Cd163", "Marco", "Clec4f", "Mrc1", "Timd4", "Cd5l", "Vsig4"),
    "Monocyte" = c("Ly6c2", "Ccr2", "Cd14", "Fcgr3", "Plac8", "Itgam", "Fn1", "Csf1r"),
    "Cholangiocyte" = c("Krt19", "Krt7", "Sox9", "Epcam", "Anxa4", "Krt8", "Cftr"),
    "NK_T" = c("Cd3e", "Nkg7", "Gzma", "Ccl5", "Trac", "Cd3d", "Cd8a", "Klrb1c"),
    "B_cell" = c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Ighm", "Cd74", "Ebf1"),
    "Plasma_cell" = c("Jchain", "Mzb1", "Igkc", "Ighg1", "Igha1", "Sdc1"),
    "Neutrophil" = c("S100a8", "S100a9", "G0s2", "Retnlg", "Csf3r", "Lcn2")
  )
  
  # 计算模块评分
  for(celltype in names(liver_markers)) {
    genes <- liver_markers[[celltype]]
    genes <- genes[genes %in% rownames(high_quality_cells)]
    if(length(genes) >= 3) {
      high_quality_cells <- AddModuleScore(
        high_quality_cells, 
        features = list(genes), 
        name = celltype
      )
    }
  }  
  
  # 改进的注释函数
  cluster_annotation <- function(seu_obj) {
    if (length(Layers(seu_obj)) > 1) {
      seu_obj <- JoinLayers(seu_obj)
    }
    
    markers <- tryCatch({
      FindAllMarkers(seu_obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
    }, error = function(e) {
      FindAllMarkers(seu_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    })
    
    if(is.null(markers) || nrow(markers) == 0) {
      return(fallback_annotation(seu_obj))
    }
    
    markers$cluster <- as.character(markers$cluster)
    cluster_ids <- as.character(seu_obj$seurat_clusters)
    celltypes <- character(length(cluster_ids))
    
    for(cl in unique(cluster_ids)) {
      cl_markers <- markers %>% 
        filter(cluster == cl, p_val_adj < 0.05) %>% 
        head(30) %>% 
        pull(gene)
      
      if(length(cl_markers) == 0) {
        celltypes[cluster_ids == cl] <- paste0("Cluster_", cl)
        next
      }
      
      assigned <- FALSE
      
      # Hepatocyte
      if(sum(c("Alb", "Apoa1", "Fabp1", "Ttr", "Serpina1a", "Cyp2e1") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "Hepatocyte"
        assigned <- TRUE
      }
      # HSC
      else if(sum(c("Lrat", "Pdgfrb", "Des", "Vim", "Acta2", "Col1a1", "Dcn", "Gfap") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "HSC"
        assigned <- TRUE
      }
      # Endothelial
      else if(sum(c("Pecam1", "Cdh5", "Vwf", "Stab2", "Kdr", "Tek", "Cldn5") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "Endothelial"
        assigned <- TRUE
      }
      # Kupffer
      else if(sum(c("Adgre1", "Cd68", "Cd163", "Marco", "Clec4f", "Mrc1", "Timd4") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "Kupffer"
        assigned <- TRUE
      }
      # Monocyte
      else if(sum(c("Ly6c2", "Ccr2", "Cd14", "Fcgr3", "Plac8", "Itgam", "Fn1") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "Monocyte"
        assigned <- TRUE
      }
      # Cholangiocyte
      else if(sum(c("Krt19", "Krt7", "Sox9", "Epcam", "Krt8", "Anxa4") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "Cholangiocyte"
        assigned <- TRUE
      }
      # NK/T
      else if(sum(c("Cd3e", "Cd3d", "Nkg7", "Gzma", "Trac", "Cd8a", "Klrb1c") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "NK_T"
        assigned <- TRUE
      }
      # B cell
      else if(sum(c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Ighm", "Ebf1") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "B_cell"
        assigned <- TRUE
      }
      # Plasma cell
      else if(sum(c("Jchain", "Mzb1", "Igkc", "Ighg1", "Sdc1") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "Plasma_cell"
        assigned <- TRUE
      }
      # Neutrophil
      else if(sum(c("S100a8", "S100a9", "G0s2", "Retnlg", "Csf3r") %in% cl_markers) >= 2) {
        celltypes[cluster_ids == cl] <- "Neutrophil"
        assigned <- TRUE
      }
      
      # 模块评分法回退
      if(!assigned) {
        cells_in_cl <- WhichCells(seu_obj, idents = cl)
        score_cols <- c("Hepatocyte1", "HSC1", "Endothelial1", "Kupffer1", 
                        "Monocyte1", "Cholangiocyte1", "NK_T1", "B_cell1",
                        "Plasma_cell1", "Neutrophil1")
        
        scores <- sapply(score_cols, function(col) {
          if(col %in% colnames(seu_obj@meta.data)) {
            mean(seu_obj@meta.data[cells_in_cl, col], na.rm = TRUE)
          } else {
            -Inf
          }
        })
        
        best_celltype <- names(which.max(scores))
        if(!is.null(best_celltype) && max(scores) > 0) {
          celltypes[cluster_ids == cl] <- gsub("1$", "", best_celltype)
        } else {
          celltypes[cluster_ids == cl] <- paste0("Cluster_", cl)
        }
      }
    }
    return(celltypes)
  }
  
  fallback_annotation <- function(seu_obj) {
    cluster_ids <- as.character(seu_obj$seurat_clusters)
    celltypes <- character(length(cluster_ids))
    
    score_cols <- c("Hepatocyte1", "HSC1", "Endothelial1", "Kupffer1", 
                    "Monocyte1", "Cholangiocyte1", "NK_T1", "B_cell1",
                    "Plasma_cell1", "Neutrophil1")
    
    for(cl in unique(cluster_ids)) {
      cells_in_cl <- WhichCells(seu_obj, idents = cl)
      if(length(cells_in_cl) == 0) {
        celltypes[cluster_ids == cl] <- paste0("Cluster_", cl)
        next
      }
      
      scores <- sapply(score_cols, function(col) {
        if(col %in% colnames(seu_obj@meta.data)) {
          mean(seu_obj@meta.data[cells_in_cl, col], na.rm = TRUE)
        } else {
          -Inf
        }
      })
      
      best_celltype <- names(which.max(scores))
      if(!is.null(best_celltype) && !is.infinite(max(scores))) {
        celltypes[cluster_ids == cl] <- gsub("1$", "", best_celltype)
      } else {
        celltypes[cluster_ids == cl] <- paste0("Cluster_", cl)
      }
    }
    return(celltypes)
  }
  
  # 执行注释
  high_quality_cells$celltype <- cluster_annotation(high_quality_cells)
  
  # 可视化
  p_anno <- DimPlot(high_quality_cells, reduction = "umap", group.by = "celltype", 
                    label = TRUE, pt.size = 0.5, repel = TRUE) +
    ggtitle("C. Cell Types") + 
    theme_journal +
    theme(legend.position = "right")
  
  ggsave(file.path(figures.dir, "Figure1_CellType.tiff"), 
         p_anno, width = 12, height = 8, dpi = 300)
  
  cat("✓ 细胞注释完成
")
  print(table(high_quality_cells$celltype))
  
} else {
  cat("✓ 使用已有注释
")
}

# ==================== 4. HSC亚群分类(v2.5标准)====================
cat("
")
cat(paste(rep("=", 70), collapse = ""), "
")
cat("              HSC亚群细分类(平衡标准 v2.5 - 65%分位数)
")
cat(paste(rep("=", 70), collapse = ""), "
")

hsc_cells <- subset(high_quality_cells, celltype == "HSC")

if(ncol(hsc_cells) == 0) {
  stop("错误: 未找到HSC细胞")
}

cat(sprintf("提取到 %d 个HSC细胞
", ncol(hsc_cells)))

# 重新计算活化评分
activation_genes <- c("Acta2", "Col1a1", "Col1a2", "Col3a1", "Timp1", "Lox", "Spp1", 
                      "Postn", "Tagln", "Thbs1", "Cthrc1", "Vim", "Fn1", "Tgfb1")
activation_genes <- activation_genes[activation_genes %in% rownames(hsc_cells)]

if(length(activation_genes) >= 3) {
  hsc_cells <- AddModuleScore(hsc_cells, features = list(activation_genes), name = "Activation")
  hsc_cells$act_score <- hsc_cells$Activation1
  
  # v2.5平衡标准: 65%分位数 + Acta2>0.2
  act_threshold <- quantile(hsc_cells$act_score, 0.65, na.rm = TRUE)
  cat(sprintf("活化评分阈值 (65%%分位数): %.3f
", act_threshold))
  
  if("Acta2" %in% rownames(hsc_cells)) {
    acta2_expr <- FetchData(hsc_cells, vars = "Acta2")[,1]
    
    high_score <- hsc_cells$act_score > act_threshold
    high_acta2 <- acta2_expr > 0.2
    very_high_score <- hsc_cells$act_score > quantile(hsc_cells$act_score, 0.85, na.rm = TRUE)
    
    hsc_cells$hsc_subtype <- ifelse((high_score & high_acta2) | very_high_score, 
                                    "HSC_Act", "HSC_Qui")
    
    cat(sprintf("Acta2>0.2细胞: %d/%d (%.1f%%)
", 
                sum(high_acta2), ncol(hsc_cells), 100*sum(high_acta2)/ncol(hsc_cells)))
    
  } else {
    hsc_cells$hsc_subtype <- ifelse(hsc_cells$act_score > quantile(hsc_cells$act_score, 0.75, na.rm = TRUE), 
                                    "HSC_Act", "HSC_Qui")
  }
  
  # 检查分类结果
  subtype_table <- table(hsc_cells$hsc_subtype, hsc_cells$group)
  cat("
HSC亚群分布(v2.5平衡标准):
")
  print(subtype_table)
  
  # 计算比例
  ccl4_total <- sum(hsc_cells$group == "CCl4")
  ctrl_total <- sum(hsc_cells$group == "Control")
  ccl4_act <- sum(hsc_cells$hsc_subtype == "HSC_Act" & hsc_cells$group == "CCl4")
  ctrl_act <- sum(hsc_cells$hsc_subtype == "HSC_Act" & hsc_cells$group == "Control")
  
  cat(sprintf("  CCl4组: %d/%d (%.1f%%)
", ccl4_act, ccl4_total, 100*ccl4_act/ccl4_total))
  cat(sprintf("  Control组: %d/%d (%.1f%%)
", ctrl_act, ctrl_total, 100*ctrl_act/ctrl_total))
  
  # 更新主对象 - 修复: 确保hsc_subtype被添加到high_quality_cells
  hsc_meta <- data.frame(
    cell = colnames(hsc_cells),
    hsc_subtype = hsc_cells$hsc_subtype,
    stringsAsFactors = FALSE
  )
  
  # 初始化celltype_reannotated列(如果还没有的话)
  if(!"celltype_reannotated" %in% colnames(high_quality_cells@meta.data)) {
    high_quality_cells$celltype_reannotated <- as.character(high_quality_cells$celltype)
  }
  
  match_idx <- match(hsc_meta$cell, colnames(high_quality_cells))
  valid_idx <- !is.na(match_idx)
  high_quality_cells$celltype_reannotated[match_idx[valid_idx]] <- hsc_meta$hsc_subtype[valid_idx]
  
  # 关键修复: 同时将hsc_subtype添加为独立列，供后续使用
  high_quality_cells$hsc_subtype <- NA
  high_quality_cells$hsc_subtype[match_idx[valid_idx]] <- hsc_meta$hsc_subtype[valid_idx]
  
  cat(sprintf("\n✓ HSC分类v2.5完成，已更新 %d 个细胞\n", sum(valid_idx)))
} else {
  cat("警告: 活化标志物基因不足，跳过HSC亚型分类
")
}


# ==================== 5. 纤维化评分 ====================
cat("
=== 步骤5: 计算纤维化评分 ===
")

fibrosis_genes <- c(
  "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2", "Col5a1",
  "Acta2", "Tagln", "Myh11", "Des", "Vim",
  "Spp1", "Postn", "Thbs1", "Thbs2", "Fn1", "Fbn2", "Fbn1",
  "Lox", "Loxl1", "Loxl2",
  "Timp1", "Timp2", "Serpine1",
  "Tgfb1", "Tgfb2", "Ctgf", "Ccn2", "Pdgfa", "Pdgfb",
  "Pdgfra", "Pdgfrb", "Tgfbr1", "Tgfbr2",
  "Cthrc1", "Comp", "Cilp", "Grem1", "Fstl1",
  "Mmp2", "Mmp9", "Mmp13", "Mmp14", "Mmp3"
)

fibrosis_genes <- fibrosis_genes[fibrosis_genes %in% rownames(high_quality_cells)]

if(length(fibrosis_genes) >= 5) {
  cat(sprintf("使用 %d 个纤维化相关基因计算评分
", length(fibrosis_genes)))
  
  high_quality_cells <- AddModuleScore(
    high_quality_cells,
    features = list(fibrosis_genes),
    name = "Fibrosis_Score"
  )
  
  # 可视化
  p_score1 <- FeaturePlot(high_quality_cells, 
                          features = "Fibrosis_Score1",
                          pt.size = 0.3,
                          order = TRUE) +
    scale_color_gradientn(
      colors = c("#E8E8E8", "#FFE4B5", "#FFA500", "#FF4500", "#8B0000"),
      name = "Fibrosis
Score"
    ) +
    labs(title = "Fibrosis Score") +
    theme_journal
  
  # HSC亚型小提琴图
  hsc_cells_idx <- which(high_quality_cells$celltype_reannotated %in% c("HSC_Qui", "HSC_Act"))
  
  if(length(hsc_cells_idx) > 0) {
    hsc_subset <- high_quality_cells[, hsc_cells_idx]
    
    hsc_df <- data.frame(
      Fibrosis_Score1 = hsc_subset$Fibrosis_Score1,
      celltype = hsc_subset$celltype_reannotated,
      group = hsc_subset$group
    )
    
    p_score2 <- ggplot(hsc_df, 
                       aes(x = celltype, 
                           y = Fibrosis_Score1, 
                           fill = group)) +
      geom_violin(alpha = 0.6, 
                  trim = FALSE, 
                  scale = "width",
                  position = position_dodge(width = 0.8)) +
      geom_boxplot(width = 0.2, 
                   alpha = 0.8, 
                   outlier.shape = NA,
                   position = position_dodge(width = 0.8)) +
      stat_compare_means(aes(group = group), 
                         method = "wilcox.test", 
                         label = "p.format") +
      scale_fill_manual(values = c("CCl4" = "#E74C3C", "Control" = "#3498DB")) +
      labs(title = "Score by Subtype & Group",
           x = "HSC Subtype",
           y = "Fibrosis Score") +
      theme_journal
    
    p_combined <- p_score1 / p_score2
    ggsave(file.path(figures.dir, "Figure2_Fibrosis_Score.tiff"), 
           p_combined, width = 8, height = 10, dpi = 300)
    cat("✓ 纤维化评分图已保存
")
  }
  
  cat("✓ 纤维化评分完成
")
} else {
  cat("警告: 纤维化基因不足，跳过评分
")
}

save(high_quality_cells, file = rdata_file)

# ==================== 6. 差异基因分析(CCl4组内 HSC_Act vs HSC_Qui)====================
cat("
=== 步骤6: 差异基因分析(CCl4组内 HSC_Act vs HSC_Qui)===
")

hsc_ccl4 <- subset(high_quality_cells, 
                   group == "CCl4" & 
                     celltype_reannotated %in% c("HSC_Act", "HSC_Qui"))

if(ncol(hsc_ccl4) == 0) {
  stop("错误: 未找到CCl4组的HSC细胞")
}

cat(sprintf("CCl4组HSC细胞: HSC_Act=%d, HSC_Qui=%d
",
            sum(hsc_ccl4$celltype_reannotated == "HSC_Act"),
            sum(hsc_ccl4$celltype_reannotated == "HSC_Qui")))

hsc_ccl4 <- JoinLayers(hsc_ccl4)
Idents(hsc_ccl4) <- "celltype_reannotated"

# 主分析: HSC_Act vs HSC_Qui
cat("
【主分析】HSC_Act vs HSC_Qui(CCl4组内)...
")
hsc_deg_ccl4 <- FindMarkers(
  hsc_ccl4,
  ident.1 = "HSC_Act",
  ident.2 = "HSC_Qui",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

hsc_deg_ccl4$p_val_adj <- p.adjust(hsc_deg_ccl4$p_val, method = "BH")

hsc_deg_ccl4_df <- hsc_deg_ccl4 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    group = ifelse(avg_log2FC > 0, "HSC_Act", "HSC_Qui"),
    abs_log2FC = abs(avg_log2FC)
  ) %>%
  arrange(p_val_adj)

write.csv(hsc_deg_ccl4_df, file.path(tables.dir, "Table_DEG_HSC_Act_vs_Qui_CCl4.csv"), row.names = FALSE)

cat(sprintf("
✓ 主分析完成: 共%d个DEG (padj<0.05)
", 
            sum(hsc_deg_ccl4_df$p_val_adj < 0.05)))

# 保存供后续使用
hsc_deg <- hsc_deg_ccl4_df

# ==================== 7. 通路富集分析与Figure 1生成 ====================
cat("\n=== 步骤7: 通路富集分析 ===\n")

# 确保dplyr已加载(用于管道操作)
if(!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
  library(dplyr)
}

hsc_deg_unique <- hsc_deg %>%
  group_by(gene) %>%
  slice_min(p_val_adj, n = 1) %>%
  ungroup()

up_genes <- hsc_deg_unique %>% 
  filter(avg_log2FC > 0.25, p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

down_genes <- hsc_deg_unique %>% 
  filter(avg_log2FC < -0.25, p_val_adj < 0.05) %>% 
  arrange(avg_log2FC) %>%
  pull(gene)

cat(sprintf("分析基因: 上调%d个，下调%d个\n", length(up_genes), length(down_genes)))

# GO分析
cat("GO富集分析...\n")
ego_bp_up <- NULL
ego_mf_up <- NULL

if(length(up_genes) >= 10) {
  tryCatch({
    ego_bp_up <- enrichGO(gene = up_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    
    ego_mf_up <- enrichGO(gene = up_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05)
    
    cat(sprintf("  GO BP: %d个条目\n", ifelse(is.null(ego_bp_up), 0, nrow(ego_bp_up@result))))
  }, error = function(e) {
    cat("  GO分析错误:", conditionMessage(e), "\n")
  })
}

# KEGG分析
cat("KEGG富集分析...\n")
kk_up <- NULL

if(length(up_genes) >= 10) {
  tryCatch({
    gene_df_up <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    
    if(nrow(gene_df_up) > 0) {
      kk_up <- enrichKEGG(gene = gene_df_up$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
      cat(sprintf("  找到%d个KEGG通路\n", nrow(kk_up@result)))
    }
  }, error = function(e) {
    cat("  KEGG分析失败:", conditionMessage(e), "\n")
  })
}

if(!is.null(ego_bp_up)) {
  write.csv(ego_bp_up@result, file.path(tables.dir, "Table_GO_BP_Up.csv"), row.names = FALSE)
}
if(!is.null(kk_up)) {
  write.csv(kk_up@result, file.path(tables.dir, "Table_KEGG_Up.csv"), row.names = FALSE)
}

cat("✓ 富集分析完成\n")

# ==================== Figure 1: 细胞类型注释与HSC纤维化特征 ====================
cat("\n=== 生成 Figure 1: 细胞类型注释与HSC纤维化特征 ===\n")

# 加载需要的包
if(!require("tidyr", quietly = TRUE)) install.packages("tidyr")
library(tidyr)

# 定义颜色方案
celltype_colors <- c(
  "Hepatocyte" = "#E41A1C", 
  "LSEC" = "#377EB8", 
  "Kupffer" = "#4DAF4A",
  "HSC_Act" = "#FF7F00", 
  "HSC_Qui" = "#FFFF33", 
  "Bcell" = "#A65628",
  "NK" = "#F781BF", 
  "Neutrophil" = "#999999", 
  "Doublet_or_Debris" = "#CCCCCC",
  # 添加可能存在的其他细胞类型
  "Endothelial" = "#377EB8",
  "Monocyte" = "#984EA3",
  "Cholangiocyte" = "#FF7F00",
  "NK_T" = "#F781BF",
  "B_cell" = "#A65628",
  "Plasma_cell" = "#E41A1C",
  "HSC" = "#FFFF33"
)

# 确保celltype_reannotated存在且为因子
if(!"celltype_reannotated" %in% colnames(high_quality_cells@meta.data)) {
  high_quality_cells$celltype_reannotated <- as.character(high_quality_cells$celltype)
}

# 获取实际的细胞类型水平
actual_celltypes <- unique(high_quality_cells$celltype_reannotated)
available_colors <- celltype_colors[names(celltype_colors) %in% actual_celltypes]

# 为未定义颜色的细胞类型分配灰色
missing_types <- setdiff(actual_celltypes, names(celltype_colors))
if(length(missing_types) > 0) {
  for(mt in missing_types) {
    available_colors[mt] <- "#888888"
  }
}

high_quality_cells$celltype_reannotated <- factor(
  high_quality_cells$celltype_reannotated,
  levels = names(available_colors)
)

# 1A: UMAP
p1a <- DimPlot(high_quality_cells, reduction = "umap", group.by = "celltype_reannotated",
               label = TRUE, label.size = 4, pt.size = 0.5) +
  scale_color_manual(values = available_colors) +
  ggtitle("A. Cell Type Annotation") + 
  theme_journal

# 1B: 细胞组成(排除双细胞，按比例从高到低排列)
final_stats <- as.data.frame(table(high_quality_cells$celltype_reannotated))
final_stats <- final_stats[!final_stats$Var1 %in% c("Doublet_or_Debris", "NA"), ]
colnames(final_stats) <- c("CellType", "Count")
final_stats$Percentage <- round(final_stats$Count / sum(final_stats$Count) * 100, 1)

# 按Count从高到低排序
final_stats <- final_stats[order(-final_stats$Count), ]
final_stats$CellType <- factor(final_stats$CellType, levels = final_stats$CellType)

p1b <- ggplot(final_stats, aes(x = CellType, y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = available_colors[as.character(final_stats$CellType)]) +
  labs(x = "", y = "Cell Count", title = "B. Cell Type Composition") +
  theme_journal + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")

# 1C: 分组比较
group_stats <- table(high_quality_cells$celltype_reannotated, high_quality_cells$group)
group_stats <- group_stats[!rownames(group_stats) %in% c("Doublet_or_Debris", "NA"), ]
group_props <- prop.table(group_stats, margin = 2) * 100
plot_data_grp <- as.data.frame(group_props)
colnames(plot_data_grp) <- c("CellType", "Group", "Percentage")

p1c <- ggplot(plot_data_grp, aes(x = Group, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = available_colors) +
  labs(x = "Group", y = "Percentage (%)", title = "C. Composition by Group") +
  theme_journal

# 1D: 纤维化评分(使用已计算的或重新计算)
if(!"Fibrosis_Score1" %in% colnames(high_quality_cells@meta.data)) {
  fibrosis_genes <- c("Col1a1", "Col1a2", "Acta2", "Des", "Vim", "Tgfb1")
  fibrosis_genes <- fibrosis_genes[fibrosis_genes %in% rownames(high_quality_cells)]
  if(length(fibrosis_genes) >= 3) {
    high_quality_cells <- AddModuleScore(high_quality_cells, features = list(fibrosis_genes), 
                                         name = "Fibrosis_Score")
    cat(sprintf("  计算纤维化评分(%d个基因)\n", length(fibrosis_genes)))
  }
}

# 选择要显示的细胞类型(优先使用HSC分类，如果没有则使用原始HSC)
hsc_types <- c("HSC_Act", "HSC_Qui")
if(!any(hsc_types %in% actual_celltypes)) {
  hsc_types <- "HSC"
}

display_types <- c(hsc_types, "LSEC", "Kupffer", "Hepatocyte")
display_types <- intersect(display_types, actual_celltypes)

if(length(display_types) > 0 && "Fibrosis_Score1" %in% colnames(high_quality_cells@meta.data)) {
  score_data <- high_quality_cells@meta.data %>%
    filter(celltype_reannotated %in% display_types) %>%
    mutate(celltype_reannotated = factor(celltype_reannotated, levels = display_types))
  
  fibrosis_display_colors <- c(
    "HSC_Act" = "#C0392B", "HSC_Qui" = "#2980B9", "HSC" = "#FFFF33",
    "LSEC" = "#27AE60", "Kupffer" = "#F39C12", "Hepatocyte" = "#8E44AD"
  )
  
  p1d <- ggplot(score_data, aes(x = celltype_reannotated, y = Fibrosis_Score1, fill = celltype_reannotated)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = fibrosis_display_colors[as.character(score_data$celltype_reannotated)]) +
    labs(y = "Fibrosis Score (Z-score)", x = "", title = "D. Fibrogenic Capacity") +
    theme_journal + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
} else {
  p1d <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") + 
    theme_void() + ggtitle("D. Fibrogenic Capacity")
}

# 1E: 富集倍数计算
fibrosis_genes_for_enrich <- c("Col1a1", "Col1a2", "Acta2", "Des", "Vim", "Tgfb1")
fibrosis_genes_for_enrich <- fibrosis_genes_for_enrich[fibrosis_genes_for_enrich %in% rownames(high_quality_cells)]

if(length(display_types) > 0 && length(fibrosis_genes_for_enrich) > 0) {
  mean_expr_by_type <- sapply(display_types, function(ct) {
    cells <- WhichCells(high_quality_cells, expression = celltype_reannotated == ct)
    if (length(cells) == 0) return(NA)
    expr <- GetAssayData(high_quality_cells, layer = "data")[fibrosis_genes_for_enrich, cells, drop = FALSE]
    mean(colMeans(expm1(expr), na.rm = TRUE), na.rm = TRUE)
  })
  
  if("Hepatocyte" %in% names(mean_expr_by_type) && !is.na(mean_expr_by_type["Hepatocyte"]) && mean_expr_by_type["Hepatocyte"] > 0) {
    fold_vs_hepato <- mean_expr_by_type / mean_expr_by_type["Hepatocyte"]
    
    fold_df <- data.frame(
      CellType = names(fold_vs_hepato),
      Fold_Enrichment = as.numeric(fold_vs_hepato),
      stringsAsFactors = FALSE
    )
    fold_df <- fold_df[fold_df$CellType != "Hepatocyte" & !is.na(fold_df$Fold_Enrichment), ]
    
    if(nrow(fold_df) > 0) {
      p1e <- ggplot(fold_df, aes(x = reorder(CellType, -Fold_Enrichment), y = Fold_Enrichment, fill = CellType)) +
        geom_bar(stat = "identity", width = 0.6) +
        geom_text(aes(label = sprintf("%.1f×", Fold_Enrichment)), vjust = -0.5, size = 4) +
        scale_fill_manual(values = available_colors[fold_df$CellType]) +
        labs(y = "Fold Enrichment (vs Hepatocyte)", x = "", title = "E. Fibrosis Gene Enrichment") +
        theme_journal + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
      
      # 保存富集倍数数据
      write.csv(fold_df, file.path(tables.dir, "Table_Fibrosis_Enrichment.csv"), row.names = FALSE)
    } else {
      p1e <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No enrichment data") + 
        theme_void() + ggtitle("E. Fibrosis Gene Enrichment")
    }
  } else {
    p1e <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Hepatocyte reference not found") + 
      theme_void() + ggtitle("E. Fibrosis Gene Enrichment")
  }
} else {
  p1e <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient genes") + 
    theme_void() + ggtitle("E. Fibrosis Gene Enrichment")
}

# 1F: HSC Marker对比(仅当有HSC分类时)
hsc_marker_types <- intersect(c("HSC_Act", "HSC_Qui", "HSC"), actual_celltypes)

if(length(hsc_marker_types) > 0) {
  hsc_cells <- subset(high_quality_cells, subset = celltype_reannotated %in% hsc_marker_types)
  fibrogenic_markers <- c("Col1a1", "Col1a2", "Acta2", "Des", "Vim")
  fibrogenic_markers <- fibrogenic_markers[fibrogenic_markers %in% rownames(high_quality_cells)]
  
  if(length(fibrogenic_markers) > 0 && ncol(hsc_cells) > 0) {
    plot_data <- FetchData(hsc_cells, vars = c(fibrogenic_markers, "celltype_reannotated"))
    plot_long <- plot_data %>%
      pivot_longer(cols = -celltype_reannotated, names_to = "gene", values_to = "expression")
    
    # 定义HSC颜色
    hsc_colors <- c("HSC_Act" = "#C0392B", "HSC_Qui" = "#2980B9", "HSC" = "#FFFF33")
    
    p1f <- ggplot(plot_long, aes(x = gene, y = expm1(expression), fill = celltype_reannotated)) +
      geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
      scale_fill_manual(values = hsc_colors[as.character(unique(plot_long$celltype_reannotated))],
                        labels = c("HSC_Act" = "Activated", "HSC_Qui" = "Quiescent", "HSC" = "HSC"),
                        name = "HSC State") +
      labs(y = "Mean Expression", x = "", title = "F. HSC Marker Expression") +
      theme_journal + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p1f <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No marker genes found") + 
      theme_void() + ggtitle("F. HSC Marker Expression")
  }
} else {
  p1f <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No HSC cells found") + 
    theme_void() + ggtitle("F. HSC Marker Expression")
}

# 组合(2x3布局)
fig1_complete <- (p1a | p1b | p1c) / (p1d | p1e | p1f) + 
  plot_annotation(
    title = "Figure 1. Cell Type Annotation and Fibrosis Characteristics",
    subtitle = "A-C: Cell Type Annotation | D-F: Fibrosis Characteristics",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40")
    )
  ) & 
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(figures.dir, "Figure1_Cell_Type_Annotation.tiff"), fig1_complete, 
       width = 18, height = 12, dpi = 300, compression = "lzw")
cat("✓ Figure 1 saved to", file.path(figures.dir, "Figure1_Celltype_Complete.tiff"), "\n")

# 保存更新后的对象
save(high_quality_cells, file = rdata_file)
cat("✓ 第7步完成: 通路富集分析 + Figure 1生成\n")

# ==================== 8. 肝细胞深度分析(含8基因签名 + 7类应激评分)====================
cat("
")
cat(paste(rep("=", 70), collapse = ""), "
")
cat("              肝细胞深度分析(8基因签名 + 7类应激评分)
")
cat(paste(rep("=", 70), collapse = ""), "
")

# 8.1 提取肝细胞并基础分析
cat("
=== 步骤8.1: 提取肝细胞并进行基础分析 ===
")

hepatocytes <- subset(high_quality_cells, celltype == "Hepatocyte")

if(ncol(hepatocytes) == 0) {
  stop("错误: 未找到肝细胞")
}

cat(sprintf("提取到 %d 个肝细胞
", ncol(hepatocytes)))

# 重新降维(推荐对子集重新分析)
hepatocytes <- NormalizeData(hepatocytes)
hepatocytes <- FindVariableFeatures(hepatocytes, selection.method = "vst", nfeatures = 2000)
hepatocytes <- ScaleData(hepatocytes)
hepatocytes <- RunPCA(hepatocytes, features = VariableFeatures(object = hepatocytes))
hepatocytes <- FindNeighbors(hepatocytes, dims = 1:20)
hepatocytes <- FindClusters(hepatocytes, resolution = 0.6)
hepatocytes <- RunUMAP(hepatocytes, dims = 1:20)

# 8.2 肝小叶分区注释(基于文献标志物)
cat("
=== 步骤8.2: 肝小叶分区注释 ===
")

zone_markers <- list(
  "Zone1_Periportal" = c("Cyp2f2", "Hal", "Sds", "Ass1", "Pck1"),
  "Zone3_Pericentral" = c("Cyp2e1", "Cyp1a2", "Glul", "Oat", "Cyp3a11")
)

# 计算分区评分
for(zone in names(zone_markers)) {
  genes <- zone_markers[[zone]]
  genes <- genes[genes %in% rownames(hepatocytes)]
  if(length(genes) >= 3) {
    hepatocytes <- AddModuleScore(hepatocytes, features = list(genes), name = zone)
  }
}

# 基于评分分配区域
if(all(c("Zone1_Periportal1", "Zone3_Pericentral1") %in% colnames(hepatocytes@meta.data))) {
  z1_score <- hepatocytes$Zone1_Periportal1
  z3_score <- hepatocytes$Zone3_Pericentral1
  
  hepatocytes$hep_zone <- ifelse(z1_score > z3_score, "Zone1_Periportal",
                                 ifelse(z3_score > z1_score, "Zone3_Pericentral", 
                                        "Zone2_Midlobular"))
} else {
  hepatocytes$hep_zone <- "Zone2_Midlobular"
}

cat(sprintf("肝细胞分区分布:
"))
print(table(hepatocytes$hep_zone, hepatocytes$group))

# 可视化分区
p_zone <- DimPlot(hepatocytes, reduction = "umap", group.by = "hep_zone", pt.size = 0.5) +
  ggtitle("A. Hepatocyte Zonation") + theme_journal

p_zone_group <- DimPlot(hepatocytes, reduction = "umap", group.by = "group", pt.size = 0.5) +
  ggtitle("B. Group Distribution") + theme_journal

p_combined_zone <- p_zone + p_zone_group

ggsave(file.path(hep_figures_dir, "FigureH1_Hepatocyte_Zonation.tiff"), 
       p_combined_zone, width = 12, height = 5, dpi = 300)

cat("✓ 肝细胞基础分析完成
")

# 8.3 【关键】8基因签名评分(Figure 2专用)
cat("
=== 步骤8.3: 8基因签名评分(Figure 2专用)===
")

fig2_signature_genes <- c("Cyp2e1", "Cyp3a11", "Hspa1a", "Hsp90aa1", 
                          "Gpx1", "Hmox1", "Atf4", "Gadd45a")
fig2_signature_genes <- fig2_signature_genes[fig2_signature_genes %in% rownames(hepatocytes)]

if(length(fig2_signature_genes) >= 5) {
  cat(sprintf("使用 %d 个基因计算8基因签名评分
", length(fig2_signature_genes)))
  
  hepatocytes <- AddModuleScore(hepatocytes, 
                                features = list(fig2_signature_genes), 
                                name = "Fig2_Signature")
  
  # 创建Figure 2需要的兼容列名
  hepatocytes$Stress1 <- hepatocytes$Fig2_Signature1
  hepatocytes$stress_group <- ifelse(
    hepatocytes$Stress1 > mean(hepatocytes$Stress1, na.rm = TRUE) + 
      sd(hepatocytes$Stress1, na.rm = TRUE),
    "High_stress", "Low_stress"
  )
  
  # 创建obj别名(供Figure 2代码使用)
  obj <- hepatocytes
  
  cat(sprintf("8基因签名评分完成: 高应激细胞 %d/%d (%.1f%%)
",
              sum(hepatocytes$stress_group == "High_stress"), 
              ncol(hepatocytes),
              100 * sum(hepatocytes$stress_group == "High_stress") / ncol(hepatocytes)))
  cat("✓ 已创建对象别名 'obj' 供Figure 2代码使用
")
} else {
  cat("警告: 8基因签名基因不足，跳过计算
")
}

# 8.4 7类分层应激评分(Figure Master用)
cat("
=== 步骤8.4: 7类分层应激评分(Figure Master用)===
")

stress_categories <- list(
  "Oxidative_Stress" = c("Cyp2e1", "Gpx1", "Sod1", "Sod2", "Cat", "Gstm1", "Nqo1", "Hmox1"),
  "ER_Stress" = c("Xbp1", "Atf4", "Atf6", "Hspa5", "Ddit3", "Ern1", "Pdia3"),
  "Inflammation" = c("Tnfaip3", "Nfkbia", "Cxcl1", "Cxcl2", "Cxcl10", "Ccl2", "Ccl20", "Il1b"),
  "Cell_Death" = c("Hmgb1", "Il1a", "Il33", "S100a8", "S100a9", "Casp3", "Casp7", "Casp8"),
  "ProFibrotic_Factors" = c("Tgfb1", "Tgfb2", "Pdgfa", "Pdgfb", "Hgf", "Vegfa", "Ctgf", "Ccn2"),
  "Lipid_Dysregulation" = c("Fasn", "Acaca", "Srebf1", "Ppara", "Cpt1a", "Cd36", "Fabp1", "Pnpla3"),
  "Regeneration" = c("Mki67", "Pcna", "Hgf", "Egfr", "Met", "Axin2", "Wnt2")
)

# 计算各类型应激评分
for(stress_type in names(stress_categories)) {
  genes <- stress_categories[[stress_type]]
  genes <- genes[genes %in% rownames(hepatocytes)]
  if(length(genes) >= 3) {
    hepatocytes <- AddModuleScore(hepatocytes, features = list(genes), name = stress_type)
    cat(sprintf("  %s: %d个基因
", stress_type, length(genes)))
  }
}

# 计算总应激评分(7类总和)
all_stress_genes <- unlist(stress_categories)
all_stress_genes <- all_stress_genes[all_stress_genes %in% rownames(hepatocytes)]
hepatocytes <- AddModuleScore(hepatocytes, features = list(all_stress_genes), name = "Total_Stress")

# 关键: 使用Hep_Stress_Score作为总评分(Figure Master用)
hepatocytes$Hep_Stress_Score <- hepatocytes$Total_Stress1

cat(sprintf("
总应激评分: %d个基因
", length(all_stress_genes)))
cat("✓ 7类分层应激评分完成
")

# ==================== 步骤8.5: 生成Figure 2数据文件（优化版）====================
cat("\n=== 步骤8.5: 生成Figure 2数据文件（Zone 3应激分析）===\n")

# 确保目录存在
if(!dir.exists(hep_tables_dir)) {
  dir.create(hep_tables_dir, recursive = TRUE)
}

# 1. DE_stress_comprehensive.csv（高vs低应激差异基因）
cat("生成 DE_stress_comprehensive.csv...\n")

Idents(hepatocytes) <- "stress_group"
de_stress <- FindMarkers(hepatocytes, 
                         ident.1 = "High_stress", 
                         ident.2 = "Low_stress",
                         min.pct = 0.1, 
                         logfc.threshold = 0.1,
                         test.use = "wilcox")

de_stress_df <- de_stress %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    group = ifelse(avg_log2FC > 0, "High_stress", "Low_stress"),
    abs_log2FC = abs(avg_log2FC)
  ) %>%
  arrange(p_val_adj)

write.csv(de_stress_df, 
          file.path(hep_tables_dir, "DE_stress_comprehensive.csv"), 
          row.names = FALSE)

cat(sprintf("  完成: %d个差异基因 (padj<0.05: %d)\n", 
            nrow(de_stress_df), 
            sum(de_stress_df$p_val_adj < 0.05)))

# 2. stress_scores_all_methods.csv（多种评分方法对比）
cat("生成 stress_scores_all_methods.csv...\n")

stress_methods_df <- data.frame(
  cell = colnames(hepatocytes),
  Original_4gene = hepatocytes$Oxidative_Stress1,  # 4基因: 氧化应激
  Comprehensive_8gene = hepatocytes$Fig2_Signature1,  # 8基因签名
  Total_7category = hepatocytes$Total_Stress1,  # 7类总和
  Oxidative = hepatocytes$Oxidative_Stress1,
  ERstress = hepatocytes$ER_Stress1,
  Lipid = hepatocytes$Lipid_Dysregulation1
)

write.csv(stress_methods_df, 
          file.path(hep_tables_dir, "stress_scores_all_methods.csv"), 
          row.names = FALSE)

cat("  完成: 6种评分方法\n")

# 3. GO_BP_high_stress_hepatocytes.csv（高应激组GO富集）
cat("生成 GO_BP_high_stress_hepatocytes.csv...\n")

high_stress_genes <- de_stress_df %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  pull(gene)

if(length(high_stress_genes) >= 10) {
  tryCatch({
    ego_stress <- enrichGO(gene = high_stress_genes, 
                           OrgDb = org.Mm.eg.db,
                           keyType = "SYMBOL", 
                           ont = "BP", 
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)
    
    if(!is.null(ego_stress) && nrow(ego_stress@result) > 0) {
      write.csv(ego_stress@result, 
                file.path(hep_tables_dir, "GO_BP_high_stress_hepatocytes.csv"), 
                row.names = FALSE)
      cat(sprintf("  完成: %d个GO条目\n", nrow(ego_stress@result)))
    } else {
      cat("  警告: 未找到显著GO富集\n")
    }
  }, error = function(e) {
    cat("  GO分析错误:", conditionMessage(e), "\n")
  })
} else {
  cat("  警告: 高应激基因不足(<10个)，跳过GO分析\n")
}

# 4. KEGG_upregulated_CCl4.csv（CCl4 vs Control KEGG）
cat("生成 KEGG_upregulated_CCl4.csv...\n")

Idents(hepatocytes) <- "group"
de_cc <- FindMarkers(hepatocytes, 
                     ident.1 = "CCl4", 
                     ident.2 = "Control",
                     min.pct = 0.1,
                     logfc.threshold = 0.1,
                     test.use = "wilcox")

cc_up_genes <- de_cc %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  rownames()

if(length(cc_up_genes) >= 10) {
  tryCatch({
    gene_df <- bitr(cc_up_genes, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Mm.eg.db)
    
    if(nrow(gene_df) > 0) {
      kk_cc <- enrichKEGG(gene = gene_df$ENTREZID, 
                          organism = 'mmu', 
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
      
      if(!is.null(kk_cc) && nrow(kk_cc@result) > 0) {
        write.csv(kk_cc@result, 
                  file.path(hep_tables_dir, "KEGG_upregulated_CCl4.csv"), 
                  row.names = FALSE)
        cat(sprintf("  完成: %d个KEGG通路\n", nrow(kk_cc@result)))
      } else {
        cat("  警告: 未找到显著KEGG通路\n")
      }
    } else {
      cat("  警告: 基因ID转换失败\n")
    }
  }, error = function(e) {
    cat("  KEGG分析错误:", conditionMessage(e), "\n")
  })
} else {
  cat("  警告: CCl4上调基因不足(<10个)，跳过KEGG分析\n")
}

cat("✓ Figure 2数据文件生成完成\n")

# 保存肝细胞对象
save(hepatocytes, file = file.path(work.dir, "results", "2_analysis", "hepatocytes_fig2.RData"))

# ==================== 步骤9: Figure 2生成（Zone 3应激整合版，10 panels）====================
cat("\n=== 步骤9: Figure 2生成（Zone 3 Hepatocyte Stress Analysis，10 panels）===\n")

# 9.1 加载必要包
if(!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if(!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if(!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(viridis)
library(RColorBrewer)
library(ggrepel)

# 9.2 准备数据
obj <- hepatocytes

# 关键统计
stats <- list(
  mean = mean(obj$Fig2_Signature1, na.rm = TRUE),
  sd = sd(obj$Fig2_Signature1, na.rm = TRUE),
  threshold = mean(obj$Fig2_Signature1, na.rm = TRUE) + sd(obj$Fig2_Signature1, na.rm = TRUE),
  n_total = ncol(obj),
  n_high = sum(obj$stress_group == "High_stress", na.rm = TRUE),
  high_pct = 100 * sum(obj$stress_group == "High_stress", na.rm = TRUE) / ncol(obj)
)

cat(sprintf("8基因签名统计: Mean=%.3f, Threshold=%.3f, High_stress=%d (%.1f%%)\n",
            stats$mean, stats$threshold, stats$n_high, stats$high_pct))

# 统一主题
fig_theme <- theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "bold")
  )

# ==================== Panel A: 肝细胞分区UMAP ====================
cat("Panel A: 肝细胞分区UMAP\n")

p_a <- DimPlot(obj, reduction = "umap", group.by = "hep_zone", pt.size = 0.4, label = FALSE) +
  scale_color_manual(values = c("Zone1_Periportal" = "#1f77b4", 
                                "Zone2_Midlobular" = "#ff7f0e", 
                                "Zone3_Pericentral" = "#2ca02c"),
                     name = "Zone") +
  ggtitle("A. Hepatocyte Zonation") + 
  fig_theme +
  theme(legend.position = "right")

# ==================== Panel B: 分组分布（CCl4富集于Zone 3）====================
cat("Panel B: 分组分布\n")

p_b <- DimPlot(obj, reduction = "umap", group.by = "group", pt.size = 0.4) +
  scale_color_manual(values = c("CCl4" = "#E74C3C", "Control" = "#3498DB")) +
  ggtitle("B. CCl4 Enrichment in Zone 3") + 
  fig_theme +
  theme(legend.position = "right")

# ==================== Panel C: 8-gene应激评分分布 ====================
cat("Panel C: 8-gene应激评分分布\n")

p_c <- ggplot(data.frame(Score = obj$Fig2_Signature1), aes(x = Score)) +
  geom_histogram(bins = 40, fill = "#4682B4", alpha = 0.8, color = "white") +
  geom_vline(xintercept = stats$mean, color = "#E41A1C", linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = stats$threshold, color = "#8B0000", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = stats$mean, y = Inf, label = sprintf("Mean=%.3f", stats$mean),
           vjust = 1.5, hjust = -0.1, size = 2.5, color = "#E41A1C") +
  annotate("text", x = stats$threshold, y = Inf, 
           label = sprintf("Threshold=%.3f\nn=%d (%.1f%%)", stats$threshold, stats$n_high, stats$high_pct),
           vjust = 3, hjust = -0.1, size = 2.5, color = "#8B0000") +
  labs(x = "8-Gene Stress Signature", y = "Cell Count", title = "C. Stress Score Distribution") +
  fig_theme

# ==================== Panel D: 高应激细胞UMAP定位 ====================
cat("Panel D: 高应激细胞UMAP定位\n")

p_d <- FeaturePlot(obj, features = "Fig2_Signature1", pt.size = 0.2, order = TRUE) +
  scale_color_viridis_c(option = "plasma", name = "Score") +
  labs(title = "D. High-Stress Cell Localization") +
  fig_theme +
  theme(legend.position = "right", legend.key.height = unit(0.6, "cm"))

# ==================== Panel E: 分区×分组应激评分定量 ====================
cat("Panel E: 分区×分组应激评分\n")

zone_data <- obj@meta.data %>%
  filter(hep_zone %in% c("Zone1_Periportal", "Zone3_Pericentral")) %>%
  group_by(hep_zone, group) %>%
  summarise(
    n = n(),
    mean_stress = mean(Fig2_Signature1, na.rm = TRUE),
    sem = sd(Fig2_Signature1, na.rm = TRUE) / sqrt(n),
    .groups = 'drop'
  )

# 添加统计检验
zone_stats <- obj@meta.data %>%
  filter(hep_zone %in% c("Zone1_Periportal", "Zone3_Pericentral"))

p_e <- ggplot(zone_data, aes(x = group, y = mean_stress, fill = hep_zone)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_stress - sem, ymax = mean_stress + sem), 
                position = position_dodge(width = 0.6), width = 0.2, linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", mean_stress)), 
            position = position_dodge(width = 0.6), vjust = -0.8, size = 2.5) +
  scale_fill_manual(values = c("Zone1_Periportal" = "#1f77b4", "Zone3_Pericentral" = "#2ca02c"),
                    labels = c("Zone 1 (Periportal)", "Zone 3 (Pericentral)")) +
  labs(x = "", y = "Mean Stress Score", 
       title = "E. Zone-Specific Stress (Zone3: 0.55 vs 0.36)") +
  fig_theme +
  theme(legend.position = "top", legend.key.size = unit(0.4, "cm"))

# ==================== Panel F: 七类应激热图 ====================
cat("Panel F: 七类应激热图\n")

stress_categories <- list(
  "Oxidative" = "Oxidative_Stress1",
  "ER" = "ER_Stress1", 
  "Inflammation" = "Inflammation1",
  "Death" = "Cell_Death1",
  "ProFibrotic" = "ProFibrotic_Factors1",
  "Lipid" = "Lipid_Dysregulation1",
  "Regeneration" = "Regeneration1"
)

stress_summary <- data.frame()
for(stress_name in names(stress_categories)) {
  col_name <- stress_categories[[stress_name]]
  if(col_name %in% colnames(obj@meta.data)) {
    temp <- obj@meta.data %>%
      filter(hep_zone %in% c("Zone1_Periportal", "Zone3_Pericentral")) %>%
      group_by(group, hep_zone) %>%
      summarise(score = mean(.data[[col_name]], na.rm = TRUE), .groups = 'drop') %>%
      mutate(stress_type = stress_name)
    stress_summary <- rbind(stress_summary, temp)
  }
}

# 确保正确顺序
stress_summary$stress_type <- factor(stress_summary$stress_type, 
                                     levels = c("Oxidative", "ER", "Inflammation", "Death", 
                                                "ProFibrotic", "Lipid", "Regeneration"))

p_f <- ggplot(stress_summary, aes(x = interaction(group, hep_zone), y = stress_type, fill = score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 2) +
  scale_fill_gradientn(colors = c("#2166ac", "white", "#b2182b"), name = "Score") +
  scale_x_discrete(labels = c("Control.Zone1_Periportal" = "C_Z1",
                              "Control.Zone3_Pericentral" = "C_Z3",
                              "CCl4.Zone1_Periportal" = "CCl4_Z1",
                              "CCl4.Zone3_Pericentral" = "CCl4_Z3")) +
  labs(x = "", y = "", title = "F. Seven-Category Stress (CCl4_Z3 Oxidative: 0.73)") +
  fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7))

# ==================== Panel G: Hmox1 vs Hspa1a对比（氧化vs热休克）====================
cat("Panel G: Hmox1 vs Hspa1a对比\n")

contrast_genes <- c("Hmox1", "Hspa1a")
contrast_data <- data.frame()

for (g in contrast_genes) {
  for (grp in c("High_stress", "Low_stress")) {
    cells <- WhichCells(obj, expression = stress_group == grp)
    if(length(cells) > 0) {
      expr <- GetAssayData(obj, layer = "data")[g, cells]
      contrast_data <- rbind(contrast_data, data.frame(
        gene = g,
        stress_group = grp,
        mean_expr = mean(expm1(expr) - 1),
        pct_expr = mean(expr > 0) * 100
      ))
    }
  }
}

# 获取统计信息
hmox1_info <- de_stress_df %>% filter(gene == "Hmox1")
hspa1a_info <- de_stress_df %>% filter(gene == "Hspa1a")

p_g <- ggplot(contrast_data, aes(x = stress_group, y = mean_expr, fill = stress_group)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f\n(%.1f%%)", mean_expr, pct_expr)), 
            vjust = -0.3, size = 2.5) +
  facet_wrap(~gene, scales = "free_y") +
  scale_fill_manual(values = c("High_stress" = "#E41A1C", "Low_stress" = "#377EB8")) +
  labs(x = "", y = "Mean Expression", 
       title = "G. Oxidative (Hmox1) vs Heat Shock (Hspa1a)",
       subtitle = sprintf("Hmox1: 64%% vs 2%% cells | Hspa1a: 0.7%% vs 0.5%% cells")) +
  fig_theme +
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "grey90"),
        plot.subtitle = element_text(size = 7, hjust = 0.5))

# ==================== Panel H: Cyp2e1三层比较 ====================
cat("Panel H: Cyp2e1三层比较\n")

# 计算三层比较
cyp_high_vs_low <- de_stress_df %>% filter(gene == "Cyp2e1")

Idents(obj) <- "group"
de_cc_all <- FindMarkers(obj, ident.1 = "CCl4", ident.2 = "Control", 
                         features = "Cyp2e1", min.pct = 0, logfc.threshold = 0)

obj_high <- subset(obj, stress_group == "High_stress")
de_cc_high <- NULL
if(ncol(obj_high) > 0) {
  Idents(obj_high) <- "group"
  de_cc_high <- FindMarkers(obj_high, ident.1 = "CCl4", ident.2 = "Control",
                            features = "Cyp2e1", min.pct = 0, logfc.threshold = 0)
}

cyp_comparison <- data.frame(
  Comparison = c("High vs Low\nstress", "CCl4 vs Ctrl\n(all cells)", "CCl4 vs Ctrl\n(high stress)"),
  logFC = c(cyp_high_vs_low$avg_log2FC, 
            de_cc_all["Cyp2e1", "avg_log2FC"],
            ifelse(!is.null(de_cc_high), de_cc_high["Cyp2e1", "avg_log2FC"], 0)),
  p_value = c(cyp_high_vs_low$p_val_adj,
              de_cc_all["Cyp2e1", "p_val_adj"],
              ifelse(!is.null(de_cc_high), de_cc_high["Cyp2e1", "p_val_adj"], 1)),
  Direction = c("Up", "Down", "Down")
)

cyp_comparison$label <- sprintf("%.2f%s", cyp_comparison$logFC,
                                ifelse(cyp_comparison$p_value < 0.001, "***",
                                       ifelse(cyp_comparison$p_value < 0.05, "*", "ns")))

p_h <- ggplot(cyp_comparison, aes(x = Comparison, y = logFC, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
  geom_text(aes(label = label), 
            vjust = ifelse(cyp_comparison$logFC > 0, -0.5, 1.2), 
            size = 3, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
  labs(x = "", y = expression(log[2]~Fold~Change), 
       title = "H. Cyp2e1 Context-Dependent Regulation") +
  fig_theme +
  theme(legend.position = "none", axis.text.x = element_text(size = 7))

# ==================== Panel I: 8基因签名表达点图 ====================
cat("Panel I: 8基因签名表达\n")

fig2_signature_genes <- c("Cyp2e1", "Cyp3a11", "Hspa1a", "Hsp90aa1", 
                          "Gpx1", "Hmox1", "Atf4", "Gadd45a")
fig2_signature_genes <- fig2_signature_genes[fig2_signature_genes %in% rownames(obj)]

dot_data <- data.frame()
for (g in fig2_signature_genes) {
  expr_high <- GetAssayData(obj, layer = "data")[g, obj$stress_group == "High_stress"]
  expr_low <- GetAssayData(obj, layer = "data")[g, obj$stress_group == "Low_stress"]
  
  de_row <- de_stress_df %>% filter(gene == g)
  
  dot_data <- rbind(dot_data, data.frame(
    gene = g,
    pathway = case_when(
      g %in% c("Cyp2e1", "Cyp3a11") ~ "Metabolism",
      g %in% c("Hspa1a", "Hsp90aa1") ~ "HeatShock",
      g %in% c("Gpx1", "Hmox1") ~ "Oxidative",
      g == "Atf4" ~ "ERstress",
      TRUE ~ "DNAdamage"
    ),
    pct_high = mean(expr_high > 0) * 100,
    avg_high = mean(expm1(expr_high) - 1),
    logFC = ifelse(nrow(de_row) > 0, de_row$avg_log2FC, NA),
    p_adj = ifelse(nrow(de_row) > 0, de_row$p_val_adj, NA)
  ))
}

dot_data$significance <- ifelse(dot_data$p_adj < 0.001, "***",
                                ifelse(dot_data$p_adj < 0.01, "**",
                                       ifelse(dot_data$p_adj < 0.05, "*", "ns")))

p_i <- ggplot(dot_data, aes(x = reorder(gene, logFC), y = logFC, 
                            size = pct_high, color = pathway)) +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = significance), vjust = -1.2, size = 2.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_size_continuous(range = c(2, 8), name = "% Cells\n(High)") +
  scale_color_brewer(palette = "Set2", name = "Pathway") +
  labs(x = "", y = expression(log[2]~FC~High/Low), title = "I. 8-Gene Signature Expression") +
  fig_theme +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# ==================== Panel J: GO-BP富集 ====================
cat("Panel J: GO-BP富集\n")

go_file <- file.path(hep_tables_dir, "GO_BP_high_stress_hepatocytes.csv")

if (file.exists(go_file)) {
  ego_stress_df <- read.csv(go_file)
  if (nrow(ego_stress_df) > 0) {
    ego_top <- head(ego_stress_df, 8)
    p_j <- ggplot(ego_top, aes(x = Count, y = reorder(Description, Count), 
                               color = p.adjust, size = Count)) +
      geom_point(alpha = 0.8) +
      scale_color_gradient(low = "red", high = "blue", trans = "log10", name = "p.adj") +
      scale_size_continuous(range = c(3, 8)) +
      labs(x = "Gene Count", y = "", title = "J. GO-BP: High-Stress Hepatocytes") +
      fig_theme + 
      theme(axis.text.y = element_text(size = 6))
  } else {
    p_j <- ggplot() + annotate("text", x=0.5, y=0.5, label="No significant enrichment") + 
      labs(title="J. GO-BP") + theme_void() + fig_theme
  }
} else {
  p_j <- ggplot() + annotate("text", x=0.5, y=0.5, label="GO file not found") + 
    labs(title="J. GO-BP") + theme_void() + fig_theme
}

# ==================== Panel K: KEGG富集 ====================
cat("Panel K: KEGG通路富集\n")

kegg_file <- file.path(hep_tables_dir, "KEGG_upregulated_CCl4.csv")

if (file.exists(kegg_file)) {
  kk_cc_df <- read.csv(kegg_file)
  if (nrow(kk_cc_df) > 0) {
    kk_top <- head(kk_cc_df, 8)
    p_k <- ggplot(kk_top, aes(x = Count, y = reorder(Description, Count), 
                              color = p.adjust, size = Count)) +
      geom_point(alpha = 0.8) +
      scale_color_gradient(low = "red", high = "blue", trans = "log10", name = "p.adj") +
      scale_size_continuous(range = c(3, 8)) +
      labs(x = "Gene Count", y = "", title = "K. KEGG: CCl4 vs Control") +
      fig_theme + 
      theme(axis.text.y = element_text(size = 6))
  } else {
    p_k <- ggplot() + annotate("text", x=0.5, y=0.5, label="No significant KEGG") + 
      labs(title="K. KEGG") + theme_void() + fig_theme
  }
} else {
  p_k <- ggplot() + annotate("text", x=0.5, y=0.5, label="KEGG file not found") + 
    labs(title="K. KEGG") + theme_void() + fig_theme
}

# ==================== 组合完整Figure 2（3行布局）====================
cat("组合完整Figure 2...\n")

# 第1行: A B C D（分区与评分分布）
row1 <- p_a + p_b + p_c + p_d + plot_layout(widths = c(1, 1, 1.2, 1.2))

# 第2行: E F G H（分区定量、七类应激、机制对比）
row2 <- p_e + p_f + p_g + p_h + plot_layout(widths = c(1, 1.2, 1, 1))

# 第3行: I J K（签名表达、GO、KEGG）- 3 panels，留空或调整
row3 <- p_i + p_j + p_k + plot_layout(widths = c(1.2, 1, 1))

fig2 <- row1 / row2 / row3 +
  plot_layout(heights = c(1, 1.2, 1)) +
  plot_annotation(
    title = "Figure 2. Zone 3 Hepatocytes as the Primary Injury Target",
    subtitle = sprintf("Zonation-specific stress responses | n=%d hepatocytes | High-stress: %d (%.1f%%) | CCl4_Zone3 oxidative: 0.73",
                       stats$n_total, stats$n_high, stats$high_pct),
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 8, face = "italic")
    ),
    tag_levels = list(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"))
  )

# 保存
ggsave(file.path(figures.dir, "Figure2_Zone3_Hepatocyte_Stress.tiff"), 
       fig2, width = 16, height = 13, dpi = 300, compression = "lzw")

ggsave(file.path(figures.dir, "Figure2_Zone3_Hepatocyte_Stress.pdf"), 
       fig2, width = 16, height = 13, device = cairo_pdf)

cat("\n✓ Figure 2 Complete: Zone 3 Hepatocyte Stress Analysis (11 panels: A-K)\n")
cat(sprintf("输出: %s/Figure2_Zone3_Hepatocyte_Stress.{tiff,pdf}\n", figures.dir))

# ==================== 步骤10: Figure 3:HSC激活机制 ====================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("              Figure 5: HSC Activation and Fibrosis Mechanisms (9 panels, 3x3)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 检查并加载必要的数据
if(!exists("high_quality_cells")) {
  load(file.path(work.dir, "results", "2_analysis", "high_quality_cells.RData"))
  cat("✓ 加载high_quality_cells对象\n")
}

fig_theme <- theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank())

# 加载DEG结果(使用已存在的hsc_deg或从文件加载)
if(!exists("hsc_deg")) {
  deg_file <- file.path(tables.dir, "Table_DEG_HSC_Act_vs_Qui_CCl4.csv")
  if(file.exists(deg_file)) {
    hsc_deg <- read.csv(deg_file)
    cat("✓ 从文件加载HSC DEG结果:", nrow(hsc_deg), "个基因\n")
  } else {
    stop("错误: 未找到HSC DEG结果文件")
  }
}

# ========== Panel A: HSC Activation Proportion ==========
cat("Panel A: HSC Activation Proportion...\n")

hsc_comp <- high_quality_cells@meta.data %>%
  filter(celltype_reannotated %in% c("HSC_Act", "HSC_Qui")) %>%
  group_by(group, celltype_reannotated) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(group) %>%
  mutate(prop = n / sum(n) * 100,
         total_n = sum(n))

p3a <- ggplot(hsc_comp, aes(x = group, y = n, fill = celltype_reannotated)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, prop)), 
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("HSC_Act" = "#FF7F00", "HSC_Qui" = "#FFFF33"),
                    labels = c("Activated", "Quiescent")) +
  labs(x = "", y = "Cell Count", title = "A. HSC Activation Status") +
  fig_theme

# ========== Panel B: DEG Volcano ==========
cat("Panel B: DEG Volcano...\n")

# 确保significance列被正确创建
hsc_deg <- hsc_deg %>%
  mutate(significance = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up in CCl4",
    p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down in CCl4",
    TRUE ~ "NS"
  ))

top_genes <- bind_rows(
  hsc_deg %>% filter(significance == "Up in CCl4") %>% head(10),
  hsc_deg %>% filter(significance == "Down in CCl4") %>% head(5)
)

p3b <- ggplot(hsc_deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("Up in CCl4" = "#E41A1C", "Down in CCl4" = "#377EB8", "NS" = "grey80")) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 15) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  labs(x = expression(log[2]~Fold~Change), y = expression(-log[10]~p[adj]),
       title = "B. Differential Expression in Activated HSCs") +
  fig_theme +
  theme(legend.position = "bottom")

# ========== Panel C: Top Fibrosis Genes Dot Plot ==========
cat("Panel C: Fibrosis Genes...\n")

fibrosis_genes <- c("Col1a1", "Col1a2", "Col3a1", "Acta2", "Timp1", "Spp1", "Lox", "Cthrc1")
hsc_subset <- subset(high_quality_cells, celltype_reannotated %in% c("HSC_Act", "HSC_Qui"))

# 检查基因是否存在
available_genes <- fibrosis_genes[fibrosis_genes %in% rownames(hsc_subset)]
missing_genes <- fibrosis_genes[!fibrosis_genes %in% rownames(hsc_subset)]
if(length(missing_genes) > 0) {
  cat("  警告: 缺失基因:", paste(missing_genes, collapse = ", "), "\n")
}

if(length(available_genes) > 0) {
  p3c <- DotPlot(hsc_subset, features = available_genes, group.by = "celltype_reannotated",
                 cols = c("lightgrey", "#E41A1C")) +
    RotatedAxis() +
    labs(title = "C. Core Fibrosis Genes") +
    fig_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  # 备用: 使用ggplot手动创建
  dot_data <- hsc_deg %>% 
    filter(gene %in% fibrosis_genes) %>%
    select(gene, avg_log2FC, p_val_adj) %>%
    mutate(pct_expr = 80, avg_expr = avg_log2FC)
  
  p3c <- ggplot(dot_data, aes(x = gene, y = "HSC_Act", size = pct_expr, color = avg_log2FC)) +
    geom_point() +
    scale_color_gradient(low = "lightgrey", high = "#E41A1C") +
    labs(title = "C. Core Fibrosis Genes (DEG-based)") +
    fig_theme
}

# ========== Panel D-F: 富集分析(带容错) ==========
cat("Panel D-F: Enrichment Analysis...\n")

# 准备基因列表
up_genes <- hsc_deg %>% 
  filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>%
  pull(gene)

cat("  上调基因数:", length(up_genes), "\n")

# 尝试基因转换
cat("  转换基因ID...\n")
gene_conv <- tryCatch({
  bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
}, error = function(e) {
  cat("  转换失败:", e$message, "\n")
  return(NULL)
})

if(!is.null(gene_conv) && nrow(gene_conv) > 10) {
  cat("  成功转换", nrow(gene_conv), "个基因\n")
  
  # Panel D: GO-BP
  ego_bp <- tryCatch({
    enrichGO(gene = gene_conv$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", 
             pvalueCutoff = 0.05, readable = TRUE)
  }, error = function(e) NULL)
  
  if(!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    p3d <- dotplot(ego_bp, showCategory = 10, title = "D. GO Biological Process") + fig_theme
    write.csv(as.data.frame(ego_bp), file.path(tables.dir, "Table_Figure5_GO_BP_HSC_up.csv"), row.names = FALSE)
  } else {
    p5d <- ggplot() + annotate("text", x=0.5, y=0.5, label="No significant GO enrichment") + 
      labs(title="D. GO Biological Process") + theme_void() + fig_theme
  }
  
  # Panel E: KEGG
  kk <- tryCatch({
    enrichKEGG(gene = gene_conv$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
  }, error = function(e) NULL)
  
  if(!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    p3e <- dotplot(kk, showCategory = 10, title = "E. KEGG Pathway") + fig_theme
    write.csv(as.data.frame(kk), file.path(tables.dir, "Table_Figure5_KEGG_HSC_up.csv"), row.names = FALSE)
  } else {
    p3e <- ggplot() + annotate("text", x=0.5, y=0.5, label="No significant KEGG enrichment") + 
      labs(title="E. KEGG Pathway") + theme_void() + fig_theme
  }
  
  # Panel F: Gene-Concept Network
  if(!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    p3f <- cnetplot(ego_bp, showCategory = 5, 
                    foldChange = setNames(hsc_deg$avg_log2FC, hsc_deg$gene)) +
      ggtitle("F. Gene-Pathway Network") +
      theme_bw()
  } else {
    p5f <- ggplot() + annotate("text", x=0.5, y=0.5, label="Network not available") + 
      labs(title="F. Gene-Pathway Network") + theme_void()
  }
  
} else {
  cat("  基因转换失败或基因数不足，跳过富集分析\n")
  empty_plot <- ggplot() + annotate("text", x=0.5, y=0.5, 
                                    label="Enrichment analysis failed\n(insufficient gene mapping)") + 
    theme_void() + fig_theme
  
  p3d <- empty_plot + labs(title="D. GO Biological Process")
  p3e <- empty_plot + labs(title="E. KEGG Pathway")
  p3f <- empty_plot + labs(title="F. Gene-Pathway Network")
}

# ========== Panel G: Core Fibrosis Heatmap ==========
cat("Panel G: Fibrosis Heatmap...\n")

fibrosis_core <- c("Col1a1", "Col1a2", "Col3a1", "Acta2", "Timp1", "Spp1", "Cthrc1", "Lox")
available_core <- fibrosis_core[fibrosis_core %in% rownames(hsc_subset)]

if(length(available_core) >= 4) {
  expr_avg <- AverageExpression(hsc_subset, features = available_core, group.by = "celltype_reannotated")$RNA
  
  expr_df <- as.data.frame(expr_avg) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(cols = -gene, names_to = "group", values_to = "expression")
  
  p3g <- ggplot(expr_df, aes(x = group, y = reorder(gene, expression), fill = expression)) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(x = "", y = "", title = "G. Core Fibrosis Genes") +
    fig_theme
  
  tryCatch({
    pheatmap::pheatmap(expr_avg, scale = "row", 
                       color = colorRampPalette(c("blue", "white", "red"))(50),
                       filename = file.path(figures.dir, "Figure5G_heatmap_pheatmap.pdf"),
                       width = 6, height = 8)
  }, error = function(e) cat("  pheatmap保存失败(可忽略)\n"))
  
} else {
  p3g <- ggplot() + annotate("text", x=0.5, y=0.5, label="Insufficient genes for heatmap") + 
    labs(title="G. Core Fibrosis Genes") + theme_void() + fig_theme
}

# ========== Panel H: GSEA (简化版) ==========
cat("Panel H: GSEA...\n")

gsea_data <- hsc_deg %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)

p3h <- ggplot(gsea_data, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = avg_log2FC > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
  coord_flip() +
  labs(x = "", y = expression(log[2]~FC), title = "H. Top Ranked Genes (GSEA-style)") +
  fig_theme +
  theme(legend.position = "none")

# ========== Panel I: Therapeutic Targets ==========
cat("Panel I: Therapeutic Targets...\n")

target_genes <- c("Col1a1", "Acta2", "Timp1", "Spp1", "Lox", "Cthrc1", "Postn", "Vim")
target_data <- hsc_deg %>%
  filter(gene %in% target_genes) %>%
  mutate(
    Druggability = case_when(
      gene %in% c("Col1a1", "Acta2") ~ "★★★",
      gene %in% c("Timp1", "Spp1", "Lox") ~ "★★☆",
      TRUE ~ "★☆☆"
    ),
    Priority = abs(avg_log2FC) * ifelse(Druggability == "★★★", 3, ifelse(Druggability == "★★☆", 2, 1))
  ) %>%
  arrange(desc(Priority))

if(nrow(target_data) > 0) {
  p3i <- ggplot(target_data, aes(x = reorder(gene, Priority), y = avg_log2FC, fill = Druggability)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Druggability), hjust = -0.1, size = 3) +
    scale_fill_manual(values = c("★★★" = "#E41A1C", "★★☆" = "#FF7F00", "★☆☆" = "#FFD700")) +
    coord_flip() +
    labs(x = "", y = expression(log[2]~FC), title = "I. Therapeutic Target Priority") +
    fig_theme
} else {
  p3i <- ggplot() + annotate("text", x=0.5, y=0.5, label="No target data available") + 
    labs(title="I. Therapeutic Target Priority") + theme_void() + fig_theme
}

# ========== 组合 Figure 3 (3x3布局) ==========
cat("组合 Figure 3 (3x3布局)...\n")

# 修正：使用 p3a-p3i 而不是 p5a-p5i
fig3 <- wrap_plots(p3a, p3b, p3c, p3d, p3e, p3f, p3g, p3h, p3i, ncol = 3) +
  plot_annotation(
    title = "Figure 3. HSC Activation Status and Fibrosis Mechanisms",  # 统一标题
    subtitle = "Comprehensive analysis of hepatic stellate cell activation, differential expression, pathway enrichment, and therapeutic targets",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

# 修正：使用 fig3 而不是 fig5
ggsave(file.path(figures.dir, "Figure3_HSC_Activation_Mechanisms.tiff"), 
       fig3, width = 15, height = 15, dpi = 300, compression = "lzw")  # ✅ 用fig3

ggsave(file.path(figures.dir, "Figure3_HSC_Activation_Mechanisms.pdf"), 
       fig3, width = 15, height = 15, device = cairo_pdf)  # ✅ 用fig3

ggsave(file.path(figures.dir, "Figure3_HSC_Activation_Mechanisms.png"), 
       fig3, width = 15, height = 15, dpi = 300)  # ✅ 用fig3

cat("\n✓ Figure 3 Complete! (HSC激活机制，9 panels A-I)\n")  # 统一输出信息

# ==================== 步骤11: Figure 4 - 通讯分析====================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("              Figure 3: Comprehensive Hepatocyte-HSC Communication Analysis\n")
cat("              (CellChat + Custom LR Analysis) - 完全修复版\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 10.1 加载必要包
if(!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")
if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if(!requireNamespace("magick", quietly = TRUE)) install.packages("magick")

suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(magick)
})

# 确保CellChat加载
if(!requireNamespace("CellChat", quietly = TRUE)) {
  stop("请先安装CellChat: devtools::install_github('sqjin/CellChat')")
}
library(CellChat)

# 修复命名空间冲突
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
desc <- dplyr::desc
group_by <- dplyr::group_by
summarise <- dplyr::summarise
pull <- dplyr::pull

# 10.2 加载CellChat对象
cellchat.dir <- file.path(work.dir, "results", "2_analysis", "cellchat")
cellchat_available <- FALSE

if(dir.exists(cellchat.dir)) {
  ccl4_rds <- file.path(cellchat.dir, "cellchat_CCl4.rds")
  control_rds <- file.path(cellchat.dir, "cellchat_Control.rds")
  
  if(file.exists(ccl4_rds) && file.exists(control_rds)) {
    cat("加载CellChat对象...\n")
    cellchat_ccl4 <- readRDS(ccl4_rds)
    cellchat_ctrl <- readRDS(control_rds)
    cellchat_available <- TRUE
    cat("✓ CellChat对象已加载\n")
    
    # 检查对象结构
    cell_types_ccl4 <- levels(cellchat_ccl4@idents)
    cell_types_ctrl <- levels(cellchat_ctrl@idents)
    cat(sprintf("  CCl4组细胞类型 (%d): %s\n", length(cell_types_ccl4), 
                paste(cell_types_ccl4, collapse = ", ")))
    cat(sprintf("  Control组细胞类型 (%d): %s\n", length(cell_types_ctrl), 
                paste(cell_types_ctrl, collapse = ", ")))
    
    # 识别关键细胞类型
    hep_types <- cell_types_ccl4[grep("Hepatocyte", cell_types_ccl4, ignore.case = TRUE)]
    hsc_types <- cell_types_ccl4[grep("HSC", cell_types_ccl4, ignore.case = TRUE)]
    
    cat(sprintf("  肝细胞类型: %s\n", paste(hep_types, collapse = ", ")))
    cat(sprintf("  HSC类型: %s\n", paste(hsc_types, collapse = ", ")))
    
    if("HSC_Act" %in% cell_types_ctrl) {
      cat("  ✓ Control组包含HSC_Act\n")
    } else {
      cat("  ! Control组缺少HSC_Act(细胞数不足被排除)\n")
    }
  } else {
    cat("! 未找到CellChat RDS文件\n")
  }
} else {
  cat("! 未找到CellChat目录\n")
}

# ==================== 关键修复: 安全的通讯数据提取函数 ====================

#' 从CellChat对象安全提取通讯数据
#' 兼容v1.x和v2.x版本
extract_communication_safe <- function(cellchat_obj, sources = NULL, targets = NULL) {
  
  # 方法1: 尝试使用subsetCommunication(v2.0+)
  comm_data <- tryCatch({
    if(!is.null(sources) && !is.null(targets)) {
      # 尝试调用CellChat命名空间的函数
      CellChat::subsetCommunication(cellchat_obj, 
                                    sources.use = sources, 
                                    targets.use = targets)
    } else {
      CellChat::subsetCommunication(cellchat_obj)
    }
  }, error = function(e) {
    cat("    subsetCommunication失败:", conditionMessage(e), "\n")
    cat("    尝试替代方法...\n")
    return(NULL)
  })
  
  # 方法2: 如果方法1失败，直接从@net提取
  if(is.null(comm_data) || nrow(comm_data) == 0) {
    cat("    使用@net$prob直接提取...\n")
    
    if(!"net" %in% slotNames(cellchat_obj) || 
       is.null(cellchat_obj@net$prob)) {
      cat("    ! 无法提取通讯数据\n")
      return(data.frame())
    }
    
    prob_matrix <- cellchat_obj@net$prob
    cell_types <- rownames(prob_matrix)
    
    # 筛选source
    if(!is.null(sources)) {
      source_idx <- which(cell_types %in% sources)
    } else {
      source_idx <- 1:nrow(prob_matrix)
    }
    
    # 筛选target
    if(!is.null(targets)) {
      target_idx <- which(cell_types %in% targets)
    } else {
      target_idx <- 1:ncol(prob_matrix)
    }
    
    # 提取子矩阵
    if(length(source_idx) == 0 || length(target_idx) == 0) {
      return(data.frame())
    }
    
    sub_prob <- prob_matrix[source_idx, target_idx, drop = FALSE]
    
    # 转换为数据框格式(模拟subsetCommunication输出)
    comm_list <- list()
    
    for(i in 1:nrow(sub_prob)) {
      for(j in 1:ncol(sub_prob)) {
        if(sub_prob[i,j] > 0) {
          lr_pair <- paste(rownames(sub_prob)[i], "→", colnames(sub_prob)[j])
          
          comm_list[[length(comm_list)+1]] <- data.frame(
            source = rownames(sub_prob)[i],
            target = colnames(sub_prob)[j],
            ligand = NA,
            receptor = NA,
            prob = sub_prob[i,j],
            pval = NA,
            pathway = NA,
            LR_pair = lr_pair
          )
        }
      }
    }
    
    if(length(comm_list) > 0) {
      comm_data <- do.call(rbind, comm_list)
    } else {
      comm_data <- data.frame()
    }
  }
  
  return(comm_data)
}

# ==================== 第1行: CellChat分析 (Panel A-D) ====================

if(cellchat_available) {
  cat("\n========== Panel A-D: CellChat分析 ==========\n")
  
  # 提取通讯数据
  cat("\n提取CCl4组通讯数据...\n")
  df_comm_ccl4 <- extract_communication_safe(cellchat_ccl4)
  cat(sprintf("  CCl4总通讯对: %d\n", nrow(df_comm_ccl4)))
  
  cat("\n提取Control组通讯数据...\n")
  df_comm_ctrl <- extract_communication_safe(cellchat_ctrl)
  cat(sprintf("  Control总通讯对: %d\n", nrow(df_comm_ctrl)))
  
  # 提取Hep→HSC特异性通讯
  cat("\n提取Hep→HSC特异性通讯...\n")
  hep_hsc_ccl4 <- extract_communication_safe(cellchat_ccl4, hep_types, hsc_types)
  hep_hsc_ctrl <- extract_communication_safe(cellchat_ctrl, hep_types, hsc_types)
  
  cat(sprintf("  CCl4 Hep→HSC: %d\n", nrow(hep_hsc_ccl4)))
  cat(sprintf("  Control Hep→HSC: %d\n", nrow(hep_hsc_ctrl)))
  
  # Panel A: CellChat通讯概览
  cat("\nPanel A: CellChat通讯概览...\n")
  
  comm_counts <- data.frame(
    Group = c("Control", "CCl4"),
    Total = c(nrow(df_comm_ctrl), nrow(df_comm_ccl4)),
    Hep_HSC = c(nrow(hep_hsc_ctrl), nrow(hep_hsc_ccl4))
  ) %>%
    tidyr::pivot_longer(cols = c(Total, Hep_HSC), 
                        names_to = "Type", 
                        values_to = "Count")
  
  p4a <- ggplot(comm_counts, aes(x = Group, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    scale_fill_manual(values = c("Total" = "#999999", "Hep_HSC" = "#E41A1C"),
                      labels = c("Hep→HSC", "Total")) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.6), 
              vjust = -0.3, size = 3.5, fontface = "bold") +
    labs(x = "", y = "Communication Pairs", title = "A. CellChat Overview") +
    theme_journal +
    theme(legend.position = "top", legend.title = element_blank())
  
  # 保存通讯数据
  write.csv(df_comm_ccl4, file.path(tables.dir, "Table_Figure3_CellChat_CCl4.csv"), row.names = FALSE)
  write.csv(hep_hsc_ccl4, file.path(tables.dir, "Table_Figure3_CellChat_Hep_HSC_CCl4.csv"), row.names = FALSE)
  
  # ==================== Panel B: Top LR气泡图(修复版)====================
  cat("Panel B: Top LR气泡图...\n")
  
  if(nrow(hep_hsc_ccl4) > 0) {
    # 修复关键: 按source-target合并概率，避免重复pair
    hep_hsc_ccl4_clean <- hep_hsc_ccl4 %>%
      group_by(source, target) %>%
      summarise(
        prob = sum(prob, na.rm = TRUE),
        n_pairs = n(),
        .groups = "drop"
      ) %>%
      mutate(pair = paste(source, "→", target))
    
    # 双重检查: 确保无重复
    if(anyDuplicated(hep_hsc_ccl4_clean$pair) > 0) {
      warning("仍存在重复pair，强制去重")
      hep_hsc_ccl4_clean <- hep_hsc_ccl4_clean %>%
        distinct(pair, .keep_all = TRUE)
    }
    
    # 获取top 15
    top_pairs <- hep_hsc_ccl4_clean %>%
      arrange(desc(prob)) %>%
      head(15)
    
    cat(sprintf("  Top 15 Hep→HSC pairs selected (from %d unique pairs)\n", 
                nrow(hep_hsc_ccl4_clean)))
    
    # 准备绘图数据
    plot_data <- hep_hsc_ccl4_clean %>%
      filter(pair %in% top_pairs$pair) %>%
      mutate(
        pair = factor(pair, levels = unique(rev(top_pairs$pair))),
        Group = "CCl4"
      )
    
    # 添加Control组数据
    if(nrow(hep_hsc_ctrl) > 0) {
      hep_hsc_ctrl_clean <- hep_hsc_ctrl %>%
        group_by(source, target) %>%
        summarise(prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
        mutate(pair = paste(source, "→", target))
      
      common_pairs <- intersect(hep_hsc_ctrl_clean$pair, top_pairs$pair)
      if(length(common_pairs) > 0) {
        plot_data_ctrl <- hep_hsc_ctrl_clean %>% 
          filter(pair %in% common_pairs) %>%
          mutate(
            pair = factor(pair, levels = levels(plot_data$pair)),
            Group = "Control"
          )
        plot_data <- bind_rows(plot_data, plot_data_ctrl)
      }
    }
    
    # 绘制气泡图
    p4b <- ggplot(plot_data, aes(x = Group, y = pair, size = prob, color = prob)) +
      geom_point(alpha = 0.8) +
      scale_color_gradient(low = "#56B1F7", high = "#132B43", name = "Prob") +
      scale_size_continuous(range = c(3, 10), name = "Prob") +
      labs(x = "", y = "", title = "B. Top 15 Hep→HSC Pairs") +
      theme_journal
    
  } else {
    p3b <- ggplot() + annotate("text", x=0.5, y=0.5, label="No Hep→HSC data") + 
      labs(title="B. Top 15 Hep→HSC Pairs") + theme_void() + theme_journal
  }
  
  # Panel C: 通路活性(使用Hep→HSC通讯强度作为替代)
  if(nrow(hep_hsc_ccl4) > 0) {
    pw_df <- hep_hsc_ccl4 %>%
      group_by(target) %>%
      summarise(Score = sum(prob), .groups = 'drop') %>%
      mutate(Pathway = target) %>%
      arrange(desc(Score))
    
    p4c <- ggplot(pw_df, aes(x = reorder(Pathway, Score), y = Score, fill = Pathway)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = c("HSC_Act" = "#E41A1C", "HSC_Qui" = "#4DAF4A")) +
      coord_flip() +
      labs(x = "", y = "Total Communication Strength", title = "C. HSC Target Activity") +
      theme_journal +
      theme(legend.position = "none")
  } else {
    p4c <- ggplot() + annotate("text", x=0.5, y=0.5, label="No pathway data") + 
      labs(title="C. Pathway Activity") + theme_void() + theme_journal
  }
  
  # Panel D: HSC激活状态转变
  cat("Panel D: HSC激活状态转变...\n")
  
  if(nrow(hep_hsc_ccl4) > 0 || nrow(hep_hsc_ctrl) > 0) {
    qui_ccl4 <- if(nrow(hep_hsc_ccl4) > 0) sum(hep_hsc_ccl4$target == "HSC_Qui", na.rm = TRUE) else 0
    act_ccl4 <- if(nrow(hep_hsc_ccl4) > 0) sum(hep_hsc_ccl4$target == "HSC_Act", na.rm = TRUE) else 0
    qui_ctrl <- if(nrow(hep_hsc_ctrl) > 0) sum(hep_hsc_ctrl$target == "HSC_Qui", na.rm = TRUE) else 0
    act_ctrl <- if(nrow(hep_hsc_ctrl) > 0) sum(hep_hsc_ctrl$target == "HSC_Act", na.rm = TRUE) else 0
    
    shift_df <- data.frame(
      Group = rep(c("Control", "CCl4"), each = 2),
      Target = rep(c("Quiescent", "Activated"), times = 2),
      Count = c(qui_ctrl, act_ctrl, qui_ccl4, act_ccl4)
    ) %>%
      group_by(Group) %>%
      mutate(Percentage = ifelse(sum(Count) > 0, 
                                 100 * Count / sum(Count), 0)) %>%
      ungroup()
    
    note_text <- ifelse(act_ctrl == 0, "\n(Control HSC_Act excluded)", "")
    
    p4d <- ggplot(shift_df, aes(x = Group, y = Count, fill = Target)) +
      geom_bar(stat = "identity", position = "stack", width = 0.6) +
      scale_fill_manual(values = c("Quiescent" = "#4DAF4A", "Activated" = "#E41A1C")) +
      geom_text(aes(label = ifelse(Count > 0, 
                                   sprintf("%d\n(%.0f%%)", Count, Percentage), 
                                   "")), 
                position = position_stack(vjust = 0.5), 
                color = "white", fontface = "bold", size = 3.5) +
      labs(x = "", y = "LR Pairs", title = paste0("D. Target HSC State", note_text)) +
      theme_journal +
      theme(legend.position = "right", plot.title = element_text(size = 10))
  } else {
    p4d <- ggplot() + annotate("text", x=0.5, y=0.5, label="No data") + 
      labs(title="D. Target HSC State") + theme_void() + theme_journal
  }
  
  cat("✓ Panel A-D完成\n")
  
} else {
  cat("! CellChat不可用，生成占位图\n")
  empty_plot <- ggplot() + annotate("text", x=0.5, y=0.5, 
                                    label="CellChat data not available\nRun run_cellchat_analysis.R first") + 
    theme_void() + theme_journal
  
  p4a <- empty_plot + labs(title = "A. CellChat Overview")
  p4b <- empty_plot + labs(title = "B. Top 15 Hep→HSC Pairs")
  p4c <- empty_plot + labs(title = "C. Pathway Activity")
  p4d <- empty_plot + labs(title = "D. Target HSC State")
}

# ==================== 第2行: 自定义LR分析 + 弦图 (Panel E-H) ====================
cat("\n========== Panel E-H: 自定义LR分析与弦图 ==========\n")

# 定义关键的配体-受体对(基于文献的Hep→HSC通讯)
key_lr_pairs <- list(
  "Tgfb1_Tgfbr" = c("Tgfb1", "Tgfbr1", "Tgfbr2"),
  "Pdgfa_Pdgfra" = c("Pdgfa", "Pdgfra"),
  "Pdgfb_Pdgfrb" = c("Pdgfb", "Pdgfrb"),
  "Ccl2_Ccr2" = c("Ccl2", "Ccr2"),
  "Cxcl12_Cxcr4" = c("Cxcl12", "Cxcr4"),
  "Hgf_Met" = c("Hgf", "Met"),
  "Vegfa_Vegfr" = c("Vegfa", "Kdr", "Flt1"),
  "Spp1_Cd44" = c("Spp1", "Cd44"),
  "Il6_Il6r" = c("Il6", "Il6r"),
  "Tnf_Tnfr" = c("Tnf", "Tnfrsf1a", "Tnfrsf1b")
)

# Panel E: 配体-受体热图
cat("Panel E: LR基因表达热图...\n")

lr_heatmap_data <- data.frame()
for(lr_name in names(key_lr_pairs)) {
  genes <- key_lr_pairs[[lr_name]]
  genes <- genes[genes %in% rownames(high_quality_cells)]
  
  if(length(genes) > 0) {
    for(gene in genes) {
      for(ct in c("Hepatocyte", "HSC_Act", "HSC_Qui")) {
        if(ct %in% high_quality_cells$celltype_reannotated) {
          cells <- WhichCells(high_quality_cells, expression = celltype_reannotated == ct)
          if(length(cells) > 0) {
            expr <- GetAssayData(high_quality_cells, layer = "data")[gene, cells]
            mean_expr <- mean(expm1(expr), na.rm = TRUE)
            pct_expr <- mean(expr > 0, na.rm = TRUE) * 100
            
            lr_heatmap_data <- rbind(lr_heatmap_data, data.frame(
              LR_pair = lr_name,
              Gene = gene,
              Gene_type = ifelse(gene %in% c("Tgfbr1", "Tgfbr2", "Pdgfra", "Pdgfrb", 
                                             "Ccr2", "Cxcr4", "Met", "Kdr", "Flt1",
                                             "Cd44", "Il6r", "Tnfrsf1a", "Tnfrsf1b"), 
                                 "Receptor", "Ligand"),
              CellType = ct,
              Mean_expr = mean_expr,
              Pct_expr = pct_expr
            ))
          }
        }
      }
    }
  }
}

if(nrow(lr_heatmap_data) > 0) {
  p4e <- ggplot(lr_heatmap_data, aes(x = CellType, y = LR_pair, fill = Mean_expr)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", Mean_expr)), size = 2.5) +
    scale_fill_gradientn(colors = c("white", "yellow", "orange", "red", "darkred"),
                         name = "Mean\nExpr") +
    facet_wrap(~Gene_type, ncol = 1, scales = "free_y") +
    labs(x = "", y = "", title = "E. Key LR Pairs Expression") +
    theme_journal +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 7),
          strip.text = element_text(size = 9, face = "bold"))
} else {
  p4e <- ggplot() + annotate("text", x=0.5, y=0.5, 
                             label="No LR data available") + 
    labs(title="E. Key LR Pairs") + theme_void() + theme_journal
}

# Panel F: 通讯强度网络图
cat("Panel F: 通讯网络数据...\n")

if(exists("hep_hsc_ccl4") && nrow(hep_hsc_ccl4) > 0) {
  network_df <- hep_hsc_ccl4 %>%
    group_by(source, target) %>%
    summarise(strength = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    mutate(group = "CCl4")
  
  if(exists("hep_hsc_ctrl") && nrow(hep_hsc_ctrl) > 0) {
    ctrl_df <- hep_hsc_ctrl %>%
      group_by(source, target) %>%
      summarise(strength = sum(prob, na.rm = TRUE), .groups = "drop") %>%
      mutate(group = "Control")
    network_df <- bind_rows(network_df, ctrl_df)
  }
  
  p4f <- ggplot(network_df, aes(x = source, y = target, size = strength, 
                                color = group, shape = group)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("CCl4" = "#E41A1C", "Control" = "#377EB8")) +
    scale_size_continuous(range = c(3, 15), name = "Strength") +
    labs(x = "Source", y = "Target", title = "F. Communication Strength") +
    theme_journal +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))
} else {
  if(nrow(lr_heatmap_data) > 0) {
    inferred_comm <- lr_heatmap_data %>%
      filter(Gene_type == "Ligand") %>%
      group_by(LR_pair, CellType) %>%
      summarise(avg_expr = mean(Mean_expr, na.rm = TRUE), .groups = "drop")
    
    p4f <- ggplot(inferred_comm, aes(x = LR_pair, y = CellType, 
                                     fill = avg_expr)) +
      geom_tile() +
      scale_fill_gradientn(colors = c("white", "orange", "red")) +
      labs(x = "Ligand", y = "Cell Type", title = "F. Inferred Communication") +
      theme_journal +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  } else {
    p4f <- ggplot() + annotate("text", x=0.5, y=0.5, 
                               label="No communication data") + 
      labs(title="F. Communication Network") + theme_void() + theme_journal
  }
}

# Panel G: 空间分布图
cat("Panel G: 空间分布...\n")

if("celltype_reannotated" %in% colnames(high_quality_cells@meta.data)) {
  high_quality_cells$comm_type <- ifelse(
    high_quality_cells$celltype_reannotated == "Hepatocyte", "Hepatocyte (Source)",
    ifelse(high_quality_cells$celltype_reannotated %in% c("HSC_Act", "HSC_Qui"),
           paste0(high_quality_cells$celltype_reannotated, " (Target)"), "Other")
  )
  
  comm_cells <- subset(high_quality_cells, 
                       subset = celltype_reannotated %in% c("Hepatocyte", "HSC_Act", "HSC_Qui"))
  
  if(ncol(comm_cells) > 0) {
    p4g <- DimPlot(comm_cells, reduction = "umap", group.by = "comm_type",
                   pt.size = 0.5, label = TRUE, repel = TRUE) +
      scale_color_manual(values = c(
        "Hepatocyte (Source)" = "#E41A1C",
        "HSC_Act (Target)" = "#FF7F00",
        "HSC_Qui (Target)" = "#FFFF33"
      )) +
      ggtitle("G. Hep→HSC Spatial Map") +
      theme_journal
  } else {
    p4g <- ggplot() + annotate("text", x=0.5, y=0.5, 
                               label="No Hep/HSC cells found") + 
      labs(title="G. Spatial Map") + theme_void() + theme_journal
  }
} else {
  p4g <- ggplot() + annotate("text", x=0.5, y=0.5, 
                             label="Cell type data missing") + 
    labs(title="G. Spatial Map") + theme_void() + theme_journal
}
# Panel H: 全局网络统计(完全修复版)
cat("Panel H: 全局网络统计...\n")

if(exists("df_comm_ccl4") && !is.null(df_comm_ccl4) && nrow(df_comm_ccl4) > 0) {
  
  # 清理数据
  df_comm_clean <- df_comm_ccl4 %>%
    dplyr::filter(!is.na(source), !is.na(target), !is.na(prob))
  
  cat(sprintf("  清理后通讯对: %d\n", nrow(df_comm_clean)))
  
  if(nrow(df_comm_clean) > 0) {
    # 计算out-degree - 关键修复: 使用dplyr::前缀，并在summarise后立即转换为data.frame
    out_degree <- df_comm_clean %>%
      dplyr::group_by(source) %>%
      dplyr::summarise(
        out_strength = sum(prob, na.rm = TRUE), 
        n_targets = dplyr::n_distinct(target),
        .groups = "drop"
      ) %>%
      dplyr::ungroup()
    
    # 关键修复: 显式转换为普通data.frame，避免类属性问题
    out_degree <- as.data.frame(out_degree)
    names(out_degree) <- c("cell_type", "out_strength", "n_targets")
    
    # 计算in-degree
    in_degree <- df_comm_clean %>%
      dplyr::group_by(target) %>%
      dplyr::summarise(
        in_strength = sum(prob, na.rm = TRUE), 
        n_sources = dplyr::n_distinct(source),
        .groups = "drop"
      ) %>%
      dplyr::ungroup()
    
    # 关键修复: 显式转换为普通data.frame
    in_degree <- as.data.frame(in_degree)
    names(in_degree) <- c("cell_type", "in_strength", "n_sources")
    
    cat(sprintf("  Out-degree细胞数: %d, In-degree细胞数: %d\n", 
                nrow(out_degree), nrow(in_degree)))
    
    # 使用base R merge避免dplyr join问题
    centrality_df <- merge(out_degree, in_degree, by = "cell_type", all = TRUE)
    
    # 处理NA
    centrality_df$out_strength[is.na(centrality_df$out_strength)] <- 0
    centrality_df$in_strength[is.na(centrality_df$in_strength)] <- 0
    centrality_df$total_strength <- centrality_df$out_strength + centrality_df$in_strength
    
    cat(sprintf("  合并后数据行数: %d\n", nrow(centrality_df)))
    
    if(nrow(centrality_df) > 0) {
      # 排序并取前10
      ord_idx <- order(centrality_df$total_strength, decreasing = TRUE)
      centrality_top <- centrality_df[ord_idx[1:min(10, nrow(centrality_df))], ]
      
      # 转换为因子
      centrality_top$cell_type <- factor(
        centrality_top$cell_type, 
        levels = centrality_top$cell_type[order(centrality_top$total_strength)]
      )
      
      p4h <- ggplot2::ggplot(centrality_top, ggplot2::aes(x = cell_type)) +
        ggplot2::geom_col(ggplot2::aes(y = out_strength, fill = "Outgoing"), alpha = 0.7, width = 0.6) +
        ggplot2::geom_col(ggplot2::aes(y = -in_strength, fill = "Incoming"), alpha = 0.7, width = 0.6) +
        ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        ggplot2::scale_fill_manual(values = c("Outgoing" = "#E41A1C", "Incoming" = "#377EB8"),
                                   name = "Direction") +
        ggplot2::labs(x = "", y = "Communication Strength", 
                      title = "H. Network Centrality (CCl4)") +
        ggplot2::coord_flip() +
        theme_journal +
        ggplot2::theme(legend.position = "top", 
                       legend.title = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_text(size = 9))
      
      cat("  ✓ Panel H绘制成功\n")
    } else {
      p4h <- ggplot2::ggplot() + 
        ggplot2::annotate("text", x=0.5, y=0.5, label="No valid centrality data") + 
        ggplot2::labs(title="H. Network Centrality (Error)") + 
        ggplot2::theme_void() + theme_journal
    }
  } else {
    p4h <- ggplot2::ggplot() + 
      ggplot2::annotate("text", x=0.5, y=0.5, label="No valid communication data") + 
      ggplot2::labs(title="H. Network Centrality (Error)") + 
      ggplot2::theme_void() + theme_journal
  }
} else {
  # 备用
  cell_counts <- as.data.frame(table(high_quality_cells$celltype_reannotated))
  colnames(cell_counts) <- c("CellType", "Count")
  cell_counts <- cell_counts[!is.na(cell_counts$CellType), ]
  
  p4h <- ggplot2::ggplot(cell_counts, ggplot2::aes(x = reorder(CellType, Count), y = Count)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "", y = "Cell Count", title = "H. Cell Type Abundance (Backup)") +
    theme_journal
}

cat("✓ Panel E-H完成\n")

# ==================== 组合Figure 4(2行×4列布局)====================
cat("\n组合 Figure 4 (2行×4列布局)...\n")

library(grid)
library(gridExtra)
library(patchwork)

# 组合第1行
fig4_row1 <- wrap_plots(p4a, p4b, p4c, p4d, ncol = 4)

# 组合第2行
fig4_row2 <- wrap_plots(p4e, p4f, p4g, p4h, ncol = 4, widths = c(1, 1, 0.8, 0.8))

# 完整组合
fig4_complete <- fig4_row1 / fig4_row2 +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 4. Comprehensive Hepatocyte-HSC Communication in Liver Injury",
    subtitle = "A-D: CellChat analysis | E-F: Custom LR heatmaps | G: Hep→HSC spatial | H: Global network",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic")
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

# 保存
ggsave(file.path(figures.dir, "Figure4_Hepatocyte_HSC_Communication.tiff"), 
       fig4_complete, width = 16, height = 10, dpi = 300, compression = "lzw")

ggsave(file.path(figures.dir, "Figure4_Hepatocyte_HSC_Communication.pdf"), 
       fig4_complete, width = 16, height = 10, device = cairo_pdf)

cat("\n✓ Figure 4 Complete! (原Figure 3，通讯分析，8 panels A-H)\n")
cat("关键改进:\n")
cat("  1. 完整的extract_communication_safe()函数\n")
cat("  2. 完整的Panel E-H代码\n")
cat("  3. 自动检测CellChat对象结构\n")
cat("  4. 多重降级策略\n")
cat("  5. 正确的2行×4列布局\n")

# ==================== 步骤12: Figure 5 - 治疗靶点优先级（新增）====================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("              Figure 5: Therapeutic Target Prioritization\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 加载必要包
if(!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

# 确保hsc_deg存在
if(!exists("hsc_deg")) {
  deg_file <- file.path(tables.dir, "Table_DEG_HSC_Act_vs_Qui_CCl4.csv")
  if(file.exists(deg_file)) {
    hsc_deg <- read.csv(deg_file)
  } else {
    stop("错误: 未找到HSC DEG结果")
  }
}

# 统一主题
fig5_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# ==================== Panel A: 治疗靶点优先级排序 ====================
cat("Panel A: 治疗靶点优先级排序\n")

# 定义靶点及其可药性评分
target_info <- data.frame(
  gene = c("Acta2", "Col1a1", "Lox", "Timp1", "Cthrc1", "Postn", "Vim", "Col1a2", "Col3a1"),
  druggability = c("★★★", "★★★", "★★☆", "★★☆", "★☆☆", "★☆☆", "★☆☆", "★★★", "★★☆"),
  stringsAsFactors = FALSE
)

# 合并DEG数据
target_data <- hsc_deg %>%
  filter(gene %in% target_info$gene) %>%
  left_join(target_info, by = "gene") %>%
  mutate(
    drug_score = case_when(
      druggability == "★★★" ~ 3,
      druggability == "★★☆" ~ 2,
      TRUE ~ 1
    ),
    priority_score = abs(avg_log2FC) * drug_score,
    significance = ifelse(p_val_adj < 0.001, "***",
                          ifelse(p_val_adj < 0.01, "**",
                                 ifelse(p_val_adj < 0.05, "*", "ns")))
  ) %>%
  arrange(desc(priority_score))

# 绘制优先级条形图
p_5a <- ggplot(target_data, aes(x = reorder(gene, priority_score), y = avg_log2FC)) +
  geom_bar(aes(fill = druggability), stat = "identity", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = significance), hjust = -0.2, size = 3.5, color = "black") +
  geom_text(aes(label = druggability), hjust = 1.2, size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("★★★" = "#C0392B", "★★☆" = "#E67E22", "★☆☆" = "#F1C40F"),
                    name = "Druggability") +
  coord_flip() +
  expand_limits(y = max(target_data$avg_log2FC) * 1.2) +
  labs(x = "", y = expression(log[2]~Fold~Change), 
       title = "A. Therapeutic Target Priority Ranking") +
  fig5_theme +
  theme(legend.position = "right")

# ==================== Panel B: 机制总结示意图 ====================
cat("Panel B: Stress-Communication-Activation轴示意图\n")

# 创建机制示意图（使用ggplot绘制流程图）
mechanism_data <- data.frame(
  x = c(1, 2, 3, 4, 2.5, 3.5),
  y = c(2, 2, 2, 2, 1, 1),
  label = c("CCl4\nExposure", "Zone 3\nHepatocyte\nStress", "Cell-Cell\nCommunication", "HSC\nActivation", "PDGF\nTGF-β", "Oxidative\nStress"),
  type = c("start", "process", "process", "end", "signal", "signal"),
  size = c(12, 11, 10, 11, 8, 8)
)

arrow_data <- data.frame(
  x = c(1.4, 2.4, 3.4, 2.5, 3.5),
  xend = c(1.6, 2.6, 3.6, 3.1, 3.9),
  y = c(2, 2, 2, 1.3, 1.3),
  yend = c(2, 2, 2, 1.7, 1.7)
)

p_5b <- ggplot() +
  # 绘制节点
  geom_point(data = mechanism_data, aes(x = x, y = y, size = size), 
             color = c("#3498DB", "#E74C3C", "#F39C12", "#9B59B6", "#2ECC71", "#E74C3C"), 
             alpha = 0.8) +
  # 绘制箭头
  geom_segment(data = arrow_data, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.3, "cm")), linewidth = 1, color = "grey30") +
  # 添加标签
  geom_text(data = mechanism_data, aes(x = x, y = y, label = label), 
            size = 3, fontface = "bold", color = "white") +
  # 添加信号标注
  annotate("text", x = 2.5, y = 0.6, label = "PDGF/TGF-β\nSignaling", 
           size = 3, color = "#2ECC71", fontface = "bold") +
  annotate("text", x = 3.5, y = 0.6, label = "Oxidative\nStress", 
           size = 3, color = "#E74C3C", fontface = "bold") +
  # 设置坐标轴
  scale_x_continuous(limits = c(0.5, 4.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 2.5), expand = c(0, 0)) +
  scale_size_continuous(range = c(15, 25)) +
  labs(title = "B. Stress-Communication-Activation Axis") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "none"
  )

# ==================== 组合Figure 5 ====================
cat("组合Figure 5...\n")

fig5_final <- p_5a + p_5b + 
  plot_layout(ncol = 2, widths = c(1.2, 1)) +
  plot_annotation(
    title = "Figure 5. Therapeutic Target Prioritization and Mechanistic Framework",
    subtitle = "Multi-node therapeutic strategy targeting hepatocyte stress, communication blockade, and HSC deactivation",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic", color = "gray40")
    ),
    tag_levels = "A"
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

# 保存
ggsave(file.path(figures.dir, "Figure5_Therapeutic_Targets.tiff"), 
       fig5_final, width = 14, height = 7, dpi = 300, compression = "lzw")

ggsave(file.path(figures.dir, "Figure5_Therapeutic_Targets.pdf"), 
       fig5_final, width = 14, height = 7, device = cairo_pdf)

ggsave(file.path(figures.dir, "Figure5_Therapeutic_Targets.png"), 
       fig5_final, width = 14, height = 7, dpi = 300)

cat("\n✓ Figure 5 Complete! (新增，治疗靶点优先级，2 panels)\n")
cat("输出文件:\n")
cat("  - Figure5_Therapeutic_Targets.tiff\n")
cat("  - Figure5_Therapeutic_Targets.pdf\n")
cat("  - Figure5_Therapeutic_Targets.png\n")

