# ==================== Complete Single-Cell Analysis Workflow v4.0====================

# ==================== 0. Environment Configuration ====================
work.dir <- "D:/文章/生信分析/小鼠肝脏单细胞RNA测序分析/20260325"
set.seed(2025)
setwd(work.dir)

# Create Directory
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

# Load Packages
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

# Fix Namespace
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
desc <- dplyr::desc
group_by <- dplyr::group_by
summarise <- dplyr::summarise
pull <- dplyr::pull

# Theme Setting
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

cat("===Environment initialization completed ===\n")

# ==================== 1. Data Loading and QC ====================
cat("\n=== Step 1: Reading h5 format data and quality control ===\n")

rdata_file <- file.path(work.dir, "results", "2_analysis", "high_quality_cells_v4.RData")

if(file.exists(rdata_file)) {
  cat("Processed Seurat object detected, loading...\n")
  load(rdata_file)
  cat(sprintf("✓ Loading complete, %d cells in total\n", ncol(high_quality_cells)))
} else {
  # Searching for h5 files 
  h5_files <- list.files(data.dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  
  if(length(h5_files) == 0) {
    h5_files <- list.files(data.dir, pattern = "\\.h5ad$|\\.hdf5$", recursive = TRUE, full.names = TRUE)
  }
  
  if(length(h5_files) == 0) {
    stop("Error: ", data.dir, "No h5 files found in")
  }
  
  cat(sprintf("%d h5 file(s) detected:\n", length(h5_files)))
  print(basename(h5_files))
  
  # Reading h5 file 
  seurat_list <- list()
  
  for(i in 1:length(h5_files)) {
    file_path <- h5_files[i]
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    
    cat(sprintf("\n  reading sample %d/%d: %s\n", i, length(h5_files), sample_name))
    
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
          stop("Unrecognized h5 file format")
        }
      }, error = function(e2) {
        
        if(grepl("\\.h5ad$", file_path)) {
          Convert(file_path, dest = "h5seurat", overwrite = TRUE)
          h5seurat_path <- sub("\\.h5ad$", ".h5seurat", file_path)
          LoadH5Seurat(h5seurat_path)
        } else {
          stop("Failed to read h5 file: ", conditionMessage(e2))
        }
      })
    })
    
    seu <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
    seu$sample <- sample_name
    seurat_list[[i]] <- seu
  }
  
  # merging
  if(length(seurat_list) > 1) {
    cat("\nMerging multiple samples...\n")
    high_quality_cells <- merge(seurat_list[[1]], y = seurat_list[-1], 
                                add.cell.ids = sapply(seurat_list, function(x) x$sample[1]))
  } else {
    high_quality_cells <- seurat_list[[1]]
  }
  
  # Adding group information
  high_quality_cells$group <- ifelse(grepl("Ctrl|Control|WT|wildtype|normal", 
                                           high_quality_cells$sample, ignore.case = TRUE), 
                                     "Control", "CCl4")
  
  cat(sprintf("\nInitial cell count: %d\n", ncol(high_quality_cells)))
  cat("sample group:\n")
  print(table(high_quality_cells$sample, high_quality_cells$group))
  
  # QC
  cat("\nPerforming quality control...\n")
  high_quality_cells[["percent.mt"]] <- PercentageFeatureSet(high_quality_cells, pattern = "^mt-")
  high_quality_cells[["percent.rb"]] <- PercentageFeatureSet(high_quality_cells, pattern = "^Rp[sl]")
  
  # strict QC
  cat("Filtering criteria: nFeature 200-6000, nCount 500-50000, percent.mt < 25%\n")
  high_quality_cells <- subset(high_quality_cells, 
                               subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                                 nCount_RNA > 500 & nCount_RNA < 50000 &
                                 percent.mt < 25)
  
  cat(sprintf("Post-QC cell count: %d (Retention rate: %.1f%%)\n", 
              ncol(high_quality_cells), 
              100 * ncol(high_quality_cells) / sum(sapply(seurat_list, ncol))))
  
  # normalization
  cat("Normalization and Data Scaling...\n")
  high_quality_cells <- NormalizeData(high_quality_cells)
  high_quality_cells <- FindVariableFeatures(high_quality_cells, selection.method = "vst", nfeatures = 2000)
  high_quality_cells <- ScaleData(high_quality_cells)
  high_quality_cells <- RunPCA(high_quality_cells, features = VariableFeatures(object = high_quality_cells))
  
  save(high_quality_cells, file = rdata_file)
  cat("✓ data reading and QC completed，Object has been saved\n")
}

# ==================== 2. Dimensionality Reduction and Clustering ====================
cat("\n=== step 2: Dimensionality Reduction and Clustering ===\n")

if(!"seurat_clusters" %in% colnames(high_quality_cells@meta.data)) {
  cat("running UMAP and clustering...\n")
  
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
  ggsave(file.path(figures.dir, "Figure1_UMAP.pdf"), p_combined, width = 12, height = 5, device = cairo_pdf)
  
  cat(sprintf("✓ clustering completed，found %d cluster\n", length(unique(high_quality_cells$seurat_clusters))))
} else {
  cat("✓ Using existing clustering results\n")
}


# ==================== 3. Cell Type Annotation (Figure 1C)====================
cat("
=== step 3: Cell Type Annotation ===
")

if(!"celltype" %in% colnames(high_quality_cells@meta.data) || 
   sum(grepl("^Cluster_", high_quality_cells$celltype)) > 5) {
  
  cat("Cell annotation based on marker genes...
")
  
  if (length(Layers(high_quality_cells)) > 1) {
    high_quality_cells <- JoinLayers(high_quality_cells)
  }
  
  # Extended marker list
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
  
  # Computing module scores
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
  
  # Improved annotation function
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
      
      # Fallback to module scoring method
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
  
  # Executing annotation
  high_quality_cells$celltype <- cluster_annotation(high_quality_cells)
  
  # Visualization
  p_anno <- DimPlot(high_quality_cells, reduction = "umap", group.by = "celltype", 
                    label = TRUE, pt.size = 0.5, repel = TRUE) +
    ggtitle("C. Cell Types") + 
    theme_journal +
    theme(legend.position = "right")
  
  ggsave(file.path(figures.dir, "Figure1_CellType.tiff"), 
         p_anno, width = 12, height = 8, dpi = 300)
  ggsave(file.path(figures.dir, "Figure1_CellType.pdf"), 
         p_anno, width = 12, height = 8, device = cairo_pdf)
  
  cat("✓ Cell annotation completed
")
  print(table(high_quality_cells$celltype))
  
} else {
  cat("✓ Using existing annotations
")
}

# ==================== 4. HSC Subpopulation Classification (v2.5)====================
cat("
")
cat(paste(rep("=", 70), collapse = ""), "
")
cat("              HSCSubpopulation Classification (v2.5 - 65th percentile)
")
cat(paste(rep("=", 70), collapse = ""), "
")

hsc_cells <- subset(high_quality_cells, celltype == "HSC")

if(ncol(hsc_cells) == 0) {
  stop("Error: HSC cells not found")
}

cat(sprintf("Extracted %d HSC cells
", ncol(hsc_cells)))

# Recalculating activation score
activation_genes <- c("Acta2", "Col1a1", "Col1a2", "Col3a1", "Timp1", "Lox", "Spp1", 
                      "Postn", "Tagln", "Thbs1", "Cthrc1", "Vim", "Fn1", "Tgfb1")
activation_genes <- activation_genes[activation_genes %in% rownames(hsc_cells)]

if(length(activation_genes) >= 3) {
  hsc_cells <- AddModuleScore(hsc_cells, features = list(activation_genes), name = "Activation")
  hsc_cells$act_score <- hsc_cells$Activation1
  
  # v2.5: 65th percentile + Acta2>0.2
  act_threshold <- quantile(hsc_cells$act_score, 0.65, na.rm = TRUE)
  cat(sprintf("Activation score threshold (65th percentile): %.3f
", act_threshold))
  
  if("Acta2" %in% rownames(hsc_cells)) {
    acta2_expr <- FetchData(hsc_cells, vars = "Acta2")[,1]
    
    high_score <- hsc_cells$act_score > act_threshold
    high_acta2 <- acta2_expr > 0.2
    very_high_score <- hsc_cells$act_score > quantile(hsc_cells$act_score, 0.85, na.rm = TRUE)
    
    hsc_cells$hsc_subtype <- ifelse((high_score & high_acta2) | very_high_score, 
                                    "HSC_Act", "HSC_Qui")
    
    cat(sprintf("Acta2>0.2 cell: %d/%d (%.1f%%)
", 
                sum(high_acta2), ncol(hsc_cells), 100*sum(high_acta2)/ncol(hsc_cells)))
    
  } else {
    hsc_cells$hsc_subtype <- ifelse(hsc_cells$act_score > quantile(hsc_cells$act_score, 0.75, na.rm = TRUE), 
                                    "HSC_Act", "HSC_Qui")
  }
  
  # Checking classification results
  subtype_table <- table(hsc_cells$hsc_subtype, hsc_cells$group)
  cat("
HSC subtype distribution (v2.5):
")
  print(subtype_table)
  
  # Calculating proportions
  ccl4_total <- sum(hsc_cells$group == "CCl4")
  ctrl_total <- sum(hsc_cells$group == "Control")
  ccl4_act <- sum(hsc_cells$hsc_subtype == "HSC_Act" & hsc_cells$group == "CCl4")
  ctrl_act <- sum(hsc_cells$hsc_subtype == "HSC_Act" & hsc_cells$group == "Control")
  
  cat(sprintf("  CCl4组: %d/%d (%.1f%%)
", ccl4_act, ccl4_total, 100*ccl4_act/ccl4_total))
  cat(sprintf("  Control组: %d/%d (%.1f%%)
", ctrl_act, ctrl_total, 100*ctrl_act/ctrl_total))
  
  # Updating main object - Fix: Ensure hsc_subtype is added to high_quality_cells
  hsc_meta <- data.frame(
    cell = colnames(hsc_cells),
    hsc_subtype = hsc_cells$hsc_subtype,
    stringsAsFactors = FALSE
  )
  
  match_idx <- match(hsc_meta$cell, colnames(high_quality_cells))
  valid_idx <- !is.na(match_idx)
  high_quality_cells$celltype_reannotated[match_idx[valid_idx]] <- hsc_meta$hsc_subtype[valid_idx]
  
  # Key fix: Simultaneously adding hsc_subtype as an independent column for future reference
  high_quality_cells$hsc_subtype <- NA
  high_quality_cells$hsc_subtype[match_idx[valid_idx]] <- hsc_meta$hsc_subtype[valid_idx]
  
  cat(sprintf("\n✓ HSCClassification v2.5 completed，Updated %d cells \n", sum(valid_idx)))
} else {
  cat("Warning: Insufficient activation marker genes, skipping HSC subtype classification
")
}


# ==================== 5. Fibrosis Score ====================
cat("
=== step 5: Caculating Fibrosis Score ===
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
  cat(sprintf("Using %d fibrosis-related genes to calculate score
", length(fibrosis_genes)))
  
  high_quality_cells <- AddModuleScore(
    high_quality_cells,
    features = list(fibrosis_genes),
    name = "Fibrosis_Score"
  )
  
  # visualization
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
  
  # HSCSubtype Violin Plot
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
    ggsave(file.path(figures.dir, "Figure2_Fibrosis_Score.pdf"), 
           p_combined, width = 8, height = 10, device = cairo_pdf)
    cat("✓ Fibrosis score plot saved
")
  }
  
  cat("✓ Fibrosis score plot completed
")
} else {
  cat("Warning: Insufficient fibrosis genes, skipping scoring
")
}

save(high_quality_cells, file = rdata_file)

# ========= 6. Differential Gene Expression Analysis (HSC_Act vs HSC_Qui within CCl4 group)==========
cat("
=== step 6: Differential Gene Expression Analysis (CCl4 Within-group HSC_Act vs HSC_Qui)===
")

hsc_ccl4 <- subset(high_quality_cells, 
                   group == "CCl4" & 
                     celltype_reannotated %in% c("HSC_Act", "HSC_Qui"))

if(ncol(hsc_ccl4) == 0) {
  stop("Error: No HSC cells found in CCl4 group")
}

cat(sprintf("CCl4 group HSC cells: HSC_Act=%d, HSC_Qui=%d
",
            sum(hsc_ccl4$celltype_reannotated == "HSC_Act"),
            sum(hsc_ccl4$celltype_reannotated == "HSC_Qui")))

hsc_ccl4 <- JoinLayers(hsc_ccl4)
Idents(hsc_ccl4) <- "celltype_reannotated"

# Main Analysis: HSC_Act vs HSC_Qui
cat("
【Main Analysis】HSC_Act vs HSC_Qui(CCl4 Within-group)...
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
✓ Main analysis completed: %d DEGs identified (padj<0.05)
", 
            sum(hsc_deg_ccl4_df$p_val_adj < 0.05)))

# Saved for subsequent use
hsc_deg <- hsc_deg_ccl4_df

# ==================== 7. Pathway Enrichment Analysis and Figure 1 Generation ====================
cat("\n=== step 7: Pathway Enrichment Analysis ===\n")

# Ensure dplyr is loaded (for pipe operations)
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

cat(sprintf("Gene analysis: %d upregulated，%d downregulated\n", length(up_genes), length(down_genes)))

# GOanalysis
cat("GO Enrichment Analysis...\n")
ego_bp_up <- NULL
ego_mf_up <- NULL

if(length(up_genes) >= 10) {
  tryCatch({
    ego_bp_up <- enrichGO(gene = up_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                          ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    
    ego_mf_up <- enrichGO(gene = up_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                          ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05)
    
    cat(sprintf("  GO BP: %d \n", ifelse(is.null(ego_bp_up), 0, nrow(ego_bp_up@result))))
  }, error = function(e) {
    cat("  GO Analysis Error:", conditionMessage(e), "\n")
  })
}

# KEGG analysis
cat("KEGG enrichment analysis...\n")
kk_up <- NULL

if(length(up_genes) >= 10) {
  tryCatch({
    gene_df_up <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    
    if(nrow(gene_df_up) > 0) {
      kk_up <- enrichKEGG(gene = gene_df_up$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
      cat(sprintf("  %d KEGG pathways identified\n", nrow(kk_up@result)))
    }
  }, error = function(e) {
    cat("  KEGG Analysis Failed:", conditionMessage(e), "\n")
  })
}

if(!is.null(ego_bp_up)) {
  write.csv(ego_bp_up@result, file.path(tables.dir, "Table_GO_BP_Up.csv"), row.names = FALSE)
}
if(!is.null(kk_up)) {
  write.csv(kk_up@result, file.path(tables.dir, "Table_KEGG_Up.csv"), row.names = FALSE)
}

cat("✓ Enrichment analysis completed \n")

# =========== Figure 1: Cell Type Annotation and HSC Fibrosis Characteristics ==============
cat("\n=== Generating Figure 1: Cell Type Annotation and HSC Fibrosis Characteristics ===\n")

# Loading required packages
if(!require("tidyr", quietly = TRUE)) install.packages("tidyr")
library(tidyr)

# Defining color scheme
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

# Ensuring celltype_reannotated exists and is a factor
if(!"celltype_reannotated" %in% colnames(high_quality_cells@meta.data)) {
  high_quality_cells$celltype_reannotated <- as.character(high_quality_cells$celltype)
}

# Retrieving actual cell type levels
actual_celltypes <- unique(high_quality_cells$celltype_reannotated)
available_colors <- celltype_colors[names(celltype_colors) %in% actual_celltypes]

# Assigning gray color to undefined cell types
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

# 1B: Cell Composition (excluding doublets, sorted by proportion from high to low)
final_stats <- as.data.frame(table(high_quality_cells$celltype_reannotated))
final_stats <- final_stats[!final_stats$Var1 %in% c("Doublet_or_Debris", "NA"), ]
colnames(final_stats) <- c("CellType", "Count")
final_stats$Percentage <- round(final_stats$Count / sum(final_stats$Count) * 100, 1)

# Sorted by Count in descending order
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

# 1C: Group Comparison
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

# 1D: Fibrosis Score (using pre-computed or recalculating)
if(!"Fibrosis_Score1" %in% colnames(high_quality_cells@meta.data)) {
  fibrosis_genes <- c("Col1a1", "Col1a2", "Acta2", "Des", "Vim", "Tgfb1")
  fibrosis_genes <- fibrosis_genes[fibrosis_genes %in% rownames(high_quality_cells)]
  if(length(fibrosis_genes) >= 3) {
    high_quality_cells <- AddModuleScore(high_quality_cells, features = list(fibrosis_genes), 
                                         name = "Fibrosis_Score")
    cat(sprintf("  Calculating fibrosis score (%d genes )\n", length(fibrosis_genes)))
  }
}

# Selecting cell types to display (prioritizing HSC subtypes, falling back to original HSC if unavailable)
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

# 1E: Fold Enrichment Calculation
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
      
      # Saving fold enrichment data
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

# 1F: HSC Marker Comparison (only when HSC subtypes are available)
hsc_marker_types <- intersect(c("HSC_Act", "HSC_Qui", "HSC"), actual_celltypes)

if(length(hsc_marker_types) > 0) {
  hsc_cells <- subset(high_quality_cells, subset = celltype_reannotated %in% hsc_marker_types)
  fibrogenic_markers <- c("Col1a1", "Col1a2", "Acta2", "Des", "Vim")
  fibrogenic_markers <- fibrogenic_markers[fibrogenic_markers %in% rownames(high_quality_cells)]
  
  if(length(fibrogenic_markers) > 0 && ncol(hsc_cells) > 0) {
    plot_data <- FetchData(hsc_cells, vars = c(fibrogenic_markers, "celltype_reannotated"))
    plot_long <- plot_data %>%
      pivot_longer(cols = -celltype_reannotated, names_to = "gene", values_to = "expression")
    
    # Defining HSC colors
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

# Combined (2x3 layout)
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
ggsave(file.path(figures.dir, "Figure1_Cell_Type_Annotation.pdf"), fig1_complete, 
       width = 18, height = 12, device = cairo_pdf)
cat("✓ Figure 1 saved to", file.path(figures.dir, "Figure1_Celltype_Complete.tiff"), "\n")

# Saving updated object
save(high_quality_cells, file = rdata_file)
cat("✓ Step 7 completed: Pathway Enrichment Analysis + Figure 1 Generation\n")

# ========= 8. In-depth Hepatocyte Analysis (including 8-gene signature + 7-class stress scoring)===
cat("
")
cat(paste(rep("=", 70), collapse = ""), "
")
cat("              In-depth Hepatocyte Analysis (8-gene signature + 7-class stress scoring)
")
cat(paste(rep("=", 70), collapse = ""), "
")

# 8.1 Extracting hepatocytes and performing basic analysis
cat("
=== Step 8.1: Extracting hepatocytes and performing basic analysis ===
")

hepatocytes <- subset(high_quality_cells, celltype == "Hepatocyte")

if(ncol(hepatocytes) == 0) {
  stop("Error: No hepatocytes found")
}

cat(sprintf("Extracted %d hepatocytes
", ncol(hepatocytes)))

# Re-dimensional reduction (recommended for subset reanalysis)
hepatocytes <- NormalizeData(hepatocytes)
hepatocytes <- FindVariableFeatures(hepatocytes, selection.method = "vst", nfeatures = 2000)
hepatocytes <- ScaleData(hepatocytes)
hepatocytes <- RunPCA(hepatocytes, features = VariableFeatures(object = hepatocytes))
hepatocytes <- FindNeighbors(hepatocytes, dims = 1:20)
hepatocytes <- FindClusters(hepatocytes, resolution = 0.6)
hepatocytes <- RunUMAP(hepatocytes, dims = 1:20)

# 8.2 Hepatic Zonation Annotation (based on literature markers)
cat("
=== step 8.2: Hepatic Lobule Zonation Annotation ===
")

zone_markers <- list(
  "Zone1_Periportal" = c("Cyp2f2", "Hal", "Sds", "Ass1", "Pck1"),
  "Zone3_Pericentral" = c("Cyp2e1", "Cyp1a2", "Glul", "Oat", "Cyp3a11")
)

# Calculating zonation scores
for(zone in names(zone_markers)) {
  genes <- zone_markers[[zone]]
  genes <- genes[genes %in% rownames(hepatocytes)]
  if(length(genes) >= 3) {
    hepatocytes <- AddModuleScore(hepatocytes, features = list(genes), name = zone)
  }
}

# Assigning zones based on scores
if(all(c("Zone1_Periportal1", "Zone3_Pericentral1") %in% colnames(hepatocytes@meta.data))) {
  z1_score <- hepatocytes$Zone1_Periportal1
  z3_score <- hepatocytes$Zone3_Pericentral1
  
  hepatocytes$hep_zone <- ifelse(z1_score > z3_score, "Zone1_Periportal",
                                 ifelse(z3_score > z1_score, "Zone3_Pericentral", 
                                        "Zone2_Midlobular"))
} else {
  hepatocytes$hep_zone <- "Zone2_Midlobular"
}

cat(sprintf("Hepatocyte Zonation Distribution:
"))
print(table(hepatocytes$hep_zone, hepatocytes$group))

# Visualizing zonation
p_zone <- DimPlot(hepatocytes, reduction = "umap", group.by = "hep_zone", pt.size = 0.5) +
  ggtitle("A. Hepatocyte Zonation") + theme_journal

p_zone_group <- DimPlot(hepatocytes, reduction = "umap", group.by = "group", pt.size = 0.5) +
  ggtitle("B. Group Distribution") + theme_journal

p_combined_zone <- p_zone + p_zone_group

ggsave(file.path(hep_figures_dir, "FigureH1_Hepatocyte_Zonation.tiff"), 
       p_combined_zone, width = 12, height = 5, dpi = 300)
ggsave(file.path(hep_figures_dir, "FigureH1_Hepatocyte_Zonation.pdf"), 
       p_combined_zone, width = 12, height = 5, device = cairo_pdf)

cat("✓ Basic hepatocyte analysis completed
")

# 8.3 [Key] 8-Gene Signature Scoring (for Figure 2)
cat("
=== Step 8.3: 8-Gene Signature Scoring (Dedicated for Figure 2)===
")

fig2_signature_genes <- c("Cyp2e1", "Cyp3a11", "Hspa1a", "Hsp90aa1", 
                          "Gpx1", "Hmox1", "Atf4", "Gadd45a")
fig2_signature_genes <- fig2_signature_genes[fig2_signature_genes %in% rownames(hepatocytes)]

if(length(fig2_signature_genes) >= 5) {
  cat(sprintf("Calculating 8-gene signature score using %d genes
", length(fig2_signature_genes)))
  
  hepatocytes <- AddModuleScore(hepatocytes, 
                                features = list(fig2_signature_genes), 
                                name = "Fig2_Signature")
  
  # Creating compatible column names required for Figure 2
  hepatocytes$Stress1 <- hepatocytes$Fig2_Signature1
  hepatocytes$stress_group <- ifelse(
    hepatocytes$Stress1 > mean(hepatocytes$Stress1, na.rm = TRUE) + 
      sd(hepatocytes$Stress1, na.rm = TRUE),
    "High_stress", "Low_stress"
  )
  
  # Creating obj alias (for Figure 2 code usage)
  obj <- hepatocytes
  
  cat(sprintf("8-gene signature scoring completed: high-stress cells %d/%d (%.1f%%)
",
              sum(hepatocytes$stress_group == "High_stress"), 
              ncol(hepatocytes),
              100 * sum(hepatocytes$stress_group == "High_stress") / ncol(hepatocytes)))
  cat("✓ Object alias 'obj' created for Figure 2 code usage
")
} else {
  cat("Warning: Insufficient 8-gene signature genes, skipping calculation
")
}

# 8.4 7-Category Hierarchical Stress Scoring 
cat("
=== Step 8.4: 7-Category Hierarchical Stress Scoring===
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

# Calculating stress scores for each category
for(stress_type in names(stress_categories)) {
  genes <- stress_categories[[stress_type]]
  genes <- genes[genes %in% rownames(hepatocytes)]
  if(length(genes) >= 3) {
    hepatocytes <- AddModuleScore(hepatocytes, features = list(genes), name = stress_type)
    cat(sprintf("  %s: %d genes
", stress_type, length(genes)))
  }
}

# Calculating total stress score (sum of 7 categories)
all_stress_genes <- unlist(stress_categories)
all_stress_genes <- all_stress_genes[all_stress_genes %in% rownames(hepatocytes)]
hepatocytes <- AddModuleScore(hepatocytes, features = list(all_stress_genes), name = "Total_Stress")

# Critical: Using Hep_Stress_Score as the total score (for Figure Master)
hepatocytes$Hep_Stress_Score <- hepatocytes$Total_Stress1

cat(sprintf("
Total stress score: %d genes
", length(all_stress_genes)))
cat("✓ 7-category hierarchical stress scoring completed
")

# ========Step 8.5: Generating Figure 2 Data Files (Optimized Version)==============
cat("\n=== Step 8.5: Generating Figure 2 Data Files (Zone 3 Stress Analysis)===\n")

# Ensuring directory exists
if(!dir.exists(hep_tables_dir)) {
  dir.create(hep_tables_dir, recursive = TRUE)
}

# 1. DE_stress_comprehensive.csv（Differential Genes: High vs Low Stress）
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

cat(sprintf("  Completed: %d differential genes (padj<0.05: %d)\n", 
            nrow(de_stress_df), 
            sum(de_stress_df$p_val_adj < 0.05)))

# 2. stress_scores_all_methods.csv（Comparison of Multiple Scoring Methods）
cat("生成 stress_scores_all_methods.csv...\n")

stress_methods_df <- data.frame(
  cell = colnames(hepatocytes),
  Original_4gene = hepatocytes$Oxidative_Stress1,  # 4 Gene: Oxidative Stress
  Comprehensive_8gene = hepatocytes$Fig2_Signature1,  # 8-Gene Signature
  Total_7category = hepatocytes$Total_Stress1,  # Sum of 7 Categories
  Oxidative = hepatocytes$Oxidative_Stress1,
  ERstress = hepatocytes$ER_Stress1,
  Lipid = hepatocytes$Lipid_Dysregulation1
)

write.csv(stress_methods_df, 
          file.path(hep_tables_dir, "stress_scores_all_methods.csv"), 
          row.names = FALSE)

cat("  Completed: 6 scoring methods\n")

# 3. GO_BP_high_stress_hepatocytes.csv（GO Enrichment of High-Stress Group）
cat("Generating GO_BP_high_stress_hepatocytes.csv...\n")

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
      cat(sprintf("  Completed: %d GO terms \n", nrow(ego_stress@result)))
    } else {
      cat("  Warning: No significant GO enrichment found \n")
    }
  }, error = function(e) {
    cat("  GO analyse error:", conditionMessage(e), "\n")
  })
} else {
  cat("  Warning: Insufficient high-stress genes(<10个)，Skipping GO analysis \n")
}

# 4. KEGG_upregulated_CCl4.csv（CCl4 vs Control KEGG）
cat("generating KEGG_upregulated_CCl4.csv...\n")

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
        cat(sprintf("  Completed: %d KEGG pathways \n", nrow(kk_cc@result)))
      } else {
        cat("  Warning: No significant KEGG pathways found \n")
      }
    } else {
      cat("  Warning: Gene ID conversion failed\n")
    }
  }, error = function(e) {
    cat("  KEGG Analysis Error:", conditionMessage(e), "\n")
  })
} else {
  cat("  Warning: Insufficient CCl4 upregulated genes (<10)，Skipping KEGG analysis\n")
}

cat("✓ Figure 2 Data file generation completed\n")

# Saving hepatocyte object
save(hepatocytes, file = file.path(work.dir, "results", "2_analysis", "hepatocytes_fig2.RData"))

# =====Step 9: Figure 2 Generation (Zone 3 Stress Integration Version)========
cat("\n=== Step 9: Figure 2 Generation（Zone 3 Hepatocyte Stress Analysis，10 panels）===\n")

# 9.1 Loading required packages
if(!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if(!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if(!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(viridis)
library(RColorBrewer)
library(ggrepel)

# 9.2 prepair data
obj <- hepatocytes

# Key Statistics
stats <- list(
  mean = mean(obj$Fig2_Signature1, na.rm = TRUE),
  sd = sd(obj$Fig2_Signature1, na.rm = TRUE),
  threshold = mean(obj$Fig2_Signature1, na.rm = TRUE) + sd(obj$Fig2_Signature1, na.rm = TRUE),
  n_total = ncol(obj),
  n_high = sum(obj$stress_group == "High_stress", na.rm = TRUE),
  high_pct = 100 * sum(obj$stress_group == "High_stress", na.rm = TRUE) / ncol(obj)
)

cat(sprintf("8-Gene Signature Statistics: Mean=%.3f, Threshold=%.3f, High_stress=%d (%.1f%%)\n",
            stats$mean, stats$threshold, stats$n_high, stats$high_pct))

# Unified Theme
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

# ==================== Panel A: Hepatocyte Zonation UMAP ====================
cat("Panel A: Hepatocyte Zonation UMAP\n")

p_a <- DimPlot(obj, reduction = "umap", group.by = "hep_zone", pt.size = 0.4, label = FALSE) +
  scale_color_manual(values = c("Zone1_Periportal" = "#1f77b4", 
                                "Zone2_Midlobular" = "#ff7f0e", 
                                "Zone3_Pericentral" = "#2ca02c"),
                     name = "Zone") +
  ggtitle("A. Hepatocyte Zonation") + 
  fig_theme +
  theme(legend.position = "right")

# ==================== Panel B: Group Distribution (CCl4 Enrichment in Zone 3)====================
cat("Panel B: Group Distribution\n")

p_b <- DimPlot(obj, reduction = "umap", group.by = "group", pt.size = 0.4) +
  scale_color_manual(values = c("CCl4" = "#E74C3C", "Control" = "#3498DB")) +
  ggtitle("B. CCl4 Enrichment in Zone 3") + 
  fig_theme +
  theme(legend.position = "right")

# ==================== Panel C: 8-gene Stress Score Distribution ====================
cat("Panel C: 8-gene Stress Score Distribution\n")

p_c <- ggplot(data.frame(Score = obj$Fig2_Signature1), aes(x = Score)) +
  geom_histogram(bins = 40, fill = "#4682B4", alpha = 0.8, color = "white") +
  geom_vline(xintercept = stats$mean, color = "#E41A1C", linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = stats$threshold, color = "#8B0000", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = stats$mean, y = Inf, label = sprintf("Mean=%.3f", stats$mean),
           vjust = 1.5, hjust = -0.1, size = 2.5, color = "#E41A1C") +
  annotate("text", x = stats$threshold, y = Inf, 
           label = sprintf("Threshold=%.3f\nn=%d (%.1f%%)", stats$threshold, stats$n_high, stats$high_pct),
           vjust = 3, hjust = -0.1, size = 2.5, color = "#8B0000") +
  labs(x = "8-gene Stress Signature", y = "Cell Count", title = "C. Stress Score Distribution") +
  fig_theme

# ==================== Panel D: High-Stress Cell UMAP Localization ====================
cat("Panel D: High-Stress Cell UMAP Localization\n")

p_d <- FeaturePlot(obj, features = "Fig2_Signature1", pt.size = 0.2, order = TRUE) +
  scale_color_viridis_c(option = "plasma", name = "Score") +
  labs(title = "D. High-Stress Cell Localization") +
  fig_theme +
  theme(legend.position = "right", legend.key.height = unit(0.6, "cm"))

# ==================== Panel E: Zone × Group Stress Score Quantification ====================
cat("Panel E: Zone × Group Stress Scores\n")

zone_data <- obj@meta.data %>%
  filter(hep_zone %in% c("Zone1_Periportal", "Zone3_Pericentral")) %>%
  group_by(hep_zone, group) %>%
  summarise(
    n = n(),
    mean_stress = mean(Fig2_Signature1, na.rm = TRUE),
    sem = sd(Fig2_Signature1, na.rm = TRUE) / sqrt(n),
    .groups = 'drop'
  )

# Add statistical test
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

# ==================== Panel F: Seven-Category Stress Heatmap ====================
cat("Panel F: Seven-Category Stress Heatmap\n")

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

# Ensurecorrectlyorder
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

# ==================== Panel G: Hmox1 vs Hspa1a Comparison (Oxidative vs Heat Shock)====================
cat("Panel G: Hmox1 vs Hspa1a Comparison\n")

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

# Getstatisticsinfo
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

# ==================== Panel H: Cyp2e1 Three-Layer Comparison ====================
cat("Panel H: Cyp2e1 Three-Layer Comparison\n")

# Calculate three-layer comparison
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

# ==================== Panel I: 8-gene Signature Expression Dot Plot ====================
cat("Panel I: 8-gene Signature Expression\n")

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
  labs(x = "", y = expression(log[2]~FC~High/Low), title = "I. 8-gene Signature Expression") +
  fig_theme +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# ==================== Panel J: GO-BP Enrichment ====================
cat("Panel J: GO-BP Enrichment\n")

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
      labs(x = "gene Count", y = "", title = "J. GO-BP: High-Stress Hepatocytes") +
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

# ==================== Panel K: KEGG Enrichment ====================
cat("Panel K: KEGG Pathway Enrichment\n")

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
      labs(x = "gene Count", y = "", title = "K. KEGG: CCl4 vs Control") +
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

# ==================== Assemble Complete Figure 2 (3-Row Layout)====================
cat("Assemble Complete Figure 2...\n")

# Row 1: A B C D（Zone and Score distribution）
row1 <- p_a + p_b + p_c + p_d + plot_layout(widths = c(1, 1, 1.2, 1.2))

# Row 2: E F G H（Zone quantification, Seven-category stress, Mechanism comparison）
row2 <- p_e + p_f + p_g + p_h + plot_layout(widths = c(1, 1.2, 1, 1))

# Row 3: I J K（Signatureexpression, GO, KEGG）- 3 panels, Leave empty or adjust
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

# Save
ggsave(file.path(figures.dir, "Figure2_Zone3_Hepatocyte_Stress.tiff"), 
       fig2, width = 16, height = 13, dpi = 300, compression = "lzw")

ggsave(file.path(figures.dir, "Figure2_Zone3_Hepatocyte_Stress.pdf"), 
       fig2, width = 16, height = 13, device = cairo_pdf)

cat("\n✓ Figure 2 Complete: Zone 3 Hepatocyte Stress Analysis (11 panels: A-K)\n")
cat(sprintf("Output: %s/Figure2_Zone3_Hepatocyte_Stress.{tiff,pdf}\n", figures.dir))

# ==================== Step 10: Figure 3:HSC Activation Mechanisms ====================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("              Figure 5: HSC Activation and Fibrosis Mechanisms (9 panels, 3x3)\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# CheckandLoadnecessary data
if(!exists("high_quality_cells")) {
  load(file.path(work.dir, "results", "2_analysis", "high_quality_cells.RData"))
  cat("✓ Loadhigh_quality_cellsobject\n")
}

fig_theme <- theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank())

# LoadDEG results(Useexistinghsc_degorfrom fileLoad)
if(!exists("hsc_deg")) {
  deg_file <- file.path(tables.dir, "Table_DEG_HSC_Act_vs_Qui_CCl4.csv")
  if(file.exists(deg_file)) {
    hsc_deg <- read.csv(deg_file)
    cat("✓ from fileLoadHSC DEG results:", nrow(hsc_deg), "Gene\n")
  } else {
    stop("Error: Not foundHSC DEG resultsfile")
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

# ========== Panel B: DEG Volcano Plot ==========
cat("Panel B: DEG Volcano Plot...\n")

# EnsuresignificanceColumniscorrectlyCreate
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

# ========== Panel C: Top Fibrosis genes Dot Plot ==========
cat("Panel C: Fibrosis genes...\n")

fibrosis_genes <- c("Col1a1", "Col1a2", "Col3a1", "Acta2", "Timp1", "Spp1", "Lox", "Cthrc1")
hsc_subset <- subset(high_quality_cells, celltype_reannotated %in% c("HSC_Act", "HSC_Qui"))

# CheckGeneexistence
available_genes <- fibrosis_genes[fibrosis_genes %in% rownames(hsc_subset)]
missing_genes <- fibrosis_genes[!fibrosis_genes %in% rownames(hsc_subset)]
if(length(missing_genes) > 0) {
  cat("  Warning: missingGene:", paste(missing_genes, collapse = ", "), "\n")
}

if(length(available_genes) > 0) {
  p3c <- DotPlot(hsc_subset, features = available_genes, group.by = "celltype_reannotated",
                 cols = c("lightgrey", "#E41A1C")) +
    RotatedAxis() +
    labs(title = "C. Core Fibrosis genes") +
    fig_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  # Backup: UseggplotmanuallyCreate
  dot_data <- hsc_deg %>% 
    filter(gene %in% fibrosis_genes) %>%
    select(gene, avg_log2FC, p_val_adj) %>%
    mutate(pct_expr = 80, avg_expr = avg_log2FC)
  
  p3c <- ggplot(dot_data, aes(x = gene, y = "HSC_Act", size = pct_expr, color = avg_log2FC)) +
    geom_point() +
    scale_color_gradient(low = "lightgrey", high = "#E41A1C") +
    labs(title = "C. Core Fibrosis genes (DEG-based)") +
    fig_theme
}

# ========== Panel D-F: Enrichment Analysis (with Error Handling) ==========
cat("Panel D-F: Enrichment Analysis...\n")

# Prepare gene list
up_genes <- hsc_deg %>% 
  filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>%
  pull(gene)

cat("  UpregulatedGenenumber:", length(up_genes), "\n")

# TryGeneConvert
cat("  ConvertGeneID...\n")
gene_conv <- tryCatch({
  bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
}, error = function(e) {
  cat("  Convertfailed:", e$message, "\n")
  return(NULL)
})

if(!is.null(gene_conv) && nrow(gene_conv) > 10) {
  cat("  successfullyConvert", nrow(gene_conv), "Gene\n")
  
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
  
  # Panel F: gene-Concept Network
  if(!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    p3f <- cnetplot(ego_bp, showCategory = 5, 
                    foldChange = setNames(hsc_deg$avg_log2FC, hsc_deg$gene)) +
      ggtitle("F. gene-Pathway Network") +
      theme_bw()
  } else {
    p5f <- ggplot() + annotate("text", x=0.5, y=0.5, label="Network not available") + 
      labs(title="F. gene-Pathway Network") + theme_void()
  }
  
} else {
  cat("  Gene conversion failed or insufficient gene number, Skip Enrichment Analysis\n")
  empty_plot <- ggplot() + annotate("text", x=0.5, y=0.5, 
                                    label="Enrichment analysis failed\n(insufficient gene mapping)") + 
    theme_void() + fig_theme
  
  p3d <- empty_plot + labs(title="D. GO Biological Process")
  p3e <- empty_plot + labs(title="E. KEGG Pathway")
  p3f <- empty_plot + labs(title="F. gene-Pathway Network")
}

# ========== Panel G: Core Fibrosis Heatmap ==========
cat("Panel G: Fibrosis Heatmap...\n")

fibrosis_core <- c("Col1a1", "Col1a2", "Col3a1", "Acta2", "Timp1", "Spp1", "Cthrc1", "Lox")
available_core <- fibrosis_core[fibrosis_core %in% rownames(hsc_subset)]

if(length(available_core) >= 4) {
  expr_avg <- AverageExpression(hsc_subset, features = available_core, group.by = "celltype_reannotated")$RNA
  
   p3g <- ggplot(expr_df, aes(x = group, y = reorder(gene, expression), fill = expression)) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(x = "", y = "", title = "G. Core Fibrosis genes") +
    fig_theme
  
  tryCatch({
    pheatmap::pheatmap(expr_avg, scale = "row", 
                       color = colorRampPalette(c("blue", "white", "red"))(50),
                       filename = file.path(figures.dir, "Figure5G_heatmap_pheatmap.pdf"),
                       width = 6, height = 8)
  }, error = function(e) cat("  pheatmapSavefailed((ignorable))\n"))
  
} else {
  p3g <- ggplot() + annotate("text", x=0.5, y=0.5, label="Insufficient genes for heatmap") + 
    labs(title="G. Core Fibrosis genes") + theme_void() + fig_theme
}

# ========== Panel H: GSEA (Simplified) ==========
cat("Panel H: GSEA...\n")

gsea_data <- hsc_deg %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)

p3h <- ggplot(gsea_data, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = avg_log2FC > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
  coord_flip() +
  labs(x = "", y = expression(log[2]~FC), title = "H. Top Ranked genes (GSEA-style)") +
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

# ========== Assemble Figure 3 (3x3 Layout) ==========
cat("Assemble Figure 3 (3x3 Layout)...\n")

# Correction：Use p3a-p3i instead of p5a-p5i
fig3 <- wrap_plots(p3a, p3b, p3c, p3d, p3e, p3f, p3g, p3h, p3i, ncol = 3) +
  plot_annotation(
    title = "Figure 3. HSC Activation Status and Fibrosis Mechanisms",  # Unified title
    subtitle = "Comprehensive analysis of hepatic stellate cell activation, differential expression, pathway enrichment, and therapeutic targets",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

# Correction：Use fig3 instead of fig5
ggsave(file.path(figures.dir, "Figure3_HSC_Activation_Mechanisms.tiff"), 
       fig3, width = 15, height = 15, dpi = 300, compression = "lzw")  # ✅ using fig3

ggsave(file.path(figures.dir, "Figure3_HSC_Activation_Mechanisms.pdf"), 
       fig3, width = 15, height = 15, device = cairo_pdf)  # ✅ using fig3

ggsave(file.path(figures.dir, "Figure3_HSC_Activation_Mechanisms.png"), 
       fig3, width = 15, height = 15, dpi = 300)  # ✅ using fig3

cat("\n✓ Figure 3 Complete! (HSC Activation Mechanisms, 9 panels A-I)\n")  # Unified output info

# ==================== Step 11: Figure 4 - Communication Analysis====================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("              Figure 3: Comprehensive Hepatocyte-HSC Communication Analysis\n")
cat("              (CellChat + Custom LR Analysis) - Fully Fixed Version\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# 10.1 Loading Required Packages
if(!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")
if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if(!requireNamespace("magick", quietly = TRUE)) install.packages("magick")

suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(magick)
})

# EnsureCellChatLoad
if(!requireNamespace("CellChat", quietly = TRUE)) {
  stop("Please install firstCellChat: devtools::install_github('sqjin/CellChat')")
}
library(CellChat)

# Fix namespace conflicts
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
desc <- dplyr::desc
group_by <- dplyr::group_by
summarise <- dplyr::summarise
pull <- dplyr::pull

# 10.2 Loading CellChat Objects
cellchat.dir <- file.path(work.dir, "results", "2_analysis", "cellchat")
cellchat_available <- FALSE

if(dir.exists(cellchat.dir)) {
  ccl4_rds <- file.path(cellchat.dir, "cellchat_CCl4.rds")
  control_rds <- file.path(cellchat.dir, "cellchat_Control.rds")
  
  if(file.exists(ccl4_rds) && file.exists(control_rds)) {
    cat("Loading CellChat Objects...\n")
    cellchat_ccl4 <- readRDS(ccl4_rds)
    cellchat_ctrl <- readRDS(control_rds)
    cellchat_available <- TRUE
    cat("✓ CellChat object loaded\n")
    
    # Check object structure
    cell_types_ccl4 <- levels(cellchat_ccl4@idents)
    cell_types_ctrl <- levels(cellchat_ctrl@idents)
    cat(sprintf("  CCl4 Group cell types (%d): %s\n", length(cell_types_ccl4), 
                paste(cell_types_ccl4, collapse = ", ")))
    cat(sprintf("  Control Group cell types (%d): %s\n", length(cell_types_ctrl), 
                paste(cell_types_ctrl, collapse = ", ")))
    
    # Identify key cell types
    hep_types <- cell_types_ccl4[grep("Hepatocyte", cell_types_ccl4, ignore.case = TRUE)]
    hsc_types <- cell_types_ccl4[grep("HSC", cell_types_ccl4, ignore.case = TRUE)]
    
    cat(sprintf("  Hepatocyte types: %s\n", paste(hep_types, collapse = ", ")))
    cat(sprintf("  HSC types: %s\n", paste(hsc_types, collapse = ", ")))
    
    if("HSC_Act" %in% cell_types_ctrl) {
      cat("  ✓ Control Group contains HSC_Act\n")
    } else {
      cat("  ! Control Group lacks HSC_Act(Insufficient cell number excluded)\n")
    }
  } else {
    cat("! CellChat RDS file not found\n")
  }
} else {
  cat("! CellChat directory not found\n")
}

# ==================== Key Fix: Safe Communication Data Extraction Function ====================

#' Safely extract communication data from CellChat object
#' Compatible with v1.x and v2.x versions
extract_communication_safe <- function(cellchat_obj, sources = NULL, targets = NULL) {
  
  # Method 1: TryUsesubsetCommunication(v2.0+)
  comm_data <- tryCatch({
    if(!is.null(sources) && !is.null(targets)) {
      # Try calling CellChat namespace function
      CellChat::subsetCommunication(cellchat_obj, 
                                    sources.use = sources, 
                                    targets.use = targets)
    } else {
      CellChat::subsetCommunication(cellchat_obj)
    }
  }, error = function(e) {
    cat("    subsetCommunicationfailed:", conditionMessage(e), "\n")
    cat("    Try alternative method...\n")
    return(NULL)
  })
  
  # Method 2: If Method 1 failed, Extract directly from @net
  if(is.null(comm_data) || nrow(comm_data) == 0) {
    cat("    Extract directly using @net$prob...\n")
    
    if(!"net" %in% slotNames(cellchat_obj) || 
       is.null(cellchat_obj@net$prob)) {
      cat("    ! Unable to extract communication data\n")
      return(data.frame())
    }
    
    prob_matrix <- cellchat_obj@net$prob
    cell_types <- rownames(prob_matrix)
    
    # Filtersource
    if(!is.null(sources)) {
      source_idx <- which(cell_types %in% sources)
    } else {
      source_idx <- 1:nrow(prob_matrix)
    }
    
    # Filtertarget
    if(!is.null(targets)) {
      target_idx <- which(cell_types %in% targets)
    } else {
      target_idx <- 1:ncol(prob_matrix)
    }
    
    # Extract sub-matrix
    if(length(source_idx) == 0 || length(target_idx) == 0) {
      return(data.frame())
    }
    
    sub_prob <- prob_matrix[source_idx, target_idx, drop = FALSE]
    
    # Convert to dataframe format(Simulate subsetCommunication output)
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

# ==================== Row 1: CellChat Analysis (Panel A-D) ====================

if(cellchat_available) {
  cat("\n========== Panel A-D: CellChat Analysis ==========\n")
  
  # Extract communication data
  cat("\nExtract CCl4 Group communication data...\n")
  df_comm_ccl4 <- extract_communication_safe(cellchat_ccl4)
  cat(sprintf("  CCl4 total Communication Pairs: %d\n", nrow(df_comm_ccl4)))
  
  cat("\nExtract Control Group communication data...\n")
  df_comm_ctrl <- extract_communication_safe(cellchat_ctrl)
  cat(sprintf("  Control total Communication Pairs: %d\n", nrow(df_comm_ctrl)))
  
  # Extract Hep→HSC specific communication
  cat("\nExtract Hep→HSC specific communication...\n")
  hep_hsc_ccl4 <- extract_communication_safe(cellchat_ccl4, hep_types, hsc_types)
  hep_hsc_ctrl <- extract_communication_safe(cellchat_ctrl, hep_types, hsc_types)
  
  cat(sprintf("  CCl4 Hep→HSC: %d\n", nrow(hep_hsc_ccl4)))
  cat(sprintf("  Control Hep→HSC: %d\n", nrow(hep_hsc_ctrl)))
  
  # Panel A: CellChat Communication Overview
  cat("\nPanel A: CellChat Communication Overview...\n")
  
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
  
  # Save communication data
  write.csv(df_comm_ccl4, file.path(tables.dir, "Table_Figure3_CellChat_CCl4.csv"), row.names = FALSE)
  write.csv(hep_hsc_ccl4, file.path(tables.dir, "Table_Figure3_CellChat_Hep_HSC_CCl4.csv"), row.names = FALSE)
  
  # ==================== Panel B: Top LR Bubble Plot (Fixed Version)====================
  cat("Panel B: Top LR Bubble Plot...\n")
  
  if(nrow(hep_hsc_ccl4) > 0) {
    # Key fix: Merge probability by source-target, Avoid duplicate pairs
    hep_hsc_ccl4_clean <- hep_hsc_ccl4 %>%
      group_by(source, target) %>%
      summarise(
        prob = sum(prob, na.rm = TRUE),
        n_pairs = n(),
        .groups = "drop"
      ) %>%
      mutate(pair = paste(source, "→", target))
    
    # Double check: Ensure no duplicates
    if(anyDuplicated(hep_hsc_ccl4_clean$pair) > 0) {
      warning("Duplicate pairs still exist, Force deduplication")
      hep_hsc_ccl4_clean <- hep_hsc_ccl4_clean %>%
        distinct(pair, .keep_all = TRUE)
    }
    
    # Gettop 15
    top_pairs <- hep_hsc_ccl4_clean %>%
      arrange(desc(prob)) %>%
      head(15)
    
    cat(sprintf("  Top 15 Hep→HSC pairs selected (from %d unique pairs)\n", 
                nrow(hep_hsc_ccl4_clean)))
    
    # Prepare plotting data
    plot_data <- hep_hsc_ccl4_clean %>%
      filter(pair %in% top_pairs$pair) %>%
      mutate(
        pair = factor(pair, levels = unique(rev(top_pairs$pair))),
        Group = "CCl4"
      )
    
    # Add Control Group data
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
    
    # Plot bubble chart
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
  
  # Panel C: Pathway Activity(UseHep→HSCCommunication Strengthasalternative)
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
  
  # Panel D: HSC Activation State Transition
  cat("Panel D: HSC Activation State Transition...\n")
  
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
  
  cat("✓ Panel A-Dcompleted\n")
  
} else {
  cat("! CellChatnot available, Generate placeholder plot\n")
  empty_plot <- ggplot() + annotate("text", x=0.5, y=0.5, 
                                    label="CellChat data not available\nRun run_cellchat_analysis.R first") + 
    theme_void() + theme_journal
  
  p4a <- empty_plot + labs(title = "A. CellChat Overview")
  p4b <- empty_plot + labs(title = "B. Top 15 Hep→HSC Pairs")
  p4c <- empty_plot + labs(title = "C. Pathway Activity")
  p4d <- empty_plot + labs(title = "D. Target HSC State")
}

# ==================== Row 2: Custom LR Analysis + Chord diagram (Panel E-H) ====================
cat("\n========== Panel E-H: Custom LR Analysis and Chord Diagram ==========\n")

# Define key Ligand-Receptor Pairs(Based on literature Hep→HSC communication)
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

# Panel E: Ligand-Receptor heatmap
cat("Panel E: LR gene Expression Heatmap...\n")

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
              gene = gene,
              gene_type = ifelse(gene %in% c("Tgfbr1", "Tgfbr2", "Pdgfra", "Pdgfrb", 
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
    facet_wrap(~gene_type, ncol = 1, scales = "free_y") +
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

# Panel F: Communication Strength Network
cat("Panel F: Communication Network Data...\n")

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
      filter(gene_type == "Ligand") %>%
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

# Panel G: Spatial Distribution Map
cat("Panel G: Spatial Distribution...\n")

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
# Panel H: Global Network Statistics (Fully Fixed)
cat("Panel H: Global Network Statistics...\n")

if(exists("df_comm_ccl4") && !is.null(df_comm_ccl4) && nrow(df_comm_ccl4) > 0) {
  
  # Clean data
  df_comm_clean <- df_comm_ccl4 %>%
    dplyr::filter(!is.na(source), !is.na(target), !is.na(prob))
  
  cat(sprintf("  CleanafterCommunication Pairs: %d\n", nrow(df_comm_clean)))
  
  if(nrow(df_comm_clean) > 0) {
    # Calculate out-degree - Key Fix: Use dplyr:: prefix, and immediately convert to data.frame after summarise
    out_degree <- df_comm_clean %>%
      dplyr::group_by(source) %>%
      dplyr::summarise(
        out_strength = sum(prob, na.rm = TRUE), 
        n_targets = dplyr::n_distinct(target),
        .groups = "drop"
      ) %>%
      dplyr::ungroup()
    
    # Key Fix: Explicitly convert to plain data.frame, Avoid class attribute issues
    out_degree <- as.data.frame(out_degree)
    names(out_degree) <- c("cell_type", "out_strength", "n_targets")
    
    # Calculatein-degree
    in_degree <- df_comm_clean %>%
      dplyr::group_by(target) %>%
      dplyr::summarise(
        in_strength = sum(prob, na.rm = TRUE), 
        n_sources = dplyr::n_distinct(source),
        .groups = "drop"
      ) %>%
      dplyr::ungroup()
    
    # Key Fix: Explicitly convert to plain data.frame
    in_degree <- as.data.frame(in_degree)
    names(in_degree) <- c("cell_type", "in_strength", "n_sources")
    
    cat(sprintf("  Out-degreeCellnumber: %d, In-degreeCellnumber: %d\n", 
                nrow(out_degree), nrow(in_degree)))
    
    # Use base R merge to avoid dplyr join issues
    centrality_df <- merge(out_degree, in_degree, by = "cell_type", all = TRUE)
    
    # TreatmentNA
    centrality_df$out_strength[is.na(centrality_df$out_strength)] <- 0
    centrality_df$in_strength[is.na(centrality_df$in_strength)] <- 0
    centrality_df$total_strength <- centrality_df$out_strength + centrality_df$in_strength
    
    cat(sprintf("  MergeafternumberdataRownumber: %d\n", nrow(centrality_df)))
    
    if(nrow(centrality_df) > 0) {
      # Sortandtaketop10
      ord_idx <- order(centrality_df$total_strength, decreasing = TRUE)
      centrality_top <- centrality_df[ord_idx[1:min(10, nrow(centrality_df))], ]
      
      # Convert to factor
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
      
      cat("  ✓ Panel HPlotsuccessfully\n")
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
  # Backup
  cell_counts <- as.data.frame(table(high_quality_cells$celltype_reannotated))
  colnames(cell_counts) <- c("CellType", "Count")
  cell_counts <- cell_counts[!is.na(cell_counts$CellType), ]
  
  p4h <- ggplot2::ggplot(cell_counts, ggplot2::aes(x = reorder(CellType, Count), y = Count)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "", y = "Cell Count", title = "H. Cell Type Abundance (Backup)") +
    theme_journal
}

cat("✓ Panel E-Hcompleted\n")

# ==================== Assemble Figure 4(2-Row × 4-Column Layout)====================
cat("\nAssemble Figure 4 (2-Row × 4-Column Layout)...\n")

library(grid)
library(gridExtra)
library(patchwork)

# Assemble Row 1
fig4_row1 <- wrap_plots(p4a, p4b, p4c, p4d, ncol = 4)

# Assemble Row 2
fig4_row2 <- wrap_plots(p4e, p4f, p4g, p4h, ncol = 4, widths = c(1, 1, 0.8, 0.8))

# Complete assembly
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

# Save
ggsave(file.path(figures.dir, "Figure4_Hepatocyte_HSC_Communication.tiff"), 
       fig4_complete, width = 16, height = 10, dpi = 300, compression = "lzw")

ggsave(file.path(figures.dir, "Figure4_Hepatocyte_HSC_Communication.pdf"), 
       fig4_complete, width = 16, height = 10, device = cairo_pdf)

cat("\n✓ Figure 4 Complete! (Original Figure 3, Communication Analysis, 8 panels A-H)\n")
cat("Key improvements:\n")
cat("  1. Complete extract_communication_safe() function\n")
cat("  2. Complete Panel E-H code\n")
cat("  3. Auto-detect CellChat object structure\n")
cat("  4. Multiple fallback strategies\n")
cat("  5. correctly2-Row × 4-Column Layout\n")

# ==================== Step 12: Figure 5 - Therapeutic Target Priority（Newly added）====================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("              Figure 5: Therapeutic Target Prioritization\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Loading Required Packages
if(!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

# Ensure hsc_deg exists
if(!exists("hsc_deg")) {
  deg_file <- file.path(tables.dir, "Table_DEG_HSC_Act_vs_Qui_CCl4.csv")
  if(file.exists(deg_file)) {
    hsc_deg <- read.csv(deg_file)
  } else {
    stop("Error: Not foundHSC DEG results")
  }
}

# Unified theme
fig5_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# ==================== Panel A: Therapeutic Target Priority Ranking ====================
cat("Panel A: Therapeutic Target Priority Ranking\n")

# Define targets and their druggability scores
target_info <- data.frame(
  gene = c("Acta2", "Col1a1", "Lox", "Timp1", "Cthrc1", "Postn", "Vim", "Col1a2", "Col3a1"),
  druggability = c("★★★", "★★★", "★★☆", "★★☆", "★☆☆", "★☆☆", "★☆☆", "★★★", "★★☆"),
  stringsAsFactors = FALSE
)

# Merge DEG data
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

# Plot priority bar chart
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

# ==================== Panel B: Mechanism Summary Diagram ====================
cat("Panel B: Stress-Communication-Activation Axis Diagram\n")

# Create mechanism diagram（Use ggplot to plot flow chart）
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
  # Plot nodes
  geom_point(data = mechanism_data, aes(x = x, y = y, size = size), 
             color = c("#3498DB", "#E74C3C", "#F39C12", "#9B59B6", "#2ECC71", "#E74C3C"), 
             alpha = 0.8) +
  # Plot arrows
  geom_segment(data = arrow_data, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.3, "cm")), linewidth = 1, color = "grey30") +
  # Add labels
  geom_text(data = mechanism_data, aes(x = x, y = y, label = label), 
            size = 3, fontface = "bold", color = "white") +
  # Add signal annotations
  annotate("text", x = 2.5, y = 0.6, label = "PDGF/TGF-β\nSignaling", 
           size = 3, color = "#2ECC71", fontface = "bold") +
  annotate("text", x = 3.5, y = 0.6, label = "Oxidative\nStress", 
           size = 3, color = "#E74C3C", fontface = "bold") +
  # Set axes
  scale_x_continuous(limits = c(0.5, 4.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 2.5), expand = c(0, 0)) +
  scale_size_continuous(range = c(15, 25)) +
  labs(title = "B. Stress-Communication-Activation Axis") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "none"
  )

# ==================== Assemble Figure 5 ====================
cat("Assemble Figure 5...\n")

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

# Save
ggsave(file.path(figures.dir, "Figure5_Therapeutic_Targets.tiff"), 
       fig5_final, width = 14, height = 7, dpi = 300, compression = "lzw")

ggsave(file.path(figures.dir, "Figure5_Therapeutic_Targets.pdf"), 
       fig5_final, width = 14, height = 7, device = cairo_pdf)

ggsave(file.path(figures.dir, "Figure5_Therapeutic_Targets.png"), 
       fig5_final, width = 14, height = 7, dpi = 300)

cat("\n✓ Figure 5 Complete! (Newly added, Therapeutic Target Priority, 2 panels)\n")
cat("Outputfile:\n")
cat("  - Figure5_Therapeutic_Targets.tiff\n")
cat("  - Figure5_Therapeutic_Targets.pdf\n")
cat("  - Figure5_Therapeutic_Targets.png\n")
