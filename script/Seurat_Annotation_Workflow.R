##########################################################
# Seurat Major and Minor Annotation Pipeline
# Author: Qing Zhong
##########################################################
# ==== Load Required Packages ====
library(Seurat)
library(tidyverse)
library(harmony)
library(aplot)

# ==== Step 1: Load Integrated Data ====
load("step03_MergeAndIntegrate.RData") 

# ==== Step 2: Define Marker Lists ====
marker_lists <- list(
  Tumor = c("cds-QBQ18417.1"),
  Stroma = c("Col1a1","Krt18"),
  MoMac = c("Cd68","Adgre1","Cd14"),
  T = c("Cd2","Cd3e","Cd8a","Klrk1"),
  B = c("Cd19"),
  Plasma = c("Mzb1"),
  Neutrophil = c("Csf3r", "S100a8", "S100a9")
)

# ==== Step 3: Major Cluster Annotation ====
obj_Merge <- obj_Merge %>%
  mutate(
    major_cluster_1 = case_when(
      RNA_snn_res.0.5 %in% c(1,2,5,11,18) ~ "Tumor",
      RNA_snn_res.0.5 %in% c(16,19) ~ "Stroma",
      RNA_snn_res.0.5 %in% c(0,4,14,9,7) ~ "Myeloid_Cell",
      RNA_snn_res.0.5 == 15 ~ "Neutrophil",
      RNA_snn_res.0.5 %in% c(3,6,10,13) ~ "T_NK_Cell",
      RNA_snn_res.0.5 == 8 ~ "Plasma_Cell",
      RNA_snn_res.0.5 == 12 ~ "B_Cell",
      RNA_snn_res.0.5 == 17 ~ "DP_CD68_CD2",
      TRUE ~ "Undefined"
    ),
    Group = recode(orig.ident,
                   "2310057_B3" = "Control",
                   "2310057_A3" = "Olaparib",
                   "2310057_B2" = "AZD1775",
                   "2310057_B1" = "Concurrent",
                   "2310057_A4" = "Sequential")
  )

# Save major annotated object
save(obj_Merge, file = "step04_Major_Annotated.RData")

# Save subclusters by major cell type except DP_CD68_CD2
unique(obj_Merge$major_cluster_1) %>%
  discard(~ .x == "DP_CD68_CD2") %>%
  walk(~ saveRDS(obj_Merge[, obj_Merge$major_cluster_1 == .x],
                 file = paste0("Step04_MajorAnnotated_subcluster_", .x, ".RDS")))

# ==== Step 4: Minor Cluster Annotation Function ====
minor_annotation_pipeline <- function(seurat_obj, dims = 1:30, resolutions = seq(0.05,1,0.05),
                                      features = c("percent.mt","percent.rb","CC.Difference"),
                                      umap_name = "umap.harmony.minor_cluster") {
  seurat_obj %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData(vars.to.regress = features) %>%
    RunPCA(npcs = max(dims)) %>%
    RunHarmony("orig.ident") %>%
    FindNeighbors(reduction = "harmony", dims = dims) %>%
    FindClusters(resolution = resolutions) %>%
    RunUMAP(dims = dims, reduction = "harmony", reduction.name = umap_name)
}

# ==== Step 5: QC Visualization ====
plot_QC <- function(obj){
  p <- VlnPlot(obj,
               features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
               group.by = "orig.ident",
               ncol = 1,
               pt.size = 0.1,
               alpha = 0.1)
  ggsave(filename = "QC.pdf", width = 6.79 * 3, height = 8.95 * 3, plot = p, device = "pdf", units = "cm")
}

# ==== Step 6: Integrate Minor Annotations and Plot Major, T cell, Myeloid Annotations ====

merge_minor_annotation <- function(obj_merge, T_Cell, Mac){
  Mac_Meta <- Mac@meta.data %>%
    mutate(sub_major = case_when(
      grepl("Mac", minor_Cluster) ~ "Mac",
      grepl("Mon", minor_Cluster) ~ "Monocyte",
      grepl("DC", minor_Cluster) ~ "DC",
      TRUE ~ "Mac"
    )) %>%
    mutate(barcodes = rownames(.)) %>%
    select(barcodes, sub_major, minor_Cluster)
  
  T_Cell_Meta <- T_Cell@meta.data %>%
    mutate(barcodes = rownames(.), sub_major = sub_major_6) %>%
    select(barcodes, sub_major, minor_Cluster)
  
  obj_merge@meta.data <- obj_merge@meta.data %>%
    mutate(barcodes = rownames(.)) %>%
    left_join(rbind(T_Cell_Meta, Mac_Meta)) %>%
    mutate(sub_major = ifelse(is.na(sub_major), major_cluster_1, sub_major),
           minor_Cluster = ifelse(is.na(minor_Cluster), major_cluster_1, minor_Cluster)) %>%
    tibble::column_to_rownames("barcodes")
  
  save(obj_merge, file = "Step06_Annotated_Minor_Integrated.RData")
  return(obj_merge)
}

# Major annotation DotPlot + UMAP
plot_major_annotation <- function(obj_merge){
  obj_merge <- obj_merge[, !grepl("DP", obj_merge$sub_major)]
  
  marker_lists2 <- list(
    Tumor = c("cds-QBQ18417.1"),
    Stroma = c("Col1a1", "Krt18"),
    MoMac = c("Cd68", "Adgre1", "Cd14"),
    T = c("Cd2", "Cd3e", "Cd8a", "Cd4", "Klrk1"),
    B = c("Cd19"),
    Plasma = c("Mzb1"),
    Neutrophil = c("Csf3r", "S100a8", "S100a9")
  )
  
  plot_colors <- c("#ffdb86", "#FAB562", "#F08A84", "#fcc5c0", "#F9CDE4", "#BE91BE",
                   "#F5FBB4", "#ccb38c", "#D8C965", "#C2D567", "#D8C965", "#D0EBBD",
                   "#8DD3C7", "#afeeee", "#D9D7C9", "#C8B291", "#88ADCC", "#C8C5D3")
  
  p <- DotPlot(obj_merge,
               assay = "RNA",
               features = unlist(marker_lists2),
               group.by = "major_cluster_1")
  
  level <- c("Tumor", "Stroma", "Myeloid_Cell", "T_NK_Cell", "B_Cell", "Plasma_Cell", "Neutrophil")
  
  p2 <- data.frame(id = factor(p[["data"]]$id, levels = rev(level)),
                   features = p[["data"]]$features.plot,
                   avg.exp = p[["data"]]$avg.exp,
                   pct.exp = p[["data"]]$pct.exp,
                   avg.exp.scaled = p[["data"]]$avg.exp.scaled) %>%
    ggplot(aes(x = id, y = features, color = avg.exp.scaled)) +
    geom_point(aes(size = pct.exp)) +
    scale_color_gradient(low = "#ffffd9", high = "#1d91c0") +
    scale_y_discrete(position = "right") +
    labs(x = "", y = "", size = "Percent Expressed", color = "Average Expression") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour = "black", size = 15),
          panel.border = element_rect(color = "black", size = 0.7, fill = NA),
          legend.text = element_text(colour = "black", size = 13),
          legend.title = element_text(colour = "black", size = 13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p1 <- data.frame(id = p[["data"]]$id) %>%
    mutate(features = "Cell Type") %>%
    ggplot(aes(x = id, y = features)) +
    geom_tile(aes(fill = id), color = "white", linewidth = 1) +
    scale_x_discrete("", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    scale_fill_manual(values = plot_colors) +
    labs(fill = "Cell Type") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(colour = "black", size = 13),
          legend.title = element_text(colour = "black", size = 13))
  
  p3 <- p2 %>% insert_top(p1, height = 0.025)
  
  ggsave(filename = "Fig_Major_Annotation.pdf", width = 6.75 * 3, height = 5.72 * 3,
         plot = p3, device = "pdf", units = "cm")
  
  p_umap <- DimPlot(obj_merge, group.by = "major_cluster_1", reduction = "umap.harmony",
                    label = TRUE, raster = FALSE) +
    scale_color_manual(values = plot_colors) +
    NoLegend()
  
  ggsave(filename = "Fig_Major_UMAP.pdf", width = 7.53 * 3, height = 7.53 * 3,
         plot = p_umap, device = "pdf", units = "cm")
}

# T cell annotation plot
plot_T_cell_annotation <- function(T_Cell){
  T_Cell <- T_Cell[, !grepl("DP", T_Cell$minor_Cluster)]
  
  marker <- rev(c("Cd2", "Cd3e", "Cd4", "Cd8a", "Foxp3", "Ifng", "Il4", "Cxcr5",
                  "Tcf7", "Ccr7", "Sell", "S1pr1", "Mki67", "Gzmk", "Hspa1a",
                  "Dnajb1", "Havcr2", "Isg15", "Xcl1", "Klrk1", "Trdv4"))
  
  plot_colors <- c("#ffdb86", "#FAB562", "#F08A84", "#fcc5c0", "#F9CDE4", "#BE91BE",
                   "#F5FBB4", "#ccb38c", "#D8C965", "#C2D567", "#D8C965", "#D0EBBD",
                   "#8DD3C7", "#afeeee", "#D9D7C9", "#C8B291", "#88ADCC", "#C8C5D3")
  
  p <- DotPlot(T_Cell,
               assay = "RNA",
               features = marker,
               group.by = "minor_Cluster")
  
  p2 <- data.frame(id = p[["data"]]$id,
                   features = p[["data"]]$features.plot,
                   avg.exp = p[["data"]]$avg.exp,
                   pct.exp = p[["data"]]$pct.exp,
                   avg.exp.scaled = p[["data"]]$avg.exp.scaled) %>%
    ggplot(aes(x = id, y = features, color = avg.exp.scaled)) +
    geom_point(aes(size = pct.exp)) +
    scale_color_gradient(low = "#ffffd9", high = "#1d91c0") +
    scale_y_discrete(position = "right") +
    labs(x = "", y = "", size = "Percent Expressed", color = "Average Expression") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour = "black", size = 15),
          panel.border = element_rect(color = "black", size = 0.7, fill = NA),
          legend.text = element_text(colour = "black", size = 13),
          legend.title = element_text(colour = "black", size = 13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p1 <- data.frame(id = p[["data"]]$id) %>%
    mutate(features = "Cell Type") %>%
    ggplot(aes(x = id, y = features)) +
    geom_tile(aes(fill = id), color = "white", linewidth = 1) +
    scale_x_discrete("", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    scale_fill_manual(values = plot_colors) +
    labs(fill = "Cell Type") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(colour = "black", size = 13),
          legend.title = element_text(colour = "black", size = 13))
  
  p3 <- p2 %>% insert_top(p1, height = 0.025)
  
  ggsave(filename = "Fig_T_NK_Annotation.pdf", width = 7.91 * 3, height = 7.65 * 3,
         plot = p3, device = "pdf", units = "cm")
  
  p_umap <- DimPlot(T_Cell, group.by = "minor_Cluster",
                    reduction = "umap.harmony.minor_cluster",
                    label = TRUE, raster = FALSE) +
    scale_color_manual(values = plot_colors) +
    NoLegend()
  
  ggsave(filename = "Fig_T_NK_UMAP.pdf", width = 7.53 * 3, height = 7.53 * 3,
         plot = p_umap, device = "pdf", units = "cm")
}

# Myeloid annotation plot
plot_myeloid_annotation <- function(Mac){
  marker <- rev(c("Adgre1", "Cd14", "Xcr1", "Slc2a1", "Spp1", "Marco", "Mrc1",
                  "Selenop", "Tnfrsf9", "Mki67", "Hspa1a", "Dnajb1"))
  
  plot_colors <- c("#ffdb86", "#FAB562", "#F08A84", "#fcc5c0", "#F9CDE4", "#BE91BE",
                   "#F5FBB4", "#ccb38c", "#D8C965", "#C2D567", "#D8C965", "#D0EBBD",
                   "#8DD3C7", "#afeeee", "#D9D7C9", "#C8B291", "#88ADCC", "#C8C5D3")
  
  p <- DotPlot(Mac,
               assay = "RNA",
               features = marker,
               group.by = "minor_Cluster")
  
  p2 <- data.frame(id = p[["data"]]$id,
                   features = p[["data"]]$features.plot,
                   avg.exp = p[["data"]]$avg.exp,
                   pct.exp = p[["data"]]$pct.exp,
                   avg.exp.scaled = p[["data"]]$avg.exp.scaled) %>%
    ggplot(aes(x = id, y = features, color = avg.exp.scaled)) +
    geom_point(aes(size = pct.exp)) +
    scale_color_gradient(low = "#ffffd9", high = "#1d91c0") +
    scale_y_discrete(position = "right") +
    labs(x = "", y = "", size = "Percent Expressed", color = "Average Expression") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour = "black", size = 15),
          panel.border = element_rect(color = "black", size = 0.7, fill = NA),
          legend.text = element_text(colour = "black", size = 13),
          legend.title = element_text(colour = "black", size = 13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p1 <- data.frame(id = p[["data"]]$id) %>%
    mutate(features = "Cell Type") %>%
    ggplot(aes(x = id, y = features)) +
    geom_tile(aes(fill = id), color = "white", linewidth = 1) +
    scale_x_discrete("", expand = c(0, 0)) +
    scale_y_discrete("", expand = c(0, 0)) +
    scale_fill_manual(values = plot_colors) +
    labs(fill = "Cell Type") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(colour = "black", size = 13),
          legend.title = element_text(colour = "black", size = 13))
  
  p3 <- p2 %>% insert_top(p1, height = 0.025)
  
  ggsave(filename = "Fig_Myeloid_Annotation.pdf", width = 6.95 * 3, height = 5.70 * 3,
         plot = p3, device = "pdf", units = "cm")
  
  p_umap <- DimPlot(Mac, group.by = "minor_Cluster",
                    reduction = "umap.harmony.minor_cluster",
                    label = TRUE, raster = FALSE) +
    scale_color_manual(values = plot_colors) +
    NoLegend()
  
  ggsave(filename = "Fig_Myeloid_UMAP.pdf", width = 7.53 * 3, height = 7.53 * 3,
         plot = p_umap, device = "pdf", units = "cm")
}

# ==== Example main workflow call ====

# QC plot
plot_QC(obj_Merge)

T_Cell <- readRDS("Step05_MajorAnnotated_subcluster_T_NK_Cell.RDS")
Mac <- readRDS("Step05_MajorAnnotated_subcluster_Myeloid_Cell.RDS")

obj_Merge <- merge_minor_annotation(obj_Merge, T_Cell, Mac)

# Major annotation plot
plot_major_annotation(obj_Merge)
# T cell and Macrophage minor annotation plot
plot_T_cell_annotation(T_Cell)

plot_myeloid_annotation(Mac)

