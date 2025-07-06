############################################################
# GSVA Analysis on T Cell Data with Multiple Gene Set Collections
# Author: Qing Zhong
############################################################

# === Load Required Libraries ===
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(GSVA)
library(ComplexHeatmap)

# === 1. Load Seurat Object ===
tcell_obj <- readRDS("step05_Annotated_Tcell_Integrated.RDS")

# === 2. Define Path to Gene Set Collections (folders containing GMT files) ===
gene_set_root <- "/src/"
gene_set_folders <- list.files(gene_set_root)

# === 3. Load Gene Sets from GMT Files into a List ===
gene_sets <- list()

for (folder in gene_set_folders) {
  # Get full paths to GMT files in the current folder
  gmt_files <- list.files(file.path(gene_set_root, folder),
                          pattern = "\\.gmt$", full.names = TRUE)
  
  # Read all GMT files and combine into one data frame
  combined_gmt <- do.call(rbind, lapply(gmt_files, clusterProfiler::read.gmt))
  
  # Split combined GMT data by gene set term to create named list of gene vectors
  gene_sets[[folder]] <- split(combined_gmt$gene, combined_gmt$term)
}

# === 4. Aggregate Expression Data by 'Group' Metadata ===
agg_expr <- AggregateExpression(tcell_obj, group.by = "Group")$RNA
agg_expr <- as.data.frame(agg_expr)

# === 5. Perform GSVA for Each Gene Set Collection ===
gsva_results <- list()

for (name in names(gene_sets)) {
  gsva_results[[name]] <- gsva(as.matrix(agg_expr),
                               gene_sets[[name]],
                               kcdf = "Poisson")  # Poisson suited for raw counts
}

# === 6. Combine All GSVA Results into One Data Frame for Visualization ===
combined_df <- bind_rows(lapply(names(gsva_results), function(name) {
  df_tmp <- as.data.frame(gsva_results[[name]])
  df_tmp$GeneSetCollection <- name
  return(df_tmp)
}), .id = NULL)

# === 7. Plot Heatmap of GSVA Scores (Example: first 5 groups) ===
Heatmap(as.matrix(combined_df[, 1:5]),
        left_annotation = rowAnnotation(GeneSetCollection = combined_df$GeneSetCollection),
        cluster_rows = FALSE,
        show_row_names = FALSE,
        name = "GSVA Score")

# === 8. Save Individual GSVA Results as CSV Files ===
output_dir <- "./Out"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(gsva_results)) {
  write.csv(as.data.frame(gsva_results[[name]]),
            file = file.path(output_dir, paste0(name, ".csv")),
            row.names = TRUE)
}
