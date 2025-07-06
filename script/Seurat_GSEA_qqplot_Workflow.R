##########################################################
# Comprehensive GSEA & Q-Q Plot Analysis Pipeline for Tumor and Myeloid Cells
# Author: Qing Zhong
# Description: This script performs differential expression, GSEA, and Q-Q plot analysis
#              for tumor epithelial cells and myeloid cells using Seurat objects.
##########################################################

# ==== Load Required Packages ====
library(Seurat)
library(harmony)
library(clusterProfiler)
library(GseaVis)
library(ggplot2)
library(dplyr)
library(ggpubr)

# ==== 1. Tumor Epithelial Cells: Sequential vs Control GSEA Analysis ====

# Load tumor epithelial Seurat object
Epithelial <- readRDS("Step04_MajorAnnotated_subcluster_Tumor.RDS")

# Load gene sets from GMT file
Pathway_List <- clusterProfiler::read.gmt("/src/msigdb.v2023.2.Mm.symbols.gmt")

# Preprocess and integrate the tumor epithelial data
Epithelial <- Epithelial %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("percent.mt", "percent.rb", "CC.Difference")) %>%
  RunPCA(npcs = 30) %>%
  RunHarmony("orig.ident") %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = seq(0.05, 1, 0.05)) %>%
  RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony.minor_cluster")

# Visualize a marker gene on UMAP (optional)
FeaturePlot(Epithelial, reduction = "umap.harmony.minor_cluster", features = "cds-QBQ18417.1")

# Differential expression between Sequential and Control groups
DEG <- FindMarkers(Epithelial, ident.1 = "Sequential", ident.2 = "Control", group.by = "Group")

# Add gene symbols as column for GSEA preparation
DEG$SYMBOL <- rownames(DEG)

# Filter genes based on expression percentages and arrange by log fold change
DEG_fil <- DEG %>%
  filter(pct.1 > 0.2 | pct.2 > 0.2) %>%
  arrange(desc(avg_log2FC))

# Prepare named vector for GSEA (log fold changes named by gene symbol)
geneList <- DEG_fil$avg_log2FC
names(geneList) <- DEG_fil$SYMBOL

# Run GSEA analysis
gsea <- clusterProfiler::GSEA(geneList, TERM2GENE = Pathway_List)

# Save GSEA results for downstream use
save(gsea, file = "./GSEA_Out_Tumor.RData")

# Extract specific pathway groups for visualization
Apoptosis_associated <- grep(pattern = "APOP", x = gsea@result$ID, value = TRUE)
Immune <- gsea@result %>%
  filter(grepl("IMMUNE", ID)) %>%
  top_n(n = 5, wt = NES) %>%
  pull(ID)

# Prepare dataframe with pathway categories
df <- gsea@result %>%
  arrange(desc(NES)) %>%
  mutate(Pathway = case_when(
    ID %in% Apoptosis_associated ~ "Apoptosis_associated",
    ID == "HALLMARK_DNA_REPAIR" ~ "HALLMARK_DNA_REPAIR",
    ID %in% Immune ~ "Immune",
    TRUE ~ "Other"
  )) %>%
  mutate(RANK = row_number(),
         Label = ifelse(abs(NES) > 2 & Pathway != "Other", ID, ""))

# Export processed GSEA dataframe
write.csv(df, file = "./GSEA_out_Tumor_Cell.csv")

# Plot NES barplot for selected pathways
p_tumor_gsea <- df %>%
  filter(Pathway != "Other") %>%
  arrange(NES) %>%
  mutate(ID = factor(ID, levels = .$ID)) %>%
  ggplot(aes(y = ID, x = NES, fill = Pathway)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b3de69", "#80b1d3", "#9e9ac8")) +
  theme_bw() +
  labs(title = "Tumor Sequential vs Control DEG GSEA Analysis")

# Print tumor GSEA plot
print(p_tumor_gsea)


# ==== 2. Myeloid Cells: Sequential vs Concurrent GSEA Analysis ====

# Load myeloid cell subset Seurat object
Myeloid <- readRDS("Step04_MajorAnnotated_subcluster_Myeloid_Cell.RDS")

# Load gene sets for myeloid (GMT file)
Pathway_List_myeloid <- clusterProfiler::read.gmt("mh.all.v2024.1.Mm.symbols.gmt")

# Differential expression: Sequential vs Concurrent groups
DEG_myeloid <- FindMarkers(Myeloid,
                           ident.1 = "Sequential",
                           ident.2 = "Concurrent",
                           group.by = "Group",
                           logfc.threshold = 0)

# Add gene symbols as column
DEG_myeloid$SYMBOL <- rownames(DEG_myeloid)

# Filter genes expressed in at least one group (>20%)
DEG_fil_myeloid <- DEG_myeloid %>%
  filter(pct.1 > 0.2 | pct.2 > 0.2) %>%
  arrange(desc(avg_log2FC))

# Prepare named vector for GSEA
geneList_myeloid <- DEG_fil_myeloid$avg_log2FC
names(geneList_myeloid) <- DEG_fil_myeloid$SYMBOL

# Run GSEA with no p-value cutoff to retain all results
gsea_myeloid <- clusterProfiler::GSEA(geneList_myeloid,
                                      TERM2GENE = Pathway_List_myeloid,
                                      pvalueCutoff = 1)

# Dotplot visualization of top 10 pathways ordered by NES
p_myeloid_dotplot <- GseaVis::dotplotGsea(data = gsea_myeloid, topn = 10,
                                          order.by = 'NES', pajust = 1)

# Customize dotplot with ggplot2
p_myeloid_custom <- p_myeloid_dotplot$df %>%
  ggplot(aes(x = NES, y = Description, size = -log10(p.adjust), color = -log10(p.adjust))) +
  geom_point() +
  facet_grid(~type, scales = "free") +
  theme_bw() +
  scale_color_gradient2(low = "grey", mid = "grey", high = "red", midpoint = -log10(0.05)) +
  labs(title = "GSEA Dotplot: Sequential vs Concurrent (Myeloid)")

# Save the dotplot figure
ggsave(plot = p_myeloid_custom, filename = "./Fig/GSEA_SequentialVsConcurrent_Myeloid.pdf", width = 9.67, height = 5.84)

# Get pathways with adjusted p-value < 0.05 for enrichment plots
significant_pathways <- p_myeloid_dotplot$df %>%
  filter(p.adjust < 0.05) %>%
  pull(Description) %>%
  as.vector()

# Generate enrichment plots for significant pathways
GseaVis::gseaNb(object = gsea_myeloid,
                geneSetID = significant_pathways,
                subPlot = 2,
                addPval = TRUE,
                pvalX = -0.00,
                pvalY = -0.1,
                rmPrefix = FALSE,
                curveCol = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                             "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                             "#cab2d6", "#6a3d9a"))


# ==== 3. Q-Q Plot Analysis for Tumor Genes ====

# Load tumor epithelial Seurat object (reuse if needed)
Tumor <- readRDS("./Step04_MajorAnnotated_subcluster_Tumor.RDS")

# Define Q-Q plot function for comparing gene expression quantiles between two groups
qq_plot <- function(seurat_data, gene) {
  # Fetch expression vector for the gene
  gene_expr <- FetchData(seurat_data, vars = gene)[[gene]]
  
  # Extract group assignments
  group_data <- seurat_data@meta.data$Group
  
  # Subset gene expression by groups
  expr_Concurrent <- gene_expr[group_data == "Concurrent"]
  expr_Sequential <- gene_expr[group_data == "Sequential"]
  
  # Compute quantiles at 1% increments
  qq_data <- data.frame(
    Concurrent = quantile(expr_Concurrent, probs = seq(0, 1, by = 0.01), na.rm = TRUE),
    Sequential = quantile(expr_Sequential, probs = seq(0, 1, by = 0.01), na.rm = TRUE)
  )
  
  # Perform paired t-test
  t_test_result <- t.test(expr_Concurrent, expr_Sequential)
  cat("T-test result for gene:", gene, "\n")
  print(t_test_result)
  
  # Create Q-Q plot with ggplot2
  p <- ggplot(qq_data, aes(x = Concurrent, y = Sequential)) +
    geom_point(color = "red", size = 1.5) +
    geom_line(color = "red", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
    labs(title = paste0("Paired Q-Q Plot for ", gene, " Expression"),
         x = "Quantiles of Concurrent Group",
         y = "Quantiles of Sequential Group") +
    theme_minimal() +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
    annotate("text", x = 1, y = 3.5, label = paste0("p = ", signif(t_test_result$p.value, 3)))
  
  return(p)
}

# Run Q-Q plot for a list of significant chemokine genes (must define `signif_chemokine` before)
# Example:
# signif_chemokine <- c("GeneA", "GeneB", "GeneC")  # Replace with your gene list

for (gene in signif_chemokine) {
  p <- qq_plot(seurat_data = Tumor, gene = gene)
  ggsave(plot = p, filename = paste0("./Result_1/q-qplot_Tumor_", gene, ".pdf"), width = 3.64, height = 3.29)
}
