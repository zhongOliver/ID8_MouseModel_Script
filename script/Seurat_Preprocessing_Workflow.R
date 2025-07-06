##############################################################
# Supplementary Code: Single-cell RNA-seq Preprocessing
# Description: Seurat v5 preprocessing, filtering, doublet removal,
# integration and annotation workflow
# Author:Qing Zhong
##############################################################

# === Part 1: Load Required Packages ===
library(Seurat)
library(DoubletFinder)
library(patchwork)
library(dplyr)
library(cowplot)
library(harmony)
library(SingleR)
library(scater)
library(future)
library(parallel)

# === Part 2: Load and Create Seurat Objects ===
folder_path <- "./Raw_Data/Out_Matrix/"
folderList <- list.files(folder_path, full.names = FALSE)
obj_01.list <- lapply(folderList, function(folder) {
  data <- Read10X(data.dir = file.path(folder_path, folder, "filtered_feature_bc_matrix"))
  seu <- CreateSeuratObject(counts = data, project = folder, min.cells = 3, min.features = 200)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^Rp[sl]")
  return(seu)
})
names(obj_01.list) <- folderList
save(obj_01.list, file = "step01_Dataimport.RData")

# === Part 3: Quality Control Filtering ===
obj_02.list <- list()
for (i in seq_along(obj_01.list)) {
  seu <- obj_01.list[[i]]
  upperThreshold <- quantile(seu$nFeature_RNA, 0.95)
  seu <- subset(seu, subset = nFeature_RNA > 300 & nFeature_RNA < upperThreshold & percent.mt < 20)
  obj_02.list[[i]] <- seu
}
names(obj_02.list) <- names(obj_01.list)
save(obj_02.list, file = "step02_DataFilter_01mt_filter.RData")

# === Part 4: Cell Cycle Correction ===
cc.genes <- readRDS("path/to/mmu_cycling_gene.list.RDS")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

CellCycle_rm.FUN <- function(seu) {
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  seu <- ScaleData(seu, vars.to.regress = c("percent.mt", "percent.rb", "CC.Difference"))
  seu <- RunPCA(seu)
  return(seu)
}

obj_03.list <- lapply(obj_02.list, CellCycle_rm.FUN)
names(obj_03.list) <- names(obj_02.list)
save(obj_03.list, file = "step02_DataFilter_02CCdifferenceCalculate.RData")

# === Part 5: Doublet Removal ===
load("ref_Mouse_all.RData")
paramSweep_v3_Seurat5<-function (seu, PCs = 1:10, sct = FALSE, num.cores = 1){ # nolint # nolint
  require(Seurat)
  require(fields)
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                 10000, replace = FALSE)]
    data <- seu@assays$RNA$counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts
    n.real.cells <- ncol(data)
  }
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct, mc.cores = num.cores)
    stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep, 
                      n.real.cells, real.cells, pK, pN, data, orig.commands, 
                      PCs, sct)
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                  sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

doubletFinder_v3_SeuratV5 <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                                       sct = FALSE, annotations = NULL){
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                 model.use = orig.commands$ScaleData.RNA$model.use, 
                                 do.scale = orig.commands$ScaleData.RNA$do.scale, 
                                 do.center = orig.commands$ScaleData.RNA$do.center, 
                                 scale.max = orig.commands$ScaleData.RNA$scale.max, 
                                 block.size = orig.commands$ScaleData.RNA$block.size, 
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                              npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                             ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbor_types[i, ] <- table(doublet_types1[neighbors - 
                                                        n_real.cells]) + table(doublet_types2[neighbors - 
                                                                                                n_real.cells])
          neighbor_types[i, ] <- neighbor_types[i, ]/sum(neighbor_types[i, 
          ])
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                    1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) = levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", 
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, 
                                                                              ct]
      }
    }
    return(seu)
  }
}

DoubleF.FUN <- function(ref_data, seu) {
  seu <- RunUMAP(seu, dims = 1:20)
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- mapping.FUN(ref_data, seu)

  sweep.res.list <- paramSweep_v3_Seurat5(seu, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- seu@meta.data$singleR_major_cell_clusters
  homotypic.prop <- modelHomotypic(annotations)

  nExp <- round(ncol(seu) * 8e-6 * nrow(seu@meta.data))
  nExp.adj <- round(nExp * (1 - homotypic.prop))

  seu <- doubletFinder_v3_SeuratV5(seu, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp)
  seu <- doubletFinder_v3_SeuratV5(seu, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp.adj, reuse.pANN = "pANN_0.25")

  seu <- subset(seu, subset = DF_hi.lo_fin == "Singlet")
  return(seu)
}

obj_04.list <- lapply(obj_03.list, function(seu) DoubleF.FUN(ref_Mouse, seu))
names(obj_04.list) <- names(obj_03.list)
save(obj_04.list, file = "step02_DataFilter_03DoubletFilter.RData")

# === Part 6: Merge and Integrate ===
obj_Merge <- do.call(merge, obj_04.list)
obj_Merge <- JoinLayers(obj_Merge)
obj_Merge <- NormalizeData(obj_Merge)
obj_Merge <- FindVariableFeatures(obj_Merge, nfeatures = 3000)
obj_Merge <- ScaleData(obj_Merge, vars.to.regress = c("percent.mt", "percent.rb", "CC.Difference"))
obj_Merge <- RunPCA(obj_Merge, npcs = 15)
obj_Merge <- RunHarmony(obj_Merge, group.by.vars = "orig.ident")
obj_Merge <- RunUMAP(obj_Merge, dims = 1:15, reduction = "harmony")
obj_Merge <- FindNeighbors(obj_Merge, reduction = "harmony", dims = 1:15)
obj_Merge <- FindClusters(obj_Merge, resolution = seq(0.05, 1, 0.05))
save(obj_Merge, file = "step03_MergeAndIntegrate.RData")

# === Part 7: Cell Type Annotation ===
Tumor_Marker <- c("cds-QBQ18417.1")
Stroma_Marker <- c("Col1a1", "Krt18")
MoMac_Marker <- c("Cd68", "Adgre1", "Cd14")
T_Marker <- c("Cd2", "Cd3e", "Cd8a", "Klrk1")
B_Marker <- c("Cd19")
Plasma_Marker <- c("Mzb1")
Neutrophil <- c("S100a8", "S100a9")

obj_Merge$major_cluster_1 <- case_when(
  obj_Merge$RNA_snn_res.0.5 %in% c(1,2,5,11,18) ~ "Tumor",
  obj_Merge$RNA_snn_res.0.5 %in% c(16,19) ~ "Stroma",
  obj_Merge$RNA_snn_res.0.5 %in% c(0,4,14,9,7) ~ "Myeloid_Cell",
  obj_Merge$RNA_snn_res.0.5 %in% c(15) ~ "Neutrophil",
  obj_Merge$RNA_snn_res.0.5 %in% c(3,6,10,13) ~ "T_NK_Cell",
  obj_Merge$RNA_snn_res.0.5 %in% c(8) ~ "Plasma_Cell",
  obj_Merge$RNA_snn_res.0.5 %in% c(12) ~ "B_Cell",
  obj_Merge$RNA_snn_res.0.5 == 17 ~ "DP_CD68_CD2",
  TRUE ~ "Undefined"
)
save(obj_Merge, file = "step04_Major_Annotated.RData")
