# Test of the overall data analysis workflow
# Not required for the actual data analysis

library(tidyverse)
library(Seurat)
library(SingleR)
library(celldex)
# library(future)
# plan(multisession) # don't do this!



# Datasets ----------------------------------------------------------------

# feature_matrix_file <- "data/R1_GNB_2014_0102_transcriptome/filtered_feature_bc_matrix.h5"
feature_matrix_file <- "data/MF220_NB_BM_Patient1_transcriptome/filtered_feature_bc_matrix.h5"

nb <-
  Read10X_h5(feature_matrix_file) %>% 
  CreateSeuratObject(
    project = "nb_R1_2014_0102",
    min.cells = 3,
    min.features = 200
  )


# QC ----------------------------------------------------------------------

nb <- PercentageFeatureSet(nb, pattern = "^MT-", col.name = "percent.mt")

nb@meta.data %>% 
  pivot_longer(!orig.ident) %>% 
  ggplot(aes(orig.ident, value)) +
  geom_violin(aes(fill = name), scale = "width") +
  geom_jitter(alpha = .2, size = 1) +
  facet_wrap(vars(name), scales = "free") +
  theme_bw() +
  NULL

# VlnPlot(nb, c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)
FeatureScatter(nb, "nCount_RNA", "nFeature_RNA") +
  FeatureScatter(nb, "nCount_RNA", "percent.mt")

nb <- subset(nb, nFeature_RNA > 200 & percent.mt < 5)  # between 200 and 2500
nb



# Normalization -----------------------------------------------------------

nb <- NormalizeData(nb)
nb



# Variable feature selection ----------------------------------------------

nb <- FindVariableFeatures(nb)
nb

VariableFeaturePlot(nb) %>%
  LabelPoints(VariableFeatures(nb)[1:10], repel = TRUE)



# Scaling -----------------------------------------------------------------

nb <- ScaleData(nb, features = rownames(nb))
nb



# PCA ---------------------------------------------------------------------

nb <- RunPCA(nb)

VizDimLoadings(nb)

DimPlot(nb, dims = c(1, 2))
# DimHeatmap(nb, dims = 1:15)

nb <- JackStraw(nb)
nb <- ScoreJackStraw(nb, dims = 1:20)
JackStrawPlot(nb, dims = 1:20)
ElbowPlot(nb, ndims = 50)

nb <-
  FindNeighbors(nb, dims = 1:20) %>% 
  FindClusters(resolution = .5)



# Non-linear dimension reduction ------------------------------------------

nb <-
  RunTSNE(nb, dims = 1:20) %>% 
  RunUMAP(dims = 1:20)
nb
nb

DimPlot(nb, reduction = "umap") + NoLegend() +
  DimPlot(nb, reduction = "tsne")



# Cell type classification ------------------------------------------------

reference_cell_types <- HumanPrimaryCellAtlasData()

predicted_cell_types <- SingleR(
  GetAssayData(nb),
  reference_cell_types,
  labels = reference_cell_types$label.main
)
predicted_cell_types

nb <- AddMetaData(nb, predicted_cell_types$labels, "SingleR_labels")

DimPlot(nb, reduction = "tsne", group.by = "SingleR_labels")

saveRDS(nb, "rds/nb.rds")
nb <- readRDS("rds/nb.rds")



# Biomarker detection -----------------------------------------------------

cluster1_markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1_markers)

cluster5_vs_1_3_markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5_vs_1_3_markers)

pbmc_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 2)

VlnPlot(pbmc, c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))


top10 <-
  pbmc_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)
top10

DoHeatmap(pbmc, top10$gene) + NoLegend()
