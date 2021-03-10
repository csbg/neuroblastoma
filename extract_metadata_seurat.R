# Export metadata from a Seurat object, as well as individual assays.
#
# Metadata is exported to metadata_seurat.csv with columns
# * cell – barcode as used by Seurat
# * sample – patient ID
# * percent_mt – percentage of mitochondrial genes
# * n_count_rna – number of molecules in RNA assay
# * n_feature_rna – number of features in RNA assay
# * n_count_sct – number of molecules in SCT assay
# * n_feature_sct – number of features in SCT assay
# * cluster_[res] – cluster IDs at resolution [res]
# * signature_[gs] – scores for gene signature [gs]
# * umap_1 ↓
# * umap_2 – UMAP coordinates
# * tsne_1 ↓
# * tsne_2 – tSNE coordinates
#
# Individual assays are exported to
# * assay_RNA_seurat.rds and
# * assay_SCT_seurat.rds
#
# @DEPI rna_integrated_seurat.rds
# @DEPO metadata_seurat.csv
# @DEPO assay_rna_seurat.rds
# @DEPO assay_sct_seurat.rds


library(Seurat)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

integrated_data <- "data_generated/rna_integrated_seurat.rds"
out_dir <- "data_generated"



# Calculations ------------------------------------------------------------

nb_signatures <- 
  read_csv("metadata/nb_signatures.csv", comment = "#") %>% 
  group_by(cell_type) %>% 
  summarise(x = list(gene)) %>% 
  deframe()

nb <- readRDS(integrated_data)

nb <-
  nb %>% 
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30) %>% 
  RunTSNE(dims = 1:30) %>% 
  FindNeighbors() %>% 
  FindClusters(resolution = 0.2) %>% 
  FindClusters(resolution = 0.5) %>% 
  FindClusters(resolution = 0.8) %>% 
  AddModuleScore(features = nb_signatures, assay = "SCT")



# Export metadata ---------------------------------------------------------

nb_metadata <- 
  nb@meta.data %>% 
  as_tibble(rownames = "cell") %>%
  select(
    cell,
    sample,
    percent_mt = percent.mt,
    n_count_rna = nCount_RNA,
    n_feature_rna = nFeature_RNA,
    n_count_sct = nCount_SCT,
    n_feature_sct = nFeature_SCT,
    cluster_0.2 = integrated_snn_res.0.2,
    cluster_0.5 = integrated_snn_res.0.5,
    cluster_0.8 = integrated_snn_res.0.8,
    starts_with("Cluster")
  ) %>% 
  rename_with(
    ~str_c("signature_", names(nb_signatures)),
    .cols = matches("Cluster\\d+")
  ) %>% 
  bind_cols(
    Embeddings(nb, "umap") %>% as_tibble() %>% rename_with(str_to_lower),
    Embeddings(nb, "tsne") %>% as_tibble() %>% rename_with(str_to_lower)
  )

nb_metadata %>%
  write_csv(path_join(c(out_dir, str_glue("metadata_seurat.csv"))))



# Export assays -----------------------------------------------------------

export_assay <- function(assay) {
  nb_assay <- 
    Assays(nb, assay) %>% 
    CreateSeuratObject(meta.data = nb@meta.data)
  nb_assay@reductions$umap <- nb@reductions$umap
  saveRDS(nb_assay, path_join(c(out_dir, str_glue("assay_{assay}_seurat.rds"))))
}
 
export_assay("RNA")
export_assay("SCT")
