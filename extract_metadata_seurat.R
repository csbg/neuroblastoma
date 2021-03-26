# Export metadata from a Seurat object, as well as individual assays.
#
# Metadata is exported to metadata_seurat.csv with columns
# * cell – barcode as used by Seurat
# * cluster_[res] – cluster IDs at resolution [res]
# * signature_[gs] – scores for gene signature [gs]
# * umap_1 ↓
# * umap_2 – UMAP coordinates
# * tsne_1 ↓
# * tsne_2 – tSNE coordinates
#
# The SCT assay is exported to assay_sct_seurat.rds.
#
# @DEPI rna_integrated_seurat.rds
# @DEPO metadata_seurat.csv


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
    Embeddings(nb, "umap") %>%
      as_tibble() %>%
      rename(umap_1_seurat = UMAP_1, umap_2_seurat = UMAP_2),
    Embeddings(nb, "tsne") %>%
      as_tibble() %>%
      rename(tsne_1_seurat = tSNE_1, tsne_2_seurat = tSNE_2)
  )

nb_metadata %>%
  write_csv(path_join(c(out_dir, str_glue("metadata_seurat.csv"))))
