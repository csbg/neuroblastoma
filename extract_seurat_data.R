# Extract various data from a Seurat object to facilitate analysis of large
# datasets when RAM is scarce.
#
# CSV files contain columns `cell` and `sample`:
# * nb_general.csv - molecule and feature counts, % mitochondrial genes
# * nb_tsne.csv - tSNE coordinates in columns `UMAP_1` and `UMAP_2`
# * nb_umap.csv - UMAP coordinates in columns `tSNE_1` and `tSNE_2`
# * nb_clusters_[res].csv - Cluster IDs in column `integrated_snn_res.[res]`,
#                           where [res] is 0.2, 0.5, or 0.8
# * nb_gene_signatures.csv - scores for NB gene signatures (see Table S2E from 
#                            Dong et al., Cancer Cell 38, 716–733 (2020))
#
# TRE files contain dendrograms in Newick tree format, which describe cluster
# similarities:
# * nb_dendrogram_[res].tre - Dendrogram for clusters at resolution [res]
#
# RDS files contain a Seurat object with a single assay:
# * nb_assay_RNA.rds - only the RNA assay
# * nb_assay_SCT.rds - only the SCT assay

library(Seurat)
library(tidyverse)
library(fs)
library(ape)



# Parameters --------------------------------------------------------------

infile <- "data_generated/all_datasets_current/nb_integrated.rds"
outdir <- "data_generated/all_datasets_current"



# Dimensional reduction and clustering ------------------------------------

nb <- readRDS(infile)

nb <-
  nb %>% 
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30) %>% 
  RunTSNE(dims = 1:30) %>% 
  FindNeighbors() %>% 
  FindClusters(resolution = 0.2) %>% 
  FindClusters(resolution = 0.5) %>% 
  FindClusters(resolution = 0.8)



# Extract general metadata ------------------------------------------------

nb@meta.data %>% 
  as_tibble(rownames = "cell") %>% 
  select(
    cell, sample,
    nCount_RNA, nFeature_RNA,
    nCount_SCT, nFeature_SCT,
    percent.mt
  ) %>% 
  write_csv(path_join(c(outdir, str_glue("nb_general.csv"))))



# Extract dimensional reductions ------------------------------------------

export_reduction <- function(reduction) {
  left_join(
    nb@meta.data %>% 
      select(sample) %>%
      as_tibble(rownames = "cell"),
    Embeddings(nb, reduction) %>%
      as_tibble(rownames = "cell"),
    by = "cell"
  ) %>%
    write_csv(path_join(c(outdir, str_glue("nb_{reduction}.csv"))))
}

export_reduction("tsne")
export_reduction("umap")



# Extract clusters and dendrograms ----------------------------------------

extract_cluster <- function(res) {
  cluster_col <- str_glue("integrated_snn_res.{res}")
  
  nb@meta.data %>%
    select(sample, {{cluster_col}}) %>%
    rownames_to_column("cell") %>% 
    write_csv(path_join(c(outdir, str_glue("nb_clusters_{res}.csv"))))
  
  Idents(nb) <- cluster_col
  nb <- BuildClusterTree(nb)
  Tool(nb, "BuildClusterTree") %>%
    write.tree(path_join(c(outdir, str_glue("nb_dendrogram_{res}.tre"))))
}

extract_cluster("0.2")
extract_cluster("0.5")
extract_cluster("0.8")



# Extract RNA and SCT assay -----------------------------------------------

extract_assay <- function(assay) {
  nb_assay <- 
    Assays(nb, assay) %>% 
    CreateSeuratObject(meta.data = nb@meta.data)
  nb_assay@reductions$umap <- nb@reductions$umap
  saveRDS(nb_assay, path_join(c(outdir, str_glue("nb_assay_{assay}.rds")))  )
}
 
extract_assay("RNA")
extract_assay("SCT")



# Extract gene signatures -------------------------------------------------

nb_signatures <- 
  read_csv("data_raw/metadata/nb_signatures.csv", comment = "#") %>% 
  group_by(cell_type) %>% 
  summarise(x = list(gene)) %>% 
  deframe()

nb <- AddModuleScore(nb, features = nb_signatures, assay = "SCT")

nb@meta.data %>% 
  select(sample, starts_with("Cluster")) %>%
  rename_with(~names(nb_signatures), .cols = starts_with("Cluster")) %>% 
  rownames_to_column("cell") %>% 
  write_csv("data_generated/all_datasets_current/nb_gene_signatures.csv")
