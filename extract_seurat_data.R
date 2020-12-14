# Extract various data from a Seurat object to facilitate analysis of large
# datasets when RAM is scarce
#
# Creates a series of CSV files with columns 'cell', 'sample', and various other
# columns, depending on the exported data
#
# Currently, the following files are generated:
# * nb_tsne.csv
# * nb_umap.csv
# * nb_clusters_0.5.csv
# * nb_clusters_0.8.csv

library(Seurat)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

infile <- "data_generated/3_datasets_sc/nb_integrated.rds"
outdir <- "data_generated/3_datasets_sc"



# Dimensional reduction and clustering ------------------------------------

nb <- readRDS(infile)

nb <-
  nb %>% 
  RunUMAP(dims = 1:30) %>% 
  RunTSNE(dims = 1:30)

nb <- 
  nb %>% 
  FindNeighbors() %>% 
  FindClusters(resolution = 0.5) %>% 
  FindClusters(resolution = 0.8)


# Extract data ------------------------------------------------------------

id_col <-
  nb@meta.data %>% 
  select(sample) %>%
  as_tibble(rownames = "cell")

export_embedding <- function(data, filename) {
  left_join(
    id_col,
    data,
    by = "cell"
  ) %>%
    write_csv(path_join(c(outdir, filename)))
}


Embeddings(nb, "tsne") %>%
  as_tibble(rownames = "cell") %>% 
  export_embedding("nb_tsne.csv")

Embeddings(nb, "umap") %>%
  as_tibble(rownames = "cell") %>% 
  export_embedding("nb_umap.csv")

nb@meta.data %>%
  select(sample, integrated_snn_res.0.5) %>%
  rownames_to_column("cell") %>% 
  write_csv(path_join(c(outdir, "nb_clusters_0.5.csv")))

nb@meta.data %>%
  select(sample, integrated_snn_res.0.8) %>%
  rownames_to_column("cell") %>% 
  write_csv(path_join(c(outdir, "nb_clusters_0.8.csv")))
