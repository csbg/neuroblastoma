# Perform subclustering on existing clusters.
# 
# For each given cluster in the complete dataset, this script performs
# dimensionality reduction (PCA, UMAP) followed by clustering (at two
# resolutions).
#
# Results are stored in a CSV file `nb_subclusters.csv` with columns
# * `cell`, `sample
# * `subcluster_0.2`, `subcluster_0.5` (subcluster assignments)
# * `subcluster_UMAP_1`, `subcluster_UMAP_2` (UMAP coordinates)

library(Seurat)
library(tidyverse)
library(fs)

rna_file <- "data_generated/all_datasets_current/nb_integrated.rds"
cluster_file <- "data_generated/all_datasets_current/nb_clusters_0.5.csv"
outdir <- "data_generated/all_datasets_current"



# Load data ---------------------------------------------------------------

nb <- readRDS(rna_file)
nb_clusters <-
  read_csv(cluster_file) %>%
  rename(cluster = integrated_snn_res.0.5)

nb@meta.data <- 
  nb@meta.data %>%
  as_tibble(rownames = "cell") %>% 
  left_join(nb_clusters, by = c("cell", "sample")) %>% 
  column_to_rownames("cell")



# Subclustering -----------------------------------------------------------

make_subcluster <- function(df, selected_cluster) {
  message("Subclustering cluster ", selected_cluster)
  nb_subset <-
    df %>%
    subset(subset = cluster == selected_cluster) %>% 
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE) %>% 
    FindNeighbors(verbose = FALSE) %>%
    FindClusters(resolution = 0.2, verbose = FALSE) %>%
    FindClusters(resolution = 0.5, verbose = FALSE)
  
  bind_cols(
    nb_subset@meta.data %>%
      as_tibble(rownames = "cell") %>% 
      select(cell, sample, starts_with("integrated")),
    Embeddings(nb_subset, "umap") %>%
      as_tibble()
  )
}

subclusters <-
  map_dfr(
    unique(nb_clusters$cluster),
    ~make_subcluster(nb, .x)
  )

subclusters %>%
  rename_with(
    ~c("subcluster_0.2", "subcluster_0.5",
       "subcluster_UMAP_1", "subcluster_UMAP_2"),
    .cols = !c(cell, sample)
  ) %>% 
  write_csv(path_join(c(outdir, str_glue("nb_subclusters.csv"))))



# Misc --------------------------------------------------------------------

nb_data %>% 
  count(cluster = integrated_snn_res.0.5, cell_type_broad) %>% 
  group_by(cluster) %>% 
  mutate(n_rel = n / sum(n) * 100) %>% 
  filter(n_rel >= 10) %>% 
  View()
