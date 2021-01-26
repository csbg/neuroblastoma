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



# Load data ---------------------------------------------------------------

metadata_superclusters <-
  left_join(
    read_csv("data_generated/all_datasets_current/nb_clusters_0.5.csv") %>% 
      select(cell, mid = integrated_snn_res.0.5),
    read_csv("data_generated/all_datasets_current/nb_singler.csv") %>% 
      extract(
        pruned.labels,
        into = "ctb",
        regex = "([^:]*)",
        remove = FALSE
      ) %>%
      mutate(
        ctb =
          as_factor(ctb) %>%
          fct_lump_prop(0.01) %>%
          fct_explicit_na()
      ) %>% 
      select(cell, sample, ctb),
    by = "cell"
  ) %>% 
  column_to_rownames("cell")

nb <- readRDS("data_generated/all_datasets_current/nb_integrated.rds")

outdir <- "data_generated/all_datasets_current"



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

make_all_subclusters <- function(supercluster) {
  nb@meta.data <- 
    metadata_superclusters %>% 
    select(sample, cluster = {{supercluster}})
  
  subclusters <-
    map_dfr(
      unique(nb@meta.data$cluster),
      ~make_subcluster(nb, .x)
    ) %>% 
    rename_with(
      ~c("subcluster_{s}_0.2",
         "subcluster_{s}_0.5",
         "subcluster_{s}_UMAP_1",
         "subcluster_{s}_UMAP_2") %>%
        map_chr(~str_glue_data(list(s = supercluster), .x)),
      .cols = !c(cell, sample)
    )
}

all_subclusters <- 
  colnames(metadata_superclusters) %>% 
  setdiff("sample") %>% 
  map(make_all_subclusters) %>% 
  reduce(left_join, by = c("cell", "sample"))

all_subclusters %>%
  write_csv(path_join(c(outdir, str_glue("nb_subclusters.csv"))))

