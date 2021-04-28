# Perform subclustering on existing clusters.
# 
# For each given cluster_50 in the complete dataset, this script performs
# PCA, alignment, UMAP, and clustering at two resolutions.
#
# Results are stored in a CSV file `subclusters.csv` with columns
# * `cell`
# * `umap_1_subcluster`
# * `umap_2_subcluster` (UMAP coordinates)
# * `subcluster_20`
# * `subcluster_50` (subcluster assignments, lowercase letters)
#
# @DEPI rna_integrated_monocle.rds
# @DEPO subclusters.csv

library(monocle3)
library(tidyverse)
source("common_functions.R")



# Parameters --------------------------------------------------------------

# integrated dataset
in_file <- "data_generated/rna_integrated_monocle.rds"

# folder where results are saved
out_file <- "data_generated/subclusters.csv"



# Load data ---------------------------------------------------------------

nb <- readRDS(in_file)



# Subclustering -----------------------------------------------------------

#' Perform subclustering.
#'
#' @param data A monocle object.
#' @param selected_cluster Cluster in `clusters(data)`
#'                         that should be subclustered.
#'
#' @return A dataframe.
make_subcluster <- function(data, selected_cluster) {
  info("Subclustering cluster {selected_cluster}")
  
  barcodes <- 
    clusters(nb) %>%
    enframe("cell", "id") %>%
    filter(id == {{selected_cluster}}) %>%
    pull(cell)
  
  subdata_k20 <-
    data[, barcodes] %>% 
    preprocess_cds() %>% 
    reduce_dimension(preprocess_method = "PCA") %>% 
    align_cds(alignment_group = "sample") %>% 
    reduce_dimension(
      reduction_method = "UMAP",
      preprocess_method = "Aligned"
    ) %>% 
    cluster_cells(k = 20, random_seed = 42)
  
  subdata_k50 <- cluster_cells(subdata_k20, k = 50, random_seed = 42)
  
  list(
    reducedDim(subdata_k20, "UMAP") %>%
      magrittr::set_colnames(c("umap_1_subcluster", "umap_2_subcluster")) %>% 
      as_tibble(rownames = "cell"),
    clusters(subdata_k20) %>%
      enframe(name = "cell", value = "subcluster_20"),
    clusters(subdata_k50) %>%
      enframe(name = "cell", value = "subcluster_50")
  ) %>%
    reduce(left_join, by = "cell") %>% 
    mutate(
      across(where(is.factor), ~fct_relabel(., ~letters[as.integer(.x)]))
    )
}

# make_subcluster(nb, "1")

subclusters <- map_dfr(
  levels(clusters(nb)),
  make_subcluster,
  data = nb
)

subclusters %>% write_csv(out_file)
