# WORK IN PROGRESS

library(tidyverse)
library(Seurat)


nb <- readRDS("data_generated/all_datasets_sc/nb_integrated.rds")
nb

nb@meta.data <- 
  nb@meta.data %>% 
  rownames_to_column("rownames") %>% 
  extract(sample, into = "sample", regex = "/[^_]*_(.*)") %>% 
  left_join(read_csv("data_raw/sample_groups.csv"), by = "sample") %>% 
  column_to_rownames("rownames")

nb <- 
  nb %>%
  FindNeighbors() %>% 
  FindClusters(resolution = 0.5)

# DimPlot(nb)

Idents(nb) <- nb@meta.data$seurat_clusters

save_conserved_markers <- function(cluster) {
  message("Finding conserved markers in cluster ", cluster)
  markers <-
    FindConservedMarkers(
      nb,
      ident.1 = cluster,
      grouping.var = "group",
      test.use = "negbinom",
      min.pct = 0.25
    ) %>% 
    as_tibble(rownames = "feature")
  
  markers %>%
    write_csv(str_glue("data_generated/all_datasets_sc/cmarkers_{cluster}.csv"))
  
  invisible(markers)
}

# save_conserved_markers("0")

Idents(nb) %>% levels() %>% walk(save_conserved_markers)
