# Find conserved markers.

library(tidyverse)
library(Seurat)
library(log4r)



# Parameters --------------------------------------------------------------

infile <- "data_generated/all_datasets_sc/nb_integrated.rds"
outdir <- "data_generated/all_datasets_sc"

# infile <- "data_generated/3_datasets_sc/nb_integrated.rds"
# outdir <- "data_generated/3_datasets_sc"



# Load data ---------------------------------------------------------------

nb <- readRDS(infile)

nb@meta.data <-
  nb@meta.data %>% 
  rownames_to_column("rownames") %>% 
  left_join(
    read_csv("data_raw/metadata/sample_groups.csv", comment = "#"),
    by = "sample"
  ) %>% 
  column_to_rownames("rownames")

nb <- 
  nb %>%
  FindNeighbors() %>% 
  FindClusters(resolution = 0.5)

Idents(nb) <- nb@meta.data$seurat_clusters



# Find conserved markers --------------------------------------------------

logger <- logger()

save_conserved_markers <- function(cluster) {
  info(logger, str_glue("Finding conserved markers in cluster {cluster}"))
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
    write_csv(str_glue("{outdir}/consmarkers_{cluster}.csv"))
}

Idents(nb) %>%
  levels() %>%
  walk(save_conserved_markers)
