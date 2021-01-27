# Analysis of differential gene expression via Seurat::FindMarkers()

library(Seurat)
library(tidyverse)
library(fs)
library(log4r)


# Parameters --------------------------------------------------------------

infile <- "data_generated/all_datasets_current/nb_integrated.rds"
outdir <- "data_generated/all_datasets_current"



# Load data ---------------------------------------------------------------

logger <- logger()
nb <- readRDS(infile)
DefaultAssay(nb) <- "RNA"
nb <- NormalizeData(nb)

nb_data <- readRDS("data_generated/all_datasets_current/nb_all_metadata.rds")
nb@meta.data <-
  nb_data %>%
  column_to_rownames("cell")

Idents(nb) <- str_c(nb@meta.data$refined_cluster, nb@meta.data$group, sep = "_")



# Find markers ------------------------------------------------------------

find_markers <- function(cluster) {
  info(logger, str_glue("Finding DE genes in cluster {cluster}"))
  
  map_dfr(
    c("II", "III", "IV"),
    function(group) {
      info(logger, str_glue("{group} vs I"))
      
      group_1 <- str_glue("{cluster}_{group}")
      group_2 <- str_glue("{cluster}_I")
      if (!(group_1 %in% Idents(nb) && group_2 %in% Idents(nb)))
        return(tibble())
        
      markers <- suppressWarnings(
        FindMarkers(
          nb,
          ident.1 = group_1,
          ident.2 = group_2,
          test.use = "negbinom",
          min.pct = .25
        )
      )
      markers %>% 
        as_tibble(rownames = "gene") %>%
        mutate(cluster = cluster, group = group) %>%
        relocate(c(cluster, group), 0)
    }
  )
}

all_markers <- map_dfr(levels(nb@meta.data$refined_cluster), find_markers)



# Save data ---------------------------------------------------------------

write_csv(all_markers, path_join(c(outdir, "dge_seurat_groups_clusterwise.csv")))