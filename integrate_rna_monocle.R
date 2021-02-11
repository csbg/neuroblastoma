# Integrate scRNA-seq datasets via alignment in monocle3 and export metadata.
#
# The intergrated dataset is exported to rna_integrated_monocle.rds.
#
# Metadata is exported to metadata_monocle.csv with columns
# * cell – barcode as used by Seurat
# * umap_1 ↓
# * umap_2 – UMAP coordinates
# * tsne_1 ↓
# * tsne_2 – tSNE coordinates
# * cluster_[k] – cluster IDs for different numbers k of nearest neighbors
# * partition_[k] – partition IDs for different k

library(monocle3) 
library(Seurat)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

# the merged dataset
merged_dataset <- "data_generated/rna_merged.rds"

# folder where results are saved
out_dir <- "data_generated"



# Integrate samples -------------------------------------------------------

nb <- readRDS(merged_dataset)

cds <-
  nb$RNA@counts %>% 
  new_cell_data_set(cell_metadata = nb@meta.data) %>% 
  preprocess_cds(verbose = TRUE) %>% 
  align_cds(alignment_group = "sample", verbose = TRUE) %>% 
  reduce_dimension(
    reduction_method = "UMAP",
    preprocess_method = "Aligned",
    verbose = TRUE
  ) %>% 
  reduce_dimension(
    reduction_method = "tSNE",
    preprocess_method = "Aligned",
    verbose = TRUE
  ) %>% 
  cluster_cells(k = 20, random_seed = 42, verbose = TRUE)

cds_50 <- cluster_cells(cds, k = 50, random_seed = 42, verbose = TRUE)

saveRDS(cds_50, path_join(c(out_dir, "rna_integrated_monocle.rds")))



# Export metadata ---------------------------------------------------------

monocle_metadata <- list(
  reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
    as_tibble(rownames = "cell"),
  reducedDim(cds, "tSNE") %>%
    magrittr::set_colnames(c("tsne_1", "tsne_2")) %>% 
    as_tibble(rownames = "cell"),
  clusters(cds) %>%
    enframe(name = "cell", value = "cluster_20"),
  partitions(cds) %>%
    enframe(name = "cell", value = "partition_20"),
  clusters(cds_50) %>% 
    enframe(name = "cell", value = "cluster_50"),
  partitions(cds_50) %>%
    enframe(name = "cell", value = "partition_50")
) %>%
  reduce(left_join, by = "cell")

monocle_metadata %>% write_csv(path_join(c(out_dir, "metadata_monocle.csv")))
