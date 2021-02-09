# Integrate scRNA-seq datasets via alignment in monocle3 and export metadata.
#
# The intergrated dataset is exported to rna_integrated_monocle.rds.
#
# Metadata is exported to metadata_monocle.csv with columns
# * cell – barcode as used by Seurat
# * cluster – cluster ID
# * partition – partition ID
# * umap_1 ↓
# * umap_2 – UMAP coordinates
# * tsne_1 ↓
# * tsne_2 – tSNE coordinates

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
  preprocess_cds() %>% 
  align_cds(alignment_group = "sample") %>% 
  reduce_dimension(reduction_method = "UMAP", preprocess_method = "Aligned") %>% 
  reduce_dimension(reduction_method = "tSNE", preprocess_method = "Aligned") %>% 
  cluster_cells()

saveRDS(cds, path_join(c(out_dir, "rna_integrated_monocle.rds")))



# Export metadata ---------------------------------------------------------

monocle_metadata <- list(
  enframe(cds@clusters$UMAP$clusters, name = "cell", value = "cluster"),
  enframe(cds@clusters$UMAP$partitions, name = "cell", value = "partition"),
  reducedDim(cds, "UMAP") %>%
    as_tibble(rownames = "cell") %>%
    rename(umap_1 = V1, umap_2 = V2),
  reducedDim(cds, "tSNE") %>%
    as_tibble(rownames = "cell") %>%
    rename(tsne_1 = V1, tsne_2 = V2)
) %>%
  reduce(left_join, by = "cell")

monocle_metadata %>% write_csv(path_join(c(out_dir, "metadata_monocle.csv")))
