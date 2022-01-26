# The intergrated dataset is exported to rna_integrated_monocle.rds.
#
# Metadata is exported to metadata_monocle.csv with columns
# * cell – barcode as used by Seurat
# * library_size – number of molecules in RNA assay
# * n_features – number of features in RNA assay
# * percent_mito – percentage of mitochondrial genes
# * sample – patient ID
# * group – parient group
# * bcds_score – doublet score
# * umap_1_unaligned ↓
# * umap_2_unaligned – UMAP coordinates before integration
# * umap_1 ↓
# * umap_2 – UMAP coordinates after integration
# * tsne_1 ↓
# * tsne_2 – tSNE coordinates
# * cluster_[k] – cluster IDs for different numbers k of nearest neighbors
# * partition_[k] – partition IDs for different k
# * signature_[gs] – scores for gene signature [gs]
#
# @DEPI rna_qcpassed.rds
# @DEPO rna_integrated_monocle.rds
# @DEPO metadata_monocle.csv

library(monocle3) 
library(Seurat)
library(UCell)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

# the QC-filtered datasets
filtered_datasets <- "data_generated/rna_qcpassed.rds"

# folder where results are saved
out_dir <- "data_generated"



# Integrate samples -------------------------------------------------------

nb <- readRDS(filtered_datasets)

cds <-
  map(
    nb,
    ~new_cell_data_set(.$RNA@counts, cell_metadata = .@meta.data)
  ) %>% 
  combine_cds(cell_names_unique = TRUE, sample_col_name = "unused")

cds <- 
  preprocess_cds(cds, verbose = TRUE) %>% 
  reduce_dimension(preprocess_method = "PCA", verbose = TRUE)

set.seed(42)
cds_aligned <- 
  align_cds(cds, alignment_group = "sample", verbose = TRUE) %>% 
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

cds_50 <- cluster_cells(cds_aligned, k = 50, random_seed = 42, verbose = TRUE)

saveRDS(cds_50, path_join(c(out_dir, "rna_integrated_monocle.rds")))



# Calculate gene signature scores -----------------------------------------

nb_signatures <- 
  read_csv("metadata/nb_signatures.csv", comment = "#") %>% 
  group_by(cell_type) %>% 
  summarise(x = list(gene)) %>% 
  deframe()

uscores <- ScoreSignatures_UCell(counts(cds), nb_signatures)


# Export metadata ---------------------------------------------------------

monocle_metadata <- list(
  colData(cds) %>%
    as_tibble(rownames = "cell") %>%
    select(!Size_Factor),
  reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1_unaligned", "umap_2_unaligned")) %>% 
    as_tibble(rownames = "cell"),
  reducedDim(cds_aligned, "UMAP") %>%
    magrittr::set_colnames(c("umap_1_monocle", "umap_2_monocle")) %>% 
    as_tibble(rownames = "cell"),
  reducedDim(cds_aligned, "tSNE") %>%
    magrittr::set_colnames(c("tsne_1_monocle", "tsne_2_monocle")) %>% 
    as_tibble(rownames = "cell"),
  clusters(cds_aligned) %>%
    enframe(name = "cell", value = "cluster_20"),
  partitions(cds_aligned) %>%
    enframe(name = "cell", value = "partition_20"),
  clusters(cds_50) %>% 
    enframe(name = "cell", value = "cluster_50"),
  partitions(cds_50) %>%
    enframe(name = "cell", value = "partition_50"),
  uscores %>% 
    as_tibble(rownames = "cell") %>% 
    rename_with(~str_replace(., "(.+)_UCell", "signature_\\1"))
) %>%
  reduce(left_join, by = "cell")

monocle_metadata %>%
  write_csv(path_join(c(out_dir, "metadata_monocle.csv")))
