# Detect doublets via scds.
#
# Use RNA counts to calculate doublet scores using
# (1) co-expression of gene pairs (`scds::cxds()`),
# (2) binary classification using artificially generated doublets
#     (`scds::bcds()`), and
# (3) a hybrid approach (`scds::cxds_bcds_hybrid()`).
#
# Generates a CSV file `doublet_scores.csv` with columns
# * `cell`
# * `cxds_score`, `bcds_score`, and `hybrid_score` (doublet scores)

library(Seurat)
library(scds)
library(SingleCellExperiment)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

# the merged dataset
merged_data <- "data_generated/rna_merged.rds"

# folder where results are saved
out_dir <- "data_generated"



# Analysis ----------------------------------------------------------------

nb <- readRDS(merged_data)

nb <-
  as.SingleCellExperiment(nb) %>%
  cxds(verb = TRUE) %>%
  bcds(verb = TRUE) %>%
  cxds_bcds_hybrid(verb = TRUE)

colData(nb) %>%
  as_tibble(rownames = "cell") %>%
  select(cell, ends_with("_score")) %>%
  write_csv(path_join(c(out_dir, "doublet_scores.csv")))
