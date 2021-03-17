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
#
# @DEPI rna_merged.rds
# @DEPO doublet_scores.csv

library(Seurat)
library(scds)
library(SingleCellExperiment)
library(tidyverse)
source("common_functions.R")



# Parameters --------------------------------------------------------------

# the merged dataset
merged_data <- "data_generated/rna_merged.rds"

# file where results are saved
out_file <- "data_generated/doublet_scores.csv"



# Analysis ----------------------------------------------------------------

detect_doublets <- function(data, sample) {
  info("Processing dataset {sample}")
  
  sce <- 
    data %>% 
    as.SingleCellExperiment() %>%
    cxds_bcds_hybrid(verb = TRUE)
  
  colData(sce) %>%
    as_tibble(rownames = "cell") %>%
    select(cell, ends_with("score"))
}

nb <- readRDS(merged_data)

set.seed(42)
doublet_scores <-
  nb %>% 
  SplitObject("sample") %>% 
  imap_dfr(detect_doublets)
  


# Save data ---------------------------------------------------------------

# ensure that order of rows is the same as in input dataset
tibble(cell = colnames(nb)) %>% 
  left_join(doublet_scores, by = "cell") %>% 
  write_csv(out_file)
