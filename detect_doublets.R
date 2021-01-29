# Detect doublets via scds.
#
# Use RNA counts to calculate doublet scores using
# (1) co-expression of gene pairs (`scds::cxds()`),
# (2) binary classification using artificially generated doublets
#     (`scds::bcds()`), and
# (3) a hybrid approach (`scds::cxds_bcds_hybrid()`).
#
# Generates a CSV file `nb_doublet.csv` with columns
# * `cell` and `sample` (as usual)
# * `cxds_score`, `bcds_score`, and `hybrid_score` (doublet scores)

library(Seurat)
library(tidyverse)
library(scds)
library(SingleCellExperiment)

nb <- readRDS("data_generated/all_datasets_current/nb_integrated.rds")
DefaultAssay(nb) <- "RNA"

nb <-
    as.SingleCellExperiment(nb) %>%
    cxds(verb = TRUE) %>%
    bcds(verb = TRUE) %>%
    cxds_bcds_hybrid(verb = TRUE)

colData(nb) %>%
    as_tibble(rownames = "cell") %>% 
    select(cell, sample, ends_with("_score")) %>%
    write_csv("data_generated/all_datasets_current/nb_doublet.csv")
