# @DEPI rna_decontaminated.rds
# @DEPI tumor_data_dong.rds
# @DEPI data_raw/adrmed
# @DEPO adrmed_class_*.csv

library(Seurat)
library(monocle3)
library(SingleR)
library(BiocParallel)
library(tidyverse)
source("common_functions.R")



# Load data ---------------------------------------------------------------

nb <- readRDS("data_generated/rna_decontaminated.rds")
nb@colData <-
  readRDS("data_generated/metadata.rds") %>%
  column_to_rownames("cell") %>%
  as("DataFrame")

nb <- nb[, colData(nb)$cellont_abbr == "NB"]

tumor_dong <- readRDS("data_generated/tumor_data_dong.rds")

ref <- readRDS("data_raw/adrmed/adrenal_medulla_Seurat.RDS")



# Predict cell types ------------------------------------------------------

predict_cell_types <- function(count_matrix, outfile) {
  info("Creating {outfile}")
  results <- SingleR(
    test = count_matrix,
    ref = ref$RNA@data,
    labels = Idents(ref),
    de.method = "wilcox",
    BPPARAM = MulticoreParam(workers = 32, progressbar = TRUE, RNGseed = 42)
  )
  
  score_colnames <- str_c(
    "score_",
    colnames(results$scores)
  )
  
  df <- as_tibble(results, rownames = "cell")
  colnames(df) <- c(
    "cell", score_colnames,
    "first_labels", "tuning_scores_first", "tuning_scores_second",
    "labels", "pruned_labels"
  )
  
  write_csv(df, outfile)
}


predict_cell_types(
  counts(nb),
  "data_generated/adrmed/adrmed_class_nb.csv"
)

# predict_cell_types(
#   tumor_jansky$RNA@counts,
#   "data_wip/adrmed_class_jansky.csv"
# )

tumor_dong %>% 
  iwalk(
    ~predict_cell_types(
      .x$RNA@counts,
      str_glue("data_generated/adrmed/adrmed_class_dong_{.y}.csv")
    )
  )
