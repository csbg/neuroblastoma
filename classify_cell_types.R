# Creates several CSV files named 'cell_types_[ref]_[labels].csv', where [ref]
# indicates the celldex reference dataset and [labels] denotes whether the broad
# or fine labels from the reference dataset were used for classification. Each
# file contains the dataframe returned by SingleR::SingleR(), with added cell
# and sample information
#
# @DEPI rna_qcpassed.rds
# @DEPO cell_types_[ref]_[labels].csv


library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)
source("common_functions.R")



# Parameters --------------------------------------------------------------

# the QC-filtered datasets
filtered_datasets <- "data_generated/rna_qcpassed.rds"

# folder where results are saved
out_dir <- "data_generated"



# Load data ---------------------------------------------------------------

nb <- readRDS(filtered_datasets)
count_matrix <- merge(nb[[1]], nb[-1])$RNA@counts


# subset cells for testing
# count_matrix <- count_matrix[, sample(colnames(count_matrix), 10)]


reference_cell_types <- list(
  # two general purpose datasets
  hpca = HumanPrimaryCellAtlasData(),
  blueprint = BlueprintEncodeData(),
  
  # comprehensive CD4+ subsets; only one B cell subset, no dendritic cells
  dice = DatabaseImmuneCellExpressionData(),
  
  # for bone marrow samples
  dmap = NovershternHematopoieticData(),
  
  # for PBMC
  monaco = MonacoImmuneData()
)



# Predict cell types ------------------------------------------------------

predicted_cell_types <- 
  list(
    ref = names(reference_cell_types),
    labels = c("label.main", "label.fine")
  ) %>% 
  cross_df() %>% 
  pmap(
    function(ref, labels) {
      info("Classifying with reference dataset '{ref}' and labels '{labels}'")
      results <- SingleR(
        test = count_matrix,
        ref = reference_cell_types[[ref]],
        labels = colData(reference_cell_types[[ref]])[, labels]
      )
      list(
        table = results,
        ref = ref,
        labels = labels
      )
    }
  )



# Save data ---------------------------------------------------------------

save_results <- function(results) {
  score_colnames <- str_c(
    "score_",
    colnames(results$table$scores)
  )
  
  df <- as_tibble(results$table, rownames = "cell")
  colnames(df) <- c(
    "cell", score_colnames,
    "first_labels", "tuning_scores_first", "tuning_scores_second",
    "labels", "pruned_labels"
  )

  write_csv(
    df,
    str_glue("{out_dir}/cell_types_{results$ref}_{results$labels}.csv")
  )
}

predicted_cell_types %>% walk(save_results)
