# Cell type assignment via SingleR for an integrated Seurat dataset
#
# Creates two CSV files:
# (1) nb_singler.csv – the dataframe returned by SingleR::SingleR(), with added
#                      cell and sample information
# (3) nb_singler_degenes.csv – @metadata$de.genes of SingleR results (wide form)


library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)
library(fs)


# small test datasets
# infile <- "data_generated/3_datasets/nb_integrated.rds"
# outdir <- "data_generated/3_datasets"
# 
# large complete datasets
infile <- "data_generated/all_datasets_current/nb_integrated.rds"
outdir <- "data_generated/all_datasets_current"

# whether to use broad or fine labels
use_fine_labels <- TRUE



# Load data ---------------------------------------------------------------

nb <- readRDS(infile)

# for testing: choose only a few cells
# nb <- subset(nb, cells = sample(ncol(nb), 50))

reference_cell_types <- HumanPrimaryCellAtlasData()



# Predict cell types ------------------------------------------------------

predicted_cell_types <- SingleR(
  GetAssayData(nb$RNA),
  reference_cell_types,
  labels = if (use_fine_labels)
    reference_cell_types$label.fine
  else
    reference_cell_types$label.main
)



# Save data ---------------------------------------------------------------

nb@meta.data %>% 
  select(sample) %>% 
  as_tibble(rownames = "cell") %>% 
  bind_cols(as_tibble(predicted_cell_types)) %>% 
  write_csv(path_join(c(outdir, "nb_singler.csv")))

tibble(
  cell_type_1 = names(predicted_cell_types@metadata$de.genes),
  lc = predicted_cell_types@metadata$de.genes
) %>%
  unnest_longer(lc, indices_to = "cell_type_2", values_to = "genes") %>% 
  rowwise() %>% 
  mutate(genes = str_c(genes, collapse = ", ")) %>% 
  relocate(cell_type_2, .after = cell_type_1) %>% 
  write_csv(path_join(c(outdir, "nb_singler_degenes.csv")))

