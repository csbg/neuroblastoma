# Cell type assignment via SingleR for an integrated Seurat dataset
#
# Creates three CSV files:
# (1) nb_singler.csv – cell barcodes and unfiltered cell types
# (2) nb_singler_details.csv – the dataframe returned by SingleR::SingleR()
# (3) nb_singler_degenes.csv – @metadata$de.genes of SingleR results (wide form)


library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)


# small test datasets
# infile <- "data_generated/3_datasets/nb_integrated.rds"
# outfile <- "data_generated/3_datasets/nb_singler"
# 
# large complete datasets
infile <- "data_generated/all_datasets_current/nb_integrated.rds"
outfile <- "data_generated/all_datasets_current/nb_singler"

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
  bind_cols(cell_type = predicted_cell_types$labels) %>% 
  write_csv(str_glue("{outfile}.csv"))

nb@meta.data %>% 
  select(sample) %>% 
  as_tibble(rownames = "cell") %>% 
  bind_cols(as_tibble(predicted_cell_types)) %>% 
  write_csv(str_glue("{outfile}_details.csv"))

tibble(
  cell_type_1 = names(predicted_cell_types@metadata$de.genes),
  lc = predicted_cell_types@metadata$de.genes
) %>%
  unnest_longer(lc, indices_to = "cell_type_2", values_to = "genes") %>% 
  rowwise() %>% 
  mutate(genes = str_c(genes, collapse = ", ")) %>% 
  relocate(cell_type_2, .after = cell_type_1) %>% 
  write_csv(str_glue("{outfile}_degenes.csv"))

