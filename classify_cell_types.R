# Cell type assignment via SingleR.
#
# Creates two CSV files:
# * cell_types_singler.csv – the dataframe returned by SingleR::SingleR(),
#                            with added cell and sample information
# * degenes_singler.csv – @metadata$de.genes of SingleR results (wide form)


library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

# the merged dataset
merged_data <- "data_generated/rna_merged.rds"

# folder where results are saved
out_dir <- "data_generated"



# Load data ---------------------------------------------------------------

nb <- readRDS(merged_data)

# for testing: choose only a few cells
# nb <- subset(nb, cells = sample(ncol(nb), 50))

reference_cell_types <- HumanPrimaryCellAtlasData()



# Predict cell types ------------------------------------------------------

predicted_cell_types <- SingleR(
  test = nb$RNA@counts,
  ref = reference_cell_types,
  labels = reference_cell_types$label.fine
)



# Save data ---------------------------------------------------------------

nb@meta.data %>% 
  select(sample) %>% 
  as_tibble(rownames = "cell") %>% 
  bind_cols(as_tibble(predicted_cell_types)) %>% 
  write_csv(path_join(c(out_dir, "cell_types_singler.csv")))

tibble(
  cell_type_1 = names(predicted_cell_types@metadata$de.genes),
  lc = predicted_cell_types@metadata$de.genes
) %>%
  unnest_longer(lc, indices_to = "cell_type_2", values_to = "genes") %>% 
  rowwise() %>% 
  mutate(genes = str_c(genes, collapse = ", ")) %>% 
  relocate(cell_type_2, .after = cell_type_1) %>% 
  write_csv(path_join(c(out_dir, "degenes_singler.csv")))

