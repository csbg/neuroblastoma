# Cell type assignment via SingleR for an integrated Seurat dataset
#
# Creates three results files:
# (1) a CSV file with cell barcodes, tSNE/UMAP coordinates, and cell types
# (2) a RDS file containing the dataframe returned by SingleR::SingleR()
# (3) a CSV file containing all data in (2)


library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)


# small test datasets
infile <- "data_generated/3_datasets_sc/nb_integrated.rds"
outfile <- "data_generated/3_datasets_sc/nb_singler"

# large complete datasets
# infile <- "data_generated/all_datasets_sc/nb_integrated.rds"
# outfile <- "data_generated/all_datasets_sc/nb_singler"

# whether to use broad or fine labels
use_fine_labels <- TRUE



# Load data ---------------------------------------------------------------

nb <- readRDS(infile)

nb <- 
  nb %>% 
  FindNeighbors() %>% 
  FindClusters(resolution = 0.5)  # or 0.8 for more clusters

# for testing: choose only a few cells
nb <- subset(nb, cells = sample(ncol(nb), 50))

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
