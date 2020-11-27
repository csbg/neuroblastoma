# extract tSNE data from a Seurat object

library(tidyverse)
library(Seurat)

nb <- readRDS("data_generated/all_datasets/nb_integrated_singler.rds")

nb_data <- inner_join(
  nb@meta.data %>% as_tibble(rownames = "cell"),
  nb$tsne@cell.embeddings %>% as_tibble(rownames = "cell"),
  by = "cell"
)
nb_data %>% head()

nb_data %>% write_csv("data_generated/all_datasets/nb_tsne_integrated.csv")