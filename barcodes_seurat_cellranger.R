# Creates a conversion table between barcode names in cellranger and Seurat.
#
# cellranger-aggregate ensures unique barcodes across samples by numbering them
# like BARCODE-1 to BARCODE-10, using the order of samples in aggregation.csv.
# By contrast, Seurat (if data is loaded with Read10X_h5()) uses the numbering
# system BARCODE-1_1 to BARCODE-1_10; sample indices correspond to the order of
# samples in the object.list argument of FindIntegrationAnchors().
#
# This script uses information scattered across several files to create a CSV
# file (nb_barcodes.csv) mapping between cellranger barcodes (column
# `cellranger_barcode`) and Seurat barcodes (column `seurat_barcode`)

library(tidyverse)

nb_data <- read_csv("data_generated/all_datasets_current/nb_clusters_0.5.csv")
nb_data

cr_cellnames <-
  read_csv("data_raw/AGGR_ALL_transcriptome/cellranger_cellnames.csv")
cr_cellnames



# Establish orders --------------------------------------------------------

cellranger_order <-
  read_csv("data_raw/AGGR_ALL_transcriptome/aggregation.csv") %>% 
  select(bsf_id = library_id) %>% 
  mutate(cr_order = row_number())
cellranger_order

sample_ids <-
  read_csv("data_raw/metadata/sample_groups.csv", comment = "#") %>% 
  {.}
sample_ids[10, 2] <- "2018_1404b"
sample_ids

seurat_order <- 
  nb_data %>% 
  mutate(se_order = str_match(cell, "_(\\d+)")[,2]) %>% 
  distinct(sample, se_order) %>%
  {.}
seurat_order[6, 1] <- "2018_1404b"
seurat_order

order_map <-
  seurat_order %>% 
  left_join(sample_ids, by = "sample") %>% 
  left_join(cellranger_order, by = "bsf_id") %>% 
  select(cr_order, se_order, group) %>%
  {.}
order_map



# Remap indices -----------------------------------------------------------

nb_data_barcodes <- 
  nb_data %>% 
  mutate(se_order = str_match(cell, "_(\\d+)")[,2]) %>% 
  left_join(order_map, by = "se_order") %>% 
  mutate(
    barcode = str_extract(cell, "\\w+"),
    cell_cr = str_c(barcode, cr_order, sep = "-")
  )
nb_data_barcodes

length(
  intersect(nb_data_barcodes$cell_cr, cr_cellnames$cell_cr)
) == nrow(nb_data) 

nb_data_barcodes %>% 
  select(cellranger_barcode = cell_cr, seurat_barcode = cell,
         sample, group, cluster = integrated_snn_res.0.5) %>% 
  write_csv("data_generated/all_datasets_current/nb_barcodes.csv")
