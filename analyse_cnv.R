library(monocle3)
library(infercnv)
library(testthat)
library(biomaRt)
library(tidyverse)


nb <- readRDS("data_generated/rna_integrated_monocle.rds")

annotation_data <-
  readRDS("data_generated/metadata.rds") %>%
  transmute(
    cell = cell,
    cell_type = case_when(
      cluster_20 %in% c("5", "32") ~ str_glue("malignant_{group}"),
      TRUE ~ as.character(cell_type_hpca_broad)
    )
  ) %>% 
  column_to_rownames("cell")
expect_length(
  intersect(colnames(counts(nb)), annotation_data$cell),
  ncol(counts(nb))
)

ref_group_names <-
  annotation_data %>% 
  filter(!cell_type %>% str_starts("malignant")) %>% 
  pull(cell_type) %>% 
  unique()

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_order <-
  getBM(
    attributes = c(
      "hgnc_symbol", "chromosome_name", "start_position", "end_position"
    ),
    filters = "hgnc_symbol",
    values = rownames(counts(nb)),
    mart = ensembl
  ) %>% 
  as_tibble() %>% 
  filter(str_length(chromosome_name) <= 2) %>% 
  mutate(chromosome_name = str_glue("chr{chromosome_name}")) %>% 
  distinct(hgnc_symbol, .keep_all = TRUE) %>% 
  column_to_rownames("hgnc_symbol")



infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts(nb),
  annotations_file = annotation_data,
  gene_order_file = gene_order,
  ref_group_names = ref_group_names,
  max_cells_per_group = 50  # for testing
)


infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "infercnv_output",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE
)
