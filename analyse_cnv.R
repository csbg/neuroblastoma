library(monocle3)
library(infercnv)
library(testthat)
library(biomaRt)
library(tidyverse)


# Load data ---------------------------------------------------------------

nb <- readRDS("data_generated/rna_integrated_monocle.rds")
nb_data <- readRDS("data_generated/metadata.rds")



# Prepare input data ------------------------------------------------------

# annotation of cells as reference or malignant
annotation_data <-
  nb_data %>%
  transmute(
    cell = cell,
    cell_type = case_when(
      cluster_50 == "8" ~ str_glue("malignant_{sample}_{group}"),
      TRUE              ~ as.character(cellont_name)
    )
  ) %>% 
  column_to_rownames("cell")

expect_length(
  intersect(colnames(counts(nb)), rownames(annotation_data)),
  ncol(counts(nb))
)

# reference names
ref_group_names <-
  annotation_data %>% 
  filter(!cell_type %>% str_starts("malignant")) %>% 
  pull(cell_type) %>% 
  unique()

# order of genes on chromosome
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = 103
)

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



# Analysis ----------------------------------------------------------------

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts(nb),            # |
  annotations_file = annotation_data,        # | inputs calculated above
  ref_group_names = ref_group_names,         # |
  gene_order_file = gene_order,              # |
  chr_exclude = c("chrX", "chrY", "chrMT"),  # excluded chromosomes
  # max_cells_per_group = 50                 # for testing
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "data_generated/infercnv_output_v2",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  plot_steps = TRUE,
  diagnostics = TRUE
)
