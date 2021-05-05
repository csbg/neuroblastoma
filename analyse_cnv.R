library(monocle3)
library(infercnv)
library(testthat)
library(biomaRt)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

# folder with input files
# root_dir <- "data_generated"
root_dir <- "/media/AGFORTELNY/PROJECTS/Neuroblastoma/analysis/wolfgang/data_generated"

# folder where infercnv results are saved
out_dir <- path_join(c(root_dir, "infercnv_output_v4"))

# posterior probabilities for filtering normal regions
post_probs <- c(0.5, 0.4, 0.3, 0.2, 0.1)



# Load data ---------------------------------------------------------------

nb <- readRDS(path_join(c(root_dir, "rna_integrated_monocle.rds")))
nb_data <- readRDS(path_join(c(root_dir, "metadata.rds")))



# Prepare input data ------------------------------------------------------

# annotation of cells as reference or malignant
annotation_data <-
  nb_data %>%
  arrange(group, sample) %>% 
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
  filter(str_detect(chromosome_name, "^(\\d\\d?|X|Y)$")) %>%
  mutate(chromosome_name = chromosome_name %>% as_factor() %>% fct_inseq()) %>%
  arrange(chromosome_name, start_position) %>%
  mutate(chromosome_name = str_glue("chr{chromosome_name}")) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  column_to_rownames("hgnc_symbol")



# Run infercnv ------------------------------------------------------------

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts(nb),            # |
  annotations_file = annotation_data,        # | inputs calculated above
  ref_group_names = ref_group_names,         # |
  gene_order_file = gene_order,              # |
  chr_exclude = NULL,                        # genes have already been filtered
  # max_cells_per_group = 50                 # uncomment for testing
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  out_dir = out_dir,
  cutoff = 0.1,                     # recommended for 10x data
  cluster_by_groups = TRUE,         # cluster observations by sample
  denoise = TRUE,                   # denoise results
  HMM = TRUE,                       # run HMM to predict CNV level
  BayesMaxPNormal = post_probs[1],  # use the first element of post_probs
  no_plot = TRUE                    # disable (slow!) plotting
)



# Filter CNVs at different posterior probabilities ------------------------

# load data
mcmc_obj <- readRDS(
  path_join(c(out_dir, "18_HMM_pred.Bayes_NetHMMi6.hmm_mode-samples.mcmc_obj"))
)
mcmc_obj@args$out_dir <- path_join(
  c(out_dir, "BayesNetOutput.HMMi6.hmm_mode-samples")
)

infercnv_obj <- readRDS(
  path_join(c(out_dir, "17_HMM_predHMMi6.hmm_mode-samples.infercnv_obj"))
)

# perform filtering for the remaining probability levels
walk(
  post_probs[-1],
  function(prob) {
    res <- filterHighPNormals(mcmc_obj, infercnv_obj@expr.data, prob)
    infercnv:::adjust_genes_regions_report(
      res[[1]], 
      input_filename_prefix = "17_HMM_predHMMi6.hmm_mode-samples", 
      output_filename_prefix =
        str_glue("HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_{prob}"),
      out_dir = out_dir)
  }
)
