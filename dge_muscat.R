library(muscat)
library(scater)
library(tidyverse)



# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_integrated_monocle.rds") %>% 
  logNormCounts()

nb@colData <-
  readRDS("data_generated/metadata.rds") %>%
  mutate(Size_Factor = colData(nb)$Size_Factor) %>% 
  column_to_rownames("cell") %>% 
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)
nb

# currently exclude group I
# there are no cells in the NB cluster (8) in group I, so DGE analysis does not work
nb <- nb[, nb$group != "I"]



# DGE analysis ------------------------------------------------------------

## Data preparation ----

nb <- prepSCE(
  nb,
  kid = "cluster_50",
  gid = "group",
  sid = "sample",
  drop = TRUE
)



## Approach 1: Pseudo-bulk ----

pb <- aggregateData(
  nb,
  assay = "counts",
  fun = "sum",
  by = c("cluster_id", "sample_id")
)
pb


exp_info <- metadata(pb)$experiment_info
model_mat <- model.matrix(~ 0 + fct_drop(exp_info$group_id))
dimnames(model_mat) <- list(exp_info$sample_id, levels(fct_drop(exp_info$group_id)))
model_mat

# dge_2_vs_1 <- pbDS(
#   pb,
#   design = model_mat,
#   contrast = limma::makeContrasts(II - I, levels = model_mat)
# )



dge_4_vs_2 <- pbDS(
  pb,
  design = model_mat,
  contrast = limma::makeContrasts(IV - II, levels = model_mat),
  min_cells = 1
)

## Approach 2: Cell-level mixed model ----

# TBD
# dge_mm <- mmDS(nb)



# Evaluation of results ---------------------------------------------------

dge_4_vs_2$table[[1]] %>% str(1)
nb_4v2 <- dge_4_vs_2$table[[1]][["8"]]

nb_4v2 %>% 
  as_tibble() %>% 
  # arrange(p_adj.loc) %>% 
  filter(p_adj.loc < 0.05) %>%
  filter(abs(logFC) > 1) %>% 
  {.}

