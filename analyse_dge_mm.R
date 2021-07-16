# Differential gene expression analysis using mixed models.
#
# @DEPI rna_decontaminated.rds
# @DEPI metadata.rds
# @DEPO dge_mm_results_[group]_[cluster].rds
# @DEPO dge_mm_results.RData

library(monocle3)
library(muscat)
library(scater)
library(tidyverse)
source("common_functions.R")



# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>%
  logNormCounts(assay.type = "soupx_counts")

nb_metadata <- readRDS("data_generated/metadata.rds")

# tumor infiltration rate
tif <-
  nb_metadata %>%
  group_by(sample) %>%
  summarise(tif = sum(cellont_abbr == "NB") / n())

nb@colData <-
  nb_metadata %>%
  mutate(
    sample =
      str_c(group, sample, sep = ".") %>%
      as_factor() %>%
      fct_relevel(str_sort) %>%
      fct_relabel(~str_extract(.x, "\\d+_\\d+")),
    Size_Factor = colData(nb)$Size_Factor
  ) %>%
  left_join(tif, by = "sample") %>%
  column_to_rownames("cell") %>%
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)

# clusters that contain more than 1% of total cells
used_clusters <-
  nb_metadata %>%
  count(cellont_cluster) %>%
  mutate(n = n / sum(n)) %>%
  filter(n > 0.01) %>%
  pull(cellont_cluster)




# Analysis ----------------------------------------------------------------

nb <- prepSCE(
  nb[, colData(nb)$cellont_cluster %in% used_clusters],
  kid = "cellont_abbr",
  gid = "group",
  sid = "sample",
  drop = FALSE
)

run_mm <- function(group, cluster) {
  info("Group {group}, cluster {cluster}")

  # subset data and ensure that empty factor levels are dropped
  nb_small <- nb[
    ,
    colData(nb)$cluster_id == cluster
    & colData(nb)$group_id %in% c("I", group)
  ]
  nb_small$cluster_id <- fct_drop(nb_small$cluster_id)
  nb_small$group_id <- fct_drop(nb_small$group_id)
  nb_small$sample_id <- fct_drop(nb_small$sample_id)
  metadata(nb_small)$experiment_info <- NULL

  # run mixed model-based DGEA
  dge <- mmDS(
    nb_small,
    coef = paste0("group_id", group),
    covs = "tif",
    n_cells = 1,
    n_threads = 1,
    verbose = TRUE
  )

  # save results
  saveRDS(
    dge,
    str_glue("data_generated/dge_mm_results_{group}_{cluster}.rds")
  )
}


params <- 
  list(
    group = c("II", "III", "IV"),
    cluster = colData(nb)$cluster_id %>% fct_drop() %>% levels()
  ) %>% 
  cross_df() %>% 
  filter(cluster != "NB")

safe_run_mm <- possibly(run_mm, NULL, FALSE)

pwalk(params, safe_run_mm)



# Save input --------------------------------------------------------------

save(
  nb,
  nb_metadata,
  used_clusters,
  file = "data_generated/dge_mm_results.RData"
)
