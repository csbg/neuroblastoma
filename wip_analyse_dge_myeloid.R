library(monocle3)
library(muscat)
library(scater)
library(nebula)
library(msigdbr)
library(enrichR)
library(latex2exp)
library(fgsea)
library(tidyverse)
source("styling.R")
source("common_functions.R")




# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_wip/cds_myeloid.rds")

nb_metadata <- readRDS("data_wip/metadata_myeloid.rds")

tif <-
  readRDS("data_generated/metadata.rds") %>%
  group_by(sample) %>%
  summarise(tif = sum(cellont_abbr == "NB") / n()) %>% 
  mutate(sample = rename_patients(sample))

nb@colData <-
  nb_metadata %>%
  mutate(
    Size_Factor = colData(nb)$Size_Factor
  ) %>%
  left_join(tif, by = "sample") %>%
  column_to_rownames("cell") %>%
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)



# Analyze data ------------------------------------------------------------

analyse_dge <- function(cell_type,
                        ref_group,
                        other_groups,
                        collapse_groups = NULL) {
  info("Analysing cell type {cell_type}, ",
       "{str_c(other_groups, collapse = '/')} vs {ref_group}")
  
  col_metadata <-
    colData(nb) %>% 
    as_tibble(rownames = "cell")
  
  # optionally, collapse group levels
  if (!is.null(collapse))
    col_metadata <-
    col_metadata %>%
    mutate(group = fct_collapse(group, !!!collapse_groups))
  
  # subset column metadata, set correct factor levels (reference must be first)
  col_metadata <-
    col_metadata %>%
    filter(
      collcluster == cell_type,
      group %in% c(ref_group, other_groups)
    ) %>% 
    mutate(
      group =
        group %>%
        fct_relevel(ref_group, other_groups) %>% 
        fct_drop()
    )
  
  # subset data
  nb_sub <- nb[, col_metadata$cell]
  colData(nb_sub) <-
    col_metadata %>% 
    column_to_rownames("cell") %>% 
    as("DataFrame")
  
  # reorder count matrix as required by nebula
  data_grouped <- group_cell(
    counts(nb_sub),
    id = colData(nb_sub)$sample,
    pred = model.matrix(~group + tif, data = colData(nb_sub)),
    offset = colData(nb_sub)$Size_Factor
  )
  
  # run analysis
  res <- nebula(
    data_grouped$count,
    id = data_grouped$id,
    pred = data_grouped$pred,
    offset = data_grouped$offset,
    verbose = TRUE
  )
  
  # format results
  res$summary %>%
    as_tibble() %>% 
    mutate(
      algorithm = res$algorithm,
      convergence = res$convergence,
      overdispersion_subject = res$overdispersion$Subject,
      overdispersion_cell = res$overdispersion$cell
    )
}

# analyse_dge("mDCs", "C", c("M", "A", "S"))

set.seed(1)

dge_results_vs_C <- map_dfr(
  nb_metadata$collcluster %>%
    levels() %>%
    set_names(),
  analyse_dge,
  ref_group = "C",
  other_groups = c("M", "A", "S"),
  .id = "cell_type"
)

dge_results_MNA_vs_other <- map_dfr(
  nb_metadata$collcluster %>%
    levels() %>%
    set_names(),
  analyse_dge,
  ref_group = "other",
  other_groups = "M",
  collapse_groups = list(other = c("A", "S")),
  .id = "cell_type"
)




calc_expression_frequency <- function(data, collapse_groups = NULL) {
  # required for calcExprFreqs()
  cds <- prepSCE(
    nb,
    kid = "collcluster",
    gid = "group",
    sid = "sample",
    drop = TRUE
  )
  
  # generate table with frequencies on sample level
  exp_frq <-
    cds %>% 
    calcExprFreqs() %>% 
    assays() %>%
    as.list() %>% 
    map_dfr(as_tibble, rownames = "gene", .id = "cell_type") %>% 
    select(!C:S)
  
  # summarise at group level
  exp_frq_groups <- 
    exp_frq %>% 
    pivot_longer(C1:S5, names_to = "sample", values_to = "frq") %>% 
    left_join(
      nb_metadata %>% distinct(sample, group),
      by = "sample"
    ) %>%
    group_by(cell_type, gene, group) %>% 
    summarise(frq = min(frq))
  
  # optionally, add frequencies in collapsed groups
  if (!is.null(collapse_groups))
    exp_frq_groups <- bind_rows(
      exp_frq_groups,
      exp_frq_groups %>% 
        ungroup() %>% 
        mutate(
          group = fct_collapse(group, !!!collapse_groups, other_level = "unused")
        ) %>% 
        filter(group != "unused") %>%
        group_by(cell_type, gene, group) %>%
        summarise(frq = min(frq))
    )
  
  # add columns to input
  data %>% 
    extract(
      comparison,
      into = c("group", "group_ref"),
      regex = "(.)(.+)",
      remove = FALSE
    ) %>% 
    left_join(
      exp_frq_groups,
      by = c("cell_type", "gene", "group")
    ) %>% 
    left_join(
      exp_frq_groups,
      by = c("cell_type", "gene", group_ref = "group")
    ) %>% 
    select(!c(group, group_ref)) %>% 
    rename(frq = frq.x, frq_ref = frq.y)
}


gather_dge_results <- function(df, suffix) {
  df %>%
    select(cell_type, gene, starts_with("logFC_g"), starts_with("p_g")) %>%
    pivot_longer(
      !c(cell_type, gene),
      names_to = c(".value", "comparison"),
      names_pattern = "(.+)_group(.+)"
    ) %>%
    mutate(comparison = str_c(comparison, suffix))
}

dge_results_wide <-
  bind_rows(
    gather_dge_results(dge_results_vs_C, "c"),
    gather_dge_results(dge_results_MNA_vs_other, "as")
  ) %>% 
  arrange(cell_type, gene) %>% 
  group_by(comparison, cell_type) %>%
  mutate(p_adj = p.adjust(p, method = "fdr")) %>%
  ungroup() %>%
  calc_expression_frequency(
    collapse_groups = list(as = c("A", "S"), c = "C")
  )




# filter_dge_results() from plot_dge_mm.R


dge_results_wide_filtered <- 
  dge_results_wide %>%
  filter_dge_results(max_p_adj = Inf, min_abs_log_fc = 0)

dge_results_wide_filtered %>% 
  filter(abs(logFC) < 10) %>% 
  ggplot(aes(logFC / log(2), -log10(p_adj))) +
  geom_point(alpha = .1)


# run section GSEA from plot_dge_mm.R



# Save data ---------------------------------------------------------------

list(
  cds = nb,
  metadata = nb_metadata,
  results_vs_C = dge_results_vs_C,
  results_MNA_vs_other = dge_results_MNA_vs_other,
  results_wide = dge_results_wide,
  results_wide_filtered = dge_results_wide_filtered,
  gsea = gsea_results,
  gene_sets = enrichr_genesets
) %>%
  saveRDS("data_wip/dge_mm_myeloid_results.rds")

