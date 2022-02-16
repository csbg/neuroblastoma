# @DEPI rna_decontaminated.rds
# @DEPI metadata.rds
# @DEPO rna_myeloid.rds
# @DEPO metadata_myeloid.rds
# @DEPO dge_results_myeloid.rds

library(monocle3)
library(muscat)
library(scuttle)
library(scater)
library(nebula)
library(fgsea)
library(msigdbr)
library(tidyverse)
library(scico)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

nb_metadata <- readRDS("data_generated/metadata.rds")

myeloid_barcodes <- 
  nb_metadata %>% 
  filter(cellont_abbr == "M") %>% 
  pull(cell)

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  magrittr::extract(, myeloid_barcodes) %>% 
  logNormCounts(assay.type = "soupx_counts")



# Alignment and clustering ------------------------------------------------

set.seed(42)
cds_my <- 
  nb %>% 
  preprocess_cds(verbose = TRUE) %>% 
  reduce_dimension(preprocess_method = "PCA", verbose = TRUE) %>% 
  align_cds(alignment_group = "sample", verbose = TRUE) %>% 
  reduce_dimension(
    reduction_method = "UMAP",
    preprocess_method = "Aligned",
    verbose = TRUE
  ) %>%
  cluster_cells(k = 20, random_seed = 42, verbose = TRUE)

my_metadata <-
  list(
    tibble(cell = rownames(colData(cds_my))),
    nb_metadata %>%
      select(cell, sample, group) %>% 
      mutate(group = rename_groups(group), sample = rename_patients(sample)),
    reducedDim(cds_my, "UMAP") %>%
      magrittr::set_colnames(c("UMAP1", "UMAP2")) %>% 
      as_tibble(rownames = "cell"),
    clusters(cds_my) %>% 
      enframe(name = "cell", value = "subcluster")
  ) %>% 
  reduce(left_join, by = "cell") %>% 
  mutate(
    collcluster =
      fct_collapse(
        subcluster,
        "classical mono" = c("2", "4", "5"),
        "mDCs" = "7",
        "nonclassical mono" = "8",
        other_level = "other"
      ) %>%
      fct_relevel("classical mono", "nonclassical mono", "mDCs")
  )



# DGE ---------------------------------------------------------------------

## Prepare data ----

nb <- cds_my

tif <-
  readRDS("data_generated/metadata.rds") %>%
  group_by(sample) %>%
  summarise(tif = sum(cellont_abbr == "NB") / n()) %>% 
  mutate(sample = rename_patients(sample))

nb@colData <-
  my_metadata %>%
  mutate(
    Size_Factor = colData(nb)$Size_Factor
  ) %>%
  left_join(tif, by = "sample") %>%
  column_to_rownames("cell") %>%
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)



## Find genes ----

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

set.seed(1)

dge_results_vs_C <- map_dfr(
  my_metadata$collcluster %>%
    levels() %>%
    set_names(),
  analyse_dge,
  ref_group = "C",
  other_groups = c("M", "A", "S"),
  .id = "cell_type"
)

dge_results_MNA_vs_other <- map_dfr(
  my_metadata$collcluster %>%
    levels() %>%
    set_names(),
  analyse_dge,
  ref_group = "other",
  other_groups = "M",
  collapse_groups = list(other = c("A", "S")),
  .id = "cell_type"
)


## Filter genes ----

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
      my_metadata %>% distinct(sample, group),
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


# identical to the function in analyse_dge.mm
filter_dge_results <- function(data,
                               max_p = Inf,
                               max_p_adj = 0.05,
                               min_abs_log_fc = 1,
                               min_freq = 0.05,
                               remove_ribosomal = TRUE) {
  res <-
    data %>%
    filter(
      abs(logFC) >= min_abs_log_fc,
      p <= max_p,
      p_adj <= max_p_adj,
      frq >= min_freq | frq_ref >= min_freq
    ) %>% 
    mutate(direction = if_else(logFC > 0, "up", "down"))
  
  if (remove_ribosomal) {
    ribo_proteins <-
      msigdbr(species = "Homo sapiens", category = "C2") %>% 
      filter(gs_name == "KEGG_RIBOSOME") %>% 
      pull(human_gene_symbol)
    res <- 
      res %>%
      filter(!gene %in% ribo_proteins)
  }
  
  res
}



dge_results_wide_filtered <- 
  dge_results_wide %>%
  filter_dge_results(max_p_adj = Inf, min_abs_log_fc = 0)

dge_results_wide_filtered %>% 
  filter(abs(logFC) < 10) %>% 
  ggplot(aes(logFC / log(2), -log10(p_adj))) +
  geom_point(alpha = .1)



# GSEA --------------------------------------------------------------------

# section copied from analyse_dge.R

#' Download enrichr databases in a format that can be used by fgsea.
#'
#' @param dbs Databases to download.
#'
#' @return A named list with names deriving from values in `dbs`. Each element
#'   is a named list. Names correspond to enrichr terms, values are character
#'   vectors that comprise all genes associated with the respective term.
get_enrichr_genesets <- function(dbs) {
  dbs %>% 
    map(
      function(db) {
        info("Downloading {db}")
        url <- paste0(
          "https://maayanlab.cloud/Enrichr/geneSetLibrary",
          "?mode=text&libraryName=",
          db
        )
        read_lines(url)
      }
    ) %>% 
    set_names(dbs) %>% 
    map(
      function(db) {
        m <- str_match(db, "(.+?)\\t\\t(.+)")
        terms <- m[, 2]
        genes <- m[, 3] %>% str_split("\\t")
        genes %>% 
          map(stringi::stri_remove_empty) %>% 
          set_names(terms)
      }
    )
}


#' Perform gene set enrichment analysis.
#'
#' @param data DGE data as returned by `filter_dge_results()`.
#' @param gene_sets Gene set as returned by `get_enrichr_genesets()`.
#'
#' @return A dataframe, comprising columns "db", "comparison", "cell_type", as
#'   well as all columns in the result of `fgseaMultilevel()`.
perform_gsea <- function(data, gene_sets) {
  data %>% 
    distinct(comparison, cell_type) %>% 
    mutate(db = list(names(gene_sets))) %>% 
    unnest_longer(db) %>% 
    pmap_dfr(
      function(comparison, cell_type, db) {
        ranked_genes <-
          data %>% 
          filter(comparison == {{comparison}}, cell_type == {{cell_type}}) %>%
          select(gene, logFC) %>%
          deframe()
        ranked_genes <- ranked_genes[!is.na(ranked_genes)]
        
        info("GSEA of comparison {comparison}, cell type {cell_type}, ",
             "db {db} ({length(ranked_genes)} genes)")
        
        fgseaMultilevel(
          gene_sets[[db]],
          ranked_genes,
          eps = 0,
          nPermSimple = 10000
        ) %>%
          as_tibble() %>%
          mutate(
            db = {{db}},
            comparison = {{comparison}},
            cell_type = {{cell_type}},
            .before = 1
          )
      }
    )
}

# already defined in secion 'enrichment analysis',
# but we repeat it here for convenience
enrichr_dbs <- c(
  "GO_Biological_Process_2018",
  "GO_Cellular_Component_2018",
  "GO_Molecular_Function_2018",
  "KEGG_2019_Human",
  "WikiPathways_2019_Human",
  "MSigDB_Hallmark_2020",
  "TRRUST_Transcription_Factors_2019"
)

enrichr_genesets <- get_enrichr_genesets(enrichr_dbs)

# remove mouse genes from TTRUST database
enrichr_genesets$TRRUST_Transcription_Factors_2019 <-
  enrichr_genesets$TRRUST_Transcription_Factors_2019 %>% 
  magrittr::extract(imap_lgl(., ~str_detect(.y, "human"))) %>% 
  set_names(str_extract, "\\w+")

gsea_results <-
  dge_results_wide_filtered %>%
  perform_gsea(enrichr_genesets)



# Save data ---------------------------------------------------------------

cds_my %>% saveRDS("data_generated/rna_myeloid.rds")
my_metadata %>% saveRDS("data_generated/metadata_myeloid.rds")

list(
  cds = nb,
  metadata = my_metadata,
  results_vs_C = dge_results_vs_C,
  results_MNA_vs_other = dge_results_MNA_vs_other,
  results_wide = dge_results_wide,
  results_wide_filtered = dge_results_wide_filtered,
  gsea = gsea_results,
  gene_sets = enrichr_genesets
) %>%
  saveRDS("data_generated/dge_results_myeloid.rds")