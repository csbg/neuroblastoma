# @DEPI rna_decontaminated.rds
# @DEPI metadata.rds
# @DEPO dge_mm_results.rds

library(monocle3)
library(muscat)
library(scater)
library(nebula)
library(msigdbr)
library(enrichR)
library(latex2exp)
library(fgsea)
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

# filter cds and update factor levels
nb <- nb[, colData(nb)$cellont_cluster %in% used_clusters]
colData(nb)$cellont_cluster <- fct_drop(colData(nb)$cellont_cluster)
colData(nb)$cellont_abbr <- fct_drop(colData(nb)$cellont_abbr)



# Analyse data ------------------------------------------------------------

#' Perform DGE analysis.
#'
#' @param cell_type Selected cell type.
#' @param ref_group Reference patient group.
#' @param other_groups Patient groups that are compared to the reference.
#' @param collapse_groups Optional list of named character vectors that allows
#'   to collapse groups prior to analysis (passed to `forcats::fct_collapse()`).
#'
#' @return A data frame with DGE results.
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
      cellont_abbr == cell_type,
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
  colData(nb)$cellont_abbr %>%
    levels() %>%
    setdiff("NB") %>%  # ignore tumor cells since there are none in the control
    set_names(),
  analyse_dge,
  ref_group = "I",
  other_groups = c("II", "III", "IV"),
  .id = "cell_type"
)

dge_results_vs_S <- map_dfr(
  colData(nb)$cellont_abbr %>%
    levels() %>%
    set_names(),
  analyse_dge,
  ref_group = "IV",
  other_groups = c("II", "III"),
  .id = "cell_type"
)

dge_results_vs_A <- map_dfr(
  colData(nb)$cellont_abbr %>%
    levels() %>%
    set_names(),
  analyse_dge,
  ref_group = "III",
  other_groups = "II",
  .id = "cell_type"
)

dge_results_MNA_vs_other <- map_dfr(
  colData(nb)$cellont_abbr %>%
    levels() %>%
    set_names(),
  analyse_dge,
  ref_group = "other",
  other_groups = "MNA",
  collapse_groups = list(MNA = "II", other = c("III", "IV")),
  .id = "cell_type"
)



#' Calculate minimum sample expression frequencies.
#' 
#' Expression frequency denotes the fraction of cells that expresses a given
#' gene. These values are summarized by group so that the minimum frequency over
#' all samples is reported.
#'
#' @param data MM DGE results.
#' @param collapse_groups see `analyse_dge()`
#'
#' @return A dataframe with additional columns
#' "frq": expression frequency in the group indicated in the respective column
#' "frq_ref": expression frequency in the reference group
calc_expression_frequency <- function(data, collapse_groups = NULL) {
  # required for calcExprFreqs()
  cds <- prepSCE(
    nb,
    kid = "cellont_abbr",
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
    select(!I:IV)
  
  # summarise at group level
  exp_frq_groups <- 
    exp_frq %>% 
    pivot_longer(starts_with("20"), names_to = "sample", values_to = "frq") %>% 
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
      regex = "(.+)_vs_(.+)",
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


#' Pivots DGE results to longer format.
#'
#' @param df Raw DGE results.
#' @param suffix String that is appended to the comparison column.
#'
#' @return A data frame with columns cell_type, gene, comparison, logFC, and p.
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
    gather_dge_results(dge_results_vs_C, "_vs_I"),
    gather_dge_results(dge_results_vs_S, "_vs_IV"),
    gather_dge_results(dge_results_vs_A, "_vs_III"),
    gather_dge_results(dge_results_MNA_vs_other, "_vs_other")
  ) %>% 
  arrange(cell_type, gene) %>% 
  group_by(comparison, cell_type) %>%
  mutate(p_adj = p.adjust(p, method = "fdr")) %>%
  ungroup() %>%
  calc_expression_frequency(
    collapse_groups = list(MNA = "II", other = c("III", "IV"))
  )



# Filter results ----------------------------------------------------------

#' Filter DGE results.
#'
#' @param data A data frame with MM DGE results.
#' @param max_p Maximum p value.
#' @param max_p_adj Maximum adjusted p value.
#' @param min_abs_log_fc Minimum absolute log fold change.
#' @param min_freq Minimum gene expression frequency.
#' @param remove_ribosomal If true, remove ribosomal proteins.
#'
#' @return A data frame with an additional column "direction".
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



# Enrichment analysis -----------------------------------------------------

#' Perform enrichment analysis via enrichr for upregulated genes in a given
#' contrast and cluster.
#'
#' @param data A DGE dataset.
#' @param comparison The selected comparison.
#' @param cell_type The selected cell type
#' @param dbs Character vector containing valid enrichr databases
#'   as returned by `enrichR::listEnrichrDbs()`.
#' @param direction Use genes that are up- or downregulated, respectively.
#'
#' @return A dataframe, which combines the dataframes returned by
#'   `enrichR::enrichr()` (empty results are removed) by adding a column "db".
enrich_genes <- function(data,
                         comparison,
                         cell_type,
                         dbs,
                         direction = c("up", "down")) {
  info("Cell type {cell_type}, comparison {comparison}, ",
       "{direction}regulated genes")
  
  direction <- match.arg(direction)
  
  top_genes <- 
    data %>% 
    filter(
      cell_type == {{cell_type}},
      comparison == {{comparison}},
      direction == {{direction}}
    ) %>% 
    pull(gene)
  
  if (length(top_genes) > 0) {
    info("genes for enrichr: {str_c(top_genes, collapse = ', ')}")
    enrichr(top_genes, dbs) %>% 
      keep(~nrow(.) > 0) %>% 
      bind_rows(.id = "db") 
  }
}

#' Perform enrichment analysis for all contrasts and cell types in a dataset.
#'
#' @param data A DGE dataset.
#' @param dbs Character vector containing valid enrichr databases
#'   as returned by `enrichR::listEnrichrDbs()`.
#'
#' @return A dataframe, which combines the dataframes returned by
#'   `enrichR::enrichr()` (empty results are removed) by adding three columns
#'   "comparison", "cell_type", and "db".
enrich_all_genes <- function(data, dbs) {
  queries <- 
    data %>%
    distinct(comparison, cell_type) %>% 
    mutate(direction = list(c("up", "down"))) %>% 
    unnest_longer(direction)
  
  enrichr_data <-
    queries %>%
    pmap(enrich_genes, data = data, dbs = dbs)
  
  queries %>%
    mutate(data = enrichr_data) %>%
    rowwise() %>%
    filter(!is.null(data)) %>%
    filter(nrow(data) > 0) %>%
    unnest(data) %>%
    mutate(comparison = as_factor(comparison)) %>% 
    separate(
      Overlap,
      into = c("overlap_size", "geneset_size"),
      convert = TRUE
    )
}

enrichr_dbs <- c(
  "GO_Biological_Process_2018",
  "GO_Cellular_Component_2018",
  "GO_Molecular_Function_2018",
  "KEGG_2019_Human",
  "WikiPathways_2019_Human",
  "MSigDB_Hallmark_2020",
  "TRRUST_Transcription_Factors_2019"
)

enrichr_results <-
  dge_results_wide %>%
  filter_dge_results() %>%
  enrich_all_genes(enrichr_dbs)



# GSEA --------------------------------------------------------------------

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



# Save results ------------------------------------------------------------

# nb <- dge$cds
# nb_metadata <- dge$metadata
# used_clusters <- dge$used_clusters
# dge_results_vs_C <- dge$results_vs_C
# dge_results_vs_S <- dge$results_vs_S
# dge_results_vs_A <- dge$results_vs_A
# dge_results_MNA_vs_other <- dge$results_MNA_vs_other
# dge_results_wide <- dge$results_wide
# dge_results_wide_filtered <- dge$results_wide_filtered
# enrichr_results <- dge$enrichr
# gsea_results <- dge$gsea
# enrichr_genesets <- dge$gene_sets

list(
  cds = nb,
  metadata = nb_metadata,
  used_clusters = used_clusters,
  results_vs_C = dge_results_vs_C,
  results_vs_S = dge_results_vs_S,
  results_vs_A = dge_results_vs_A,
  results_MNA_vs_other = dge_results_MNA_vs_other,
  results_wide = dge_results_wide,
  results_wide_filtered = dge_results_wide_filtered,
  enrichr = enrichr_results,
  gsea = gsea_results,
  gene_sets = enrichr_genesets
) %>%
  saveRDS("data_generated/dge_mm_results.rds")
