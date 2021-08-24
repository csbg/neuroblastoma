# Differential gene expression analysis using mixed models.
#
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

nb <- nb[, colData(nb)$cellont_cluster %in% used_clusters]



# Analyse data ------------------------------------------------------------

analyse_dge <- function(cell_type) {
  info("Analysing cell type {cell_type}")
  
  nb_sub <- nb[, colData(nb)$cellont_abbr == cell_type]
  model_mat <- model.matrix(~group + tif, data = colData(nb_sub))
  
  data_grouped <- group_cell(
    counts(nb_sub),
    id = colData(nb_sub)$sample,
    pred = model_mat,
    offset = colData(nb_sub)$Size_Factor
  )
  
  res <- nebula(
    data_grouped$count,
    id = data_grouped$id,
    pred = data_grouped$pred,
    offset = data_grouped$offset,
    verbose = TRUE
  )
  
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
dge_results <- map_dfr(
  colData(nb)$cellont_abbr %>%
    unique() %>%
    setdiff("NB") %>%
    set_names(),
  analyse_dge,
  .id = "cell_type"
)


#' Calculate minimum sample expression frequencies.
#' 
#' Expression frequency denotes the fraction of cells that expresses a given
#' gene. These values are summarized by group so that the minimum frequency over
#' all samples is reported.
#'
#' @param data MM DGE results.
#'
#' @return A dataframe with additional columns
#' "frq": expression frequency in the group indicated in the respective column
#' "frq_ref": expression frequency in the reference group
calc_expression_frequency <- function(data) {
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
    summarise(frq = min(frq)) %>% 
    mutate(frq_ref = if_else(group == "I", frq, NA_real_)) %>% 
    fill(frq_ref) %>% 
    filter(group != "I")
  
  # add columns to input
  data %>%
    left_join(exp_frq_groups, by = c("cell_type", "gene", "group"))
}

dge_results_wide <- 
  dge_results %>% 
  select(cell_type, gene, starts_with("logFC_g"), starts_with("p_g")) %>% 
  pivot_longer(
    !c(cell_type, gene),
    names_to = c(".value", "group"),
    names_pattern = "(.+)_group(.+)"
  ) %>% 
  group_by(group, cell_type) %>% 
  mutate(p_adj = p.adjust(p, method = "fdr")) %>% 
  ungroup() %>% 
  calc_expression_frequency()



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
#' @param group The selected contrast.
#' @param cell_type The selected cell type
#' @param dbs Character vector containing valid enrichr databases
#'   as returned by `enrichR::listEnrichrDbs()`.
#' @param direction Use genes that are up- or downregulated, respectively.
#'
#' @return A dataframe, which combines the dataframes returned by
#'   `enrichR::enrichr()` (empty results are removed) by adding a column "db".
enrich_genes <- function(data,
                         group,
                         cell_type,
                         dbs,
                         direction = c("up", "down")) {
  info("Cell type {cell_type}, group {group}, {direction}regulated genes")
  
  direction <- match.arg(direction)
  
  top_genes <- 
    data %>% 
    filter(
      cell_type == {{cell_type}},
      group == {{group}},
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
#'   "group", "cell_type", and "db".
enrich_all_genes <- function(data, dbs) {
  queries <- 
    data %>%
    distinct(group, cell_type) %>% 
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
    mutate(group = as_factor(group)) %>% 
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
#' @return A dataframe, comprising columns "db", "group", "cell_type", as well
#'   as all columns in the result of `fgseaMultilevel()`.
perform_gsea <- function(data, gene_sets) {
  data %>% 
    distinct(group, cell_type) %>% 
    mutate(db = list(names(gene_sets))) %>% 
    unnest_longer(db) %>% 
    pmap_dfr(
      function(group, cell_type, db) {
        ranked_genes <-
          data %>% 
          filter(group == {{group}}, cell_type == {{cell_type}}) %>%
          select(gene, logFC) %>%
          deframe()
        ranked_genes <- ranked_genes[!is.na(ranked_genes)]
        
        info("GSEA of group {group}, cell type {cell_type}, ",
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
            group = {{group}},
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



# DGE in tumor cluster ----------------------------------------------------

nb_tumor <-
  nb[
    ,
    colData(nb)$group %in% c("II", "IV")
    & colData(nb)$cellont_abbr == "NB"
  ]
colData(nb_tumor)$tumor_subcluster <- 
  read_csv("metadata/nb_subclusters.csv") %>% 
  deframe() %>% 
  magrittr::extract(rownames(colData(nb_tumor)))
colData(nb_tumor)$group <- fct_drop(colData(nb_tumor)$group)
colData(nb_tumor)$sample <- fct_drop(colData(nb_tumor)$sample)


analyse_dge_tumor <- function(subclusters) {
  info("Analysing tumor subclusters {str_c(subclusters, collapse = '+')}")
  
  nb_sub <- nb_tumor[, colData(nb_tumor)$tumor_subcluster %in% subclusters]
  
  model_mat <- model.matrix(
    ~group + tif,
    data =
      colData(nb_sub) %>%
      as_tibble(rownames = "cell") %>% 
      mutate(group = group %>% fct_drop() %>% fct_rev()) %>%
      column_to_rownames("cell")
  )
  colnames(model_mat) <-
    colnames(model_mat) %>%
    str_replace("group(.+)", "\\1_vs_IV")
  
  data_grouped <- group_cell(
    counts(nb_sub),
    id = colData(nb_sub)$sample,
    pred = model_mat,
    offset = colData(nb_sub)$Size_Factor
  )
  
  res <- nebula(
    data_grouped$count,
    id = data_grouped$id,
    pred = data_grouped$pred,
    offset = data_grouped$offset,
    verbose = TRUE
  )
  
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
dge_results_tumor <- map_dfr(
  list("all" = 1:4, "1" = 1, "2" = 2, "3" = 3, "4" = 4),
  analyse_dge_tumor,  
  .id = "subclusters"
)


# similar to `calc_expression_frequency()`
calc_expression_frequency_tumor <- function(data) {
  cds <- prepSCE(
    nb_tumor,
    kid = "tumor_subcluster",
    gid = "group",
    sid = "sample",
    drop = TRUE
  )
  
  # generate table with frequencies on sample level
  exp_frq_subclusters <-
    cds %>% 
    calcExprFreqs() %>% 
    assays() %>%
    as.list() %>%
    map_dfr(as_tibble, rownames = "gene", .id = "subclusters")
  
  colData(cds)$cluster_id <- factor("all")
  
  exp_frq_all <-
    cds %>% 
    calcExprFreqs() %>% 
    assay("all") %>%
    as_tibble(rownames = "gene") %>% 
    mutate(subclusters = "all")
  
  exp_frq <- 
    bind_rows(exp_frq_subclusters, exp_frq_all) %>% 
    select(!c(II, IV))
  exp_frq
  
  # summarise at group level
  exp_frq_groups <-
    exp_frq %>% 
    pivot_longer(starts_with("20"), names_to = "sample", values_to = "frq") %>% 
    left_join(
      nb_metadata %>% distinct(sample, group),
      by = "sample"
    ) %>%
    group_by(subclusters, gene, group) %>%
    summarise(frq = min(frq)) %>%
    pivot_wider(
      names_from = group,
      names_prefix = "frq_",
      values_from = frq
    ) %>% 
    rename(frq = frq_II, frq_ref = frq_IV)
  
  # add columns to input
  data %>%
    left_join(exp_frq_groups, by = c("subclusters", "gene"))
}

dge_results_tumor_wide <- 
  dge_results_tumor %>% 
  select(subclusters, gene, logFC = logFC_II_vs_IV, p = p_II_vs_IV) %>% 
  group_by(subclusters) %>%
  mutate(p_adj = p.adjust(p, method = "fdr")) %>%
  ungroup() %>%
  calc_expression_frequency_tumor()



# Save results ------------------------------------------------------------

# nb <- dge$cds
# nb_metadata <- dge$metadata
# used_clusters <- dge$used_clusters
# dge_results <- dge$results
# dge_results_wide <- dge$results_wide
# dge_results_wide_filtered <- dge$results_wide_filtered
# enrichr_results <- dge$enrichr
# gsea_results <- dge$gsea
# enrichr_genesets <- dge$gene_sets
# nb_tumor <- dge$cds_tumor
# dge_results_tumor <- dge$results_tumor
# dge_results_tumor_wide <- dge$results_tumor_wide

list(
  cds = nb,
  metadata = nb_metadata,
  used_clusters = used_clusters,
  results = dge_results,
  results_wide = dge_results_wide,
  results_wide_filtered = dge_results_wide_filtered,
  enrichr = enrichr_results,
  gsea = gsea_results,
  gene_sets = enrichr_genesets,
  cds_tumor = nb_tumor,
  results_tumor = dge_results_tumor,
  results_tumor_wide = dge_results_tumor_wide
) %>%
  saveRDS("data_generated/dge_mm_results.rds")
