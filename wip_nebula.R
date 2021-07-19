# Differential gene expression analysis using mixed models.
#
# @DEPI rna_decontaminated.rds
# @DEPI metadata.rds
# @DEPO dge_mm_results_[group]_[cluster].rds
# @DEPO dge_mm_results.RData

library(monocle3)
library(muscat)
library(scater)
library(nebula)
library(tidyverse)
source("common_functions.R")

library(msigdbr)
library(enrichR)
library(latex2exp)
library(fgsea)



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
dge_mm_results <- map_dfr(
  colData(nb)$cellont_abbr %>%
    unique() %>%
    setdiff("NB") %>%
    set_names(),
  analyse_dge,
  .id = "cell_type"
)

saveRDS(dge_mm_results, "data_generated/dge_mm_results.rds")



# Filter results ----------------------------------------------------------

#' Filter DGE results.
#'
#' @param data A data frame with MM DGE results.
#' @param max_p Maximum p value.
#' @param max_p_adj Maximum adjusted p value.
#' @param min_abs_log_fc Minimum absolute log fold change.
#' @param max_abs_log_fc Maximum absolute log fold change.
#' @param remove_ribosomal If true, remove ribosomal proteins.
#'
#' @return A data frame with an additional column "direction".
filter_dge_results <- function(data,
                               max_p = Inf,
                               max_p_adj = 0.05,
                               min_abs_log_fc = 1,
                               max_abs_log_fc = 10,
                               remove_ribosomal = TRUE) {
  res <-
    data %>%
    filter(
      p <= max_p,
      p_adj <= max_p_adj,
      abs(logFC) %>% between(min_abs_log_fc, max_abs_log_fc)
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


# Plot results ------------------------------------------------------------

## Prepare data ----

load("data_generated/dge_pb_results.RData")

# are there problematic enes with convergence <= -20
any(dge_mm_results$convergence <= -20)

dge_mm_results_wide <- 
  dge_mm_results %>% 
  select(cell_type, gene, starts_with("logFC_g"), starts_with("p_g")) %>% 
  pivot_longer(
    !c(cell_type, gene),
    names_to = c(".value", "group"),
    names_pattern = "(.+)_group(.+)"
  ) %>% 
  group_by(group, cell_type) %>% 
  mutate(p_adj = p.adjust(p, method = "fdr")) %>% 
  ungroup()


## Volcano plots ----

ggplot(dge_mm_results_wide, aes(logFC, -log10(p))) +
  geom_point(size = 0.1) +
  facet_grid(vars(group), vars(cell_type)) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 25)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("dge_mm/volcano")


## Comparison to pseudobulk ----

dge_results %>% 
  filter(contrast != "tif") %>% 
  select(gene, cell_type = cluster_id,
         contrast, logFC_pb = logFC, p_pb = p_val) %>% 
  extract(contrast, into = "group", regex = "(.+)_vs_I") %>% 
  left_join(
    dge_mm_results_wide %>% rename(logFC_mm = logFC, p_mm = p),
    by = c("gene", "cell_type", "group")
  ) %>% 
  ggplot(aes(logFC_pb, logFC_mm)) +
  geom_point(size = 0.1) +
  facet_grid(vars(group), vars(cell_type)) +
  # coord_fixed() +
  coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("dge_mm/comparison_pb_mm")
# ggsave_default("dge_mm/comparison_pb_mm_full")




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

#' Perform enrichment anlalysis for all contrasts and cell types in a dataset.
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
  dge_mm_results_wide %>%
  filter_dge_results() %>%
  enrich_all_genes(enrichr_dbs)



# Plot enrichment ---------------------------------------------------------

#' Generate a dotplot for enrichment terms.
#' 
#' Dots that meet given limits of adjusted p value, overlap size, and odds ratio
#' are bordered.
#'
#' @param data Results from `enrich_all_genes()`.
#' @param db Enrichr database for which results should be plotted.
#' @param direction Plot results for `"up"`- or `"down"`-regulated genes.
#' @param top_n Only plot the terms with the highest odds ratio per cluster.
#' @param max_p_adj Maximum adjusted p value, …
#' @param min_odds_ratio … minimum odds ratio, and …
#' @param min_overlap_size … minimum overlap size required for bordered dots.
#' @param log_odds_cap Upper boundary of the color scale.
#' @param filename Name of output file. If `"auto"`, derive from database name.
#' @param ... Additional parameters passed to `ggsave_default()`.
#'
#' @return A ggplot object.
plot_enrichr_dots <- function(data,
                              db,
                              direction = "up",
                              top_n = 5,
                              max_p_adj = 0.05,
                              min_odds_ratio = 5,
                              min_overlap_size = 1,
                              log_odds_cap = 2,
                              filename = "auto",
                              ...) {
  data_selected <-
    data %>% 
    filter(
      db == {{db}},
      direction == {{direction}},
      Odds.Ratio >= 1
    )
  
  top_terms <- 
    data_selected %>% 
    group_by(group, cell_type) %>%
    slice_max(n = top_n, order_by = Odds.Ratio, with_ties = FALSE) %>%
    pull(Term) %>% 
    unique()
  
  data_vis <- 
    data_selected %>%
    filter(Term %in% top_terms) %>% 
    mutate(
      is_significant =
        Adjusted.P.value <= max_p_adj &
        Odds.Ratio >= min_odds_ratio &
        overlap_size >= min_overlap_size,
      Odds.Ratio = pmin(Odds.Ratio, 10^log_odds_cap),
      Term =
        as_factor(Term) %>%
        fct_reorder(str_c(cell_type, group), n_distinct)
    )
  
  p <- 
    ggplot(data_vis, aes(group, Term, size = -log10(Adjusted.P.value))) +
    geom_point(aes(color = log10(Odds.Ratio))) +
    geom_point(data = data_vis %>% filter(is_significant), shape = 1) +
    scale_color_distiller(
      palette = "Reds",
      direction = 1,
      limits = c(0, log_odds_cap)
    ) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(cell_type), drop = FALSE, nrow = 1) +
    labs(
      y = "",
      color = TeX("log_{10} (odds ratio)"),
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue( "Enrichr results ({db})"),
      caption = str_glue(
        "top {top_n} terms per cell type; ",
        "bordered circles: adjusted p value <= {max_p_adj}, ",
        "odds ratio >= {min_odds_ratio}, ",
        "overlap size >= {min_overlap_size}; ",
        "color scale capped at {log_odds_cap}"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  if (filename == "auto") {
    filename <- str_glue("dge_mm/enrichr_{db}")
  }
  ggsave_default(filename, ...)
  p
}

plot_enrichr_dots(enrichr_results,
                  db = "MSigDB_Hallmark_2020",
                  width = 420)



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
  dge_mm_results_wide %>%
  filter_dge_results(max_p_adj = Inf, min_abs_log_fc = 0) %>% 
  perform_gsea(enrichr_genesets)



# Plot GSEA ---------------------------------------------------------------

## GSEA results ----

gsea_results %>%
  ggplot(aes(NES, -log10(padj))) +
  geom_point(alpha = .25)


#' Generate a dotplot for terms enriched by GSEA.
#' 
#' Select the top "important" terms (i.e., terms that meet given limits of
#' adjusted p-value and normalized enrichment score, NES) per contrast and
#' cluster. Plot dots for these terms in all conditions, border the important
#' ones.
#'
#' @param data Results from `perform_gsea()`.
#' @param db Enrichr database for which results should be plotted.
#' @param top_n_positive Only plot n terms with the highest NES per cluster.
#' @param top_n_negative Only plot n terms with the lowest NES per cluster.
#' @param max_p_adj Maximum adjusted p value and …
#' @param min_NES minimum absolute NES required for bordered dots.
#' @param filename Name of output file. If `"auto"`, derive from database name.
#' @param ... Additional parameters passed to `ggsave_default()`.
#'
#' @return A ggplot object.
plot_gsea_dots <- function(data,
                           db,
                           top_n_positive = 5L,
                           top_n_negative = 5L,
                           max_p_adj = 0.05,
                           min_abs_NES = 1,
                           filename = "auto",
                           ...) {
  data_top_terms <-
    data %>% 
    filter(db == {{db}}, padj <= max_p_adj, abs(NES) >= min_abs_NES) %>% 
    group_by(group, cell_type)
  
  top_terms_pos <- 
    data_top_terms %>% 
    slice_max(n = top_n_positive, order_by = NES, with_ties = FALSE) %>%
    pull(pathway) %>%
    unique()
  
  top_terms_neg <- 
    data_top_terms %>% 
    slice_min(n = top_n_negative, order_by = NES, with_ties = FALSE) %>%
    pull(pathway) %>%
    unique()
  
  data_vis <- 
    data %>% 
    filter(db == {{db}}, pathway %in% c(top_terms_pos, top_terms_neg)) %>% 
    mutate(
      is_significant =
        padj <= max_p_adj &
        abs(NES) >= min_abs_NES,
      pathway =
        as_factor(pathway) %>%
        # fct_reorder(NES * is_significant, sum, na.rm = TRUE)
        fct_reorder(NES * -log10(padj), sum, na.rm = TRUE)
    )
  
  if (nlevels(data_vis$pathway) > 5) {
    horizontal_grid <-
      geom_hline(
        yintercept = seq(5, nlevels(data_vis$pathway), 5),
        size = 0.5,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  color_limit <- max(abs(data_vis$NES))
  
  p <- 
    ggplot(data_vis, aes(group, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    geom_point(data = data_vis %>% filter(is_significant), shape = 1) +
    scale_color_distiller(
      palette = "PiYG",
      direction = 1,
      limits = c(-color_limit, color_limit)
    ) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(cell_type), nrow = 1) +
    labs(
      y = "",
      color = "normalized\nenrichment\nscore",
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue("GSEA results ({db})"),
      caption = str_glue(
        "top {top_n_positive} positively and top {top_n_negative} negatively ",
        "enriched terms per cell_type; ",
        "bordered circles: adjusted p value <= {max_p_adj}, ",
        "absolute NES >= {min_abs_NES}"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  if (filename == "auto") {
    filename <- str_glue("dge_mm/gsea_{db}")
  }
  ggsave_default(filename, width = 400, ...)
  p
}

plot_gsea_dots(gsea_results,
               db = "MSigDB_Hallmark_2020")
