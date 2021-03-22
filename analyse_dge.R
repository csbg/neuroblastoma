# Differential gene expression analysis.
#
# @DEPI rna_decontaminated.rds
# @DEPI metadata.rds

library(muscat)
library(scater)
library(tidyverse)
library(viridis)
library(patchwork)
library(ggpmisc)
library(ComplexHeatmap)
library(enrichR)
library(latex2exp)
library(fgsea)
library(msigdbr)
source("common_functions.R")



# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  logNormCounts(assay.type = "soupx_counts")

nb@colData <-
  readRDS("data_generated/metadata.rds") %>%
  mutate(
    sample =
      str_c(group, sample, sep = ".") %>% 
      as_factor() %>%
      fct_relevel(str_sort) %>% 
      fct_relabel(~str_extract(.x, "\\d+_\\d+")),
    Size_Factor = colData(nb)$Size_Factor
  ) %>% 
  column_to_rownames("cell") %>% 
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)

# TODO: remove this! find a way that cluster 9 only comprises nb cells
ignored_barcodes <-
  readRDS("data_generated/metadata.rds") %>%
  filter(cluster_50 == "9", group == "I") %>%
  pull(cell)
nb <- nb[, !colnames(nb) %in% ignored_barcodes]



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

exp_info <-
  metadata(pb)$experiment_info %>%
  column_to_rownames("sample_id")
model_mat <- model.matrix(~ 0 + group_id, data = exp_info)
colnames(model_mat) <- levels(exp_info$group_id)
model_mat


dge <- pbDS(
  pb,
  design = model_mat,
  contrast = limma::makeContrasts(
    II_vs_I = II - I,
    III_vs_I = III - I,
    IV_vs_I = IV - I,
    levels = model_mat
  ),
  min_cells = 1
)

dge_results <-
  resDS(nb, dge, frq = TRUE) %>%
  as_tibble() %>%
  mutate(
    cluster_id =
      as_factor(cluster_id) %>%
      fct_expand(levels(colData(nb)$cluster_id)) %>%
      fct_inseq()
  )


## Approach 2: Cell-level mixed model ----

# TBD
# dge_mm <- mmDS(nb)



# Filter DGE results ------------------------------------------------------

#' Filter DGE results
#'
#' @param data A dataframe as returned by `muscat::resDS()`.
#' @param max_p Maximum p value.
#' @param max_p_adj Maximum adjusted p value.
#' @param min_abs_log_fc Minimum absolute log fold change.
#' @param min_freq Minimum gene expression frequency (fracions of cells that
#'   expresses a given gene.)
#' @param remove_ribosomal If true, remove ribosomal proteins.
#'
#' @return A data frame.
filter_dge_results <- function(data,
                               max_p = Inf,
                               max_p_adj = 0.05,
                               min_abs_log_fc = 1,
                               min_freq = 0.1,
                               remove_ribosomal = TRUE) {
  res <-
    data %>%
    filter(
      I.frq > min_freq |
      (contrast == "II_vs_I" & II.frq > min_freq) |
      (contrast == "III_vs_I" & III.frq > min_freq) |
      (contrast == "IV_vs_I" & IV.frq > min_freq)
    ) %>%
    filter(
      p_val < max_p,
      p_adj.loc < max_p_adj,
      abs(logFC) > min_abs_log_fc
    )
  
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

dge_results_filtered <- filter_dge_results(dge_results)



# Plot results ------------------------------------------------------------

## Volcano plot ----

#' Make one volcano plot per cluster for a given contrast.
#'
#' @param data A DGE dataset.
#' @param contrast A contrast present in the dataset.
#' @param max_p_adj see `filter_dge_results()`.
#' @param min_abs_log_fc see `filter_dge_results()`.
#' @param min_freq see `filter_dge_results()`.
#' @param point_colors A named character vector specifying colors for different
#'   classes of genes.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_volcano <- function(data,
                         contrast,
                         max_p_adj = 0.05,
                         min_abs_log_fc = 1,
                         min_freq = 0.1,
                         filename = NULL) {
  group_left <- str_extract(contrast, "^[^_]+")
  group_right <- str_extract(contrast, "[^_]$")
  
  data_filtered <- 
    data %>%
    filter(contrast == {{contrast}}) %>% 
    mutate(
      is_significant =
        p_adj.loc < max_p_adj &
        abs(logFC) > min_abs_log_fc,
      is_frequent =
        !!sym(str_glue("{group_left}.frq")) > min_freq |
        !!sym(str_glue("{group_right}.frq")) > min_freq,
      diff_expressed = is_significant & is_frequent
    )
  
  genes_count <- 
    data_filtered %>% 
    count(cluster_id, diff_expressed) %>% 
    pivot_wider(names_from = "diff_expressed", values_from = "n")

  p <-
    ggplot(data_filtered, aes(logFC, -log10(p_val))) +
    geom_point(
      aes(color = diff_expressed),
      alpha = .5,
      size = 0.5
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = `TRUE`),
      npcx = 0.05,
      npcy = 0.9,
      size = 3,
      color = "black",
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = `FALSE`),
      npcx = 0.05,
      npcy = 0.8,
      size = 3,
      color = "gray80",
      na.rm = TRUE
    ) +
    scale_color_manual(
      "differentially\nexpressed",
      values = c("gray80", "black")
    ) +
    facet_wrap(vars(cluster_id), ncol = 6, drop = FALSE) +
    labs(
      title = str_glue("DGE using contrast {contrast}"),
      caption = str_glue(
        "gene expression frequency > {min_freq}, ",
        "adjusted p value < {max_p_adj}, ",
        "absolute log fold change > {min_abs_log_fc}"
      )
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  ggsave_default(filename)
  p
}

plot_volcano(dge_results, "II_vs_I", filename = "dge/volcano_II")
plot_volcano(dge_results, "III_vs_I", filename = "dge/volcano_III")
plot_volcano(dge_results, "IV_vs_I", filename = "dge/volcano_IV")



## Upset plot ----

#' Plot an upset plot showing how many DGE genes are shared between the
#' clusters.
#'
#' @param data A DGE dataset.
#' @param contrast A contrast present in the dataset.
#' @param filename Name of output file.
#'
#' @return An object of class Heatmap.
plot_upset <- function(data, contrast, filename = NULL) {
  p <-
    data %>%
    filter(contrast == {{contrast}}) %>%
    select(cluster_id, gene) %>%
    chop(gene) %>%
    deframe() %>%
    make_comb_mat() %>% 
    UpSet(
      comb_order = order(comb_size(.), decreasing = TRUE),
      comb_col = viridis(nrow(.))[comb_degree(.)],
      column_title = "Number of DE genes shared between clusters",
      row_title = "cluster"
    )
  
  ggsave_default(filename, plot = p, width = 297, height = 100)
  p
}

plot_upset(dge_results_filtered, "II_vs_I", filename = "dge/upset_II")
plot_upset(dge_results_filtered, "III_vs_I", filename = "dge/upset_III")
plot_upset(dge_results_filtered, "IV_vs_I", filename = "dge/upset_IV")



## Violin plot ----

#' Title
#'
#' @param data A DGE dataset.
#' @param contrast A contrast and ...
#' @param cluster ... cluster for which all samples should be plotted.
#' @param top_n Plot the top differentially expressed genes.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_violin <- function(data, contrast, cluster, top_n = 12, filename = NULL) {
  group_left <- str_extract(contrast, "^[^_]+")
  group_right <- str_extract(contrast, "[^_]$")
  
  top_genes <- 
    data %>%
    filter(contrast == {{contrast}}, cluster_id == {{cluster}}) %>% 
    slice_max(n = top_n, logFC, with_ties = FALSE) %>% 
    pull(gene)
  
  p <-
    plotExpression(
      nb[, nb$cluster_id == cluster &
           nb$group_id %in% c(group_left, group_right)],
      top_genes,
      x = "sample_id",
      colour_by = "group_id",
      ncol = 4,
      point_size = 1
    ) +
    ggtitle(str_glue("Top {top_n} upregulated genes in cluster {cluster} ",
                     "of contrast {contrast}")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave_default(filename)
  p
}

plot_violin(dge_results_filtered, "II_vs_I", cluster = "1",
            filename = "dge/violin_II_1")
plot_violin(dge_results_filtered, "II_vs_I", cluster = "13",
            filename = "dge/violin_II_13")
plot_violin(dge_results_filtered, "IV_vs_I", cluster = "13",
            filename = "dge/violin_IV_13")



## Pseudo-bulk heatmap ----

#' Plot a pseudobulk heatmap.
#'
#' @param data A DGE dataset.
#' @param pb_data Pseudobulk data as generated by `muscat::aggregateData()`.
#' @param contrast The selected contrast.
#' @param top_n Plot the top differentially expressed genes.
#' @param clusters A character vector denoting clusters that should be included.
#' @param filename Name of output file.
#'
#' @return An object of class Heatmap.
plot_pbheatmap <- function(data,
                           pb_data,
                           contrast,
                           top_n = Inf,
                           clusters = NULL,
                           filename = NULL) {
  group_left <- str_extract(contrast, "^[^_]+")
  group_right <- str_extract(contrast, "[^_]$")
  
  top_genes <- 
    data %>% 
    filter(contrast == {{contrast}}) %>% 
    {
      if (!is.null(clusters))
        filter(., cluster_id %in% {{clusters}})
      else
        .
    } %>% 
    group_by(cluster_id) %>% 
    slice_max(logFC, n = top_n, with_ties = FALSE)
  
  sample_data <- 
    colData(pb_data) %>%
    as_tibble(rownames = "sample") %>%
    filter(group_id %in% c(group_left, group_right))
  
  mat <- 
    top_genes %>% 
    group_map(~assays(pb_data)[[.y$cluster_id]][.x$gene, sample_data$sample]) %>%
    reduce(rbind) %>% 
    t() %>% 
    scale() %>% 
    t()
  
  p <-
    Heatmap(
      mat,
      col = viridis::viridis(10),
      
      cluster_rows = FALSE,
      show_row_names = FALSE,
      row_split = top_genes$cluster_id,
      row_title_rot = 0,
      
      cluster_columns = FALSE,
      column_split = sample_data$group_id,
      column_title = "group",
      
      left_annotation = rowAnnotation(
        logFC = top_genes$logFC,
        col = list(
          logFC = circlize::colorRamp2(
            breaks = c(min(top_genes$logFC), 0, max(top_genes$logFC)),
            colors = c("blue", "white", "red")
          )
        )
      ),
      
      top_annotation = HeatmapAnnotation(
        foo = anno_block(
          gp = gpar(fill = 3:4),
          labels = unique(sample_data$group_id),
        )
      ),

      heatmap_legend_param = list(title = "scaled\nmean\nlogcounts")
    ) %>% 
    draw(
      row_title = "cluster",
      column_title = "Pseudobulk heatmap",
      column_title_gp = gpar(fontface = "bold", fontsize = 14)
    )
  
  ggsave_default(filename, plot = p, width = 150, height = 297)
  p
}

pb_log <- aggregateData(
  nb,
  assay = "logcounts",
  fun = "mean",
  by = c("cluster_id", "sample_id")
)

plot_pbheatmap(dge_results_filtered, pb_log,
               contrast = "II_vs_I", filename = "dge/heatmap_II")
plot_pbheatmap(dge_results_filtered, pb_log,
               contrast = "III_vs_I", filename = "dge/heatmap_III")
plot_pbheatmap(dge_results_filtered, pb_log,
               contrast = "IV_vs_I", filename = "dge/heatmap_IV")



# Enrichment analysis -----------------------------------------------------

#' Perform enrichment analysis via enrichr for upregulated genes in a given
#' contrast and cluster.
#'
#' @param data A DGE dataset.
#' @param contrast The selected contrast.
#' @param cluster The selected cluster.
#' @param dbs Character vector containing valid enrichr databases
#'   as returned by `enrichR::listEnrichrDbs()`.
#' @param direction Use genes that are up- or downregulated, respectively.
#'
#' @return A dataframe, which combines the dataframes returned by
#'   `enrichR::enrichr()` (empty results are removed) by adding a column "db".
enrich_genes <- function(data,
                         contrast,
                         cluster,
                         dbs,
                         direction = c("up", "down")) {
  message("Cluster ", cluster, ", contrast ", contrast, ", ",
          direction, "regulated genes")
  
  direction <- match.arg(direction)
  
  top_genes <- 
    data %>% 
    filter(
      cluster_id == {{cluster}},
      contrast == {{contrast}},
      if (direction == "up") logFC > 0 else logFC < 0
    ) %>% 
    pull(gene)
  
  if (length(top_genes) > 0) {
    message("  genes for enrichr: ", str_c(top_genes, collapse = ", "))
    enrichr(top_genes, dbs) %>% 
      keep(~nrow(.) > 0) %>% 
      bind_rows(.id = "db") 
  }
}

#' Perform enrichment anlalysis for all contrasts and clusters in a dataset.
#'
#' @param data A DGE dataset.
#' @param dbs Character vector containing valid enrichr databases
#'   as returned by `enrichR::listEnrichrDbs()`.
#'
#' @return A dataframe, which combines the dataframes returned by
#'   `enrichR::enrichr()` (empty results are removed) by adding three columns
#'   "contrast", "cluster", and "db".
enrich_all_genes <- function(data,
                             dbs = c(
                               "GO_Biological_Process_2018",
                               "GO_Cellular_Component_2018",
                               "GO_Molecular_Function_2018",
                               "KEGG_2019_Human",
                               "WikiPathways_2019_Human"
                             )) {
  queries <- 
    data %>%
    distinct(contrast, cluster = cluster_id) %>% 
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
    mutate(contrast = as_factor(contrast) %>% fct_relevel(str_sort)) %>% 
    separate(
      Overlap,
      into = c("overlap_size", "geneset_size"),
      convert = TRUE
    )
}

enrichr_results <- enrich_all_genes(dge_results_filtered)



#' Generate a dotplot for enrichment terms.
#' 
#' Only enrichment terms with a minimum overlap size are included in the plot.
#' Dots that meet given p value and odds ratio requirements are bordered.
#'
#' @param data Results from `enrich_all_genes()`.
#' @param db Enrichr database for which results should be plotted.
#' @param max_p_adj Maximum adjusted p value, and …
#' @param min_odds_ratio … minimum odds ratio required for bordered dots.
#' @param min_overlap_size Minimum overlap size of an included term.
#' @param log_odds_cap Upper boundary of the color scale.
#' @param direction Plot results for `"up"`- or `"down"`-regulated genes.
#' @param filename Name of output file.
#' @param ... Additional parameters passed to `ggsave_default()`.
#'
#' @return A ggplot object.
plot_enrichr_dots <- function(data,
                              db,
                              max_p_adj = 0.05,
                              min_odds_ratio = 10,
                              min_overlap_size = 2,
                              log_odds_cap = 2,
                              direction = "up",
                              filename = NULL,
                              ...) {
  data_filtered <-
    data %>% 
    filter(
      db == {{db}},
      direction == {{direction}},
      overlap_size >= min_overlap_size
    ) %>% 
    group_by(Term) %>%
    filter(
      min(Adjusted.P.value) < max_p_adj,
      max(Odds.Ratio) > min_odds_ratio
    ) %>%
    ungroup() %>%
    select(contrast, cluster, Term, Odds.Ratio, Adjusted.P.value) %>% 
    mutate(
      Term = as_factor(Term) %>% fct_reorder(Odds.Ratio, max),
      is_significant =
        Adjusted.P.value < max_p_adj & Odds.Ratio > min_odds_ratio
    ) %>% 
    mutate(Odds.Ratio = pmin(Odds.Ratio, 10^log_odds_cap))
    
  p <- 
    ggplot(data_filtered, aes(cluster, Term, size = -log10(Adjusted.P.value))) +
    geom_point(aes(color = log10(Odds.Ratio))) +
    geom_point(data = data_filtered %>% filter(is_significant), shape = 1) +
    scale_color_distiller(
      palette = "Reds",
      direction = 1,
      limits = c(0, log_odds_cap)
    ) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(contrast), drop = FALSE) +
    labs(
      y = "",
      color = TeX("log_{10} (odds ratio)"),
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue( "Enrichr results ({db})"),
      caption = str_glue(
        "minimum overlap size: {min_overlap_size}; ",
        "bordered circles: adjusted p value < {max_p_adj}, ",
        "odds ratio > {min_odds_ratio}; ",
        "color scale capped at {log_odds_cap}"
      )
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  ggsave_default(filename, ...)
  p
}


plot_enrichr_dots(enrichr_results,
                  db = "GO_Biological_Process_2018",
                  min_odds_ratio = 50,
                  filename = "dge/dots_bioproc",
                  height = 840)

plot_enrichr_dots(enrichr_results,
                  min_odds_ratio = 10,
                  db = "GO_Cellular_Component_2018",
                  filename = "dge/dots_cellcomp")

plot_enrichr_dots(enrichr_results,
                  db = "GO_Molecular_Function_2018",
                  min_odds_ratio = 30,
                  filename = "dge/dots_molfun",
                  width = 420)

plot_enrichr_dots(enrichr_results,
                  db = "KEGG_2019_Human",
                  min_odds_ratio = 10,
                  filename = "dge/dots_kegg")

plot_enrichr_dots(enrichr_results,
                  db = "WikiPathways_2019_Human",
                  min_odds_ratio = 30,
                  filename = "dge/dots_wiki",
                  width = 420,
                  height = 250)



# TODO: select interesting genes from enrichment analysis
enrichr_results %>% 
  filter(
    contrast == "II_vs_I",
    db == "KEGG_2019_Human",
    direction == "up",
    # cluster == "1",
    Adjusted.P.value < 0.05,
    # Odds.Ratio > 10,
    # overlap_size >= 2
  ) %>% 
  # pull(Genes) %>%
  # str_split(";") %>%
  # flatten_chr() %>%
  # as_factor() %>%
  # fct_count(sort = TRUE) %>%
  View() %>% 
  {.}






# GSEA --------------------------------------------------------------------

get_enrichr_genesets <- function(dbs) {
  dbs %>% 
    map(
      function(db) {
        message("Downloading ", db)
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

enrichr_genesets <- get_enrichr_genesets(c("GO_Molecular_Function_2018",
                                           "WikiPathways_2019_Human"))



do_gsea <- function(data, gene_sets) {
  data %>% 
    distinct(contrast, cluster_id) %>% 
    mutate(gene_set = list(names(gene_sets))) %>% 
    unnest_longer(gene_set) %>% 
    pmap_dfr(
      function(contrast, cluster_id, gene_set) {
        ranked_genes <-
          data %>% 
          filter(contrast == {{contrast}}, cluster_id == {{cluster_id}}) %>%
          select(gene, logFC) %>%
          deframe()
        ranked_genes <- ranked_genes[!is.na(ranked_genes)]
        
        message("GSEA of contrast ", contrast, ", cluster ", cluster_id,
                ", gene set ", gene_set, " (", length(ranked_genes), " genes)")
        
        fgseaMultilevel(
          gene_sets[[gene_set]],
          ranked_genes,
          eps = 0
        ) %>%
          as_tibble() %>%
          mutate(
            gene_set = {{gene_set}},
            contrast = {{contrast}},
            cluster = {{cluster_id}},
            .before = 1
          )
      }
    )
}


gsea_results <-
  dge_results %>%
  filter_dge_results(min_abs_log_fc = 0, min_freq = 0.05) %>%
  do_gsea(enrichr_genesets)
  



gsea_results %>% 
  ggplot(aes(NES, -log10(padj))) +
  geom_point(alpha = .1)




plot_gsea_dots <- function(data,
                           db,
                           max_p_adj = 0.05,
                           min_NES = 1,
                           filename = NULL,
                           plot_params = list()) {
  data_filtered <-
    data %>% 
    filter(db == {{db}}) %>% 
    group_by(pathway) %>%
    filter(
      min(padj) < max_p_adj,
      max(NES) > min_NES
    ) %>%
    ungroup() %>%
    select(contrast, cluster, pathway, NES, padj) %>% 
    mutate(
      pathway = as_factor(pathway) %>% fct_reorder(NES, max),
      is_significant =
        padj < max_p_adj & NES > min_NES
    )
  
  p <- 
    ggplot(data_filtered, aes(cluster, pathway, size = -log10(padj))) +
    geom_point(aes(color = NES)) +
    geom_point(data = data_filtered %>% filter(is_significant), shape = 1) +
    scale_color_distiller(palette = "Reds", direction = 1) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(contrast), drop = FALSE) +
    labs(
      y = "",
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue( "GSEA results ({db})"),
      caption = str_glue(
        "bordered circles: adjusted p value < {max_p_adj}, ",
        "normalized enrichment score > {min_NES}"
      )
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  rlang::exec(ggsave_default, filename, !!!plot_params)
  p
}

gsea_results %>% 
  plot_gsea_dots("GO_Molecular_Function_2018",
                 filename = "dge/gsea_molfun")






top_terms <- 
  gsea_molfun_II_1 %>% 
  arrange(padj) %>% 
  head(20) %>% 
  pull(pathway)

plotEnrichment(
  pathways$GO_Molecular_Function_2018[[
    gsea_molfun_II_1 %>% 
      arrange(padj) %>% 
      magrittr::extract2(1, 1)
  ]],
  dge_results %>%
    filter_dge_results(min_abs_log_fc = 0, min_freq = 0.05, max_p_adj = 0.05) %>%
    filter(contrast == "II_vs_I", cluster_id == "1") %>%
    select(gene, logFC) %>% 
    deframe(),
)

plotGseaTable(
  pathways$GO_Molecular_Function_2018[top_terms],
  dge_results %>%
    filter_dge_results(min_abs_log_fc = 0, min_freq = 0.05, max_p_adj = 0.05) %>%
    filter(contrast == "II_vs_I", cluster_id == "1") %>%
    select(gene, logFC) %>% 
    deframe(),
  gsea_molfun_II_1,
  gseaParam = 0.5
)


read_lines("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Molecular_Function_2018")



# Unused ------------------------------------------------------------------


#' Shorten enrichment terms.
#' 
#' Enrichment terms are shortened depending on their type:
#' * For GO terms, the part preceding the GO identifier is truncated
#'   to at most `max_length` characters.
#' * For other terms, the complete term is truncated to at most `max_length`
#'   characters. The function ensures that truncated terms are unique.
#'
#' @param s A character vector containing enrichment terms.
#' @param max_length Maximum length of the shortened term.
#'
#' @return A character vector containing truncated terms.
shorten_go_term <- function(s, max_length = 30) {
  str_match(s, "(.*)( \\(GO:\\d+\\))") %>% 
    magrittr::set_colnames(c("match", "description", "goid")) %>% 
    as_tibble() %>% 
    mutate(
      description = coalesce(description, s),
      short = case_when(
        str_length(description) > max_length ~
          description %>% str_sub(end = max_length - 1) %>% paste0("…"),
        TRUE ~ description
      ),
      goid = coalesce(goid, ""),
      result = str_c(short, goid) %>% make.unique(sep = "")
    ) %>% 
    pull(result)
}


#' Plot enrichment results
#'
#' @param data Results from enrichtment analysis.
#' @param contrast The selected contrast.
#' @param direction Select results for up- or downregulated genes, respectively.
#' @param n_bars Number of bars to plot per subchart.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_enrichr_bars <- function(data,
                              contrast,
                              direction = "up",
                              n_bars = 10,
                              filename = NULL) {
  data <- 
    data %>%
    filter(contrast == {{contrast}}, direction == {{direction}}) %>%
    mutate(Term = shorten_go_term(Term, 30)) %>% 
    group_by(cluster, db) %>% 
    slice_max(Combined.Score, n = n_bars, with_ties = FALSE) %>% 
    ungroup()
  
  n_row <- n_distinct(data$db)
  dbs <- unique(data$db)
  clusters <-
    unique(data$cluster) %>%
    str_sort(numeric = TRUE)
  
  p <-
    list(
      cluster = clusters,
      db = dbs
    ) %>%
    cross_df() %>%
    mutate(
      show_cluster = db == dbs[1],
      show_db = cluster == clusters[1]
    ) %>%
    pmap(
      function(cluster, db, show_db, show_cluster) {
        data <-
          data %>%
          filter(cluster == {{cluster}}, db == {{db}})
        
        if (nrow(data) == 0L)
          return(textGrob("no data", gp = gpar(fontsize = 20)))
        
        data %>%
          mutate(
            Term =
              as_factor(Term) %>%
              fct_reorder(Combined.Score)
          ) %>%
          ggplot(aes(Term, Combined.Score)) +
          geom_col(aes(fill = Combined.Score), show.legend = FALSE) +
          xlab(if (show_db) db else "") +
          ylab(NULL) +
          coord_flip() +
          ggtitle(if (show_cluster) paste("cluster", cluster) else "") +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title.y = element_text(hjust = 0.5, face = "bold"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
          )
      }
    ) %>%
    wrap_plots(nrow = n_row) +
    plot_annotation(str_glue("Enrichment analysis for contrast {contrast}"))
  ggsave_default(filename, width = 1000, height = 420)
  p
}

plot_enrichr_bars(enrichr_results, "II_vs_I",
                  filename = "dge/enrichr_bars_II")  
plot_enrichr_bars(enrichr_results, "III_vs_I",
                  filename = "dge/enrichr_bars_III")  
plot_enrichr_bars(enrichr_results, "IV_vs_I",
                  filename = "dge/enrichr_bars_IV")  





plot_enrichr_halfdots <- function(data,
                                  db,
                                  max_p_adj = 0.05,
                                  min_odds_ratio = 2,
                                  filename = NULL,
                                  plot_params = list()) {
  data_filtered <-
    data %>% 
    filter(db == {{db}}) %>% 
    group_by(Term) %>%
    filter(
      min(Adjusted.P.value) < max_p_adj,
      max(Odds.Ratio) > min_odds_ratio
    ) %>%
    ungroup() %>%
    select(contrast, cluster, Term, Odds.Ratio, Adjusted.P.value) %>% 
    mutate(
      cluster = 
        cluster %>% 
        fct_drop() %>%
        {.},
      Term =
        as_factor(Term) %>%
        fct_reorder(Odds.Ratio, max) %>% 
        fct_drop() %>%
        {.}
    )
  
  scale_factor <- max(-log10(data_filtered$Adjusted.P.value)) / 2.5
  
  p <- 
    ggplot(data_filtered) +
    geom_arc_bar(
      aes(
        x0 = as.integer(cluster),
        y0 = as.integer(Term),
        fill = log10(Odds.Ratio),
        r0 = 0,
        r = sqrt(-log10(Adjusted.P.value)) / scale_factor,
        start = 0,
        end = pi
      ),
      color = NA
    ) +
    geom_arc_bar(
      aes(
        x0 = as.integer(cluster),
        y0 = as.integer(Term),
        fill = log10(Odds.Ratio),
        r0 = 0,
        r = sqrt(-log10(Adjusted.P.value)) / scale_factor,
        start = pi,
        end = 2*pi
      ),
      color = NA
    ) +
    scale_x_continuous(
      breaks = seq(nlevels(data_filtered$cluster)),
      labels = levels(data_filtered$cluster),
      minor_breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(nlevels(data_filtered$Term)),
      labels = levels(data_filtered$Term),
      minor_breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    coord_fixed(
      xlim = c(0.5, nlevels(data_filtered$cluster) + 0.5),
      ylim = c(0.5, nlevels(data_filtered$Term) + 0.5)
    ) +
    facet_wrap(vars(contrast), drop = FALSE) +
    labs(
      y = "",
      fill = TeX("log_{10} (odds ratio)"),
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue( "Enrichr results ({db})"),
      caption = str_glue(
        "adjusted p value < {max_p_adj}, ",
        "odds ratio > {min_odds_ratio}"
      )
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  # {.}
  
  rlang::exec(ggsave_default, filename, !!!plot_params)
  p
}

plot_enrichr_halfdots(enrichr_results,
                      db = "GO_Cellular_Component_2018",
                      filename = "dge/dots_cellcomp_half")
plot_enrichr_dots(enrichr_results,
                  db = "GO_Cellular_Component_2018",
                  filename = "dge/dots_cellcomp")




enrichr_results %>% 
  filter(
    db == "GO_Cellular_Component_2018",
  ) %>%
  group_by(Term) %>%
  filter(min(Adjusted.P.value) < 0.05, max(Odds.Ratio) > 2) %>%
  ungroup() %>% 
  # filter(
  #   Adjusted.P.value < 0.05,
  #   Odds.Ratio > 2,
  # ) %>%
  select(contrast, cluster, Term, Odds.Ratio, Adjusted.P.value)