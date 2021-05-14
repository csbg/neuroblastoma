# Differential gene expression analysis.
#
# @DEPI rna_decontaminated.rds
# @DEPI metadata.rds

library(monocle3)
library(muscat)
library(scater)
library(tidyverse)
library(viridis)
library(patchwork)
library(ggpmisc)
library(ComplexHeatmap)
library(enrichR)
library(latex2exp)
library(scico)
library(fgsea)
library(msigdbr)
library(openxlsx)
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



# DGE analysis ------------------------------------------------------------

nb <- prepSCE(
  nb,
  kid = "cellont_cluster",
  gid = "group",
  sid = "sample",
  drop = TRUE
)

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
      fct_relevel(levels(colData(nb)$cluster_id))
  )



# Filter DGE results ------------------------------------------------------

#' Filter DGE results and add a column "direction".
#'
#' @param data A dataframe as returned by `muscat::resDS()`.
#' @param max_p Maximum p value.
#' @param max_p_adj Maximum adjusted p value. Note: The FDR is calculated by
#'   muscat/edgeR via p.adjust(method = "BH").
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
  # lookup table of gene frequencies in groups
  gene_frq <-
    data %>% 
    select(gene, cluster_id, ends_with("frq")) %>%
    distinct() %>% 
    pivot_longer(
      ends_with("frq"),
      names_to = "frq_col",
      names_pattern = "(.+)\\.",
      values_to = "frq"
    )
  
  res <-
    data %>%
    extract(
      contrast,
      into = c("contrast_left", "contrast_right"),
      regex = "(.+)_vs_(.+)",
      remove = FALSE
    ) %>% 
    left_join(
      gene_frq,
      by = c("gene", "cluster_id", contrast_left = "frq_col")
    ) %>%
    left_join(
      gene_frq,
      by = c("gene", "cluster_id", contrast_right = "frq_col")
    ) %>%
    filter(
      frq.x >= min_freq | frq.y >= min_freq,
      p_val <= max_p,
      p_adj.loc <= max_p_adj,
      abs(logFC) >= min_abs_log_fc
    ) %>% 
    mutate(direction = if_else(logFC > 0, "up", "down")) %>% 
    select(!c(matches("frq"), contrast_left, contrast_right))
  
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

# for visualization and EnrichR analysis, filter with default settings
dge_results_filtered <- filter_dge_results(dge_results)

# for GSEA, only remove low-frequent and ribosomal genes
dge_results_filtered_fgsea <- filter_dge_results(
  dge_results,
  max_p_adj = Inf,
  min_abs_log_fc = 0
)



# Plot results ------------------------------------------------------------

## Volcano plot ----

#' Make one volcano plot per cluster for a given contrast.
#'
#' @param data A DGE dataset.
#' @param contrast A contrast present in the dataset.
#' @param max_p_adj see `filter_dge_results()`.
#' @param min_abs_log_fc see `filter_dge_results()`.
#' @param min_freq see `filter_dge_results()`.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_volcano <- function(data,
                         contrast,
                         max_p_adj = 0.05,
                         min_abs_log_fc = 1,
                         min_freq = 0.1,
                         filename = NULL) {
  contrast_groups <- str_match(contrast, "(.+)_vs_(.+)")[1,]
  group_left <- contrast_groups[2]
  group_right <- contrast_groups[3]
  
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
      diff_expressed =
        case_when(
          is_significant & is_frequent & logFC > 0 ~ "up",  
          is_significant & is_frequent & logFC < 0 ~ "down",
          TRUE ~ "no change"
        ) %>% 
        factor(levels = c("up", "no change", "down"))
    )

  genes_count <- 
    data_filtered %>% 
    count(cluster_id, diff_expressed, .drop = FALSE) %>% 
    pivot_wider(names_from = "diff_expressed", values_from = "n")

  dot_colors <- c(
    `no change` = "gray80",
    up = "red",
    down = "blue"
  )

  p <-
    ggplot(data_filtered, aes(logFC, -log10(p_val))) +
    geom_point(
      aes(color = diff_expressed),
      alpha = .5,
      size = 0.5
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = down),
      npcx = 0.2,
      npcy = 0.9,
      size = 3,
      color = dot_colors["down"],
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = `no change`),
      npcx = 0.5,
      npcy = 0.9,
      size = 3,
      color = dot_colors["no change"],
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = up),
      npcx = 0.8,
      npcy = 0.9,
      size = 3,
      color = dot_colors["up"],
      na.rm = TRUE
    ) +
    scale_color_manual(
      "differentially\nexpressed",
      values = dot_colors,
      guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
      drop = FALSE
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



## Violin plot ----

#' Title
#'
#' @param data A DGE dataset.
#' @param contrast A contrast and ...
#' @param cluster ... cluster for which all samples should be plotted.
#' @param top_n Plot the top differentially expressed genes.
#' @param direction Plot up- or downregulated genes.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_violin <- function(data,
                        contrast,
                        cluster,
                        top_n = 10,
                        direction = c("up", "down"),
                        filename = NULL) {
  contrast_groups <- str_match(contrast, "(.+)_vs_(.+)")[1,]
  group_left <- contrast_groups[2]
  group_right <- contrast_groups[3]
  direction <- match.arg(direction)
  
  top_genes <- 
    data %>%
    filter(
      contrast == {{contrast}},
      cluster_id == {{cluster}},
      direction == {{direction}}
    ) %>% 
    slice_max(n = top_n, logFC, with_ties = FALSE) %>% 
    pull(gene)
  
  p <-
    plotExpression(
      nb[, nb$cluster_id == cluster &
           nb$group_id %in% c(group_left, group_right)],
      top_genes,
      x = "sample_id",
      colour_by = "group_id",
      point_size = 1,
      ncol = NULL
    ) +
    labs(
      caption = str_glue(
        "{direction}regulated genes, cluster {cluster}, contrast {contrast}"
      )
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # calculate facet grid dimensions to ensure that all plots have the same size
  n_facet <- wrap_dims(length(top_genes))
  ggsave_default(
    filename,
    height = 30 + 60 * n_facet[1],
    width = 20 + 70 * n_facet[2]
  )
  p
}

plot_violin(dge_results_filtered, "II_vs_I", cluster = "1",
            filename = "dge/violin_II_vs_I_c1_up")

dge_results_filtered %>% 
  distinct(contrast, cluster = as.character(cluster_id), direction) %>% 
  pwalk(
    function(contrast, cluster, direction) {
      info("Plotting {direction}regulated genes ",
           "in cluster {cluster} in contrast {contrast}")
      plot_violin(
        dge_results_filtered,
        top_n = Inf,
        contrast = contrast,
        cluster = cluster,
        direction = direction,
        filename = str_glue("dge/violin_{contrast}_c{cluster}_{direction}")
      )
    }
  )



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
  contrast_groups <- str_match(contrast, "(.+)_vs_(.+)")[1,]
  group_left <- contrast_groups[2]
  group_right <- contrast_groups[3]
  
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
      col = scico(9, palette = "roma", direction = -1),
      
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
            breaks = c(min(data$logFC), 0, max(data$logFC)),
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
  info("Cluster {cluster}, contrast {contrast}, {direction}regulated genes")
  
  direction <- match.arg(direction)
  
  top_genes <- 
    data %>% 
    filter(
      cluster_id == {{cluster}},
      contrast == {{contrast}},
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

#' Perform enrichment anlalysis for all contrasts and clusters in a dataset.
#'
#' @param data A DGE dataset.
#' @param dbs Character vector containing valid enrichr databases
#'   as returned by `enrichR::listEnrichrDbs()`.
#'
#' @return A dataframe, which combines the dataframes returned by
#'   `enrichR::enrichr()` (empty results are removed) by adding three columns
#'   "contrast", "cluster", and "db".
enrich_all_genes <- function(data, dbs) {
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

enrichr_dbs <- c(
  "GO_Biological_Process_2018",
  "GO_Cellular_Component_2018",
  "GO_Molecular_Function_2018",
  "KEGG_2019_Human",
  "WikiPathways_2019_Human",
  "MSigDB_Hallmark_2020"
)
enrichr_results <- enrich_all_genes(dge_results_filtered, enrichr_dbs)



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
    group_by(contrast, cluster) %>%
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
        fct_reorder(str_c(cluster, contrast), n_distinct)
    )
  
  p <- 
    ggplot(data_vis, aes(contrast, Term, size = -log10(Adjusted.P.value))) +
    geom_point(aes(color = log10(Odds.Ratio))) +
    geom_point(data = data_vis %>% filter(is_significant), shape = 1) +
    scale_color_distiller(
      palette = "Reds",
      direction = 1,
      limits = c(0, log_odds_cap)
    ) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(cluster), drop = FALSE, nrow = 1) +
    labs(
      y = "",
      color = TeX("log_{10} (odds ratio)"),
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue( "Enrichr results ({db})"),
      caption = str_glue(
        "top {top_n} terms per cluster; ",
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
    filename <- str_glue("dge/enrichr_{db}")
  }
  ggsave_default(filename, ...)
  p
}


plot_enrichr_dots(enrichr_results,
                  db = "GO_Biological_Process_2018",
                  height = 400,
                  width = 400)

plot_enrichr_dots(enrichr_results,
                  db = "GO_Cellular_Component_2018")

plot_enrichr_dots(enrichr_results,
                  db = "GO_Molecular_Function_2018",
                  width = 420)

plot_enrichr_dots(enrichr_results,
                  db = "KEGG_2019_Human")

plot_enrichr_dots(enrichr_results,
                  db = "WikiPathways_2019_Human",
                  width = 420,
                  height = 250)

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
#' @return A dataframe, comprising columns "db", "contrast", "cluster", as well
#'   as all columns in the result of `fgseaMultilevel()`.
perform_gsea <- function(data, gene_sets) {
  data %>% 
    distinct(contrast, cluster_id) %>% 
    mutate(db = list(names(gene_sets))) %>% 
    unnest_longer(db) %>% 
    pmap_dfr(
      function(contrast, cluster_id, db) {
        ranked_genes <-
          data %>% 
          filter(contrast == {{contrast}}, cluster_id == {{cluster_id}}) %>%
          select(gene, logFC) %>%
          deframe()
        ranked_genes <- ranked_genes[!is.na(ranked_genes)]
        
        info("GSEA of contrast {contrast}, cluster {cluster_id}, ",
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
            contrast = {{contrast}},
            cluster = {{cluster_id}},
            .before = 1
          )
      }
    )
}


enrichr_genesets <- get_enrichr_genesets(enrichr_dbs)
gsea_results <- perform_gsea(dge_results_filtered_fgsea, enrichr_genesets)
# gsea_results %>% saveRDS("~/Desktop/gsea_results.rds")
# gsea_results <- readRDS("~/Desktop/gsea_results.rds")

gsea_results %>%
  ggplot(aes(NES, -log10(padj))) +
  geom_point(alpha = .1)



#' Generate a dotplot for terms enriched by GSEA.
#' 
#' Dots that meet given limits of adjusted p value and normalized enrichment
#' score (NES) are bordered.
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
                           top_n_positive = 5,
                           top_n_negative = 5,
                           max_p_adj = 0.05,
                           min_abs_NES = 1,
                           filename = "auto",
                           ...) {
  data_selected <-
    data %>% 
    filter(db == {{db}}, padj <= max_p_adj & abs(NES) >= min_abs_NES) %>% 
    group_by(contrast, cluster)
  
  top_terms_pos <- 
    data_selected %>% 
    slice_max(n = top_n_positive, order_by = NES, with_ties = FALSE) %>%
    pull(pathway) %>%
    unique()
  
  top_terms_neg <- 
    data_selected %>% 
    slice_min(n = top_n_negative, order_by = NES, with_ties = FALSE) %>%
    pull(pathway) %>%
    unique()
    
  data_vis <- 
    data_selected %>%
    ungroup() %>% 
    filter(pathway %in% c(top_terms_pos, top_terms_neg)) %>% 
    mutate(
      is_significant =
        padj <= max_p_adj &
        abs(NES) >= min_abs_NES,
      pathway =
        as_factor(pathway) %>%
        fct_reorder(NES, sum, na.rm = TRUE)
      )
  
  color_limit <- max(abs(data_vis$NES))
  
  p <- 
    ggplot(data_vis, aes(contrast, pathway, size = -log10(padj))) +
    geom_point(aes(color = NES)) +
    geom_point(data = data_vis %>% filter(is_significant), shape = 1) +
    scale_color_distiller(
      palette = "PiYG",
      direction = 1,
      limits = c(-color_limit, color_limit)
    ) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(cluster), drop = FALSE, nrow = 1) +
    labs(
      y = "",
      color = "normalized\nenrichment\nscore",
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue("FGSEA results ({db})"),
      caption = str_glue(
        "top {top_n_positive} positively and top {top_n_negative} negatively ",
        "enriched terms per cluster; ",
        "bordered circles: adjusted p value <= {max_p_adj}, ",
        "absolute NES >= {min_abs_NES}"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),
      panel.spacing = unit(0, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  if (filename == "auto") {
    filename <- str_glue("dge/gsea_{db}")
  }
  ggsave_default(filename, ...)
  p
}


plot_gsea_dots(gsea_results,
               db = "MSigDB_Hallmark_2020",
               width = 400)

plot_gsea_dots(gsea_results,
               db = "GO_Biological_Process_2018",
               width = 600,
               height = 600)

plot_gsea_dots(gsea_results,
               db = "KEGG_2019_Human",
               width = 400,
               height = 350)

plot_gsea_dots(gsea_results,
               db = "WikiPathways_2019_Human",
               width = 800,
               height = 400)



#' List genes that influence GSEA the most, for all databases and separately for
#' positively and negatively enriched terms.
#'
#' @param data Results from `perform_gsea()`.
#' @param top_n Only select the terms with the highest NES per cluster.
#' @param max_p_adj Maximum adjusted p value and …
#' @param min_abs_NES minimum absolute NES required for selecting a term.
#'
#' @return A dataframe with columns "db", "direction", "gene", and "count".
get_gsea_genes <- function(data,
                           top_n = 5,
                           max_p_adj = 0.05,
                           min_abs_NES = 1) {
  data <- 
    data %>%
    mutate(direction = if_else(NES > 0, "positive", "negative"))
  
  params <-
    data %>% 
    distinct(db, direction)
  
  pmap_dfr(
    params,
    function(db, direction) {
      data %>% 
        filter(db == {{db}}, direction == {{direction}}) %>% 
        filter(padj <= max_p_adj & abs(NES) >= min_abs_NES) %>%
        group_by(contrast, cluster) %>%
        slice_max(n = top_n, order_by = abs(NES), with_ties = FALSE) %>%
        pull(leadingEdge) %>%
        flatten_chr() %>%
        fct_count() %>%
        arrange(desc(n)) %>% 
        rename(gene = f, count = n) %>% 
        mutate(db = db, direction = direction, .before = 1)
    }
  )
}

gsea_genes <- get_gsea_genes(gsea_results)
gsea_genes



# Export tables -----------------------------------------------------------

# create Excel sheet
header_style <- createStyle(textDecoration = "bold")
wb <- createWorkbook()

# workbook with DE genes
table_data <- 
  dge_results_filtered %>% 
  select(contrast, cluster = cluster_id, gene, logFC,
         p_value = p_val, FDR = p_adj.loc, I.frq:IV.frq)

addWorksheet(wb, "de_genes")
writeData(wb, "de_genes", table_data, headerStyle = header_style)
freezePane(wb, "de_genes", firstRow = TRUE)

# workbook with enrichr results
table_data <- 
  enrichr_results %>% 
  filter(direction == "up", Odds.Ratio >= 1) %>%
  select(
    db, contrast, cluster, term = Term, FDR = Adjusted.P.value,
    odds_ratio = Odds.Ratio, genes = Genes, overlap_size, geneset_size
  ) %>%
  arrange(db, contrast, cluster, desc(odds_ratio)) %>% 
  mutate(odds_ratio = pmin(odds_ratio, 1e99))

addWorksheet(wb, "enrichr")
writeData(wb, "enrichr", table_data, headerStyle = header_style)
freezePane(wb, "enrichr", firstRow = TRUE)


# workbook with fgsea results
table_data <- 
  gsea_results %>% 
  select(
    db, contrast, cluster, term = pathway, FDR = padj, NES,
    genes = leadingEdge
  ) %>% 
  rowwise() %>% 
  mutate(genes = str_c(genes, collapse = ";")) %>%
  ungroup() %>% 
  filter(abs(NES) >= 1) %>%
  arrange(db, contrast, cluster, desc(NES))
addWorksheet(wb, "fgsea")
writeData(wb, "fgsea", table_data, headerStyle = header_style)
freezePane(wb, "fgsea", firstRow = TRUE)


# workbook with leading edge genes from FGSEA
# addWorksheet(wb, "gsea_genes")
# writeData(wb, "gsea_genes", gsea_genes, headerStyle = header_style)
# freezePane(wb, "gsea_genes", firstRow = TRUE)

# save file
saveWorkbook(wb, "plots/dge/dge_data.xlsx", overwrite = TRUE)
