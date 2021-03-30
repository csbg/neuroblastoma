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
library(scico)
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

# for k=50, cluster 8 comprises a few cells in group I
ignored_barcodes <-
  colData(nb) %>%
  as_tibble(rownames = "cell") %>% 
  filter(cluster_50 == "8", group == "I") %>%
  pull(cell)
nb <- nb[, !colnames(nb) %in% ignored_barcodes]



# DGE analysis ------------------------------------------------------------

nb <- prepSCE(
  nb,
  kid = "cluster_50",
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
      fct_inseq()
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
  res <-
    data %>%
    filter(
      I.frq >= min_freq |
      (contrast == "II_vs_I" & II.frq >= min_freq) |
      (contrast == "III_vs_I" & III.frq >= min_freq) |
      (contrast == "IV_vs_I" & IV.frq >= min_freq)
    ) %>%
    filter(
      p_val <= max_p,
      p_adj.loc <= max_p_adj,
      abs(logFC) >= min_abs_log_fc
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
  group_left <- str_extract(contrast, "^[^_]+")
  group_right <- str_extract(contrast, "[^_]$")
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
  message("Cluster ", cluster, ", contrast ", contrast, ", ",
          direction, "regulated genes")
  
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
                               "WikiPathways_2019_Human",
                               "MSigDB_Hallmark_2020"
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
#' @param direction Plot results for `"up"`- or `"down"`-regulated genes.
#' @param top_n Only plot the terms with the highest odds ratio per cluster.
#' @param max_p_adj Maximum adjusted p value, …
#' @param min_odds_ratio … minimum odds ratio, and …
#' @param min_overlap_size … minimum overlap size required for bordered dots.
#' @param log_odds_cap Upper boundary of the color scale.
#' @param filename Name of output file.
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
                              filename = NULL,
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
        Adjusted.P.value < max_p_adj &
        Odds.Ratio > min_odds_ratio &
        overlap_size >= min_overlap_size,
      Odds.Ratio = pmin(Odds.Ratio, 10^log_odds_cap),
      Term =
        as_factor(Term) %>%
        fct_reorder(str_c(cluster, contrast), n_distinct)
    )
  
  # make plot
  p <- 
    ggplot(data_vis, aes(cluster, Term, size = -log10(Adjusted.P.value))) +
    geom_point(aes(color = log10(Odds.Ratio))) +
    geom_point(data = data_vis %>% filter(is_significant), shape = 1) +
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
        "top {top_n} terms per cluster; ",
        "bordered circles: adjusted p value <= {max_p_adj}, ",
        "odds ratio >= {min_odds_ratio}, ",
        "overlap size >= {min_overlap_size}; ",
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
                  filename = "dge/enrichr_GO_Biological_Process_2018",
                  height = 400,
                  width = 400)

plot_enrichr_dots(enrichr_results,
                  db = "GO_Cellular_Component_2018",
                  filename = "dge/enrichr_GO_Cellular_Component_2018")

plot_enrichr_dots(enrichr_results,
                  db = "GO_Molecular_Function_2018",
                  filename = "dge/enrichr_GO_Molecular_Function_2018",
                  width = 420)

plot_enrichr_dots(enrichr_results,
                  db = "KEGG_2019_Human",
                  filename = "dge/enrichr_KEGG_2019_Human")

plot_enrichr_dots(enrichr_results,
                  db = "WikiPathways_2019_Human",
                  filename = "dge/enrichr_WikiPathways_2019_Human",
                  width = 420,
                  height = 250)

plot_enrichr_dots(enrichr_results,
                  db = "MSigDB_Hallmark_2020",
                  filename = "dge/enrichr_MSigDB_Hallmark_2020")



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
