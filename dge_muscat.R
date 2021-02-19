library(muscat)
library(scater)
library(tidyverse)
library(viridis)
library(patchwork)
library(ggpmisc)
library(ComplexHeatmap)
library(enrichR)
source("common_functions.R")



# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_integrated_monocle.rds") %>% 
  logNormCounts()

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
#' @param max_p_adj Maximum adjusted p value.
#' @param min_abs_log_fc Minimum absolute log fold change.
#' @param min_freq Minimum gene expression frequency (fracions of cells that
#'   expresses a given gene.)
#'
#' @return A data frame.
filter_dge_results <- function(data,
                               max_p_adj = 0.05,
                               min_abs_log_fc = 1,
                               min_freq = 0.1) {
  data %>%
    filter(
      I.frq > min_freq |
      (contrast == "II_vs_I" & II.frq > min_freq) |
      (contrast == "III_vs_I" & III.frq > min_freq) |
      (contrast == "IV_vs_I" & IV.frq > min_freq)
    ) %>%
    filter(
      p_adj.loc < max_p_adj,
      abs(logFC) > min_abs_log_fc
    )
}

dge_results_filtered <- 
  filter_dge_results(dge_results)



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
                         point_colors = c(
                           "frequent, DE" = "red",
                           "DE" = "black",
                           "not DE" =  "gray80"
                         ),
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
      type = case_when(
        is_significant & is_frequent ~ "frequent, DE",
        is_significant ~ "DE",
        TRUE ~ "not DE"
      )
    ) 
  
  genes_count <- 
    data_filtered %>% 
    count(cluster_id, type) %>% 
    pivot_wider(names_from = "type", values_from = "n")
  
  p <-
    ggplot(data_filtered, aes(logFC, -log10(p_adj.loc))) +
    geom_point(
      aes(color = type),
      alpha = .5,
      size = 1
    ) +
    geom_hline(yintercept = -log10(max_p_adj), size = .25) +
    geom_vline(xintercept = c(-min_abs_log_fc, min_abs_log_fc), size = .25) +
    geom_text_npc(
      data = genes_count,
      aes(label = `frequent, DE`),
      npcx = 0.05,
      npcy = 0.9,
      size = 3,
      color = point_colors["frequent, DE"],
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = DE),
      npcx = 0.05,
      npcy = 0.8,
      size = 3,
      color = point_colors["DE"],
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = `not DE`),
      npcx = 0.05,
      npcy = 0.7,
      size = 3,
      color = point_colors["not DE"],
      na.rm = TRUE
    ) +
    scale_color_manual(values = point_colors) +
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

plot_violin(dge_results_filtered, "II_vs_I",
            cluster = "1", filename = "dge/violin_II_1")
plot_violin(dge_results_filtered, "II_vs_I",
            cluster = "13", filename = "dge/violin_II_13")
plot_violin(dge_results_filtered, "IV_vs_I",
            cluster = "13", filename = "dge/violin_IV_13")



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
  
  ggsave_default(filename, plot = p)
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
#'
#' @return A dataframe, which combines the dataframes returned by
#'   `enrichR::enrichr()` (empty results are removed) by adding a column "db".
enrich_genes <- function(data, contrast, cluster, dbs) {
  message("Cluster ", cluster, " in contrast ", contrast)
  top_genes <- 
    data %>% 
    filter(
      cluster_id == {{cluster}},
      contrast == {{contrast}},
      logFC > 0
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
    distinct(contrast, cluster = cluster_id)
  
  enrichr_data <-
    queries %>% 
    pmap(enrich_genes, data = data, dbs = dbs)
  
  queries %>% 
    mutate(data = enrichr_data) %>% 
    rowwise() %>% 
    filter(!is.null(data)) %>%
    filter(nrow(data) > 0) %>% 
    unnest(data)
}

enrichr_results <- enrich_all_genes(dge_results_filtered)




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
#' @param n_bars Number of bars to plot per subchart.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_enrichr_bars <- function(data, contrast, n_bars = 10, filename = NULL) {
  data <- 
    data %>%
    filter(contrast == {{contrast}}) %>%
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

enrichr_results %>% 
  filter(contrast == "II_vs_I") %>% 
  mutate(Term = str_sub(Term, end = 40)) %>% 
  group_by(cluster, db) %>% 
  slice_max(Combined.Score, n = 10, with_ties = FALSE) %>% 
  filter(cluster == "1", db == "GO_Cellular_Component_2018") %>% 
  mutate(Term = as_factor(Term) %>% fct_reorder(Combined.Score)) %>% 
  View()

x[[8]] %>% View()







# from here, work in progress ...
plot_enrichr_heatmap <- function(data,
                                 min_combined_score = 3000,
                                 db = "GO_Biological_Process_2018") {
  data <-
    data %>% 
    filter(db == {{db}}) %>% 
    select(contrast, Term, cluster, Combined.Score) %>%
    arrange(contrast, cluster) %>% 
    pivot_wider(
      names_from = c(contrast, cluster),
      names_sep = ".",
      values_from = Combined.Score
    ) %>% 
    column_to_rownames("Term")
  
  col_data <-
    tibble(cols = colnames(data)) %>%
    separate(cols, into = c("contrast", "cluster"), sep = "\\.")

  mat <- 
    data %>%
    rowwise() %>% 
    mutate(row_sum = mean(c_across(everything()), na.rm = TRUE)) %>%
    arrange(desc(row_sum)) %>%
    select(!row_sum) %>% 
    filter(if_any(everything(), ~. > min_combined_score)) %>% 
    as.matrix() %>%
    magrittr::set_colnames(col_data$cluster)
  
  Heatmap(
    mat,
    col = circlize::colorRamp2(
      seq(0, quantile(mat, 0.95, na.rm = TRUE), length.out = 10),
      plasma(10, direction = -1)
    ),
    # col = plasma(10, direction = -1),
    na_col = "gray95",
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = HeatmapAnnotation(
      contrast = col_data$contrast,
      col = list(
        contrast = c(
          "II_vs_I" = "red",
          "III_vs_I" = "green",
          "IV_vs_I" = "blue"
        )
      ),
      show_legend = FALSE
    ),
    column_split = col_data$contrast
  )
}

plot_enrichr_heatmap(enrichr_results)
ggsave_default("dge/heatmap")

  


enrichr_results %>% 
  ggplot(aes(Combined.Score)) +
  geom_histogram(binwidth = 100) +
  scale_y_log10()






# Unused ------------------------------------------------------------------

# cross table
table(colData(nb)$cluster_id, colData(nb)$group_id) %>% sum()

# UMAP with expression colored
top_12 <- 
  dge_results_filtered %>%
  filter(contrast == "II_vs_I", cluster_id == "9") %>% 
  slice_max(n = 4, logFC, with_ties = FALSE) %>%
  pull(gene) %>%
  {.}

top_12 %>% 
  map(
    ~plotReducedDim(
      nb,
      "UMAP",
      colour_by = .x,
      point_size = 0.1,
      point_alpha = 1,
      add_legend = FALSE
    ) +
      coord_fixed()
  ) %>% 
  wrap_plots()
