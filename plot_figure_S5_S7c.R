# @DEPI rna_decontaminated.rds
# @DEPI rna_myeloid.rds
# @DEPI metadata.rds
# @DEPI metadata_myeloid.rds
# @DEPI dge_results.rds
# @DEPI dge_results_myeloid.rds

library(monocle3)
library(scuttle)
library(tidyverse)
library(latex2exp)
library(ComplexHeatmap)
library(CellChat)
library(scico)
library(RColorBrewer)
library(patchwork)
library(ggtext)
source("common_functions.R")
source("styling.R")

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(2, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(2, "pt"),
  TITLE_PADDING = unit(1, "mm")
)



# Load data ---------------------------------------------------------------

cds_nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  logNormCounts(assay.type = "soupx_counts")
cds_my <- readRDS("data_generated/rna_myeloid.rds")

my_metadata <- readRDS("data_generated/metadata_myeloid.rds")
nb_metadata <- readRDS("data_generated/metadata.rds")

dge_nb <- readRDS("data_generated/dge_results.rds")
dge_my <- readRDS("data_generated/dge_results_myeloid.rds")

# are there problematic genes with convergence <= -20 ?
stopifnot(
  all(dge_my$results_vs_C$convergence > -20),
  all(dge_my$results_MNA_vs_other$convergence > -20)
)

markers <- read_csv("metadata/myeloid_markers.csv")



# Figures -----------------------------------------------------------------

## S5a ----

plot_celltype_heatmap <- function(cluster_col = collcluster,
                                  clusters = NULL,
                                  label = c("fine", "broad"),
                                  lump_prop = 0.01) {
  # prepare data
  label <- match.arg(label)
  
  # generate matrix of cell type abundances
  make_matrix <- function(ref) {
    cell_type_column <- rlang::sym(str_glue("cell_type_{ref}_{label}"))
    
    my_metadata %>% 
      rename(cluster = {{cluster_col}}) %>% 
      {if (!is.null(clusters)) filter(., cluster %in% {{clusters}}) else .} %>% 
      left_join(
        nb_metadata %>% select(cell, starts_with("cell_type")),
        by = "cell"
      ) %>% 
      mutate(
        cell_type =
          as_factor(!!cell_type_column) %>%
          fct_infreq() %>% 
          fct_explicit_na("Unknown") %>% 
          fct_relabel(~str_c(ref, .x, sep = "_")) %>% 
          fct_lump_prop(lump_prop, other_level = str_glue("{ref}_other"))
      ) %>% 
      count(cluster, cell_type) %>%
      group_by(cluster) %>%
      mutate(n_rel = n / sum(n)) %>%
      select(!n) %>%
      ungroup() %>%
      arrange(cell_type) %>% 
      pivot_wider(names_from = "cell_type", values_from = "n_rel") %>%
      arrange(cluster) %>% 
      column_to_rownames("cluster") %>%
      as.matrix() %>%
      replace_na(0)
  }
  
  mat <- 
    map(
      c("blueprint", "hpca", "dice", "dmap", "monaco"),
      make_matrix
    ) %>% 
    reduce(cbind)
  
  # set up column metadata
  col_metadata <-
    tibble(colname = colnames(mat)) %>% 
    left_join(
      read_csv("metadata/celldex_celltypes.csv", comment = "#"),
      by = "colname"
    ) %>% 
    separate(
      colname,
      into = c("ref", "cell_type"),
      extra = "merge",
      remove = FALSE
    ) %>% 
    mutate(
      ref =
        as_factor(ref) %>% 
        fct_recode(
          "Human Primary Cell Atlas" = "hpca",
          "Blueprint/ENCODE" = "blueprint",
          "DICE" = "dice",
          "Novershtern" = "dmap",
          "Monaco" = "monaco"
        ),
      abbr = factor(abbr, levels = names(CELL_TYPE_COLORS))
    ) %>% 
    group_by(ref) %>% 
    arrange(abbr, .by_group = TRUE) %>% 
    ungroup()
  
  mat <- mat[, col_metadata$colname]
  colnames(mat) <- col_metadata$cell_type
  
  # export source data
  mat %>% 
    magrittr::set_colnames(
      col_metadata %>% 
        unite(ref, cell_type, col = "ref_type") %>% 
        pull(ref_type)
    ) %>% 
    as_tibble(rownames = "cluster") %>% 
    save_table("source_data/figure_S5a", "Figure S5a")
  
  # make plot
  set.seed(2)
  Heatmap(
    mat %>% magrittr::set_rownames(rename_myeloid(rownames(.))),
    col = colorRampPalette(brewer.pal(9, "YlOrBr"))(100),
    
    heatmap_legend_param = list(
      at = c(0, 1),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    name = "relative\nabundance",
    
    row_title = "assigned cell type",
    row_title_side = "right",
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_dend_gp = gpar(lwd = 0.5),
    row_dend_width = unit(3, "mm"),
    
    column_split = col_metadata$ref,
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    cluster_columns = FALSE,
    column_gap = unit(0.5, "mm"),  # 50 x 4 cells, 4 gaps
    
    width = unit(93, "mm"),
    height = unit(7.28, "mm"),
    
    right_annotation = rowAnnotation(
      cell_type = rownames(mat),
      col = list(cell_type = MYELOID_COLORS),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
  ) %>%
    draw(
      gap = unit(50, "mm"),
      column_title = "cell type in reference dataset",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    )
}

(p <- plot_celltype_heatmap(collcluster))
ggsave_publication("S5a_myeloid_types", plot = p, width = 14, height = 6)



## S5b ----

plot_umap <- function() {
  # prepare data
  adjust_positions <- tribble(
    ~label, ~dx, ~dy,
    "classical mono", 0, 0,
    "nonclassical mono", 4.5, -1,
    "mDCs", 0, 2,
    "other", 0, 0
  )
  
  cluster_labels <-
    my_metadata %>% 
    group_by(label = collcluster) %>% 
    summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)) %>% 
    left_join(adjust_positions, by = "label") %>%
    mutate(
      label = as_factor(label) %>% rename_myeloid(),
      UMAP1 = UMAP1 + replace_na(dx, 0),
      UMAP2 = UMAP2 + replace_na(dy, 0)
    )
  
  # export source data
  my_metadata %>% 
    select(cell, umap_1 = UMAP1, umap_2 = UMAP2, cell_type = collcluster) %>% 
    save_table("source_data/figure_S5b", "Figure S5b")
  
  # make plot
  ggplot(my_metadata, aes(UMAP1, UMAP2)) +
    geom_point(
      aes(color = collcluster),
      size = .001,
      shape = 16,
      show.legend = FALSE
    ) +
    geom_text(
      data = cluster_labels,
      aes(label = label),
      size = BASE_TEXT_SIZE_MM
    ) +
    scale_color_manual(
      name = "cell type",
      values = MYELOID_COLORS,
      guide = guide_legend(override.aes = list(size = 1))
    ) +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(),
      plot.background = element_blank()
    )
}

plot_umap()
ggsave_publication("S5b_umap_myeloid", type = "png",
                   height = 5, width = 6, bg = "transparent")



## S5c ----

plot_gsea <- function(db = "MSigDB_Hallmark_2020",
                      comparisons = NULL,
                      cell_types = NULL,
                      top_n_positive = 5L,
                      top_n_negative = 5L,
                      max_p_adj = 0.05,
                      min_abs_NES = 1) {
  # prepare data
  data <- 
    dge_my$gsea %>%
    filter(comparison != "Mas") %>% 
    mutate(comparison = factor(comparison)) %>% 
    filter(db == {{db}})
  
  if (!is.null(comparisons))
    data <- 
      data %>% 
      filter(comparison %in% {{comparisons}})
  
  if (!is.null(cell_types))
    data <- 
      data %>% 
      filter(cell_type %in% {{cell_types}})
  
  data_top_terms <-
    data %>% 
    filter(
      padj <= max_p_adj,
      abs(NES) >= min_abs_NES
    ) %>% 
    group_by(comparison, cell_type)
  
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
    filter(pathway %in% c(top_terms_pos, top_terms_neg)) %>% 
    mutate(
      pathway =
        as_factor(pathway) %>%
        fct_reorder(NES * -log10(padj), sum, na.rm = TRUE),
      cell_type =
        cell_type %>%
        factor(levels = names(MYELOID_COLORS)) %>%  
        rename_myeloid_newline(),
      comparison = 
        comparison %>% 
        str_sub(1, 1) %>% 
        fct_relevel("M", "A", "S")
    )
  
  if (nlevels(data_vis$pathway) > 5) {
    horizontal_grid <-
      geom_hline(
        yintercept = seq(6, nlevels(data_vis$pathway), 5) - 0.5,
        size = BASE_LINE_SIZE,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  
  color_limit <- max(abs(dge_my$gsea$NES), na.rm = TRUE)
  
  # export source data
  data_vis %>% 
    select(cell_type, comparison, pathway, padj, NES) %>% 
    save_table("source_data/figure_S5c", "Figure S5c")
  
  # make plot  
  ggplot(data_vis, aes(comparison, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    xlab("comparison (NB subtype vs control)") +
    ylab(str_glue("{str_replace_all(db, '_', ' ')} gene set")) +
    scale_color_gsea(
      "normalized\nenrichment\nscore",
      limits = c(-color_limit, color_limit),
      breaks = c(-color_limit, 0, color_limit),
      labels = function(x) round(x, 2),
      guide = guide_colorbar(
        barheight = unit(15, "mm"),
        barwidth = unit(2, "mm"),
        ticks = FALSE
      )
    ) +
    scale_size_area(
      name = TeX("-log_{10} p_{adj}"),
      max_size = 3,
      limits = c(0, 40),
      breaks = c(0, 20, 40),
      labels = c("0", "20", "40 or higher"),
      oob = scales::oob_squish
    )  +
    coord_fixed() +
    facet_wrap(vars(cell_type), nrow = 1) +
    theme_nb(grid = FALSE) +
    theme(
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      legend.margin = margin(),
      panel.spacing = unit(-.5, "pt"),
      strip.text = element_text(angle = 90, hjust = 0),
    )
}

plot_gsea()
ggsave_publication("S5c_gsea", width = 8, height = 7.8)



## S5d (UMAP) ----

plot_genes <- function(...,
                       dataset = c("myeloid", "all"),
                       colorbar_quantiles = c(0, .99),
                       source_data_filename = NULL,
                       source_data_sheet = NULL) {
  # prepare data
  dataset <- match.arg(dataset)
  if (dataset == "myeloid") {
    metadata <- my_metadata
    cds <- cds_my
  } else {
    metadata <-
      nb_metadata %>% 
      transmute(
        cell,
        sample,
        UMAP1 = umap_1_monocle,
        UMAP2 = umap_2_monocle,
        group = rename_groups(group)
      )
    cds <- cds_nb
  }
  
  genes <- c(...)
  unknown_genes <- setdiff(genes, rownames(cds))
  if (length(unknown_genes) > 0) {
    info("Ignoring genes: {unknown_genes}")
    genes <- setdiff(genes, unknown_genes)
  }
  
  plot_data <- 
    metadata %>% 
    left_join(
      logcounts(cds)[genes, , drop = FALSE] %>% 
        as.matrix() %>% 
        t() %>% 
        scale() %>% 
        as_tibble(rownames = "cell") %>% 
        pivot_longer(!cell, names_to = "gene", values_to = "expression"),
      by = "cell"
    ) %>% 
    mutate(gene = factor(gene, levels = genes)) %>% 
    arrange(expression)
  
  colorbar_limits <- quantile(plot_data$expression, colorbar_quantiles)
  
  # export source data
  plot_data %>% 
    select(cell, group, UMAP1, UMAP2, gene, expression) %>% 
    mutate(group = if_else(group == "C", "control group", "NB patients")) %>% 
    pivot_wider(names_from = gene, values_from = expression) %>% 
    save_table(source_data_filename, source_data_sheet)
  
  # make plot
  ggplot(plot_data, aes(UMAP1, UMAP2)) +
    geom_point(
      aes(color = expression),
      size = 0.001,
      shape = 16,
    ) +
    scico::scale_color_scico(
      "scaled\nexpression",
      palette = "acton",
      direction = -1,
      limits = colorbar_limits,
      breaks = colorbar_limits,
      labels = c("low", "high"),
      oob = scales::oob_squish_any,
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(15, "mm"),
        direction = "horizontal",
        ticks = FALSE,
        title.hjust = 1,
        title.vjust = 1
      )
    ) +
    scale_x_continuous("gene", position = "top") +
    coord_fixed() +
    facet_grid(
      vars(if_else(group == "C", "control group", "NB patients")),
      vars(gene)
    ) +
    theme_nb(grid = FALSE) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.spacing.x = unit(2, "mm"),
      panel.border = element_blank(),
      plot.background = element_blank()
    )
}

plot_genes(
  "CD163", "CXCL2", "TIMP1", "EREG",
  dataset = "myeloid",
  source_data_filename = "source_data/figure_S5d_UMAP",
  source_data_sheet = "Figure S5d (UMAP)"
)
ggsave_publication("S5d_genes", type = "png",
                   width = 10, height = 8, bg = "transparent")


## S5d (Dotplot) ----

plot_significance <- function(...) {
  # prepare data
  selected_genes <- factor(c(...))
  selected_comparisons <- c("II_vs_I", "III_vs_I", "IV_vs_I")
  
  plot_data <- 
    dge_nb$results_wide %>% 
    filter(
      cell_type == "M",
      gene %in% selected_genes,
      comparison %in% selected_comparisons
    ) %>% 
    mutate(
      logFC = logFC / log(2),
      comparison =
        comparison %>%
        factor(levels = selected_comparisons) %>%
        fct_rev() %>% 
        rename_contrast_long(with_format = TRUE),
      gene = factor(gene, levels = selected_genes)
    )
  
  color_limit <- max(abs(plot_data$logFC))
  
  # export source data
  plot_data %>% 
    select(cell_type:logFC, p_adj) %>% 
    mutate(comparison = str_replace_all(comparison, "<.+?>", "")) %>% 
    save_table("source_data/figure_S5d_dotplot", "Figure S5d (Dotplot)")
  
  # make plot
  ggplot(plot_data, aes(gene, comparison, size = -log10(p_adj))) +
    geom_point(aes(color = logFC)) +
    geom_point(data = plot_data %>% filter(frq >= 0.05), shape = 1) +
    scale_color_distiller(
      palette = "RdBu",
      direction = -1,
      limits = c(-color_limit, color_limit),
      labels = function(x) round(x, 2),
      breaks = c(-color_limit, 0, color_limit),
      guide = guide_colorbar(
        title = "log fold change",
        barheight = unit(10, "mm"),
        barwidth = unit(2, "mm"),
        ticks = FALSE
      )
    ) +
    scale_size_area(
      name = TeX("-log_{10} p_{adj}"),
      max_size = 3,
      limits = c(0, 4),
      breaks = c(0, 2, 4),
      labels = c("0", "2", "4 or higher"),
      oob = scales::oob_squish
    ) +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_markdown(),
      legend.box = "horizontal",
      legend.box.just = "bottom",
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm")
    )
}

plot_significance("CD163", "CXCL2", "TIMP1", "EREG")
ggsave_publication("S5d_significance", width = 8, height = 2.5)



## S7c ----

wrap_plots(
  plot_genes(
    "IL10",
    dataset = "all",
    colorbar_quantiles = c(0, 0.995),
    source_data_filename = "source_data/figure_S7c_left",
    source_data_sheet = "Figure S7c"
  ) +
    ggtitle("all cells") +
    theme(strip.text.y = element_blank()),
  plot_genes(
    "IL10",
    dataset = "myeloid",
    source_data_filename = "source_data/figure_S7c_right",
    source_data_sheet = "Figure S7c"
  ) +
    ggtitle("myeloid cells"),
  nrow = 1,
  guides = "collect"
) &
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(
      size = BASE_TEXT_SIZE_PT,
      hjust = 0.5,
      margin = margin(b = -1, unit = "mm")
    ),
    strip.text.x = element_blank(),
  )
ggsave_publication("S7c_IL10", type = "png",
                   width = 11, height = 11, bg = "transparent")

# combine the two exported source data tables
bind_rows(
  all = read.xlsx("tables/source_data/figure_S7c_left.xlsx"),
  myeloid = read.xlsx("tables/source_data/figure_S7c_right.xlsx"),
  .id = "cells"
) %>% 
  save_table("source_data/figure_S7c", "Figure S7c")
file_delete("tables/source_data/figure_S7c_left.xlsx")
file_delete("tables/source_data/figure_S7c_right.xlsx")



# Data --------------------------------------------------------------------

gene_pathways <- 
  CellChatDB.human$interaction %>% 
  as_tibble() %>% 
  select(interaction_name, pathway_name) %>% 
  mutate(genes = str_split(interaction_name, "_")) %>% 
  unnest_longer(genes) %>% 
  distinct(pathway_name, genes) %>% 
  group_by(Gene = genes) %>% 
  summarise(
    Pathways =
      unique(pathway_name) %>%
      str_sort() %>%
      str_c(collapse = ", ")
  )

rename_columns <- function(s) {
  comparison <- str_sub(s, 7)
  case_when(
    str_starts(s, "log") ~ paste0("Log fold change (", comparison, ")"),
    str_starts(s, "p_adj") ~ paste0("Adjusted p-value (", comparison, ")"),
    TRUE ~ s
  )
}



## S5 ----

dge_my$gsea %>% 
  arrange(db, comparison, cell_type, desc(NES)) %>%
  mutate(
    comparison = rename_contrast_long_alt(comparison),
    cell_type = rename_myeloid(cell_type),
    leadingEdge = map_chr(leadingEdge, str_c, collapse = ", ")
  ) %>%
  filter(
    db %in% c("MSigDB_Hallmark_2020", "TRRUST_Transcription_Factors_2019")
  ) %>% 
  select(
    "Database" = db,
    "Comparison" = comparison,
    "Cell type" = cell_type,
    "Pathway" = pathway,
    "Normalized Enrichment Score" = NES,
    "Adjusted p-value" = padj,
    "Leading edge genes" = leadingEdge
  ) %>%
  split(.$`Cell type`) %>%
  map(select, !`Cell type`) %>%
  save_table("data_S5_gsea_myeloid")



## S6 ----

dge_my$results_wide_filtered %>% 
  arrange(comparison, cell_type, desc(logFC)) %>%
  mutate(
    logFC = logFC / log(2),  # nebula returns natural log fold changes
    comparison = rename_contrast_long_alt(comparison),
  ) %>%
  select(!c(p, frq, frq_ref, direction)) %>%
  pivot_wider(names_from = comparison, values_from = c(logFC, p_adj)) %>%
  relocate(1, 2, 3, 7, 4, 8, 5, 9, 6, 10) %>%
  rename_with(rename_columns) %>%
  left_join(gene_pathways, by = c(gene = "Gene")) %>%
  split(.$cell_type) %>%
  map(select, !cell_type) %>%
  save_table("data_S6_dge_myeloid")