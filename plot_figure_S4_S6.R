# @DEPI rna_myeloid.rds
# @DEPI metadata.rds
# @DEPI metadata_myeloid.rds
# @DEPI dge_results_myeloid.rds

library(monocle3)
library(scuttle)
library(tidyverse)
library(latex2exp)
library(ComplexHeatmap)
library(CellChat)
library(scico)
library(RColorBrewer)
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
dge <- readRDS("data_generated/dge_results_myeloid.rds")

# are there problematic genes with convergence <= -20 ?
stopifnot(
  all(dge$results_vs_C$convergence > -20),
  all(dge$results_MNA_vs_other$convergence > -20)
)

markers <- read_csv("metadata/myeloid_markers.csv")



# Figures -----------------------------------------------------------------

## S4a ----

plot_celltype_heatmap <- function(cluster_col = collcluster,
                                  clusters = NULL,
                                  label = c("fine", "broad"),
                                  lump_prop = 0.01) {
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
  rownames(mat) <- rename_myeloid(rownames(mat))
  
  # draw heatmap  
  set.seed(2)
  Heatmap(
    mat,
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
    column_gap = unit(1, "mm"),  # 50 x 4 cells, 4 gaps
    width = unit(95, "mm"),
    height = unit(7.28, "mm")
  ) %>%
    draw(
      gap = unit(50, "mm"),
      column_title = "cell type in reference dataset",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    )
}

(p <- plot_celltype_heatmap(collcluster))
ggsave_publication("S4a_myeloid_types", plot = p, width = 14, height = 6)



## S4b ----

plot_umap <- function() {
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
ggsave_publication("S4b_umap_myeloid", type = "png",
                   height = 5, width = 6, bg = "transparent")



## S4c ----

plot_gsea <- function(db = "MSigDB_Hallmark_2020",
                      comparisons = NULL,
                      cell_types = NULL,
                      top_n_positive = 5L,
                      top_n_negative = 5L,
                      max_p_adj = 0.05,
                      min_abs_NES = 1) {
  data <- 
    dge$gsea %>%
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
        yintercept = seq(5, nlevels(data_vis$pathway), 5),
        size = BASE_LINE_SIZE,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  
  color_limit <- max(abs(dge$gsea$NES), na.rm = TRUE)
  
  ggplot(data_vis, aes(comparison, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    xlab("comparison (genetic subtype vs control)") +
    ylab(NULL) +
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
ggsave_publication("S4c_gsea", width = 8, height = 8)



## S4d ----

plot_genes <- function(...,
                       dataset = c("myeloid", "all"),
                       colorbar_quantiles = c(0, .99)) {
  dataset <- match.arg(dataset)
  if (dataset == "myeloid") {
    metadata <- my_metadata
    cds <- cds_my
  } else {
    metadata <-
      nb_metadata %>% 
      transmute(
        cell,
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
  
  ggplot(plot_data, aes(UMAP1, UMAP2)) +
    geom_point(
      aes(color = expression),
      size = 0.001,
      shape = 16,
    ) +
    scale_color_viridis_c(
      "scaled\nexpression",
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
    coord_fixed() +
    facet_grid(vars(if_else(group == "C", "C", "M/A/S")), vars(gene)) +
    theme_nb(grid = FALSE) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.spacing.x = unit(2, "mm"),
      panel.border = element_blank(),
      plot.background = element_blank()
    )
}

plot_genes("CD163", "CXCL2", "VEGFA", "EREG")
ggsave_publication("S4d_genes", type = "png",
                   width = 10, height = 8, bg = "transparent")


## S6c ----

plot_genes("IL10", dataset = "all", colorbar_quantiles = c(0, 0.995)) +
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  )
ggsave_publication("S6c_IL10_all", type = "png",
                   width = 5, height = 11, bg = "transparent")



## S6d ----

plot_genes("IL10") +
  theme(strip.text.x = element_blank())
ggsave_publication("S6d_IL10_myeloid", type = "png",
                   width = 5, height = 11, bg = "transparent")



# Tables ------------------------------------------------------------------

## S4 ----

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

dge$results_wide_filtered %>% 
  arrange(comparison, cell_type, desc(logFC)) %>%
  mutate(
    logFC = logFC / log(2),  # nebula returns natural log fold changes
    comparison = recode(comparison, Ac = "A vs C", Mc = "M vs C",
                        Sc = "S vs C", Mas = "M vs A+S"),
  ) %>%
  select(!c(p, frq, frq_ref, direction)) %>%
  pivot_wider(names_from = comparison, values_from = c(logFC, p_adj)) %>%
  relocate(1, 2, 3, 7, 4, 8, 5, 9, 6, 10) %>%
  rename_with(rename_columns) %>%
  left_join(gene_pathways, by = c(gene = "Gene")) %>%
  split(.$cell_type) %>%
  map(select, !cell_type) %>%
  save_table("Sxxx_dge_myeloid")
