library(monocle3)
library(scuttle)
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

markers <- read_csv("metadata/myeloid_markers.csv")



# Analysis ----------------------------------------------------------------

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

cds_my %>% saveRDS("data_wip/cds_myeloid.rds")
my_metadata %>% saveRDS("data_wip/metadata_myeloid.rds")



# Plots -------------------------------------------------------------------

my_metadata %>% 
  count(subcluster) %>% 
  mutate(n_rel = n / sum(n) * 100)



## UMAPs ----

plot_umap <- function(cluster_col) {
  cluster_labels <-
    my_metadata %>% 
    group_by(cluster = {{cluster_col}}) %>% 
    summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
  
  ggplot(my_metadata, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = {{cluster_col}}), size = 0.1, show.legend = FALSE) +
    geom_text(data = cluster_labels, aes(label = cluster)) +
    coord_fixed() +
    theme_bw() +
    theme(panel.grid = element_blank())
}

plot_umap(subcluster)
ggsave_default("myeloid/umap_subcluster")

plot_umap(collcluster)
ggsave_default("myeloid/umap_collcluster")



## Groupwise abundances ----

plot_abundance_pies <- function(cluster_col) {
  my_metadata %>% 
    group_by(cluster = {{cluster_col}}) %>% 
    count(group) %>% 
    mutate(n = n / sum(n) * 100) %>% 
    filter(as.numeric(cluster) <= 11) %>% 
    ggplot(aes("", n, fill = group)) +
    geom_col() +
    scale_x_discrete(NULL) +
    scale_fill_manual(values = GROUP_COLORS) +
    coord_polar("y") +
    facet_wrap(vars(cluster), nrow = 1) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank()
    )
}

plot_abundance_pies(subcluster)
ggsave_default("myeloid/abundances_groups_subcluster")

plot_abundance_pies(collcluster)
ggsave_default("myeloid/abundances_groups_collcluster", width = 150)


plot_abundance_bars <- function(cluster_col) {
  my_metadata %>% 
    group_by(cluster = {{cluster_col}}) %>% 
    count(sample) %>% 
    mutate(n = n / sum(n) * 100) %>% 
    filter(as.numeric(cluster) <= 11) %>% 
    ggplot(aes("", n, fill = sample)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = PATIENT_COLORS) +
    xlab(NULL) +
    ylab("relative abundance (%)") +
    facet_wrap(vars(cluster), nrow = 2) +
    theme_nb() +
    theme(
      legend.key.width = unit(2, "mm"),
      legend.key.height = unit(2, "mm"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_blank()
    )
}

plot_abundance_bars(subcluster)
ggsave_default("myeloid/abundances_patients_subcluster", width = 150, height = 50)

plot_abundance_bars(collcluster)
ggsave_default("myeloid/abundances_patients_collcluster", width = 70, height = 50)



## Markers ----

plot_cell_type <- function(t) {
  plot_dots(
    logcounts(cds_my),
    markers %>% 
      filter(cell_type == t) %>% 
      pull(gene),
    my_metadata$collcluster
  ) +
    scale_y_discrete(t) +
    theme(legend.position = "none")
}



unique(markers$cell_type) %>% 
  map(plot_cell_type) %>% 
  wrap_plots(ncol = 1)
ggsave_default("myeloid/markers", height = 700)


walk(
  unique(markers$cell_type),
  function(t) {
    plot_cell_type(t)
    ggsave_default(str_glue("myeloid/markers_{t}"), width = 100)
  }
)


## Cell types ----

plot_celltype_heatmap <- function(cluster_col,
                                  clusters = NULL,
                                  label = c("broad", "fine"),
                                  lump_prop = 0) {
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
  
  # draw heatmap  
  ht_opt(
    simple_anno_size = unit(1.5, "mm"),
    COLUMN_ANNO_PADDING = unit(1, "pt"),
    DENDROGRAM_PADDING = unit(1, "pt"),
    HEATMAP_LEGEND_PADDING = unit(1, "mm"),
    ROW_ANNO_PADDING = unit(1, "pt"),
    TITLE_PADDING = unit(1, "mm")
  )
  
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
    
    row_title = "subcluster",
    row_title_side = "right",
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_dend_gp = gpar(lwd = 0.5),
    row_dend_width = unit(3, "mm"),
    
    column_split = col_metadata$ref,
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    cluster_columns = FALSE,
  ) %>%
    draw(
      gap = unit(50, "mm"),
      column_title = "cell type in reference dataset",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    )
}

(p <- plot_celltype_heatmap(subcluster, 1:10))
ggsave_default("myeloid/celltype_heatmap_broad_subcluster",
               plot = p, width = 150, height = 50)

(p <- plot_celltype_heatmap(subcluster, 1:10, label = "fine", lump_prop = .01))
ggsave_default("myeloid/celltype_heatmap_fine_subcluster",
               plot = p, width = 150, height = 80)

(p <- plot_celltype_heatmap(collcluster))
ggsave_default("myeloid/celltype_heatmap_broad_collcluster",
               plot = p, width = 150, height = 35)

(p <- plot_celltype_heatmap(collcluster, label = "fine", lump_prop = .01))
ggsave_default("myeloid/celltype_heatmap_fine_collcluster",
               plot = p, width = 150, height = 60)

