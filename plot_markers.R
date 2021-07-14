# Plot canonical cell type and neuroblastoma markers.
# Exports several plots to plots/markers.
#
# @DEPI metadata.rds
# @DEPI rna_decontaminated.rds

library(monocle3)
library(scuttle)
library(tidyverse)
library(scico)
library(ggnewscale)
library(patchwork)
source("common_functions.R")
source("styling.R")



# Functions ---------------------------------------------------------------

#' Add subclusters to the NB cluster.
#'
#' @param metadata Cell metadata.
#' @param cluster_col Column with supercluster IDs.
#' @param nb_cluster_name Name of the NB cluster in this column.
#' @param subcluster_file File with subcluster information.
#'
#' @return The dataframe passed as `metadata`, with the NB cluster in the column
#'   given by `cluster_col` subdivided into subclusters.
subdivide_tumor_cluster <- function(metadata,
                                    cluster_col,
                                    nb_cluster_name,
                                    subcluster_file) {
  subclusters <- 
    read_csv(subcluster_file, col_types = "ci") %>% 
    mutate(
      tumor_subcluster =
        letters[tumor_subcluster] %>% 
        as_factor() %>%
        fct_infreq() %>% 
        fct_relabel(~str_glue("{str_sub(nb_cluster_name, end = -2L)}{.})"))
    )
  
  old_levels <-
    metadata %>% 
    pull({{cluster_col}}) %>% 
    levels()
  
  tumor_cluster_pos <- str_which(old_levels, fixed(nb_cluster_name))
  
  new_levels <- c(
    old_levels[1:(tumor_cluster_pos - 1)],
    levels(subclusters$tumor_subcluster),
    old_levels[(tumor_cluster_pos + 1):length(old_levels)]
  )
  
  metadata %>% 
    left_join(subclusters, by = "cell") %>%
    mutate(
      {{cluster_col}} :=
        case_when(
          {{cluster_col}} == nb_cluster_name ~ as.character(tumor_subcluster),
          TRUE ~ as.character({{cluster_col}})
        ) %>%
        as_factor() %>%
        fct_relevel(new_levels)
    ) %>%
    select(!tumor_subcluster)
}


# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  logNormCounts(assay.type = "soupx_counts")

nb_metadata <- readRDS("data_generated/metadata.rds")

colData(nb) <-
  nb_metadata %>%
  mutate(Size_Factor = colData(nb)$Size_Factor) %>% 
  # subdivide_tumor_cluster(
  #   cluster_col = cellont_cluster,
  #   nb_cluster_name = "NB (8)",
  #   subcluster_file = "metadata/nb_subclusters.csv"
  # ) %>% 
  column_to_rownames("cell") %>% 
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)

markers <- read_csv("metadata/cell_markers.csv", comment = "#")



# Canonical cell type markers ---------------------------------------------

#' Plot manually selected canonical cell type markers for each cluster.
#' Oder cluster by manually assigned cell types.
#'
#' @param cluster_col Column of the cell metadata that contains cluster ids.
#' @param filename Name of the output file.
#'
#' @return A ggplot2 object
plot_canonical_markers <- function(cluster_col, filename = NULL) {
  y_annotation_data <-
    markers %>%
    arrange(desc(row_number())) %>%
    mutate(
      super_type = case_when(
        str_starts(cell_type, "l") ~ "leukocyte",
        str_starts(cell_type, "T") ~ "T cell",
        str_starts(cell_type, "m") ~ "myeloid",
        TRUE ~ cell_type
      ),
      r = row_number()
    ) %>%
    group_by(label = super_type) %>%
    summarise(
      yintercept = first(r) - 0.5,
      label_y = mean(r)
    )

  x_annotation_data <-
    tibble(level = levels(colData(nb)$cellont_cluster)) %>% 
    extract(level, into = "label", regex = "(\\w+)") %>%
    mutate(label = as_factor(label), r = row_number()) %>%
    group_by(label) %>%
    summarize(
      xmin = first(r) - 0.5,
      xmax = last(r) + 0.5,
      label_x = mean(r)
    ) %>%
    mutate(
      fill = case_when(
        row_number() %% 2 == 0 ~ "white",
        TRUE                   ~ "gray90"
      )
    )
  
  p <-
    plot_dots(
      logcounts(nb),
      rev(markers$gene),
      colData(nb) %>% as_tibble() %>% pull({{cluster_col}}),
      panel_annotation = x_annotation_data
    ) +
    scale_x_discrete(labels = function(x) str_replace(x, " ", "\n")) +
    geom_hline(
      yintercept = y_annotation_data$yintercept,
      linetype = "dashed",
      size = 0.25
    ) +
    geom_text(
      data = y_annotation_data,
      aes(
        x = nlevels(colData(nb)$cellont_cluster) + 2,
        y = label_y,
        label = label
      ),
      size = 4,
      hjust = 0
    ) +
    theme(legend.position = "left") +
    NULL

  ggsave_default(filename, height = 297, width = 350)
  p
}

plot_canonical_markers(cellont_cluster, filename = "markers/overview")



# NB markers --------------------------------------------------------------

#' Normalize expression values.
#' 
#' Values are first limited to the given lower and upper quantile (which are
#' calculated after exclusing zero values) and then normalized to the range [0,
#' 1].
#'
#' @param x Numeric vector.
#' @param lower_quantile Quantile for calculating the lower limit.
#' @param upper_quantile Quantile for calculating the upper limit.
#'
#' @return The quantile-normalized numeric vector.
normalize_expression <- function(x, lower_quantile, upper_quantile) {
  q <- quantile(x[x > 0], c(lower_quantile, upper_quantile))
  x_limited <- pmax(pmin(x, q[2]), q[1])
  (x_limited - min(x_limited)) / (max(x_limited) - min(x_limited))
}


#' Color single cells on a dimensional reduction plot by normalized expression.
#'
#' @param counts Count matrix whose rownames correspond to features.
#' @param x Vector of x-coordinates for each cell.
#' @param y Vector of y-coordinates for each cell.
#' @param features Genes to plot.
#' @param lower_quantile Cut low expression values at this quantile.
#' @param upper_quantile Cut high expression values at this quantile.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_features <- function(counts, x, y, features,
                          lower_quantile = 0.05, upper_quantile = 0.95,
                          filename = NULL) {
  expr_data <-
    counts[features, , drop = FALSE] %>%
    Matrix::t() %>%
    apply(
      2,
      normalize_expression,
      lower_quantile = lower_quantile,
      upper_quantile = upper_quantile
    ) %>%
    as_tibble()
  
  p <-
    tibble(x = {{x}}, y = {{y}}) %>% 
    bind_cols(expr_data) %>% 
    pivot_longer(!1:2, names_to = "feature", values_to = "expression") %>%
    arrange(expression) %>%
    ggplot(aes(x, y)) +
    geom_point(aes(color = expression), size = 0.1) +
    scale_color_scico("relative\nnormalized\nexpression", palette = "bamako") +
    coord_fixed() +
    facet_wrap(vars(feature)) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )

  ggsave_default(filename, width = 420, height = 297)
  p
}

nb_markers <-
  markers %>%
  filter(cell_type == "neuroblastoma") %>%
  pull(gene)

plot_features(
  logcounts(nb),
  x = colData(nb)$umap_1_monocle,
  y = colData(nb)$umap_2_monocle,
  nb_markers,
  filename = "markers/neuroblastoma"
)



#' Plot markers for NB cells and cells with highest scores in a selected gene
#' signature. Group cells by cluster.
#'
#' @param signature_col Column with gene signature scores.
#' @param top_prop Fraction of cells with highest scores to select.
#' @param filename Name of the output file.
#'
#' @return A ggplot2 object
plot_nb_markers_of_high_signature_cells <- function(signature_col,
                                                    top_prop = 0.05,
                                                    filename = NULL) {
  selected_cells <- union(
    colData(nb) %>% 
      as_tibble(rownames = "cell") %>% 
      slice_max(prop = top_prop, order_by = {{signature_col}}) %>% 
      pull(cell),
    colData(nb) %>% 
      as_tibble(rownames = "cell") %>% 
      filter(cellont_abbr == "NB") %>% 
      pull(cell)
  )
  
  nb_subset <- nb[, selected_cells]
  
  nb_markers <-
    markers %>% 
    filter(cell_type == "neuroblastoma") %>% 
    pull(gene)
  
  p <- plot_dots(
    logcounts(nb_subset),
    nb_markers,
    colData(nb_subset)$cluster_50
  )
  
  ggsave_default(filename, height = 100)
  p
}

plot_nb_markers_of_high_signature_cells(signature_mesenchymal,
                                        filename = "markers/nb_mesenchymal")
plot_nb_markers_of_high_signature_cells(signature_ncc_like,
                                        filename = "markers/nb_ncc_like")



# PanglaoDB ---------------------------------------------------------------

panglaodb <- 
  read_tsv("metadata/PanglaoDB_markers_27_Mar_2020.tsv.gz") %>% 
  set_names(colnames(.) %>% str_replace_all(" ", "_")) %>%
  filter(
    str_detect(species, "Hs"),
    canonical_marker == 1
  ) %>% 
  mutate(cell_type = str_replace(cell_type, "/", " or "))

panglaodb_cell_types <- 
  panglaodb %>%
  filter(
    organ %in% c("Immune system", "Blood", "Brain") |
    cell_type == "Hematopoietic stem cells"
  ) %>% 
  count(organ, cell_type)


#' Plot canonical markers as defined by the PanglaoDB database in a dotplot
#' and highlight manual cluster classifications.
#'
#' @param cell_type Cell type available in the database.
#' @param cluster_col Column of the cell metadata that contains cluster ids.
#' @param save_plot If `TRUE`, save the plot in an automatically derived folder.
#'
#' @return A ggplot object
plot_panglaodb_markers <- function(cell_type, cluster_col,
                                   folder = "panglaodb", save_plot = TRUE) {
  info("Plotting markers for {cell_type}")
  
  counts <- logcounts(nb)
  
  markers <-
    panglaodb %>%
    filter(cell_type == {{cell_type}}) %>%
    pull(official_gene_symbol) %>%
    unique() %>% 
    str_sort(decreasing = TRUE)
  
  organ <- 
    panglaodb %>%
    filter(cell_type == {{cell_type}}) %>%
    slice(1) %>% 
    pull(organ)
  
  known_markers <- intersect(markers, rownames(counts))
  missing_markers <- setdiff(markers, rownames(counts))
  if (length(missing_markers) > 0)
    missing_markers <- str_glue(
      "Not observed: {str_c(missing_markers, collapse = ', ')}"
    )
  else
    missing_markers <- waiver()

  figure_height <- max(6 * length(known_markers), 70)
  
  x_annotation_data <-
    tibble(level = levels(colData(nb)$cellont_cluster)) %>% 
    extract(level, into = "label", regex = "(\\w+)") %>%
    mutate(label = as_factor(label), r = row_number()) %>%
    group_by(label) %>%
    summarize(
      xmin = first(r) - 0.5,
      xmax = last(r) + 0.5,
      label_x = mean(r)
    ) %>%
    mutate(
      fill = case_when(
        row_number() %% 2 == 0 ~ "white",
        TRUE                   ~ "gray90"
      )
    )
  
  p <-
    plot_dots(
      counts,
      known_markers,
      colData(nb) %>% as_tibble() %>% pull({{cluster_col}}),
      panel_annotation = x_annotation_data
    ) +
    scale_x_discrete(labels = function(x) str_replace(x, " ", "\n")) +
    labs(
      title = cell_type,
      subtitle = "PanglaoDB canonical markers",
      caption = missing_markers
    )

  if (save_plot)
    ggsave_default(str_glue("markers/{folder}/{organ}/{cell_type}"),
                   height = figure_height, width = 250)
  p
}


plot_panglaodb_markers("B cells", cellont_cluster, folder = "panglaodb")

walk(
  panglaodb_cell_types$cell_type,
  plot_panglaodb_markers,
  cluster_50,
  folder = "panglaodb"
)



# Receptors/ligands -------------------------------------------------------

plot_pathway_genes <- function(pathways = NULL,
                               interactions = NULL,
                               y_group = c("pathway", "interaction"),
                               filename = NULL,
                               ...) {
  y_group <- match.arg(y_group)
  
  complexes <-
    CellChat::CellChatDB.human$complex %>% 
    as_tibble(rownames = "complex") %>% 
    mutate(across(subunit_1:subunit_4, na_if, "")) %>% 
    unite(!complex, col = "genes", na.rm = TRUE) %>%
    deframe()
  
  interaction_genes <-
    CellChat::CellChatDB.human$interaction %>%
    as_tibble() %>%
    {
      if (!is.null(pathways))
        filter(., pathway_name %in% pathways)
      else
        .
    } %>% 
    {
      if (!is.null(interactions))
        filter(., interaction_name %in% interactions)
      else
        .
    } %>% 
    select(interaction_name:receptor) %>%
    pivot_longer(c(ligand, receptor), names_to = "type", values_to = "genes") %>%
    mutate(genes = genes %>% recode(!!!complexes) %>% str_split("_")) %>%
    unchop(genes) %>%
    distinct(genes, .keep_all = TRUE) %>%
    arrange(pathway_name, type)
  
  x_annotation_data <-
    tibble(level = levels(colData(nb)$cellont_cluster)) %>% 
    extract(level, into = "label", regex = "(\\w+)") %>%
    mutate(label = as_factor(label), r = row_number()) %>%
    group_by(label) %>%
    summarize(
      xmin = first(r) - 0.5,
      xmax = last(r) + 0.5,
      label_x = mean(r)
    ) %>%
    mutate(
      fill = case_when(
        row_number() %% 2 == 0 ~ "white",
        TRUE                   ~ "gray90"
      )
    )
  
  y_annotation_data <-
    interaction_genes %>%
    {
      if (y_group == "pathway")
        mutate(., label = as_factor(pathway_name), r = row_number())
      else
        mutate(., label = as_factor(interaction_name), r = row_number())
    } %>% 
    group_by(label) %>%
    summarise(
      yintercept = first(r) - 0.5,
      label_y = mean(r)
    )
  
  p <-
    plot_dots(
      logcounts(nb),
      interaction_genes$genes,
      colData(nb)$cellont_cluster,
      panel_annotation = x_annotation_data
    ) +
    scale_x_discrete(labels = function(x) str_replace(x, " ", "\n")) +
    geom_hline(
      yintercept = y_annotation_data$yintercept[-1],
      linetype = "dashed",
      size = 0.25
    ) +
    geom_text(
      data = y_annotation_data,
      aes(
        x = nlevels(colData(nb)$cellont_cluster) + 2,
        y = label_y,
        label = label
      ),
      size = 4,
      hjust = 0
    ) +
    new_scale_fill() +
    geom_tile(
      data = interaction_genes %>% mutate(y = row_number()),
      aes(x = .25, y = y, width = .5, height = 1, fill = type)
    ) +
    scale_fill_manual(values = c("#e66101", "#5e3c99")) +
    theme(
      axis.line = element_blank(),
      legend.position = "left"
    ) +
    NULL
  
  ggsave_default(filename, ...)
  p
}

plot_pathway_genes(
  pathways = c("TNF", "IFN-I", "IFN-II", "CD6", "CD46", "CD99"),
  filename = "markers/receptors_ligands_pathways"  
)

plot_pathway_genes(
  interactions = c("PTN_NCL", "NAMPT_INSR", "MIF_CD74_CXCR4", "MIF_CD74_CD44",
                   "MDK_NCL", "MDK_ITGA4_ITGB1", "HLA-A_CD8A", "COL6A1_CD44",
                   "CD99_PILRA", "CD99_CD99", "APP_CD74", "ALCAM_CD6"),
  y_group = "interaction",
  filename = "markers/receptors_ligands_interactions",
  width = 420
)


# CellChat::CellChatDB.human$interaction %>%
#   as_tibble() %>% 
#   filter(receptor == "CD8A")



# Publication figures -----------------------------------------------------

## Figure 1d ----

selected_markers <-
  read_csv("metadata/cell_markers_publication.csv", comment = "#")

plot_canonical_markers <- function() {
  # remove other cells
  nb <- nb[, colData(nb)$cellont_abbr != "other"]
  colData(nb)$cellont_abbr <- fct_drop(colData(nb)$cellont_abbr)
  colData(nb)$cellont_cluster <- fct_drop(colData(nb)$cellont_cluster)
  
  y_annotation_data <-
    selected_markers %>%
    arrange(desc(row_number())) %>%
    mutate(r = row_number()) %>%
    group_by(label = cell_type) %>%
    summarise(
      yintercept = first(r) - 0.5,
      label_y = mean(r)
    )
  
  x_annotation_data <-
    tibble(level = levels(colData(nb)$cellont_cluster)) %>% 
    extract(level, into = "label", regex = "(\\w+)") %>%
    mutate(label = as_factor(label), r = row_number()) %>%
    group_by(label) %>%
    summarize(
      xmin = first(r) - 0.5,
      xmax = last(r) + 0.5,
      label_x = mean(r)
    ) %>%
    mutate(
      fill = case_when(
        row_number() %% 2 == 0 ~ "white",
        TRUE                   ~ "gray90"
      )
    )
  
  n_col <- nlevels(colData(nb)$cellont_cluster)
  
  p <-
    plot_dots(
      logcounts(nb),
      rev(selected_markers$gene),
      colData(nb)$cellont_cluster,
      panel_annotation = x_annotation_data
    ) +
    scale_x_discrete(
      name = "cluster / cell type",
      labels = function(x) str_replace(x, "(.+) \\((.+)\\)", "\\2\n\\1"),
      expand = expansion(add = 0.5)
    ) +
    scale_y_discrete(NULL, expand = expansion(add = 0.5)) +
    scale_color_dotplot(
      "scaled average expression",
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(15, "mm"),
        label.position = "top",
        title.vjust = 0.1
      )
    ) +
    scale_radius("% expressed", range = c(0, 2.5)) +
    coord_fixed(
      xlim = c(1, n_col),
      clip = "off"
    ) +
    geom_hline(
      yintercept = y_annotation_data$yintercept,
      linetype = "dashed",
      size = 0.25
    ) +
    geom_text(
      data = y_annotation_data,
      aes(
        x = n_col + 1,
        y = label_y,
        label = label
      ),
      size = BASE_TEXT_SIZE_MM,
      hjust = 0
    ) +
    theme_nb(grid = FALSE) +
    theme(
      legend.box.just = "bottom",
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.position = "top",
      legend.spacing = unit(0, "mm"),
      legend.margin = margin(0, 1, -2, 1, "mm"),
      panel.border = element_rect(color = "black", size = .25)
    )
  p
}

plot_canonical_markers()
ggsave_publication("1d_markers", height = 12, width = 7)


## Figure S1c ----
  
plot_mesenchymal <- function(top_prop = 0.05) {
  data_highlight <-
    bind_rows(
      mesenchymal =
        nb_metadata %>%
        slice_max(signature_mesenchymal, prop = top_prop) %>% 
        select(
          cell, umap_1_monocle, umap_2_monocle,
          value = signature_mesenchymal
        ),
      `NCC-like` =
        nb_metadata %>%
        slice_max(signature_ncc_like, prop = top_prop) %>% 
        select(
          cell, umap_1_monocle, umap_2_monocle,
          value = signature_ncc_like
        ),
      .id = "signature"
    ) %>% 
    arrange(value)
  
  nb_metadata %>%
    ggplot(aes(umap_1_monocle, umap_2_monocle)) +
    geom_point(color = "gray90", size = 0.001, shape = 16) +
    geom_point(
      data = data_highlight,
      aes(color = value),
      size = 0.001,
      shape = 16
    ) +
    scale_x_continuous("UMAP1", breaks = c(-10, 0, 10)) +
    scale_y_continuous("UMAP2", breaks = c(-10, 0, 10)) +
    scale_color_distiller(
      name = "gene\nsign.\nscore",
      palette = "RdPu",
      guide = guide_colorbar(
        barwidth = unit(2, "mm"),
        barheight = unit(15, "mm")
      )
    ) +
    coord_fixed() +
    facet_wrap(vars(signature)) +
    theme_nb(grid = FALSE) +
    theme(
      legend.position = c(.95, .32)
    )
}

plot_mesenchymal()
ggsave_publication("S1c_mesenchymal", type = "png", width = 9, height = 5)


## Figure S1d ----

subplot_nb_dots <- function(signature_col, title = NULL, top_prop = 0.05) {
  selected_cells <- union(
    colData(nb) %>% 
      as_tibble(rownames = "cell") %>% 
      slice_max(prop = top_prop, order_by = {{signature_col}}) %>% 
      pull(cell),
    colData(nb) %>% 
      as_tibble(rownames = "cell") %>% 
      filter(cellont_abbr == "NB") %>% 
      pull(cell)
  )
  
  nb_subset <- nb[, selected_cells]
  
  nb_markers <-
    markers %>% 
    filter(cell_type == "neuroblastoma") %>% 
    pull(gene)
  
  plot_dots(
    logcounts(nb_subset),
    nb_markers,
    colData(nb_subset)$cluster_50
  ) +
    scale_x_discrete("cluster") +
    scale_y_discrete(NULL) +
    scale_color_dotplot(
      "scaled average expression",
      limits = c(-0.5, 2.5),
      breaks = c(0, 1, 2),
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(15, "mm"),
        label.position = "top",
        title.vjust = 0.1
      )
    ) +
    scale_radius(
      "% expressed",
      range = c(0, 2.5),
      limits = c(0, 100)
    ) +
    ggtitle(title) +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      plot.title = element_text(
        size = BASE_TEXT_SIZE_PT,
        hjust = 0.5,
        margin = margin(b = 1, unit = "mm")
      )
    )
}

wrap_plots(
  guide_area(),
  subplot_nb_dots(signature_mesenchymal, "mesenchymal"),
  subplot_nb_dots(signature_ncc_like, "NCC-like") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ),
  guides = "collect",
  design = "AA\nBC",
  heights = c(0.2, 0.8),
  widths = c(12, 19)
) &
  theme(
    legend.box.just = "bottom",
    legend.key.height = unit(1, "mm"),
    legend.key.width = unit(1, "mm"),
    legend.position = "top",
    legend.spacing = unit(0, "mm"),
    legend.margin = margin(0, 1, 0, 1, "mm"),
    plot.margin = margin(0, 1, 0, 1, "mm"),
    
  )

ggsave_publication("s1d_dots", width = 9, height = 4)
