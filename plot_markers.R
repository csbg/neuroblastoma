# Plot canonical cell type and neuroblastoma markers.
# Exports several plots to plots/markers.
#
# @DEPI metadata.rds
# @DEPI assay_sct_seurat.rds

library(Seurat)
library(tidyverse)
library(scico)
library(gginnards)
source("common_functions.R")



# Load data ---------------------------------------------------------------

nb <- readRDS("data_generated/assay_sct_seurat.rds")

nb@meta.data <-
  readRDS("data_generated/metadata.rds") %>%
  column_to_rownames("cell")



# Canonical cell type markers ---------------------------------------------

#' Make a Seurat-style dotplot of features vs clusters.
#'
#' @param counts Count matrix whose rownames correspond to features.
#' @param features Vector of features, plotted in the given order.
#' @param groups Factor of groups to which each cell belongs, plotted in the
#'   order of its levels.
#' @param min_exp Lower limit for the scaled average expression.
#' @param max_exp Upper limit for the scaled average expression.
#' @param panel_annotation A dataframe used for drawing rectangles on the panel
#'   background. Must contain three columns "xmin", "xmax" (aesthetics for
#'   `geom_rect()`), and "fill" (color name).
#'
#' @return A ggplot object.
plot_dots <- function(counts, features, groups,
                      min_exp = -2.5, max_exp = 2.5, panel_annotation = NULL) {
  scale_and_limit <- function(x) {
    scale(x)[,1] %>% 
      pmax(min_exp) %>% 
      pmin(max_exp)
  }
  
  vis_data <- 
    counts[features, , drop = FALSE] %>% 
    Matrix::t() %>% 
    as_tibble(rownames = "cell") %>% 
    group_by(id = groups) %>% 
    summarise(
      across(
        where(is.numeric),
        list(
          avg_exp = ~mean(expm1(.)),
          pct_exp = ~length(.[. > 0]) / length(.) * 100
        ),
        .names = "{.col}__{.fn}"
      )
    ) %>% 
    mutate(across(ends_with("avg_exp"), scale_and_limit)) %>%
    pivot_longer(
      !id,
      names_to = c("feature", ".value"),
      names_pattern = "(.+)__(.+)"
    ) %>% 
    mutate(feature = factor(feature, levels = features))
  
  if (!is.null(panel_annotation)) {
    panel_bg <- list(
      geom_point(color = "white"),  # initialize discrete coordinate system
      geom_rect(
        data = panel_annotation,
        aes(xmin = xmin, xmax = xmax,
            ymin = 0.5, ymax = nlevels(vis_data$feature) + 0.5,
            fill = fill),
        show.legend = FALSE,
        inherit.aes = FALSE,
      ),
      scale_fill_identity()
    )
  } else {
    panel_bg <- NULL
  }
    
  ggplot(vis_data, aes(id, feature)) +
    panel_bg +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    scale_x_discrete("cluster", expand = c(0, 0)) +
    scale_y_discrete("canonical marker", expand = c(0, 0)) +
    scale_color_scico(
      "scaled\naverage\nexpression",
      palette = "oslo",
      direction = -1,
      aesthetics = "color"
    ) +
    scale_radius("% expressed", range = c(0, 6)) +
    coord_fixed(
      # xlim = c(0.5, nlevels(vis_data$id) + 0.5),
      # ylim = c(0.5, nlevels(vis_data$feature) + 0.5),
      clip = "off"
    ) +
    theme_classic() +
    theme(panel.grid = element_blank())
}

markers <- read_csv("metadata/cell_markers.csv", comment = "#")

cluster_info <-
  read_csv("metadata/clusters.csv") %>%
  mutate(
    cell_type =
      as_factor(cell_type) %>%
      fct_relevel("T", "NK", "B", "M", "D", "E", "NB", "other")
  ) %>%
  arrange(cell_type, cluster)

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
  cluster_info %>%
  mutate(r = row_number()) %>%
  group_by(label = cell_type) %>%
  summarize(
    xmin = first(r) - 0.5,
    xmax = last(r) + 0.5,
    label_x = mean(r)
  ) %>%
  mutate(fill = case_when(row_number() %% 2 == 0 ~ "white", TRUE ~ "gray90"))

plot_dots(
  nb$SCT@data,
  rev(markers$gene),
  fct_relevel(
    nb@meta.data$cluster_50,
    as.character(cluster_info$cluster)
  ),
  panel_annotation = x_annotation_data
) +
  geom_text(
    data = x_annotation_data,
    aes(x = label_x, y = length(markers$gene) + 2, label = label),
    size = 4
  ) +
  geom_hline(
    yintercept = y_annotation_data$yintercept,
    linetype = "dashed",
    size = 0.25
  ) +
  geom_text(
    data = y_annotation_data,
    aes(x = 22, y = label_y, label = label),
    size = 4,
    hjust = 0
  ) +
  theme(legend.position = "left") +
  NULL


ggsave_default("markers/canonical_markers_clusters",
               height = 297, width = 250)



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
  nb$SCT@counts,
  x = nb@meta.data$umap_1_monocle,
  y = nb@meta.data$umap_2_monocle,
  nb_markers,
  filename = "markers/nb_markers"
)
