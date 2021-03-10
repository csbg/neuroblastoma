# Plot canonical cell type and neuroblastoma markers.
# Exports several plots to plots/markers.

library(Seurat)
library(tidyverse)
library(scico)
source("common_functions.R")



# Load data ---------------------------------------------------------------

nb <- readRDS("data_generated/assay_SCT_seurat.rds")

nb@meta.data <-
  readRDS("data_generated/metadata.rds") %>%
  column_to_rownames("cell")



# Canonical cell type markers ---------------------------------------------

markers <- read_csv("metadata/cell_markers.csv", comment = "#")

Idents(nb) <- 
  fct_relevel(
    nb@meta.data$cluster_50,
    # T cells
    "3", "5", "18",
    
    # NK cells
    "4", "6",
    
    # B cells
    "2", "9", "12", "17", "19", "21",
    
    # monocyte
    "1", "15", "16", "22",
    
    # pDC
    "14",
    
    # erythroblast
    "13",
    
    # bone marrow
    "10", "11", 
    
    # other
    "7", "20",
    
    # NB
    "8"
  )

DotPlot(nb, features = rev(markers$gene)) +
  coord_flip() +
  scale_color_scico(palette = "oslo", direction = -1)
ggsave_default("markers/canonical_markers_clusters",
               height = 297, width = 210)
ggsave_default("markers/canonical_markers_clusters",
               type = "pdf", height = 297, width = 210, crop = FALSE)



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
#' @param data A Seurat object.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param features Genes to plot.
#' @param lower_quantile Cut low expression values at this quantile.
#' @param upper_quantile Cut high expression values at this quantile.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_features <- function(data, x, y, features,
                          lower_quantile = 0.05, upper_quantile = 0.95,
                          filename = NULL) {
  expr_data <-
    data$RNA@data[features, , drop = FALSE] %>%
    as.matrix() %>% 
    t() %>%
    apply(
      2,
      normalize_expression,
      lower_quantile = lower_quantile,
      upper_quantile = upper_quantile
    ) %>%
    as_tibble()
  
  p <-
    bind_cols(
      data@meta.data %>%
        as_tibble() %>%
        select({{x}}, {{y}}),
      expr_data
    ) %>%
    pivot_longer(!1:2, names_to = "feature", values_to = "expression") %>%
    arrange(expression) %>%
    ggplot(aes({{x}}, {{y}})) +
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

plot_features(nb, umap_1_monocle, umap_2_monocle, nb_markers,
              filename = "markers/nb_markers_new")
