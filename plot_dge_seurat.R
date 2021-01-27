library(Seurat)
library(tidyverse)
library(fs)
library(patchwork)
library(viridis)
library(pheatmap)
library(grid)

ggsave_default <- function(filename, width = 297, height = 210,
                           crop = TRUE, ...) {
  if (is.null(filename))
    return()
  
  filename <- str_glue("plots/{filename}.png")
  filename %>% path_dir() %>% dir_create()
  
  ggsave(filename, dpi = 300, units = "mm", limitsize = FALSE,
         width = width, height = height, ...)
  
  if (crop) knitr::plot_crop(filename)
  
  invisible(filename)
}



# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/all_datasets_current/nb_assay_RNA.rds") %>% 
  NormalizeData()
nb@meta.data <-
  readRDS("data_generated/all_datasets_current/nb_metadata.rds") %>%
  column_to_rownames("cell")



# DE genes ----------------------------------------------------------------

de_genes_cluster <-
  read_csv("data_generated/all_datasets_current/dge_seurat_groups_clusterwise.csv")
nb <- ScaleData(nb, features = de_genes_cluster$gene)


plot_de_gcwise <- function(data, de_genes, group_1, group_2, cluster,
                           top_n = 10, col_min = NULL, col_max = NULL,
                           silent = FALSE) {
  message("Plotting DE genes in cluster ", cluster, ", ",
          group_1, " vs ", group_2)
  
  top_genes <- 
    de_genes %>% 
    filter(group == {{group_1}}, cluster == {{cluster}}) %>% 
    slice_max(n = top_n, order_by = avg_logFC, with_ties = FALSE) %>% 
    pull(gene)
  if (length(top_genes) == 0)
    return(list(gtable = textGrob("(no data)")))
  
  cell_data <- 
    data@meta.data %>% 
    filter(
      group %in% c({{group_1}}, {{group_2}}),
      refined_cluster == {{cluster}}
    ) %>% 
    arrange(desc(group))
  
  mat <- data$RNA@scale.data[top_genes, rownames(cell_data)]
  col_min <- col_min %||% quantile(mat, 0.05)
  col_max <- col_max %||% quantile(mat, 0.95)
  
  pheatmap(
    mat,
    color = viridis(100),
    breaks = seq(col_min, col_max, length.out = 101),
    legend = FALSE,
    border_color = NA,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    annotation_col = cell_data %>% select(sample, group),
    annotation_legend = FALSE,
    annotation_colors = list(
      group = c(I = "black", II = "#e41a1c", III = "#377eb8", IV = "#4daf4a")
    ),
    silent = silent
  )
}

plot_de_gcwise(nb, de_genes_cluster, group_1 = "II", group_2 = "I", cluster = "5a")


plot_de <- function(data, de_genes, groups, clusters) {
  col_min <- quantile(data$RNA@scale.data, 0.05)
  col_max <- quantile(data$RNA@scale.data, 0.95)
  label_groups <- c("label", groups)
  
  map(
    clusters,
    function(c) {
      plots <- map(
        label_groups,
        function(g) {
          if (g == "label") {
            textGrob(c, gp = gpar(fontsize = 20))
          } else {
            p <- plot_de_gcwise(data, de_genes, group_1 = g, group_2 = "I",
                                cluster = c, col_min = col_min, col_max = col_max,
                                silent = TRUE)
            p$gtable  
          }
        }
      )
    } 
  ) %>%
    flatten() %>% 
    wrap_plots(nrow = length(clusters), widths = c(1, rep(10, length(groups))))
}

plot_de(nb, de_genes_cluster,
        groups = c("II", "III", "IV"), clusters = c("5a", "5b"))
ggsave_default("de/de_test", height = 100)

plot_de(nb, de_genes_cluster,
        groups = c("II", "III", "IV"), clusters = levels(nb@meta.data$refined_cluster))
ggsave_default("de/de_clusters", height = 1400)