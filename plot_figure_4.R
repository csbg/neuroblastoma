# @DEPI ccc_cellchat_object.rds
# @DEPI ccc_signaling_data.rds

library(CellChat)
library(tidyverse)
library(ComplexHeatmap)
library(latex2exp)
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

cellchat <- readRDS("data_generated/ccc_cellchat_object.rds") 
sig_data <- readRDS("data_generated/ccc_signaling_data.rds")



# Figures -----------------------------------------------------------------

## 4a ----

plot_n_interactions <- function() {
  mat <- cellchat@net$count
  bar_colors <- CELL_TYPE_COLORS[colnames(mat)]
  total_target <- colSums(mat)
  total_source <- rowSums(mat)
  bar_ylim <- c(0, max(total_source, total_target))
  
  Heatmap(
    mat,
    col = RColorBrewer::brewer.pal(9, "Reds"),
    name = "Number of interactions",
    
    cluster_rows = FALSE, 
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    row_names_side = "left",
    row_title = "Source (ligand)", 
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    
    cluster_columns = FALSE, 
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    column_title = "Target (receptor)",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    
    width = unit(20, "mm"),
    height = unit(20, "mm"),
    
    top_annotation = HeatmapAnnotation(
      count_bar = anno_barplot(
        total_target,
        ylim = bar_ylim,
        border = FALSE, 
        axis = FALSE,
        gp = gpar(fill = bar_colors, col = bar_colors)
      ),
      count_text = anno_text(
        total_target, 
        gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      ), 
      simple_anno_size_adjust = TRUE,
      show_annotation_name = FALSE
    ), 
    
    right_annotation = rowAnnotation(
      count_text = anno_text(
        total_source, 
        gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      ), 
      count_bar = anno_barplot(
        total_source,
        ylim = bar_ylim,
        border = FALSE, 
        axis = FALSE,
        gp = gpar(fill = bar_colors, col = bar_colors)
      ),
      simple_anno_size_adjust = TRUE,
      show_annotation_name = FALSE
    ), 
    
    heatmap_legend_param = list(
      at = c(min(mat), max(mat)),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      title_position = "leftcenter-rot", 
      border = NA, 
      legend_height = unit(20, "mm"), 
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      grid_width = unit(2, "mm")
    )
  )
}

(p <- plot_n_interactions())
ggsave_publication("4a_n_interactions", plot = p, width = 6, height = 5)



## 4b ----

plot_contribution_celltype <- function(cell_type = "NB", signif = 0.05) {
   sig_data %>%
    filter(
      pathway_name %in% cellchat@netP$pathways, 
      source == {{cell_type}}, 
      target != {{cell_type}},
      pval <= signif
    ) %>%
    group_by(interaction = interaction_name_2) %>%
    summarise(
      prob.norm = sum(prob.norm),
      count = n(),
      pathway_name = first(pathway_name)
    ) %>%
    ungroup() %>%
    mutate(prob.avg = prob.norm / count) %>% 
    mutate(interaction = fct_reorder(interaction, prob.norm)) %>% 
    ggplot(aes(interaction, prob.norm)) +
    geom_col(fill = "#23cfc3", width = 0.85) +
    geom_text(
      aes(y = 2.68, label = pathway_name),
      size = BASE_TEXT_SIZE_MM
    ) +
    xlab("Ligand-receptor pair") +
    scale_y_continuous(
      "Total outgoing communication score",
      expand = expansion(0)
    ) +
    coord_flip() +
    theme_nb(grid = FALSE) +
    theme(
      axis.ticks.y = element_blank(),
      panel.border = element_blank()
    )
}

plot_contribution_celltype()
ggsave_publication("4b_outgoing_comm_score", width = 5.5, height = 5)


## 4c ----

plot_selected_dots <- function(pathways, source_type = "NB") {
  vis_data <- 
    sig_data %>% 
    filter(
      pathway_name %in% {{pathways}}, 
      source == {{source_type}}
    ) %>%
    mutate(pathway_name = factor(pathway_name, levels = pathways))
  
  ggplot(vis_data, aes(target, interaction_name_2)) +
    geom_point(
      aes(color = prob.norm, size = pmin(5, -log10(pval)))
    ) + 
    scale_radius(
      TeX("-log_{10} p-Value_{cap.}"),
      range = c(1, 5),
      guide = "none"
    ) +
    scale_color_distiller(
      "relative communication score",
      type = "seq",
      palette = "Reds", 
      direction = 1,
      breaks = round(range(vis_data$prob.norm), 2),
      guide = guide_colorbar(
        barheight = unit(15, "mm"),
        barwidth = unit(2, "mm"),
        ticks = FALSE,
        title.vjust = .9
      )
    ) +
    xlab("Cell types (source -> target)") +
    ylab("Ligand-receptor pairs and pathways") +
    facet_grid(
      vars(pathway_name), 
      vars(source),
      scales = "free",
      space = "free",
      switch = "both"
    ) + 
    theme_nb(grid = FALSE) + 
    theme(
      axis.text.x = element_text(
        size = BASE_TEXT_SIZE_PT, 
        margin = margin(t = 1, b = 1, r = 1, l = 1, unit = "mm")
      ),
      axis.text.y = element_text(size = BASE_TEXT_SIZE_PT),
      axis.ticks.length.x = unit(2, "mm"),
      axis.ticks.x = element_line(
        arrow = arrow(length = unit(1,"mm"), ends = "first", type = "closed")
      ),
      axis.title = element_text(size = BASE_TEXT_SIZE_PT),
      
      panel.border = element_rect(color = "black", size = BASE_LINE_SIZE),
      panel.spacing = unit(-BASE_LINE_SIZE / 2, "pt"),
      
      strip.background.y = element_rect(
        fill = "#FBDE96",
        size = BASE_LINE_SIZE
      ),
      strip.background.x = element_rect(
        fill = "lightblue", 
        size = BASE_LINE_SIZE
      ),
      strip.text.x =  element_text(
        size = BASE_TEXT_SIZE_PT, 
        margin = margin(t = 1, b = 1, r = 4, l = 4, unit = "mm"),
        angle = 0
      ),
      strip.text.y.left = element_text(
        size = BASE_TEXT_SIZE_PT, 
        margin = margin(t = 0.5, b = 0.5, r = 1, l = 1, unit = "mm"),
        angle = 0
      ),
      strip.switch.pad.grid = unit(0, "pt"),
      
      legend.box.spacing = unit(1, "mm"),
      legend.key.size = unit(5, "mm"),
      legend.margin = margin(t = 2, b = 0, r = 0, l = 0, unit = "mm"),
      legend.position = "right",
      legend.text = element_text(size = BASE_TEXT_SIZE_PT),
      legend.title = element_text(size = BASE_TEXT_SIZE_PT, angle = 90),
    )
}

plot_selected_dots(c("MK", "MIF", "COLLAGEN", "PTN", "APP", "ALCAM", "THY1"))
ggsave_publication("4c_dots", height = 5, width = 7.5)



## 4d ----

plot_centrality <- function(pathway) {
  centralities <- cellchat@netP$centr[[pathway]]
  
  mat <- 
    centralities %>% 
    as_tibble() %>% 
    select(
      Sender = outdeg,
      Receiver = indeg,
      Mediator = flowbet,
      Influencer = info
    ) %>%
    mutate(
      across(everything(), ~. / max(.)),
      cell_type = names(centralities$outdeg)
    ) %>%
    column_to_rownames("cell_type") %>% 
    as.matrix() %>% 
    t()
  
  Heatmap(
    mat, 
    col = RColorBrewer::brewer.pal(9, "YlOrRd"),
    name = "Importance",
    
    cluster_rows = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    
    cluster_columns = FALSE, 
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    column_title = str_glue("{pathway} signaling pathway network"),
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    
    width = unit(40, "mm"),
    height = unit(20, "mm"), 
    
    bottom_annotation = HeatmapAnnotation(
      group = colnames(mat), 
      col = list(group = CELL_TYPE_COLORS[colnames(mat)]), 
      which = "column", 
      show_legend = FALSE, 
      show_annotation_name = FALSE, 
      simple_anno_size = grid::unit(0.2, "cm")
    ), 
    
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT, fontface = "plain"), 
      title_position = "leftcenter-rot", 
      border = NA, 
      at = c(0, 1), 
      legend_height = unit(20, "mm"), 
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      grid_width = unit(2, "mm")
    )
  )
}

(p <- plot_centrality("MIF"))
ggsave_publication("4d_importance_MIF", plot = p, width = 6, height = 3.5)



## 4e ----

(p <- plot_centrality("MK"))
ggsave_publication("4e_importance_MK", plot = p, width = 6, height = 3.5)



## 4f ----

make_matrix <- function(gene, cell_type) {
  metadata <- 
    cellchat@meta %>% 
    as_tibble(rownames = "cell")
  
  barcodes <- 
    metadata %>% 
    filter(cellont_abbr == {{cell_type}}) %>% 
    pull(cell)
  
  cellchat@data.signaling[gene, barcodes, drop = FALSE] %>%
    t() %>%
    as.matrix() %>%
    magrittr::set_colnames("logexp") %>%
    as_tibble(rownames = "cell") %>%
    left_join(metadata, by = "cell") %>%
    transmute(
      gene = gene,
      cell_type = cell_type,
      sample = rename_patients(sample),
      group = rename_groups(group),
      logexp = logexp / max(logexp)
    )
}

plot_violin <- function(genes, cell_types) {
  plot_data <-
    list(gene = genes, cell_type = cell_types) %>%
    cross_df() %>% 
    pmap_dfr(make_matrix) %>% 
    mutate(gene = as_factor(gene), cell_type = as_factor(cell_type))
    
  ggplot(plot_data, aes(sample, logexp)) +
    geom_violin(
      aes(color = group, fill = group),
      size = BASE_LINE_SIZE,
      scale = "width",
      width = 0.8,
      show.legend = FALSE
    ) +
    stat_summary(geom = "point", fun = mean, size = .2) +
    xlab("patient") +
    scale_y_continuous(
      "log-normalized expression",
      limits = c(0, 1)
    ) +
    scale_fill_manual(values = GROUP_COLORS, aesthetics = c("color", "fill")) +
    facet_grid(
      vars(gene),
      vars(cell_type),
      scales = "free_x",
      space = "free_x",
      switch = "y"
    ) +
    theme_nb(grid = FALSE) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.placement = "outside"
    )
}

plot_violin(
  genes = c("MIF", "CD74", "CXCR4", "CD44"),
  cell_types = c("NB", "M")
)
ggsave_publication("4f_violins_MIF", width = 8.5, height = 9)



## 4g ----

plot_violin(
  genes = c("MDK", "LRP1", "NCL", "ALK"),
  cell_types = c("NB", "M")
)
ggsave_publication("4g_violins_MK", width = 8.5, height = 9)



# Tables ------------------------------------------------------------------

## S6 ----

sig_data %>% 
  select(
    "Cell type expressing ligand" = source,
    "Cell type expressing receptor" = source,
    Ligand = ligand,
    Receptor = receptor,
    Score = prob,
    "Normalized score" = prob.norm,
    "P value" = pval,
    Pathway = pathway_name,
    "Type of interaction" = annotation
  ) %>%
  save_table("S6_ccc", sheet_name = "Interactions")
