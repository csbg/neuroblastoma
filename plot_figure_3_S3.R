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

## 3a ----

plot_n_interactions <- function() {
  # prepare data
  rc_order <- names(CELL_TYPE_ABBREVIATIONS)
  mat <- cellchat@net$count[rc_order, rc_order]
  total_target <- colSums(mat)
  total_source <- rowSums(mat)
  bar_ylim <- c(0, max(total_source, total_target))
  
  # export source data
  mat %>% 
    as_tibble(rownames = "source") %>% 
    save_table("source_data/figure_3a", "Figure 3a")
  
  # make plot
  Heatmap(
    mat,
    col = RColorBrewer::brewer.pal(9, "YlOrBr"),
    name = "number of\ninteractions",
    
    cluster_rows = FALSE, 
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    row_names_side = "left",
    row_title = "source (ligand)", 
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    
    cluster_columns = FALSE, 
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    column_title = "target (receptor)",
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
        gp = gpar(fill = "gray70", col = "gray70")
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
        gp = gpar(fill = "gray70", col = "gray70")
      ),
      simple_anno_size_adjust = TRUE,
      show_annotation_name = FALSE
    ), 
    
    heatmap_legend_param = list(
      at = c(min(mat), max(mat)),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      title_position = "topleft", 
      border = NA, 
      legend_height = unit(20, "mm"), 
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      grid_width = unit(2, "mm")
    )
  )
}

(p <- plot_n_interactions())
ggsave_publication("3a_n_interactions", plot = p, width = 6, height = 5)



## 3b ----

plot_contribution_celltype <- function(cell_type = "NB", signif = 0.05) {
  # prepare data
  vis_data <- 
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
    mutate(interaction = fct_reorder(interaction, prob.norm))
  
  # export source data
  vis_data %>% 
    select(interaction, pathway_name, total_outgoing_comm_score = prob.norm) %>% 
    arrange(desc(total_outgoing_comm_score)) %>% 
    save_table("source_data/figure_3b", "Figure 3b")
  
  # make plot
  ggplot(vis_data, aes(interaction, prob.norm)) +
    geom_col(fill = "#ffb03e", width = 0.85) +
    geom_hline(yintercept = 0, size = BASE_LINE_SIZE) +
    geom_text(
      aes(y = 2.68, label = pathway_name),
      size = BASE_TEXT_SIZE_MM,
      fontface = "italic"
    ) +
    xlab("ligand-receptor pair") +
    scale_y_continuous(
      "total outgoing communication score from NB",
      expand = expansion(0)
    ) +
    coord_flip() +
    theme_nb(grid = FALSE) +
    theme(
      axis.title.x = element_text(hjust = 1),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey92",
        size = BASE_LINE_SIZE
      )
    )
}

plot_contribution_celltype()
ggsave_publication("3b_outgoing_comm_score", width = 5.5, height = 5)


## 3c ----

plot_selected_dots <- function(interactions, source_type = "NB") {
  #prepare data
  vis_data <- 
    sig_data %>% 
    filter(
      interaction_name %in% {{interactions}}, 
      source == {{source_type}}
    ) %>%
    mutate(
      pathway_name = factor(pathway_name),
      source = fct_recode(source, "from NB (source) to" = "NB"),
      target = factor(target, levels = names(CELL_TYPE_ABBREVIATIONS))
    )
  
  # export source data
  vis_data %>% 
    select(source, target, pathway = pathway_name,
           ligand_receptor_pair = interaction_name_2,
           relative_comm_score = prob.norm) %>% 
    save_table("source_data/figure_3c", "Figure 3c")
  
  # make plot  
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
      palette = "YlOrBr",
      direction = 1,
      breaks = round(range(vis_data$prob.norm), 2),
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(15, "mm"),
        ticks = FALSE,
        title.position = "top"
      )
    ) +
    xlab("cell type (target)") +
    ylab("ligand-receptor pair\nand pathway") +
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
        fill = "gray95",
        size = BASE_LINE_SIZE
      ),
      strip.background.x = element_rect(
        fill = "gray95",
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
      
      legend.position = c(-.4, -.2),
      legend.direction = "horizontal",
      legend.text = element_text(size = BASE_TEXT_SIZE_PT),
      legend.title = element_text(size = BASE_TEXT_SIZE_PT),
    )
}

plot_selected_dots(c("MDK_NCL", "MDK_LRP1", "MIF_CD74_CXCR4", "MIF_CD74_CD44"))
ggsave_publication("3c_dots", height = 3.4, width = 7)



## 3d ----

plot_centrality <- function(pathway,
                            row_names_side = "left",
                            source_data_filename = NULL,
                            source_data_sheet = NULL) {
  #prepare data
  centralities <- cellchat@netP$centr[[pathway]]
  
  mat <- 
    centralities %>% 
    as_tibble() %>% 
    select(
      sender = outdeg,
      receiver = indeg,
      mediator = flowbet,
      influencer = info
    ) %>%
    mutate(
      across(everything(), ~. / max(.)),
      cell_type = names(centralities$outdeg) %>% 
        factor(levels = names(CELL_TYPE_ABBREVIATIONS))
    ) %>%
    arrange(cell_type) %>% 
    column_to_rownames("cell_type") %>% 
    as.matrix() %>% 
    t()
  
  # export source data
  mat %>% 
    as_tibble(rownames = "centrality_measure") %>% 
    save_table(source_data_filename, source_data_sheet)
  
  # make plot
  Heatmap(
    mat, 
    col = RColorBrewer::brewer.pal(9, "PuBu"),
    name = "importance",
    
    cluster_rows = FALSE,
    row_names_side = row_names_side,
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    
    cluster_columns = FALSE, 
    column_names_side = "top",
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    column_title = str_glue("cell type"),
    column_title_side = "top",
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    
    width = unit(20, "mm"),
    height = unit(10, "mm"), 
    
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT, fontface = "plain"), 
      title_position = "leftcenter-rot", 
      border = NA, 
      at = c(0, 1), 
      legend_height = unit(10, "mm"), 
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      grid_width = unit(2, "mm")
    )
  )
}

(p <- plot_centrality(
  "MIF",
  source_data_filename = "source_data/figure_3d",
  source_data_sheet = "Figure 3d"
))
ggsave_publication("3d_importance_MIF", plot = p, width = 5, height = 3)



## 3e ----

(p <- plot_centrality(
  "MK",
  row_names_side = "right",
  source_data_filename = "source_data/figure_3e",
  source_data_sheet = "Figure 3e"
))
ggsave_publication("3e_importance_MK", plot = p, width = 6, height = 3.5)



## 3f ----

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
      cell,
      gene,
      cell_type,
      sample = rename_patients(sample),
      group = rename_groups(group),
      logexp = logexp / max(logexp)
    )
}

plot_violin <- function(genes,
                        cell_types,
                        y_axis_pos = "left",
                        source_data_filename = NULL,
                        source_data_sheet = NULL) {
  # prepare data
  plot_data <-
    list(gene = genes, cell_type = cell_types) %>%
    cross_df() %>% 
    pmap_dfr(make_matrix) %>% 
    mutate(gene = as_factor(gene), cell_type = as_factor(cell_type))
  
  if (y_axis_pos == "left") {
    switch <- "y"
    position <- "left"
  } else {
    switch <- NULL
    position <- "right"
  }
  
  # export source data
  plot_data %>% 
    save_table(source_data_filename, source_data_sheet)
  
  # make plot
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
      "log-normalized expression of gene",
      limits = c(0, 1),
      position = position
    ) +
    scale_fill_manual(values = GROUP_COLORS, aesthetics = c("color", "fill")) +
    facet_grid(
      vars(gene),
      vars(cell_type),
      scales = "free_x",
      space = "free_x",
      switch = switch
    ) +
    theme_nb(grid = FALSE) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y.right = element_text(margin = margin(l = 6.75, unit = "pt")),
      strip.placement = "outside",
      strip.switch.pad.grid = unit(0, "mm")
    )
}

plot_violin(
  genes = c("MIF", "CD74", "CXCR4", "CD44"),
  cell_types = c("NB", "M"),
  source_data_filename = "source_data/figure_3f",
  source_data_sheet = "Figure 3f"
)
ggsave_publication("3f_violins_MIF", width = 8.5, height = 9)



## 3g ----

plot_violin(
  genes = c("MDK", "NCL", "LRP1"),
  cell_types = c("NB", "M"),
  y_axis_pos = "right",
  source_data_filename = "source_data/figure_3g",
  source_data_sheet = "Figure 3g"
)
ggsave_publication("3g_violins_MK", width = 8.5, height = 6.97)



## S3 ----

# modifies CellChat:::netVisual_hierarchy1, which is called by
# netVisual_individual, so that it exports node and edge weights
# as source data for figure S3
save_source_data <- function(prob, vertex_weight, interaction) {
  list(
    edges = as_tibble(prob, rownames = "source"),
    nodes = tibble(cell_type = colnames(prob), weight = vertex_weight)
  ) %>% 
    save_table(str_glue("source_data/figure_S3_{interaction}"))
}

trace(
  netVisual_hierarchy1,
  tracer = substitute(save_source_data(net, vertex.weight, title.name)),
  where = netVisual_individual
)
# untrace(netVisual_hierarchy1, where = netVisual_individual)


plot_netvisual <- function(layout = "hierarchy", width = 16, height = 16) {
  colors_ccc_network <-
    CELL_TYPE_COLORS[c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")]
  
  target_cells <- 2:6
  
  pdf("plots/final/S3_interactions.pdf",
      width = width / 2.54,
      height = height / 2.54)
  
  netVisual_individual(
    cellchat,
    signaling = "MIF",
    pairLR.use = "MIF_CD74_CXCR4",
    vertex.receiver = target_cells,
    color.use = colors_ccc_network,
    layout = layout
  )
  
  netVisual_individual(
    cellchat,
    signaling = "MIF",
    pairLR.use = "MIF_CD74_CD44",
    vertex.receiver = target_cells,
    color.use = colors_ccc_network,
    layout = layout
  )
  
  netVisual_individual(
    cellchat,
    signaling = "MK",
    pairLR.use = "MDK_NCL",
    vertex.receiver = target_cells,
    color.use = colors_ccc_network,
    layout = layout
  )
  
  netVisual_individual(
    cellchat,
    signaling = "MK",
    pairLR.use = "MDK_LRP1",
    vertex.receiver = target_cells,
    color.use = colors_ccc_network,
    layout = layout
  )
  
  dev.off()
}

plot_netvisual()


# combine the resulting source data tables into a single sheet
merge_source_data <- function() {
  write_data <- function(data, i, start_row) {
    
    writeData(
      wb,
      "Figure S3",
      data,
      startRow = start_row + (i - 1) * 26,
      headerStyle = createStyle(textDecoration = "bold")
    )
  }
  
  files <- c(
    "tables/source_data/figure_S3_MIF - (CD74+CXCR4).xlsx",
    "tables/source_data/figure_S3_MIF - (CD74+CD44).xlsx",
    "tables/source_data/figure_S3_MDK - NCL.xlsx",
    "tables/source_data/figure_S3_MDK - LRP1.xlsx"
  )
  
  wb <- createWorkbook()
  addWorksheet(wb, "Figure S3")
  
  iwalk(
    files,
    function(file, i) {
      interaction_name <- str_match(file, "S3_(.+?)\\.xlsx")[,2]
      str_glue("{i}: {interaction_name}") %>% 
        as.character() %>% 
        write_data(i, 1)
      
      str_glue("{i}.1: edge weights") %>% 
        as.character() %>% 
        write_data(i, 3)
      
      data_edges <- read.xlsx(file, sheet = "edges")
      write_data(data_edges, i, 5)
      
      str_glue("{i}.2: node weights") %>% 
        as.character() %>% 
        write_data(i, 15)
      
      data_nodes <- read.xlsx(file, sheet = "nodes")
      write_data(data_nodes, i, 17)
      
      file_delete(file)
    }
  )
  
  saveWorkbook(wb, "tables/source_data/figure_S3.xlsx", overwrite = TRUE)
}

merge_source_data()



# Data --------------------------------------------------------------------

## S2 ----

sig_data %>% 
  select(
    "Cell type expressing ligand" = source,
    "Cell type expressing receptor" = target,
    Ligand = ligand,
    Receptor = receptor,
    Score = prob,
    "Normalized score" = prob.norm,
    "P value" = pval,
    Pathway = pathway_name,
    "Type of interaction" = annotation
  ) %>%
  save_table("data_S2_ccc", sheet_name = "Interactions")
