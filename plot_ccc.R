library(CellChat)
library(tidyverse)
library(ComplexHeatmap)
library(latex2exp)
source("styling.R")



# Load data ---------------------------------------------------------------

cellchat <- readRDS("data_generated/ccc_cellchat_object.rds") 
sig_data <- readRDS("data_generated/ccc_signaling_data.rds")



# Figure 3a ---------------------------------------------------------------

# a modified version of CellChat::netVisual_heatmap()

plot_n_interactions <- function(object, 
                                comparison = c(1, 2), 
                                measure = "count", 
                                signaling = NULL, 
                                slot.name = c("netP", "net"), 
                                color.use = NULL, 
                                color.heatmap = "Reds", 
                                title.name = NULL, 
                                width = NULL, 
                                height = NULL, 
                                cluster.rows = F, 
                                cluster.cols = F, 
                                sources.use = NULL, 
                                targets.use = NULL, 
                                remove.isolate = FALSE, 
                                row.show = NULL, 
                                col.show = NULL) {
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    }
    else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  }
  else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[, , signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    }
    else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      }
      else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), 
                                   c(color.heatmap[1], "#f7f7f7", 
                                     color.heatmap[2]))
    colorbar.break <- 
      c(
        round(
          min(mat, na.rm = T), 
          digits = nchar(
            sub(
              ".*\\.(0*).*","\\1", min(mat, na.rm = T)
            )
          ) + 1
        ), 
        0, 
        round(
          max(mat, na.rm = T), 
          digits = nchar(
            sub(
              ".*\\.(0*).*", "\\1", max(mat, na.rm = T)
            )
          ) + 1
        )
      )                                                                                                                                                      
  }
  else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = 
        (grDevices::colorRampPalette(
          (RColorBrewer::brewer.pal(n = 9, name = color.heatmap))
        )
        )(100)
    }
    colorbar.break <- c(
      round(
        min(mat, na.rm = T), 
        digits = nchar(sub(".*\\.(0*).*","\\1", min(mat, na.rm = T))) + 1
      ), 
      round(
        max(mat, na.rm = T), 
        digits = nchar(sub(".*\\.(0*).*", "\\1", max(mat, na.rm = T))) + 1
      )
    )
  }
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- 
    HeatmapAnnotation(
      df = df, 
      col = list(group = color.use),
      which = "column", 
      show_legend = FALSE, 
      show_annotation_name = FALSE, 
      simple_anno_size = grid::unit(2, "mm")
    )
  
  row_annotation <- 
    HeatmapAnnotation(
      df = df, 
      col = list(group = color.use), 
      which = "row",
      show_legend = FALSE, 
      show_annotation_name = FALSE, 
      simple_anno_size = grid::unit(2, "mm")
    )
  
  ha1 = 
    rowAnnotation(
      StrengthText = anno_text(
        colSums(abs(mat)), 
        gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      ), 
      Strength = anno_barplot(
        rowSums(abs(mat)),
        border = FALSE, 
        gp = gpar(fill = color.use, col = color.use),
        axis = FALSE
      ),
      simple_anno_size_adjust = T,
      show_annotation_name = FALSE
    )
  
  ha2 = 
    HeatmapAnnotation(
      Strength = anno_barplot(
        colSums(abs(mat)), 
        border = FALSE, 
        gp = gpar(fill = color.use, col = color.use),
        axis = FALSE
      ), 
      StrengthText = anno_text(
        colSums(abs(mat)), 
        gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      ),
      show_annotation_name = FALSE
    )
  
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
  
  Heatmap(
    mat,
    col = color.heatmap.use,
    na_col = "white", 
    name = legend.name, 
    bottom_annotation = col_annotation, 
    left_annotation = row_annotation, 
    top_annotation = ha2, 
    right_annotation = ha1, 
    cluster_rows = cluster.rows, 
    cluster_columns = cluster.rows, 
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    column_title = "Targets (Receptor)", row_title = "Sources (Ligand)", 
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
    column_names_rot = 90, row_title_rot = 90,  row_names_rot = 0, 
    column_title_side = "bottom",
    row_names_side = "left",
    width = unit(20, "mm"),
    height = unit(20, "mm"),
    heatmap_legend_param = 
      list(
        title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
        title_position = "leftcenter-rot", 
        border = NA, 
        legend_height = unit(20, "mm"), 
        labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
        grid_width = unit(2, "mm")
      )
  )
}

p1 <- plot_n_interactions(
  cellchat, 
  color.use = CELL_TYPE_COLORS[c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")],
)
p1
ggsave_publication("3a_n_interactions", plot = p1, width = 6, height = 5)



# Figure 3b ---------------------------------------------------------------

# cumulative outgoing score for NB barplot

plot_contribution_celltype <- function(data, pathways, celltype,
                                       signif = 0.05) {
  df <- 
    data %>%
    filter(
      pathway_name %in% pathways, 
      source == celltype, 
      target != celltype,
      pval <= signif
    ) %>%
    group_by(interaction = interaction_name_2) %>%
    summarise(
      prob.norm = sum(prob.norm),
      count = n()
    ) %>%
    ungroup() %>%
    mutate(prob.avg = prob.norm / count) 
  
  label <- 
    tibble(
      interaction = data$interaction_name_2, 
      pathway_name = data$pathway_name
    ) %>%
    filter(interaction %in% df$interaction) %>%
    unique()
  
  df <-
    df %>%
    left_join(label, by = "interaction") %>%
    arrange(prob.norm) %>%
    mutate(interaction = factor(interaction, levels = interaction))
  
  ggplot(df, aes(prob.norm, interaction)) +
    geom_col(fill = "#23cfc3", width = 0.85) +
    geom_text(
      aes(x = 2.68, y = interaction, label = pathway_name),
      size = BASE_TEXT_SIZE_MM
    ) +
    scale_x_continuous(breaks = c(1, 3, 5)) +
    xlab("sum(outgoing communication score)") +
    ylab("Ligand-Receptor pairs") +
    theme_nb(grid = FALSE)
}

p2 <- plot_contribution_celltype(
  sig_data, 
  cellchat@netP$pathways,
  celltype = "NB",
  signif = 0.05
)
p2

ggsave_publication("3b_outgoing_comm_score", width = 5.5, height = 5)



# Figure 3c ---------------------------------------------------------------

# LRs dotplot  

plot_selected_dots <- function(data, pathways, source_celltype) {
  # exclude all interactions not within pathways vector
  # and order pathways by input vector
  data <- 
    data %>% 
    filter(
      pathway_name %in% pathways, 
      source == source_celltype
    ) %>%
    mutate(pathway_name = factor(pathway_name, levels = pathways)) 
  
  p <- 
    ggplot(
      data,
      aes(
        target,
        interaction_name_2,
        color = prob.norm,
        size = pmin(5, -log10(pval))
      )
    ) +
    geom_point() + 
    scale_radius(
      TeX("-log_{10} p-Value_{cap.}"),
      range = c(1,5),
      guide = "none"
    ) +
    scale_color_distiller(
      TeX("relative communication score"),
      type = "seq",
      palette = "Reds", 
      direction = 1
    ) +
    xlab("Celltypes: Source -> Target") +
    ylab("Ligand-Receptor pairs & Pathways") +
    facet_grid(
      vars(pathway_name), 
      vars(source),
      scales = "free",
      space = "free",
      switch = "both"
    ) + 
    theme_nb(grid = TRUE) + 
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
      panel.spacing = unit(0, "lines"),
      
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
        face = "bold",
        margin = margin(t = 1, b = 1, r = 4, l = 4, unit = "mm"),
        angle = 0
      ),
      strip.text.y.left = element_text(
        size = BASE_TEXT_SIZE_PT, 
        face = "bold",
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
      legend.title.align = 0.5,
    )
}

p3 <- plot_selected_dots(
  sig_data, 
  pathways = c("MK", "MIF", "COLLAGEN", "PTN", "APP", "ALCAM", "THY1"), 
  source_celltype = "NB"
)
p3

ggsave_publication("3c_dots", height = 5, width = 7.5)



# Figure 3d, e ------------------------------------------------------------

# network centrality plots for MIF & MK
# a modified version of CellChat::netAnalysis_signalingRole_network()

plot_centrality <- function(object,
                            signaling,
                            slot.name = "netP", 
                            measure = c("outdeg", "indeg", "flowbet", "info"), 
                            measure.name = c("Sender", "Receiver",
                                             "Mediator", "Influencer"), 
                            color.use = NULL,
                            color.heatmap = "YlOrRd",
                            font.size = BASE_TEXT_SIZE_PT,
                            font.size.title = BASE_TEXT_SIZE_PT, 
                            cluster.rows = FALSE,
                            cluster.cols = FALSE) {
    
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality`!")
  }
  centr <- slot(object, slot.name)$centr[signaling]
  for (i in 1:length(centr)) {
    centr0 <- centr[[i]]
    mat <- matrix(unlist(centr0), ncol = length(centr0), 
                  byrow = FALSE)
    mat <- t(mat)
    rownames(mat) <- names(centr0)
    colnames(mat) <- names(centr0$outdeg)
    if (!is.null(measure)) {
      mat <- mat[measure, ]
      if (!is.null(measure.name)) {
        rownames(mat) <- measure.name
      }
    }
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    if (is.null(color.use)) {
      color.use <- scPalette(length(colnames(mat)))
    }
    color.heatmap.use = 
      (grDevices::colorRampPalette(
        (RColorBrewer::brewer.pal(n = 9, name = color.heatmap))
      )
      )(100)
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    col_annotation <- 
      HeatmapAnnotation(
        df = df, 
        col = list(group = cell.cols.assigned), 
        which = "column", 
        show_legend = FALSE, 
        show_annotation_name = FALSE, 
        simple_anno_size = grid::unit(0.2, "cm")
      )
    
    ht1 <- Heatmap(
      mat, 
      col = color.heatmap.use,
      na_col = "white", 
      name = "Importance",
      bottom_annotation = col_annotation, 
      cluster_rows = cluster.rows,
      cluster_columns = cluster.rows, 
      row_names_side = "left",
      row_names_rot = 0, 
      row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      column_title = paste0(names(centr[i])," signaling pathway network"), 
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
      width = unit(40, "mm"),
      height = unit(20, "mm"), 
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT, fontface = "plain"), 
        title_position = "leftcenter-rot", 
        border = NA, 
        at = c(
          round(min(mat, na.rm = TRUE), digits = 1), 
          round(max(mat, na.rm = TRUE), digits = 1)
        ), 
        legend_height = unit(20, "mm"), 
        labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT), 
        grid_width = unit(2, "mm")
      )
    )
  }
  ht1
}

p4 <- plot_centrality(
  cellchat,
  signaling = "MIF", 
  color.use = CELL_TYPE_COLORS[c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")]
)
p4
ggsave_publication("3d_importance_MIF", plot = p4, width = 6, height = 3.5)

p5 <- plot_centrality(
  cellchat,
  signaling = "MK", 
  color.use = CELL_TYPE_COLORS[c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")]
)
p5
ggsave_publication("3e_importance_MK", plot = p5, width = 6, height = 3.5)



# Figure 3f, g ------------------------------------------------------------

# violin plots for MK & MIF

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
ggsave_publication("3f_violins_MIF", width = 8.5, height = 9)

plot_violin(
  genes = c("MDK", "LRP1", "NCL", "ALK"),
  cell_types = c("NB", "M")
)
ggsave_publication("3g_violins_MK", width = 8.5, height = 9)

