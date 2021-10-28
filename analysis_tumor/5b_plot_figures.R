source("analysis_tumor/packages.R")
source("styling.R")



# Load data ---------------------------------------------------------------

# unaligned dataset
unaligned_cds <- readRDS("analysis_tumor/data_generated/cds_unaligned.rds")
colData(unaligned_cds)$Group <- rename_groups(colData(unaligned_cds)$group)
colData(unaligned_cds)$Patient <- rename_patients(colData(unaligned_cds)$sample)

# aligned, clustered dataset
cds <- readRDS("analysis_tumor/data_generated/cds.rds")
colData(cds)$Group <- rename_groups(colData(cds)$group)
colData(cds)$Patient <- rename_patients(colData(cds)$sample)


# GSEA results
EnrData <- readRDS("analysis_tumor/data_generated/gsea.rds")

# trajectory 
trajectory_cds <- readRDS("analysis_tumor/data_generated/cds_trajectory.rds")

# patient-cluster enrichment
Enr_Patients <- readRDS("analysis_tumor/data_generated/enrichment_patient_cluster.rds")

# literature marker enrichment
Enr_Literature <- readRDS("analysis_tumor/data_generated/enrichment_jansky.rds")

# Differentially Expressed genes Cl4
DEG <- readRDS("analysis_tumor/data_generated/deg.rds")



# Plot settings -----------------------------------------------------------

Side_quadratic_panel <- 3.6
Side_long_panel_x <- 6
Side_long_panel_y <- 7.2
Side_Heatmap_x <- 3
Side_Violin_x <- 3.8



# Figure 2a ---------------------------------------------------------------

plot_UMAP_patient <- function(cds) {
  UMAP <-
    reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1", "umap_2")) %>%
    as_tibble(rownames = "cell") %>%
    inner_join(
      colData(cds) %>%
        as_tibble() %>%
        mutate(cell = colnames(cds)),
      by = "cell"
    ) %>%
    select(1:9) %>%
    mutate(patient = rename_patients(sample))
  
  ggplot(UMAP, aes(umap_1, umap_2)) +
    geom_point(
      aes(color = patient),
      size = BASE_LINE_SIZE,
      shape = 16
    ) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_color_manual(
      "Patient",
      values = PATIENT_COLORS,
      drop = TRUE,
      guide = guide_legend(
        direction = "vertical",
        title.position = "top",
        override.aes = list(
          shape = 16,
          size = BASE_TEXT_SIZE_MM
        ),
        keywidth = 0.1,
        keyheight = 0.2,
        label.position = "left"
      )
    ) +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      legend.key.height = unit(0.3, "cm"),
      legend.box = "vertical",
      legend.box.just = "left",
      legend.direction = "vertical",
      legend.justification = "top",
      plot.margin = margin(r = 0.5, t = 0.5, unit = "mm")
    )
}

plot_UMAP_patient(unaligned_cds) 

ggsave_publication(
  "2a_unaligned_patients",
  width = Side_quadratic_panel * 1.5,
  height = Side_quadratic_panel
)



# Figure 2b, c ------------------------------------------------------------

plot_UMAP_group_patient <- function(cds, group_ID) {
  cds <- cds[, colData(cds)$group == group_ID]
  
  UMAP <-
    reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
    as_tibble(rownames = "cell") %>%
    inner_join(
      colData(cds) %>% as_tibble(rownames = "cell"),
      by = "cell"
    ) %>%
    select(1:9) %>%
    mutate(
      group = rename_groups(group),
      patient = rename_patients(sample)
    )
  
  ggplot(UMAP, aes(umap_1,umap_2)) +
    geom_point(
      size = BASE_LINE_SIZE,
      shape = 16,
      aes(color = patient
      )) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    annotate(
      "text",
      label = str_glue("Group {rename_groups(group_ID)}"),
      x = -4,
      y = 4.25,
      size = BASE_TEXT_SIZE_MM
    ) +
    scale_x_continuous(limits = c(-5, 5)) +
    scale_y_continuous(limits = c(-4.5, 4.5)) +
    scale_color_manual(values = PATIENT_COLORS, drop = TRUE) +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      plot.margin = margin(r = 0.5, t = 0.5, unit = "mm"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
}

plot_UMAP_group_patient(cds, "II")

ggsave_publication(
  "2b_UMAP_M",
  width = Side_quadratic_panel,
  height = Side_quadratic_panel
)


plot_UMAP_group_patient(cds, "IV")

ggsave_publication(
  "2c_UMAP_S",
  width = Side_quadratic_panel,
  height = Side_quadratic_panel,
  legends = FALSE
)



# Figure 2d ---------------------------------------------------------------

plot_UMAP_cluster <- function(cds) {
  UMAP <-
    reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
    as_tibble(rownames = "cell") %>%
    inner_join(
      colData(cds) %>% as_tibble(rownames = "cell"),
      by = "cell"
    ) %>%
    select(1:9 | "cluster") %>%
    mutate(
      group = rename_groups(group),
      patient = rename_patients(sample)
    )
  
  ggplot(UMAP, aes(umap_1,umap_2)) +
    geom_point(aes(color = cluster), size = .01, shape = 16) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    annotate(
      geom = "text",
      label = c("1", "2", "3", "4"),
      x = c(-2, 2.5, 1, -2.7),
      y = c(2, -1.5, 1, -2.7),
      size = BASE_TEXT_SIZE_PT
    ) +
    scale_x_continuous(limits = c(-5, 5)) +
    scale_y_continuous(limits = c(-4.5, 4.5)) +
    scale_color_manual(values = SUBCLUSTER_COLORS) +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(r = 0.5, t = 0.5, unit = "mm"),
    )
}

plot_UMAP_cluster(cds)

ggsave_publication(
  "2d_UMAP_subclusters",
  width = Side_quadratic_panel,
  height = Side_quadratic_panel
)



# Figure 2e ---------------------------------------------------------------

plot_enrich_cluster_patient <- function(data,
                                        circle_significant = FALSE,
                                        sig_p = 0.05) {
  
  data <- 
    data %>% 
    mutate(is_significant = Adjusted.P.value <= sig_p) 
  
  if (circle_significant) {
    circle_significant <- geom_point(
      data = data %>% filter(is_significant),
      shape = 1,
      stroke = BASE_LINE_SIZE * 2,
      alpha = 0.5
    )
  } else {
    circle_significant <- NULL
  }
  
  color_limit <- log2(10)
  
  ggplot(data, aes(Cluster, Patient, size = -log10(Adjusted.P.value))) +
    geom_point(aes(color = log2(Odds.Ratio))) +
    circle_significant +
    xlab("NB subcluster") +
    ylab(NULL) +
    scale_color_gsea(
      name = TeX("log_2 odds ratio"),
      limits = c(-3, 3),
      breaks = c(-2, 0, 2),
      guide = guide_colorbar(
        barheight = unit(1.2, "cm"),
        barwidth = unit(2, "mm"),
        label.position = "left",
        title.vjust = 0.2,
        order = 2
      )
    ) +
    scale_size(
      name = TeX("-log_{10} p_{adj}"),
      range = c(0, 2),
      breaks = c(0, 1, 3, 10),
      guide = guide_legend(
        order = 1,
        label.position = "left"
      )
    )  +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      aspect.ratio = 3,
      legend.box = "vertical",
      legend.box.just = "left",
      legend.direction = "vertical",
      legend.justification = "right",
      legend.title.align = 1,
      legend.title = element_text(angle = 0),
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(0.1, "mm"),
      legend.position = "left",
      legend.spacing = unit(0.5, "mm"),
      axis.text.x = element_text(),
      legend.margin = margin(0, 0.1, 0, -3, "mm"),
      plot.margin = margin(t = 0.5, r = 0.1, unit = "mm")
    )
}

Enr_Patients %>% 
  mutate(
    Patient = 
      Patient %>% 
      factor(levels = names(PATIENT_COLORS)) %>%
      fct_rev()
  ) %>%
  plot_enrich_cluster_patient()

ggsave_publication(
  "2e_patient_cluster_enrichment",
  width = Side_quadratic_panel,
  height = Side_quadratic_panel, 
)



# Figure 2f ---------------------------------------------------------------

plot_DEG_heatmap <- function(data_DEG, data_CDS, legend = FALSE) {
  genes <- 
    data_DEG$gene_short_name %>% 
    unique()
  
  matrix <-
    exprs(data_CDS)
  
  cell_group <- 
    colData(cds) %>% 
    as_tibble(rownames = "cell") %>% 
    select(cell, cluster)
  
  gene_group <- tibble(
    id = genes,
    gene_group = genes
  )  
  
  EXPR <- 
    aggregate_gene_expression(
      data_CDS, 
      gene_group,
      cell_group,
      norm_method = "log",
      scale_agg_values = FALSE
    ) %>% 
    as.matrix()
  
  EXPR_scale <- t(scale(t(EXPR)))
  
  colnames(EXPR_scale) <- colnames(EXPR)
  
  ht_opt(
    legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    DENDROGRAM_PADDING = unit(0.5, "mm"),
    TITLE_PADDING = unit(c(0.5, 0.5), "mm"),
    DIMNAME_PADDING = unit(0.5, "mm"),
    HEATMAP_LEGEND_PADDING = unit(4, "mm")
  )
  
  H <- Heatmap(
    EXPR_scale,
    col = rev(brewer.pal(8, "RdBu")),
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    show_column_dend = TRUE,
    show_column_names = TRUE,
    show_row_names = FALSE,
    #column_order = c("1", "2", "3", "4"),
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    column_names_rot = 0,
    column_dend_height = unit(3, "pt"),
    column_dend_gp = gpar(lwd = BASE_LINE_SIZE),
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    column_title = "NB subcluster",
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    column_title_side = "bottom",
    column_dend_side = "bottom",
    show_heatmap_legend = legend, 
    heatmap_legend_param = list(
      title = latex2exp::TeX("log_{10} expr_{scaled}"),
      direction = "vertical",
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      grid_height = unit(1, "mm"),
      grid_width = unit(1, "mm"),
      legend_height = unit(1.2, "cm"),
      legend_width = unit(2, "mm"),
      title_position = "topcenter"
    )
  ) %>% 
    draw(
      padding = unit(c(0.1, 4, 0.5, 0.5), "mm"),
      heatmap_legend_side = "left"
    )
}

p <- plot_DEG_heatmap(DEG, cds, legend = TRUE)

ggsave_publication(
  "2f_heatmap",
  plot = p,
  width = Side_Heatmap_x,
  height = Side_long_panel_y
)



# REVISE: Figure 2g ---------------------------------------------------------------

plot_gsea(
  data = EnrData[["Selected"]],
  min_OR = 1,
  circle_significant = FALSE
)

ggsave_publication(
  "NB_GSEA", 
  width = Side_long_panel_x, 
  height = Side_long_panel_y,
  legends = FALSE
)

ggsave_publication(
  "GSEA_guides", 
  width = Side_long_panel_x * 1.25, 
  height = Side_long_panel_y,
  legends = TRUE
)



# REVISE: Figure 2h ---------------------------------------------------------------

top_markers <- list(
  Cl1 = tibble(gene = c("STAT3", "JUN")),
  Cl2 = tibble(gene = c("E2F1", "PCNA")),
  Cl3 = tibble(gene = c("SOX11", "KIF5B"))
)

plot_violin2(cds, top_markers, 1)

ggsave_publication(
  "nb_violin_CL1",
  width = Side_Violin_x,
  height = Side_long_panel_y / 3, 
  legend = FALSE
)

plot_violin2(cds, top_markers, 2)

ggsave_publication(
  "nb_violin_CL2",
  width = Side_Violin_x,
  height = Side_long_panel_y / 3, 
  legend = FALSE
)

plot_violin2(cds, top_markers, 3)

ggsave_publication(
  "nb_violin_CL3",
  width = Side_Violin_x,
  height = Side_long_panel_y / 3, 
  legend = FALSE
)



# Figure 2i ---------------------------------------------------------------

plot_enrich_literature_markers <- function(data,
                                           circle_significant = FALSE,
                                           sig_p = 0.05) {
  data <- 
    data %>% 
    mutate(is_significant = Adjusted.P.value <= sig_p) 
  
  if (circle_significant) {
    circle_significant <- geom_point(
      data = data %>% filter(is_significant),
      shape = 1,
      stroke = BASE_LINE_SIZE * 2,
      alpha = 0.5
    )
  } else {
    circle_significant <- NULL
  }
  
  color_limit <- log2(10)
  
  ggplot(
    data, 
    aes(
      Cluster, 
      rename_adrenal_cells(Signature), 
      size = -log10(Adjusted.P.value)
    )
  ) +
    geom_point(
      aes(color = log2(Odds.Ratio)
      )
    ) +
    circle_significant +
    scale_x_discrete("NB subcluster", limits = factor(1:4)) +
    scale_y_discrete(NULL) +
    scale_color_gsea(
      name = TeX("log_2 odds ratio"),
      limits = c(-3,3),
      breaks = c(-2,0,2),
      guide = guide_colorbar(
        barheight = unit(12, "mm"),
        barwidth = unit(2, "mm"),
        label.position = "left",
        title.vjust = 0.2,
        order = 2
      )
    ) +
    scale_size(
      name = TeX("-log_{10} p_{adj}"),
      range = c(0, 2),
      breaks = c(0, 1, 3, 10),
      guide = guide_legend(
        order = 1,
        label.position = "left"
      )
    ) +
    coord_fixed() +
    expand_limits(x = 0.5) +
    theme_nb(grid = FALSE) +
    theme(
      aspect.ratio = 3,
      legend.box = "vertical",
      legend.box.just = "left",
      legend.direction = "vertical",
      legend.justification = "center",
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.spacing = unit(1, "mm"),
      legend.margin = margin(0, 1, 0, -3, "mm"),
      plot.margin = margin(0.5, 0, 0, 0, "mm")
    )
}

plot_enrich_literature_markers(Enr_Literature)

ggsave_publication(
  "2i_literature_marker_enrichment",
  width = Side_quadratic_panel * 1.3,
  height = Side_quadratic_panel,
)



# Figure 2j ---------------------------------------------------------------

plot_trajectory <- function(cds) {
  # get UMAP coordinates of trajectory nodes
  ica_space_df <-
    cds@principal_graph_aux$UMAP$dp_mst %>%
    t() %>% 
    magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
    as_tibble(rownames = "node_id")
  
  # get edges of trajectory graph and add UMAP coordinates of their nodes
  edge_df <-
    cds@principal_graph$UMAP %>%
    igraph::as_data_frame() %>%
    left_join(
      ica_space_df %>%
        rename(
          source_umap_1 = umap_1,
          source_umap_2 = umap_2
        ),
      by = c(from = "node_id")
    ) %>%
    left_join(
      ica_space_df %>%
        rename(
          target_umap_1 = umap_1,
          target_umap_2 = umap_2
        ),
      by = c(to = "node_id")
    )
  
  colData(cds) %>% 
    as_tibble(rownames = "cell") %>% 
    left_join(
      reducedDim(cds, "UMAP") %>%
        magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
        as_tibble(rownames = "cell"),
      by = "cell"
    ) %>% 
    ggplot(aes(umap_1, umap_2)) +
    geom_point(
      aes(color = cytoTRACE),
      shape = 16,
      size = BASE_LINE_SIZE / 2,
      show.legend = TRUE
    ) +
    geom_segment(
      data = edge_df,
      aes(
        x = source_umap_1, 
        y = source_umap_2,
        xend = target_umap_1, 
        yend = target_umap_2
      ),
      size = BASE_LINE_SIZE,
      color = "#525252",
      alpha = 0.75
    ) +
    scale_x_continuous("UMAP1", limits = c(-5, 5)) +
    scale_y_continuous("UMAP2", limits = c(-4.5, 4.5)) +
    scale_color_cytotrace(
      "CytoTRACE", 
      guide = guide_colorbar(
        barheight = unit(12, "mm"),
        barwidth = unit(2, "mm")
      )
    )  +
    coord_fixed() +
    theme_nb(grid = FALSE)  +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.key.height = unit(0.45,"cm"),
      legend.key.width = unit(0.45,"cm"),
      plot.margin = margin(t = 0.5, r = 0.5, unit = "mm")
    )
}

plot_trajectory(trajectory_cds)

ggsave_publication(
  "2j_trajectory",
  width = Side_quadratic_panel * 1.5,
  height = Side_quadratic_panel
)




# Supplementary figure 2 --------------------------------------------------

# unaligned UMAP colored by aligned cluster
plot_UMAP_alignedcluster(unaligned_cds)

ggsave_publication(
  "S2_UMAP_aligned_CL", 
  width = 4, 
  height = 4,
  legends = FALSE
)

ggsave_publication(
  "S2_UMAP_aligned_CL_guides", 
  width = 4, 
  height = 4,
  legends = TRUE
)

# cell cycle bargraph 

plot_CC_patient(cds, labs = FALSE)

ggsave_publication(
  "S2_CellCycle_Patients", 
  width = 4, 
  height = 4,
  legends = FALSE
)

ggsave_publication(
  "S2_CellCycle_Patients_guides", 
  width = 4,
  height = 8,
  legends = TRUE
)


# violin plot CytoTRACE ccore 
plot_Violin_cyto(cds, labs = TRUE)

ggsave_publication(
  "S2_Violin_CytoTRACE", 
  width = 4, 
  height = 4,
  legends = FALSE
)


# violin plot library size
plot_Violin_cluster(cds)

ggsave_publication(
  "S2_Violin_library", 
  width = 4, 
  height = 4,
  legends = FALSE
)


# violin plot number of features
plot_Violin_cluster(cds)

ggsave_publication(
  "S2_Violin_features", 
  width = 4, 
  height = 4,
  legends = FALSE
)

# violin plot percent mitochondrial genes
plot_Violin_cluster(cds)

ggsave_publication(
  "S2_Violin_mitopercent", 
  width = 4, 
  height = 4,
  legends = FALSE
)
