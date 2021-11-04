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



# Figure 2g ---------------------------------------------------------------

plot_gsea <- function(data,
                      circle_significant = FALSE,
                      max_p_adj = 0.05,
                      min_OR = 5,
                      min_overlap_size = 2) {
  data <- 
    data %>%
    separate(
      Overlap,
      into = c("overlap_size", "geneset_size"),
      convert = TRUE
    )
  
  top_terms <-
    data %>% 
    filter(
      Adjusted.P.value <= max_p_adj,
      Odds.Ratio >= min_OR,
      overlap_size >= min_overlap_size
    ) %>% 
    pull(Term) %>%
    unique()
  
  data_vis <- 
    data %>% 
    filter(Term %in% top_terms) %>% 
    mutate(
      is_significant =
        Adjusted.P.value <= max_p_adj &
        Odds.Ratio >= min_OR &
        overlap_size >= min_overlap_size,
      Term = as_factor(Term) 
    ) %>%
    mutate(
      Term =
        Term %>% 
        fct_reorder(log2(Odds.Ratio), max) %>%
        fct_reorder(cluster, max)
    )
  
  if (nlevels(data_vis$Term) > 5) {
    horizontal_grid <-
      geom_hline(
        yintercept = seq(5, nlevels(data_vis$Term), 5),
        size = BASE_LINE_SIZE,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  
  color_limit <- log2(64)
  
  if (circle_significant) {
    circle_significant <- geom_point(
      data = data_vis %>% filter(is_significant),
      shape = 1,
      stroke = BASE_LINE_SIZE * 2,
      alpha = 0.5
    )
  } else {
    circle_significant <- NULL
  }
  
  p <- 
    ggplot(data_vis, aes(cluster, Term, size = -log10(Adjusted.P.value))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = log2(Odds.Ratio))) +
    circle_significant +
    xlab("NB subcluster") +
    ylab(NULL) +
    scale_color_gsea_2(
      name = TeX("log_2 odds ratio"),
      breaks = c(0, 2, 4, 6),
      limit = 
        c(0, color_limit),
      guide = guide_colorbar(
        barheight = unit(1.2, "cm"),
        barwidth = unit(2, "mm"),
        title.vjust = 0.1,
        label.position = "left",
        order = 2
      )
    ) +
    scale_size(
      name = TeX("-log_{10} p_{adj}"),
      range = c(0, 2),
      breaks = c(0.1, 1.3, 3, 10),
      guide = guide_legend(
        order = 1,
        label.position = "left"
      )
    )  +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      legend.box = "vertical",
      legend.box.just = "left",
      legend.direction = "vertical",
      legend.justification = "left",
      legend.title = element_text(size = BASE_TEXT_SIZE_PT, angle = 0),
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(0.5, "mm"),
      legend.position = "left",
      legend.spacing = unit(2, "mm"),
      legend.margin = margin(0, 0, 0, -3, "mm"),
      panel.spacing = unit(-.1, "pt"),
      plot.margin = margin(t = .5, r = 0.1, unit = "mm")
    )
  
  p
}

selected_terms <- c(
  # "TNF-alpha signalling via NF-kB",
  # "dopaminergic neuron differentiation",
  # "IL2 effects mediated by PI3K",
  "IL3-mediated signaling events",
  "IL6-mediated signaling events",
  "Trk receptor signaling mediated by PI3K and PLC-gamma",
  "PDGFR-beta signaling pathway",
  "ERG-JUN-FOS DNA-protein complex",
  "JUND-FOSB-SMAD3-SMAD4 complex",
  # "CALM1-FKBP38-BCL2 complex",
  
  "Cell cycle",
  # "DNA replication",
  "DNA repair",
  "PCNA-DNA polymerase delta complex",
  "nucleoside metabolic process",
  "DNA synthesome complex (17 subunits)",
  "6S methyltransferase and RG-containing Sm proteins complex",
  "BRD4-RFC complex",
  "histone mRNA metabolic process",
  # "DNA synthesome complex (17 subunits)",
  
  "intermediate filament bundle assembly",
  "negative regulation of glial cell proliferation",
  "peripheral nervous system neuron development",
  "noradrenergic neuron differentiation",
  "myelin maintenance",
  "BRAF-RAF1-14-3-3 complex",
  "Kinase maturation complex 1",
  "IKBKB-CDC37-KIAA1967-HSP90AB1-HSP90AA1 complex",
  # "ITGA3-ITB1-BSG complex",
  
  "Regulation of Telomerase",
  "p53-SP1 complex",
  "CAND1-CUL4A-RBX1 complex",
  "Er-alpha-p53-hdm2 complex",
  "NUMB-TP53-MDM2 complex",
  
  # new
  "Epithelial Mesenchymal Transition",
  "glial cell differentiation"
)

EnrData %>% 
  mutate(Term = rename_enrichment_terms(Term)) %>%
  filter(
    db %in% c(
      "NCI-Nature_2016",
      "CORUM",
      "GO_Biological_Process_2018",
      "KEGG_2019_Human",
      "MSigDB_Hallmark_2020"
    ),
    Term %in% selected_terms
  ) %>% 
  plot_gsea(
    min_OR = 1,
    circle_significant = FALSE
  )

ggsave_publication(
  "2g_enrichr", 
  width = Side_long_panel_x * 1.5, 
  height = Side_long_panel_y,
)



# Figure 2h ---------------------------------------------------------------

plot_violin2 <- function(data_cds, genes) {
  data_cds <- scuttle::logNormCounts(
    data_cds,
    size.factors = colData(data_cds)$Size_Factor,
  )
  
  data <- 
    logcounts(data_cds)[genes, , drop = FALSE] %>%
    t() %>%
    as.matrix() %>%
    as_tibble(rownames = "cell") %>%
    pivot_longer(
      !cell, 
      names_to = "gene", 
      values_to = "logexp"
    ) %>%
    left_join(
      data_cds@colData %>%
        as_tibble(rownames = "cell") %>% 
        select(cell, group, sample, cluster, Size_Factor),
      by = "cell"
    ) %>% 
    group_by(gene) %>% 
    mutate(logexp = logexp / max(logexp)) %>% 
    ungroup()
  
  ggplot(data, aes(cluster, logexp)) +
    geom_violin(
      aes(fill = cluster, color = NULL),
      scale = "width",
      size = BASE_LINE_SIZE,
      alpha = 0.65,
      trim = TRUE
    ) +
    geom_boxplot(
      aes(fill = NA, size = 1),
      notch = TRUE,
      size = BASE_LINE_SIZE,
      outlier.shape = NA,
      width = 1
    ) +
    xlab("NB subcluster") +
    ylab(NULL) +
    scale_y_continuous(
      name = TeX("log_2 expression_{norm.}"),
      limits = c(0, 1),
      breaks = c(0, 1),
      expand = expansion(add = .025)
    ) +
    scale_fill_manual(values = SUBCLUSTER_COLORS) +
    coord_flip() +
    facet_grid(~gene,scales = "free") +
    theme_nb(grid = FALSE) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
}


plot_violin2(cds, c("STAT3", "JUN"))

ggsave_publication(
  "2h_1_violins",
  width = Side_Violin_x,
  height = Side_long_panel_y / 3, 
)

plot_violin2(cds, c("E2F1", "PCNA"))

ggsave_publication(
  "2h_2_violins",
  width = Side_Violin_x,
  height = Side_long_panel_y / 3
)

plot_violin2(cds, c("SOX11", "KIF5B"))

ggsave_publication(
  "2h_3_violins",
  width = Side_Violin_x,
  height = Side_long_panel_y / 3,
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



# Figure S2b --------------------------------------------------------------

plot_CC_patient <- function(cds) {
  meta <-
    colData(cds) %>% 
    as_tibble() %>%
    select(group, phase, sample) %>%
    mutate(
      patient = as_factor(sample) %>% rename_patients(),
      group = rename_groups(group)
    )
  
  ggplot(meta, aes(patient, fill = patient)) +
    geom_bar() +
    facet_grid(vars(phase), scales = "free", space = "free_x") +
    xlab(NULL) +
    ylab(NULL) +
    scale_fill_manual(
      "Patients",
      values = PATIENT_COLORS,
      guide = guide_legend(
        title.position = "top",
        override.aes = list(shape = 16, size = BASE_LINE_SIZE),
        keywidth = 0.1,
        keyheight = 0.2,
        label.position = "top"
      )
    ) +
    theme_nb(grid = FALSE) +
    theme(
      axis.text.y = element_blank(),
      legend.position = "none",
      panel.spacing = unit(-.5, "pt"),
      plot.margin = margin(t = .5, r = 0.5, l = -0.5, unit = "mm")
    )
}

plot_CC_patient(cds)

ggsave_publication(
  "S2b_cell_cycle_patients", 
  width = 4, 
  height = 4
)



# Figure S2c --------------------------------------------------------------

colData(cds) %>% 
  as_tibble() %>% 
  ggplot(aes(cluster, fill = cluster)) +
  geom_bar(show.legend = FALSE) +
  xlab(NULL) +
  ylab("Number of Cells") +
  scale_color_manual(values = SUBCLUSTER_COLORS) +
  theme_nb(grid = FALSE)

ggsave_publication(
  "S2c_cells_per_cluster", 
  width = 4, 
  height = 4,
)



# Figure S2d --------------------------------------------------------------

plot_UMAP_cluster_aligned <- function(cds) {
  UMAP <-
    reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
    as_tibble(rownames = "cell") %>%
    inner_join(
      colData(cds) %>%
        as_tibble(rownames = "cell"),
      by = "cell"
    ) %>%
    select(1:9 | "aligned_cluster") %>%
    mutate(
      group = rename_groups(group),
      patient = rename_patients(sample)
    )
  
  ggplot(UMAP, aes(umap_1, umap_2)) +
    geom_point(aes(color = aligned_cluster), size = .01, shape = 16) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_color_manual(
      values = SUBCLUSTER_COLORS,
      guide = guide_legend(
        override.aes = list(shape = 16, size = 1)
      )
    ) +
    theme_nb(grid = FALSE) +
    theme(
      legend.position = "none",
      plot.margin = margin(r = 0.5, t = 0.5, unit = "mm"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
}

plot_UMAP_cluster_aligned(unaligned_cds)

ggsave_publication(
  "S2d_UMAP_unaligned", 
  width = 4, 
  height = 4,
)



# Figure S2e–h ------------------------------------------------------------

plot_violin_cluster <- function(cds, var, y_lab) {
  meta <-
    colData(cds) %>% 
    as_tibble() %>%
    select(cluster, var = {{var}}) 
  
  ggplot(meta, aes(cluster, var)) +
    geom_violin(
      aes(fill = cluster),
      scale = "width"
    ) +
    geom_boxplot(
      aes(fill = cluster),
      notch = TRUE,
      alpha = 0.1,
      size = 0.5,
      varwidth = TRUE,
    ) +
    xlab("NB subcluster") +
    ylab(y_lab) +
    scale_fill_manual(values = SUBCLUSTER_COLORS) +
    theme_nb(grid = FALSE) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = .5, r = 0.5, l = 0.1, unit = "mm")
    )
}

plot_violin_cluster(cds, library_size, "Library Size")

ggsave_publication(
  "S2e_violin_library", 
  width = 4, 
  height = 4
)


plot_violin_cluster(cds, n_features, "Number of Features")

ggsave_publication(
  "S2f_violin_nfeat", 
  width = 4, 
  height = 4
)


plot_violin_cluster(cds, percent_mito, "% Mitochondrial Reads")

ggsave_publication(
  "S2g_violin_percentmt", 
  width = 4, 
  height = 4
)

plot_violin_cluster(trajectory_cds, cytoTRACE, "CytoTRACE score")

ggsave_publication(
  "S2h_violin_cytotrace", 
  width = 4, 
  height = 4
)
