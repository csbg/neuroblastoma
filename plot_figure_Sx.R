# @DEPI rna_myeloid.rds
# @DEPI metadata.rds
# @DEPI metadata_myeloid.rds
# @DEPI dge_results_myeloid.rds

library(monocle3)
library(tidyverse)
library(latex2exp)
library(ComplexHeatmap)
library(scico)
library(RColorBrewer)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

cds_my <- readRDS("data_generated/rna_myeloid.rds")
my_metadata <- readRDS("data_generated/metadata_myeloid.rds")
nb_metadata <- readRDS("data_generated/metadata.rds")
dge <- readRDS("data_generated/dge_results_myeloid.rds")

# are there problematic genes with convergence <= -20 ?
stopifnot(
  all(dge$results_vs_C$convergence > -20),
  all(dge$results_MNA_vs_other$convergence > -20)
)



# UMAPs -------------------------------------------------------------------

plot_umap <- function(cluster_col = collcluster) {
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

plot_umap()
ggsave_default("myeloid/umap")



# Abundances --------------------------------------------------------------

plot_abundance <- function(cluster_col = collcluster) {
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

plot_abundance()
ggsave_default("myeloid/abundances_patients", width = 70, height = 50)



# Markers -----------------------------------------------------------------

plot_dots(
  logcounts(cds_my),
  c("CD14", "CD36", "HLA-DRA", "HLA-DRB1", "HLA-DRB2", "HLA-DRB3",
    "HLA-DPA1", "HLA-DPB1", "CFP", "TGFB1", "PILRA",
    "CD163", "CD74", "CD86", "FCGR3A", "C1QA") %>% rev(),
  my_metadata$collcluster
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave_default("myeloid/markers", height = 100)



# Cell types --------------------------------------------------------------

plot_celltype_heatmap <- function(cluster_col = collcluster,
                                  clusters = NULL,
                                  label = c("fine", "broad"),
                                  lump_prop = 0.01) {
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

(p <- plot_celltype_heatmap(collcluster))
ggsave_default("myeloid/cell_types", plot = p, width = 150, height = 60)



# Consistent genes --------------------------------------------------------

consistent_genes <- 
  dge$results_wide_filtered %>%
  mutate(comparison = rename_contrast(comparison)) %>% 
  filter(
    abs(logFC) > log(2),
    comparison %>% str_ends("c")
  ) %>%
  group_by(cell_type, direction, comparison) %>%
  summarise(genes = list(unique(gene))) %>%
  ungroup() %>% 
  pivot_wider(names_from = comparison, values_from = genes) %>% 
  rowwise() %>%
  mutate(
    shared_AM = Ac %>% intersect(Mc) %>% length(),
    shared_AS = Ac %>% intersect(Sc) %>% length(),
    shared_MS = Mc %>% intersect(Sc) %>% length(),
    shared_AMS = Ac %>% intersect(Mc) %>% intersect(Sc) %>% length()
  ) %>%
  ungroup() %>%
  mutate(
    across(
      starts_with("shared"),
      ~if_else(direction == "up", ., -.)
    )
  )

consistent_genes %>% 
  select(!c(Ac, Mc, Sc)) %>% 
  pivot_longer(
    starts_with("shared"),
    names_to = "shared",
    names_prefix = "shared_",
    values_to = "n"
  ) %>% 
  mutate(shared = fct_relevel(shared, "AM", "AS", "MS")) %>% 
  ggplot(aes(cell_type, n, fill = shared)) +
  geom_col() +
  geom_hline(yintercept = 0, size = BASE_LINE_SIZE) +  
  xlab("cell type") +
  ylab("number of consistently down- and upregulated genes") +
  theme_nb()
ggsave_default("myeloid/number_of_consistent_genes", width = 60, height = 60)



# GSEA results ------------------------------------------------------------

plot_gsea_dots <- function(data,
                           db,
                           top_n_positive = 5L,
                           top_n_negative = 5L,
                           max_p_adj = 0.05,
                           min_abs_NES = 1,
                           filename = "auto",
                           ...) {
  data_top_terms <-
    data %>% 
    filter(db == {{db}}, padj <= max_p_adj, abs(NES) >= min_abs_NES) %>% 
    group_by(comparison, cell_type)
  
  top_terms_pos <- 
    data_top_terms %>% 
    slice_max(n = top_n_positive, order_by = NES, with_ties = FALSE) %>%
    pull(pathway) %>%
    unique()
  
  top_terms_neg <- 
    data_top_terms %>% 
    slice_min(n = top_n_negative, order_by = NES, with_ties = FALSE) %>%
    pull(pathway) %>%
    unique()
  
  data_vis <- 
    data %>% 
    filter(db == {{db}}, pathway %in% c(top_terms_pos, top_terms_neg)) %>% 
    mutate(
      is_significant =
        padj <= max_p_adj &
        abs(NES) >= min_abs_NES,
      pathway =
        as_factor(pathway) %>%
        fct_reorder(NES * -log10(padj), sum, na.rm = TRUE),
      cell_type = factor(cell_type)
    )
  
  if (nlevels(data_vis$pathway) > 5) {
    horizontal_grid <-
      geom_hline(
        yintercept = seq(5, nlevels(data_vis$pathway), 5),
        size = 0.5,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  color_limit <- max(abs(data_vis$NES))
  
  p <- 
    ggplot(data_vis, aes(comparison, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    geom_point(data = data_vis %>% filter(is_significant), shape = 1) +
    scale_color_distiller(
      palette = "PiYG",
      direction = 1,
      limits = c(-color_limit, color_limit)
    ) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(cell_type), nrow = 1) +
    labs(
      y = "",
      color = "normalized\nenrichment\nscore",
      size = TeX("-log_{10} (p_{adj})"),
      title = db
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  if (filename == "auto") {
    filename <- str_glue("myeloid/gsea_{db}")
  }
  ggsave_default(filename, width = 400, ...)
  p
}

dge$gsea %>% 
  filter(comparison != "Mas") %>% 
  mutate(
    cell_type = fct_recode(
      cell_type,
      "class" = "classical mono",
      "nonclass" = "nonclassical mono"
    )
  ) %>% 
  plot_gsea_dots(db = "TRRUST_Transcription_Factors_2019", height = 250)

dge$gsea %>% 
  filter(comparison == "Mas") %>% 
  mutate(
    cell_type = fct_recode(
      cell_type,
      "class" = "classical mono",
      "nonclass" = "nonclassical mono"
    )
  ) %>% 
  plot_gsea_dots(db = "MSigDB_Hallmark_2020",
                 height = 150, filename = "myeloid/gsea_Mas")



# Tables ------------------------------------------------------------------

## S4 ----

gene_pathways <- 
  CellChatDB.human$interaction %>% 
  as_tibble() %>% 
  select(interaction_name, pathway_name) %>% 
  mutate(genes = str_split(interaction_name, "_")) %>% 
  unnest_longer(genes) %>% 
  distinct(pathway_name, genes) %>% 
  group_by(Gene = genes) %>% 
  summarise(
    Pathways =
      unique(pathway_name) %>%
      str_sort() %>%
      str_c(collapse = ", ")
  )

rename_columns <- function(s) {
  comparison <- str_sub(s, 7)
  case_when(
    str_starts(s, "log") ~ paste0("Log fold change (", comparison, ")"),
    str_starts(s, "p_adj") ~ paste0("Adjusted p-value (", comparison, ")"),
    TRUE ~ s
  )
}

dge$results_wide_filtered %>% 
  arrange(comparison, cell_type, desc(logFC)) %>%
  mutate(
    logFC = logFC / log(2),  # nebula returns natural log fold changes
    comparison = recode(comparison, Ac = "A vs C", Mc = "M vs C",
                        Sc = "S vs C", Mas = "M vs A+S"),
  ) %>%
  select(!c(p, frq, frq_ref, direction)) %>%
  pivot_wider(names_from = comparison, values_from = c(logFC, p_adj)) %>%
  relocate(1, 2, 3, 7, 4, 8, 5, 9, 6, 10) %>%
  rename_with(rename_columns) %>%
  left_join(gene_pathways, by = c(gene = "Gene")) %>%
  split(.$cell_type) %>%
  map(select, !cell_type) %>%
  save_table("Sxxx_dge_myeloid")




# Files for Loupe Browser -------------------------------------------------

reducedDim(cds_my, "UMAP") %>% 
  magrittr::set_colnames(c("X Coordinate", "Y Coordinate")) %>% 
  as_tibble(rownames = "Barcode") %>% 
  write_csv("~/Desktop/my_umap.csv")

my_metadata %>% 
  select(Barcode = cell, MyeloidCluster = collcluster) %>% 
  write_csv("~/Desktop/my_types.csv")





# Gene plots --------------------------------------------------------------

markers <- read_csv("metadata/myeloid_markers.csv")
markers


## Violins ----

diff_genes <-
  dge$results_wide_filtered %>% 
  filter(comparison != "Mas", gene %in% markers$gene, p_adj < 0.05) %>% 
  pull(gene) %>% 
  unique()
  {.}



make_matrix <- function(gene, cell_type) {
  barcodes <- 
    my_metadata %>% 
    filter(collcluster == {{cell_type}}) %>% 
    pull(cell)
  
  logcounts(cds_my)[gene, barcodes, drop = FALSE] %>%
    t() %>%
    as.matrix() %>%
    magrittr::set_colnames("logexp") %>%
    as_tibble(rownames = "cell") %>%
    left_join(my_metadata, by = "cell") %>%
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
    mutate(
      gene = as_factor(gene),
      cell_type = as_factor(cell_type),
      cancer_state = fct_collapse(group, C = "C", other_level = "MAS")
    )
  
  ggplot(plot_data, aes(cancer_state, logexp)) +
    geom_violin(
      aes(color = cancer_state, fill = cancer_state),
      size = BASE_LINE_SIZE,
      scale = "width",
      width = 0.8,
      show.legend = FALSE
    ) +
    stat_summary(geom = "point", fun = mean, size = .2) +
    scale_y_continuous(
      "log-normalized expression",
      limits = c(0, 1)
    ) +
    scale_fill_manual(
      values = c(GROUP_COLORS, MAS = "purple"),
      aesthetics = c("color", "fill")
    ) +
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
  genes = diff_genes,
  cell_types = c("classical mono", "nonclassical mono", "mDCs", "other")
)
ggsave_default("myeloid/genes_violin_sample", height = 400, width = 150)
ggsave_default("myeloid/genes_violin_group", height = 400, width = 150)
ggsave_default("myeloid/genes_violin_state", height = 400, width = 150)



## Pathway genes ----

selected_genes <- list(
  "TNF-alpha Signaling via NF-kB" = c("IL1B", "FOS", "NFKB1", "MYC", "IFNGR2",
                                      "CD44", "RELB", "SOD2", "MAP3K8"),
  "Interferon Gamma Response" = c("IFI44", "NFKB1", "MYD88", "HLA-A", "CD86",
                                  "HLA-DQA1", "HLA-DQB1", "TNFSF10"),
  "Myc Targets V1" = c("TRIM28", "PCNA", "CDK4", "MCM5", "HDAC2"),
  "E2F Targets" = c("DNMT1", "EZH2", "MKI67", "PCNA")
)

plot_pathway_genes <- function(db = "MSigDB_Hallmark_2020",
                               pathways = selected_genes,
                               level = c("group", "sample"),
                               cell_types = c("classical mono",
                                              "nonclassical mono",
                                              "mDCs",
                                              "other")) {
  level <- match.arg(level)
  
  row_metadata <-
    pathways %>%
    enframe("pathway", "gene") %>% 
    mutate(pathway = as_factor(pathway)) %>% 
    unnest_longer(gene) %>% 
    filter(gene %in% rownames(dge$cds))
  
  # collect barcodes per column
  col_metadata <-
    dge$metadata %>% 
    filter(
      collcluster %in% {{cell_types}},
      cell %in% colnames(dge$cds)
    ) %>% 
    mutate(
      cell_type = factor(collcluster),
    ) %>% 
    arrange(group) %>% 
    {
      if (level == "sample")
        group_by(., cell_type, group, sample)
      else
        group_by(., cell_type, group)
    } %>% 
    summarise(cells = list(cell))
  
  # make expression matrix
  mat <-
    dge$cds %>% 
    logcounts() %>% 
    magrittr::extract(row_metadata$gene, ) %>% 
    as.matrix()
  
  mat <- 
    col_metadata %>% 
    pull(cells) %>% 
    map(~mat[, .] %>% rowMeans()) %>%
    purrr::reduce(cbind)
  
  # scale within cell types
  if (level == "sample") {
    mat <- cbind(
      mat[, 1:16] %>% t() %>% scale() %>% t(),
      mat[, 17:32] %>% t() %>% scale() %>% t()
    )
  } else {
    mat <- cbind(
      mat[, 1:4] %>% t() %>% scale() %>% t(),
      mat[, 5:8] %>% t() %>% scale() %>% t(),
      mat[, 9:12] %>% t() %>% scale() %>% t(),
      mat[, 13:16] %>% t() %>% scale() %>% t()
    )
  }
  
  mat <- replace_na(mat, min(mat, na.rm = TRUE))
  
  # plot heatmap
  Heatmap(
    mat,
    name = "expression",
    col = circlize::colorRamp2(
      seq(
        min(mat),
        max(mat),
        length.out = 9
      ),
      viridisLite::cividis(9)
    ),
    border = FALSE,
    heatmap_legend_param = list(
      at = c(min(mat), max(mat)),
      labels = c("low", "high"),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    
    show_row_names = TRUE,
    row_split = row_metadata$pathway,
    row_title_rot = 0,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    show_row_dend = TRUE,
    row_gap = unit(.5, "mm"),
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    cluster_columns = FALSE,
    column_split = col_metadata$cell_type,
    show_column_names = FALSE,
    column_gap = unit(.5, "mm"),
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    width = unit(.5 * (length(cell_types) - 1) + ncol(mat) * 2, "mm"),
    height = unit(.5 * (length(pathways) - 1) + nrow(mat) * 2, "mm"),
    
    top_annotation = HeatmapAnnotation(
      group = col_metadata$group,
      col = list(
        group = GROUP_COLORS
      ),
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        group = list(
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    )
  )
}

(p <- plot_pathway_genes())
ggsave_default("myeloid/genes_pathway_heatmap",
               plot = p, width = 120, height = 180)


selected_genes_my <- 
  markers %>% 
  group_by(cell_type) %>% 
  summarise(gene = list(gene)) %>% 
  deframe()

(p <- plot_pathway_genes(pathways = selected_genes_my))
ggsave_default("myeloid/genes_my_heatmap",
               plot = p, width = 120, height = 500)



## Hexbin UMAP ----

plot_umap_hexbin <- function(...) {
  genes <- c(...)
  unknown_genes <- setdiff(genes, rownames(cds_my))
  if (length(unknown_genes) > 0) {
    info("Ignoring genes: {unknown_genes}")
    genes <- setdiff(genes, unknown_genes)
  }
  my_metadata %>% 
    left_join(
      logcounts(cds_my)[genes, , drop = FALSE] %>% 
        as.matrix() %>% 
        t() %>% 
        as_tibble(rownames = "cell") %>% 
        pivot_longer(!cell, names_to = "gene", values_to = "expression"),
      by = "cell"
    ) %>% 
    mutate(gene = factor(gene, levels = genes)) %>% 
    arrange(expression) %>%
    ggplot(aes(UMAP1, UMAP2)) +
    geom_hex(aes(weight = expression), bins = 50) +
    scale_fill_viridis_c(
      "sum of expression",
      limits = c(0, 30),
      breaks = c(0, 30),
      labels = c("0", "30 or higher"),
      oob = scales::oob_squish_any
    ) +
    coord_fixed() +
    facet_grid(vars(gene), vars(if_else(group == "C", "C", "M/A/S"))) +
    theme_nb(grid = FALSE)
}

plot_umap_hexbin("CD14", "FCGR3A", "EREG", "CD163", "CXCL2", "G0S2",
                 "IL10", "MSR1", "VEGFA", "PILRA", "TGFB1",
                 "CD86", "CD58", "CXCR4")
ggsave_default("myeloid/genes_umap_hexbin", height = 500)

plot_umap_hexbin(diff_genes)
ggsave_default("myeloid/genes_umap_hexbin_diffgenes", height = 500)
