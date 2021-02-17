library(muscat)
library(scater)
library(tidyverse)
library(viridis)
library(patchwork)
library(ggpmisc)
library(fs)
library(ComplexHeatmap)
library(enrichR)


ggsave_default <- function(filename, plot = NULL, width = 297, height = 210,
                           crop = TRUE, ...) {
  if (is.null(filename))
    return()
  
  filename <- str_glue("plots/{filename}.png")
  filename %>% path_dir() %>% dir_create()
  if (is.null(plot)) {
    ggsave(filename, dpi = 300, units = "mm", limitsize = FALSE,
           width = width, height = height, ...)  
  } else {
    png(filename, res = 300, width = width, height = height, units = "mm")
    print(plot)
    dev.off()
  }
  
  if (crop) knitr::plot_crop(filename)
  
  invisible(filename)
}



# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_integrated_monocle.rds") %>% 
  logNormCounts()

nb@colData <-
  readRDS("data_generated/metadata.rds") %>%
  mutate(
    sample =
      str_c(group, sample, sep = ".") %>% 
      as_factor() %>%
      fct_relevel(str_sort) %>% 
      fct_relabel(~str_extract(.x, "\\d+_\\d+")),
    Size_Factor = colData(nb)$Size_Factor
  ) %>% 
  column_to_rownames("cell") %>% 
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)
nb

# exclude group I
# there are no cells in the NB cluster (8) in group I, so DGE analysis does not work
# nb <- nb[, nb$group != "I"]



# DGE analysis ------------------------------------------------------------

## Data preparation ----

nb <- prepSCE(
  nb,
  kid = "cluster_50",
  gid = "group",
  sid = "sample",
  drop = TRUE
)



## Approach 1: Pseudo-bulk ----

pb <- aggregateData(
  nb,
  assay = "counts",
  fun = "sum",
  by = c("cluster_id", "sample_id")
)
pb

exp_info <-
  metadata(pb)$experiment_info %>%
  column_to_rownames("sample_id")
model_mat <- model.matrix(~ 0 + group_id, data = exp_info)
colnames(model_mat) <- levels(exp_info$group_id)
model_mat


dge <- pbDS(
  pb,
  design = model_mat,
  contrast = limma::makeContrasts(
    II_vs_I = II - I,
    III_vs_I = III - I,
    IV_vs_I = IV - I,
    levels = model_mat
  ),
  min_cells = 1
)

dge_results <-
  resDS(nb, dge, frq = TRUE) %>%
  as_tibble() %>%
  mutate(
    cluster_id =
      as_factor(cluster_id) %>%
      fct_expand(levels(colData(nb)$cluster_id)) %>%
      fct_inseq()
  )


## Approach 2: Cell-level mixed model ----

# TBD
# dge_mm <- mmDS(nb)



# Evaluation of results ---------------------------------------------------

filter_dge_results <- function(data,
                               max_p_adj = 0.05,
                               min_abs_log_fc = 1,
                               min_freq = 0.1) {
  data %>%
    filter(
      I.frq > min_freq |
      (contrast == "II_vs_I" & II.frq > min_freq) |
      (contrast == "III_vs_I" & III.frq > min_freq) |
      (contrast == "IV_vs_I" & IV.frq > min_freq)
    ) %>%
    filter(
      p_adj.loc < max_p_adj,
      abs(logFC) > min_abs_log_fc
    )
}

dge_results_filtered <- 
  filter_dge_results(dge_results)



## Volcano plot ----

plot_volcano <- function(data,
                         contrast,
                         max_p_adj = 0.05,
                         min_abs_log_fc = 1,
                         min_freq = 0.1,
                         point_colors = c(
                           "frequent, DE" = "red",
                           "DE" = "black",
                           "not DE" =  "gray80"
                         ),
                         filename = NULL) {
  group_left <- str_extract(contrast, "^[^_]+")
  group_right <- str_extract(contrast, "[^_]$")
  
  data_filtered <- 
    data %>%
    filter(contrast == {{contrast}}) %>% 
    mutate(
      is_significant =
        p_adj.loc < max_p_adj &
        abs(logFC) > min_abs_log_fc,
      is_frequent =
        !!sym(str_glue("{group_left}.frq")) > min_freq |
        !!sym(str_glue("{group_right}.frq")) > min_freq,
      type = case_when(
        is_significant & is_frequent ~ "frequent, DE",
        is_significant ~ "DE",
        TRUE ~ "not DE"
      )
    ) 
  
  genes_count <- 
    data_filtered %>% 
    count(cluster_id, type) %>% 
    pivot_wider(names_from = "type", values_from = "n")
  
  p <-
    ggplot(data_filtered, aes(logFC, -log10(p_adj.loc))) +
    geom_point(
      aes(color = type),
      alpha = .5,
      size = 1
    ) +
    geom_hline(yintercept = -log10(max_p_adj), size = .25) +
    geom_vline(xintercept = c(-min_abs_log_fc, min_abs_log_fc), size = .25) +
    geom_text_npc(
      data = genes_count,
      aes(label = `frequent, DE`),
      npcx = 0.05,
      npcy = 0.9,
      size = 3,
      color = point_colors["frequent, DE"],
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = DE),
      npcx = 0.05,
      npcy = 0.8,
      size = 3,
      color = point_colors["DE"],
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = `not DE`),
      npcx = 0.05,
      npcy = 0.7,
      size = 3,
      color = point_colors["not DE"],
      na.rm = TRUE
    ) +
    scale_color_manual(values = point_colors) +
    facet_wrap(vars(cluster_id), ncol = 6, drop = FALSE) +
    labs(
      title = str_glue("DGE using contrast {contrast}"),
      caption = str_glue(
        "gene expression frequency > {min_freq}, ",
        "adjusted p value < {max_p_adj}, ",
        "absolute log fold change > {min_abs_log_fc}"
      )
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  ggsave_default(filename)
  p
}

plot_volcano(dge_results, "II_vs_I", filename = "dge/volcano_II")
plot_volcano(dge_results, "III_vs_I", filename = "dge/volcano_III")
plot_volcano(dge_results, "IV_vs_I", filename = "dge/volcano_IV")



## Upset plot ----

plot_upset <- function(data, contrast, filename = NULL) {
  p <-
    data %>%
    filter(contrast == {{contrast}}) %>%
    select(cluster_id, gene) %>%
    chop(gene) %>%
    deframe() %>%
    make_comb_mat() %>% 
    UpSet(
      comb_order = order(comb_size(.), decreasing = TRUE),
      comb_col = viridis(nrow(.))[comb_degree(.)]
    )
  
  ggsave_default(filename, plot = p, width = 297, height = 100)
  p
}

plot_upset(dge_results_filtered, "II_vs_I", filename = "dge/upset_II")
plot_upset(dge_results_filtered, "III_vs_I", filename = "dge/upset_III")
plot_upset(dge_results_filtered, "IV_vs_I", filename = "dge/upset_IV")


## Violin plot ----

plot_violin <- function(data, contrast, cluster, top_n = 12, filename = NULL) {
  group_left <- str_extract(contrast, "^[^_]+")
  group_right <- str_extract(contrast, "[^_]$")
  
  top_genes <- 
    dge_results_filtered %>%
    filter(contrast == {{contrast}}, cluster_id == {{cluster}}) %>% 
    slice_max(n = top_n, logFC, with_ties = FALSE) %>% 
    pull(gene)
  
  p <-
    plotExpression(
      nb[, nb$cluster_id == cluster &
           nb$group_id %in% c(group_left, group_right)],
      top_genes,
      x = "sample_id",
      colour_by = "group_id",
      ncol = 4,
      point_size = 1
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave_default(filename)
  p
}

plot_violin(dge_results_filtered, "II_vs_I", "1", filename = "dge/violin_II_1")
plot_violin(dge_results_filtered, "II_vs_I", "13", filename = "dge/violin_II_13")
plot_violin(dge_results_filtered, "IV_vs_I", "13", filename = "dge/violin_IV_13")



## Pseudo-bulk heatmap ----

pb_log <- aggregateData(
  nb,
  assay = "logcounts",
  fun = "mean",
  by = c("cluster_id", "sample_id")
)

plot_heatmap <- function(contrast,
                         top_n = 5,
                         clusters = NULL,
                         filename = NULL) {
  group_left <- str_extract(contrast, "^[^_]+")
  group_right <- str_extract(contrast, "[^_]$")
  
  top_genes <- 
    dge_results_filtered %>% 
    filter(contrast == {{contrast}}) %>% 
    {
      if (!is.null(clusters))
        filter(., cluster_id %in% {{clusters}})
      else
        .
    } %>% 
    group_by(cluster_id) %>% 
    slice_max(logFC, n = top_n, with_ties = FALSE)
  
  sample_data <- 
    colData(pb_log) %>%
    as_tibble(rownames = "sample") %>%
    filter(group_id %in% c(group_left, group_right))
  
  mat <- 
    top_genes %>% 
    group_map(~assays(pb_log)[[.y$cluster_id]][.x$gene, sample_data$sample]) %>%
    reduce(rbind) %>% 
    t() %>% 
    scale() %>% 
    t()
  
  p <-
    Heatmap(
      mat,
      col = viridis::viridis(10),
      
      cluster_rows = FALSE,
      row_split = top_genes$cluster_id,
      show_row_names = FALSE,
      
      cluster_columns = FALSE,
      column_split = sample_data$group_id,
      
      left_annotation = rowAnnotation(
        logFC = top_genes$logFC,
        col = list(
            logFC = circlize::colorRamp2(
            breaks = c(min(top_genes$logFC), 0, max(top_genes$logFC)),
            colors = c("blue", "white", "red")
          )
        )
      ),
      
      heatmap_legend_param = list(title = "scaled\nmean\nlogcounts")
    ) %>% 
    draw(row_title = "cluster", column_title = "group")
  
  ggsave_default(filename, plot = p)
  p
}

plot_heatmap("II_vs_I", top_n = Inf, filename = "dge/heatmap_II")
plot_heatmap("III_vs_I", top_n = Inf, filename = "dge/heatmap_III")
plot_heatmap("IV_vs_I", top_n = Inf, filename = "dge/heatmap_IV")



## Enrichment analysis ----

get_enrichr_results <- function(data, contrast, dbs) {
  map_dfr(
    set_names(unique(data$cluster_id)),
    function(cluster) {
      message("Enriching cluster ", cluster)
      top_genes <- 
        data %>% 
        filter(
          cluster_id == {{cluster}},
          contrast == {{contrast}},
          logFC > 0
        ) %>% 
        pull(gene)
      message("  genes: ", str_c(top_genes, collapse = ", "))
      
      if (length(top_genes) > 0) {
        res <-
          enrichr(top_genes, dbs) %>% 
          bind_rows(.id = "db")
        if (nrow(res) > 0)
          res
      }
        
    },
    .id = "cluster"
  ) %>% 
    as_tibble()
}

enrichr_dbs <- c(
  "GO_Biological_Process_2018",
  "GO_Cellular_Component_2018",
  "GO_Molecular_Function_2018",
  "KEGG_2019_Human",
  "WikiPathways_2019_Human"
)

enrichr_II <- get_enrichr_results(dge_results_filtered, "II_vs_I", enrichr_dbs)



plot_enrichr <- function(data, n_bars = 10, filename = NULL) {
  data <- 
    data %>%
    mutate(
      n_genes = str_extract(Overlap, "\\d+") %>% as.integer(),
      Term = str_sub(Term, end = 40)
    ) %>% 
    group_by(cluster, db) %>% 
    slice_max(Combined.Score, n = n_bars, with_ties = FALSE)
  
  # max_x <- max(data$Combined.Score)
  n_row <- n_distinct(data$db)
  clusters <- unique(data$cluster) %>% str_sort(numeric = TRUE)
  dbs <- unique(data$db)
  
  p <-
    list(
      cluster = clusters,
      db = dbs
    ) %>%
    cross_df() %>%
    mutate(
      show_cluster = db == dbs[1],
      show_db = cluster == clusters[1]
    ) %>% 
    pmap(
      function(cluster, db, show_db, show_cluster) {
        data <-
          data %>%
          filter(cluster == {{cluster}}, db == {{db}})
        if (nrow(data) == 0L)
          return(textGrob("no data", gp = gpar(fontsize = 20)))
  
        data %>%
          mutate(Term = as_factor(Term) %>% fct_reorder(Combined.Score)) %>%
          ggplot(aes(Term, Combined.Score)) +
          geom_col(aes(fill = Combined.Score), show.legend = FALSE) +
          xlab(if (show_db) db else "") +
          ylab(NULL) +
          coord_flip() +
          ggtitle(if (show_cluster) paste("cluster", cluster) else "") +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title.y = element_text(hjust = 0.5, face = "bold"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
          )
      }
    ) %>%
    wrap_plots(nrow = n_row)
  ggsave_default(filename, width = 840, height = 420)
  p
}

plot_enrichr(enrichr_II, filename = "dge/enrich_II")  





"ELAVL4, IFI44L, GSTM3, MEX3A, SOX11, MYCN, NRXN1, KIF5C, ZNF804A, RBMS3, GAP43, RBP1, NSG1, UCHL1, PHOX2B, LINC00682, MAB21L2, CPE, HAND2, HAND2-AS1, ISL1, MAP1B, TUBB2A, TUBB2B, FSCN1, ICA1, PCOLCE, PCSK1N, SYP, PAGE2, TCEAL2, BEX1, TCEAL7, PNMA5, STMN4, CLU, STMN2, NFIB, CNTFR, GADD45G, TMOD1, TEAD1, GAL, CCND1, CTTN, AP002387.2, PHOX2A, NCAM1, FXYD6, ANK3, H2AFY2, KIF21A, SYT1, MAB21L1, RTN1, CHGA, MEG3, CKB, CDK5R1, MAPT, TBX2, RAC3, SEC11C, CHGB, ID1, EEF1A2, ATCAY, ELAVL3, TMEM59L, TTC9B, PLD3, IGLC2, RBFOX2" %>% 
  str_split(", ") %>% 
  magrittr::extract2(1) %>% 
  cat(sep = "\n")



# Unused ------------------------------------------------------------------

# cross table
table(colData(nb)$cluster_id, colData(nb)$group_id) %>% sum()

# UMAP with expression colored
top_12 <- 
  dge_results_filtered %>%
  filter(contrast == "II_vs_I", cluster_id == "9") %>% 
  slice_max(n = 4, logFC, with_ties = FALSE) %>%
  pull(gene) %>%
  {.}

top_12 %>% 
  map(
    ~plotReducedDim(
      nb,
      "UMAP",
      colour_by = .x,
      point_size = 0.1,
      point_alpha = 1,
      add_legend = FALSE
    ) +
      coord_fixed()
  ) %>% 
  wrap_plots()
