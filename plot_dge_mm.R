# Plot mixed model DGE results.
#
# @DEPI dge_mm_results_[group]_[cluster].rds
# @DEPI dge_mm_results.RData

library(monocle3)
library(scater)
library(muscat)
library(msigdbr)
library(viridis)
library(patchwork)
library(ggpmisc)
library(ComplexHeatmap)
library(latex2exp)
library(scico)
library(openxlsx)
library(ggrepel)
library(ggbeeswarm)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

load("data_generated/dge_mm_results.RData")
load("data_generated/dge_pb_results.RData")

dge_mm_results <-
  list(
    cell_type = c("T", "NK", "B", "M", "SC", "pDC", "E"),
    group = c("II", "III", "IV")
  ) %>% 
  cross_df() %>%
  pmap_dfr(
    function(group, cell_type) {
      str_glue("data_generated/dge_mm_results_{group}_{cell_type}.rds") %>% 
        readRDS() %>% 
        magrittr::extract2(1) %>% 
        as_tibble() %>% 
        mutate(contrast = paste0(group, "_vs_I"))
    }
  )



# Filter results ----------------------------------------------------------

#' Filter DGE results and add a column "direction".
#'
#' @param data Mixed model DGE results.
#' @param contrast_frq Named list of character vectors indicating which groups
#'   should be included for frequency filtering in each contrast.
#' @param max_p Maximum p value.
#' @param max_p_adj Maximum adjusted p value. Note: The FDR is calculated by
#'   muscat/edgeR via p.adjust(method = "BH").
#' @param min_abs_log_fc Minimum absolute log fold change.
#' @param min_freq Minimum gene expression frequency (fracions of cells that
#'   expresses a given gene.)
#' @param remove_ribosomal If true, remove ribosomal proteins.
#'
#' @return A data frame.
filter_dge_mm_results <- function(data,
                                  contrast_frq,
                                  max_p = Inf,
                                  max_p_adj = 0.05,
                                  min_abs_log_fc = 1,
                                  min_freq = 0.1,
                                  remove_ribosomal = TRUE) {
  # lookup table of gene frequencies in groups
  gene_frq <-
    calcExprFreqs(nb) %>% 
    assays() %>%
    as.list() %>% 
    map_dfr(as_tibble, rownames = "gene", .id = "cluster_id") %>%
    select(cluster_id, gene, I:IV) %>% 
    pivot_longer(I:IV, names_to = "frq_col", values_to = "frq")
  
  contrast_frq <- map_dfr(
    contrast_frq,
    ~gene_frq %>% 
      filter(frq_col %in% .x) %>% 
      group_by(gene, cluster_id) %>% 
      summarise(frq = max(frq)),
    .id = "contrast"
  )
  
  res <-
    data %>%
    left_join(
      contrast_frq,
      by = c("gene", "cluster_id", "contrast")
    ) %>% 
    filter(
      frq >= min_freq,
      p_val <= max_p,
      p_adj.loc <= max_p_adj,
      abs(logFC) >= min_abs_log_fc
    ) %>% 
    mutate(direction = if_else(logFC > 0, "up", "down")) %>% 
    select(!matches(".frq"))
  
  if (remove_ribosomal) {
    ribo_proteins <-
      msigdbr(species = "Homo sapiens", category = "C2") %>% 
      filter(gs_name == "KEGG_RIBOSOME") %>% 
      pull(human_gene_symbol)
    res <- 
      res %>%
      filter(!gene %in% ribo_proteins)
  }
  
  res
}


contrast_frq <- list(
  II_vs_I = c("I", "II"),
  III_vs_I = c("I", "III"),
  IV_vs_I = c("I", "IV")
)

dge_mm_results_filtered <- filter_dge_mm_results(dge_mm_results, contrast_frq)



# Volcano plots -----------------------------------------------------------

plot_volcano <- function(data,
                         data_filtered,
                         contrast,
                         max_p_adj = 0.05,
                         min_abs_log_fc = 1,
                         min_freq = 0.1,
                         filename = NULL) {
  genes_count <-
    data_filtered %>%
    count(cluster_id, direction, .drop = FALSE) %>%
    pivot_wider(names_from = "direction", values_from = "n")
  
  dot_colors <- c(
    up = "red",
    down = "blue"
  )
  
  p <-
    ggplot(data, aes(logFC, -log10(p_val))) +
    geom_point(
      color = "gray80",
      alpha = .5,
      size = 0.5
    ) +
    geom_point(
      data = data_filtered,
      aes(color = direction),
      alpha = .5,
      size = 0.5
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = down),
      npcx = 0.2,
      npcy = 0.9,
      size = 3,
      color = dot_colors["down"],
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = up),
      npcx = 0.8,
      npcy = 0.9,
      size = 3,
      color = dot_colors["up"],
      na.rm = TRUE
    ) +
    scale_y_continuous(limits = c(NA, 40)) +
    scale_color_manual(
      "differentially\nexpressed",
      values = dot_colors,
      guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
      drop = FALSE
    ) +
    facet_wrap(vars(cluster_id), drop = FALSE) +
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

plot_volcano(dge_mm_results, dge_mm_results_filtered, "II_vs_I",
             filename = "dge_mm/volcano_II")



# Comparison to pseudobulk ------------------------------------------------

dge_results %>%
  filter(p_val < 0.1) %>% 
  left_join(
    dge_mm_results %>% filter(p_val < 0.1),
    by = c("gene", "cluster_id", "contrast"),
    suffix = c("_pb", "_mm")
  ) %>% 
  filter(contrast == "II_vs_I") %>% 
  ggplot(aes(logFC_pb, logFC_mm)) +
  geom_point(alpha = .25, size = .25) +
  coord_fixed() +
  facet_wrap(vars(cluster_id)) +
  theme_bw()
ggsave_default("dge_mm/comparison_pb_II")




# Violin plots ------------------------------------------------------------

dge_mm_results_filtered %>% 
  group_by(contrast, direction) %>% 
  arrange(desc(logFC), .by_group = TRUE) %>% 
  slice(1:3) %>% 
  pwalk(
    function(gene, contrast, cluster_id, direction, ...) {
      plotExpression(nb, gene, x = "sample_id",
                     colour_by = "group_id")
      filename <- str_glue("dge_mm/violin_{contrast}_{cluster_id}_{gene}_{direction}")
      ggsave_default(filename)
    }
  )
