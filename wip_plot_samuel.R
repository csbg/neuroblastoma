library(monocle3)
library(metR)
library(tidyverse)
library(latex2exp)
source("common_functions.R")
source("styling.R")



# Enrichments -------------------------------------------------------------

gsea_data <- readRDS("~/Desktop/NB_Enrichment_Cl4.rds")
gsea_data

plot_gsea <- function(data,
                      db,
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
      db == {{db}},
      Adjusted.P.value <= max_p_adj,
      Odds.Ratio >= min_OR,
      overlap_size >= min_overlap_size
    ) %>% 
    pull(Term) %>%
    unique()
  
  data_vis <- 
    data %>% 
    filter(
      db == {{db}},
      Term %in% top_terms
    ) %>% 
    mutate(
      is_significant =
        Adjusted.P.value <= max_p_adj &
        Odds.Ratio >= min_OR &
        overlap_size >= min_overlap_size,
      Term =
        Term %>% 
        str_match("(.+) Homo") %>% 
        magrittr::extract(, 2) %>% 
        as_factor() %>%
        fct_reorder(log2(Odds.Ratio), sum, na.rm = TRUE),
    )
  # return(data_vis)
  
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
  
  color_limit <- log2(10)  # max(abs(data_vis$NES))
  
  if (circle_significant) {
    circle_significant <- geom_point(
      data = data_vis %>% filter(is_significant),
      shape = 1
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
    scale_color_gsea(
      name = TeX("log_2 odds ratio"),
      limits = c(0, color_limit),
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(10, "mm"),
        label.position = "top",
        title.vjust = 0.1,
        order = 2
      )
    ) +
    scale_size_area(
      name = TeX("-log_{10} p_{adj}"),
      max_size = 2.5,
      guide = guide_legend(order = 1)
    )  +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      legend.box = "horizontal",
      legend.box.just = "bottom",
      legend.direction = "horizontal",
      legend.justification = "right",
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.position = c(1, 1.07),
      legend.spacing = unit(0, "mm"),
      legend.margin = margin(0, 1, -3, 1, "mm"),
      panel.spacing = unit(-.5, "pt"),
      plot.margin = margin(t = .7, unit = "cm")
    )
  
  p
}

plot_gsea(gsea_data, "Panther_2016")
ggsave_publication("2_enrichment", width = 6, height = 8)



# Trajectory --------------------------------------------------------------

cds <- readRDS("~/Desktop/NB_CDS_trajectory.rds")
cds

tumor_data <- 
  read_csv("rna_velocity/tumor_velocity_umap.csv") %>%
  mutate(tumor_subcluster = factor(tumor_subcluster) %>% fct_inseq())
tumor_data


make_vector_grid <- function(df, bins = 50) {
  get_breaks <- function(col) {
    seq(
      df %>% pull({{col}}) %>% min(),
      df %>% pull({{col}}) %>% max(),
      length.out = bins + 1
    )[-bins - 1]
  }  
  
  breaks_x <- get_breaks(umap_1) 
  breaks_y <- get_breaks(umap_2)
  
  bin_centers_x <- breaks_x + (breaks_x[2] - breaks_x[1]) / 2
  bin_centers_y <- breaks_y + (breaks_y[2] - breaks_y[1]) / 2
  
  vecs <- 
    df %>% 
    mutate(
      x = findInterval(umap_1, breaks_x, all.inside = TRUE),
      y = findInterval(umap_2, breaks_y, all.inside = TRUE),
    ) %>% 
    group_by(x, y) %>% 
    summarise(dx = mean(dx), dy = mean(dy))
  
  cross_df(list(x = 1:bins, y = 1:bins)) %>% 
    mutate(
      center_x = bin_centers_x[x],
      center_y = bin_centers_y[y],
    ) %>% 
    left_join(vecs, by = c("x", "y")) %>%
    replace_na(list(dx = 0, dy = 0))
}

vg <- make_vector_grid(tumor_data)




plot_direction <- function(plot_velocity = TRUE, plot_trajectory = TRUE) {
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
  
  if (plot_velocity) {
    plot_velocity <- list(
      geom_streamline(
        data = vg,
        aes(
          x = center_x,
          y = center_y,
          dx = dx,
          dy = dy,
          alpha = ..step..,
        ),
        # min.L = 2.5,
        res = 5,
        n = 40,
        arrow = NULL,
        show.legend = FALSE,
        size = BASE_LINE_SIZE / 2
      ),
      scale_alpha(limits = c(0, NA), oob = scales::oob_squish)
    )
  } else {
    plot_velocity <- NULL
  }
  
  if (plot_trajectory) {
    plot_trajectory <- geom_segment(
      data = edge_df,
      aes(
        x = source_umap_1, 
        y = source_umap_2,
        xend = target_umap_1, 
        yend = target_umap_2
      ),
      size = BASE_LINE_SIZE,
      lineend = "round",
      color = "#4393c3"
    )
  } else {
    plot_trajectory <- NULL
  }
  
  colData(cds) %>% 
    as_tibble(rownames = "cell") %>% 
    left_join(
      reducedDim(cds, "UMAP") %>%
        magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
        as_tibble(rownames = "cell"),
      by = "cell"
    ) %>% 
    ggplot(aes(umap_1, umap_2)) +
    # geom_point(color = "white") +
    geom_point(
      aes(color = cytoTRACE),
      shape = 16,
      size = 0.1,
      show.legend = FALSE
    ) +
    plot_velocity +
    plot_trajectory +
    scale_color_distiller(palette = "PiYG") +
    coord_fixed() +
    theme_nb(grid = FALSE)  
}


plot_direction()
ggsave_publication("2_trajectory", width = 4, height = 4)




# alterlatively, plot the vector field
ggplot(tumor_data, aes(umap_1, umap_2)) +
  geom_point(aes(color = tumor_subcluster)) +
  geom_segment(
    data = slice_sample(tumor_data, prop = 0.25),
    aes(xend = umap_1 + dx * 10, yend = umap_2 + dy * 10),
    arrow = arrow(length = unit(1, "mm"))
  ) +
  coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())



# Violin plots ------------------------------------------------------------
# 
# plot_violin2 <- function(data_cds,
#                          genes,
#                          cluster_number) {
  
genes <- genes[[cluster_number]]$gene

data_cds <- scuttle::logNormCounts(
  data_cds,
  size.factors = colData(data_cds)$Size_Factor,
  # subset.row = genes,
  # normalize.all = T
)

data <- 
  logcounts(data_cds)[genes, , drop = FALSE] %>%
  t() %>%
  as.matrix() %>%
  as_tibble(rownames = "cell") %>%
  pivot_longer(!cell, 
               names_to = "gene", 
               values_to = "logexp") %>%
  left_join(
    data_cds@colData %>%
      as_tibble(rownames = "cell") %>% 
      select(cell, group, sample, cluster, Size_Factor),
    by = "cell"
  )

data %>% 
  # filter(gene == "STAT3") %>% 
  group_by(cluster, gene) %>% 
  mutate(logexp = logexp / max(logexp)) %>% 
  ungroup() %>% 
  ggplot(aes(cluster, logexp)) +
  geom_violin(
    aes(fill = cluster),
    scale = "width",
    size = BASE_LINE_SIZE,
    draw_quantiles = T,
    na.rm = T
  ) +
  geom_jitter(
    width = 0.25,
    alpha = 0.1,
    size = 0.1
  ) +
  stat_summary_bin(geom = "point", fun = mean, size = 0.1) +
  xlab(NULL) +
  scale_y_continuous(
    name = TeX("log_2 expr_{}"),
    # limits = c(0, 1),
    # breaks = c(0, 0.5, 1),
    # expand = expansion(add = .025)
    oob = scales::oob_keep
  ) +
  scale_fill_manual(values = NBCLUSTER_COLORS) +
  scale_size_area(max_size = 0.01) +
  coord_flip() +
  facet_wrap(vars(gene),ncol = 4) +
  theme_nb(grid = FALSE) +
  theme(legend.position = "none")
# }

data_cds <- readRDS("~/Desktop/NB_CDS_aligned.rds")
genes <- readRDS("~/Desktop/NB_cluster_markers.rds")
cluster_number <- "CL1"

tumor_subclusters <- 
  read_csv("rna_velocity/tumor_velocity_umap.csv") %>%
  transmute(
    cell = cell,
    tumor_subcluster = factor(tumor_subcluster) %>% fct_inseq()
  ) %>% 
  deframe()

colData(data_cds)$cluster <- tumor_subclusters[rownames(data_cds@colData)]




cluster_palette <- colorspace::qualitative_hcl(4)

NBCLUSTER_COLORS <- c(
  "1" = cluster_palette[3],
  "2" = cluster_palette[1],
  "3" = cluster_palette[2],
  "4" = cluster_palette[4]
)


plot_violin2(data_cds, genes, "CL1")


data %>% 
  group_by(cluster, gene) %>% 
  mutate(logexp = logexp / max(logexp))