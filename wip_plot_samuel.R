library(tidyverse)
source("common_functions.R")
source("styling.R")

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
