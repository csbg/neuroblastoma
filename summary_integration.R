library(tidyverse)



# Load data ---------------------------------------------------------------

nb_data <-
  read_csv(
    "data_generated/all_datasets/nb_tsne_integrated.csv",
    col_types = cols(
      sample = "f",
      seurat_clusters = "f",
      SingleR_cells = "f",
      SingleR_clusters = "f"
    )
  ) %>% 
  mutate(
    seurat_clusters = fct_inseq(seurat_clusters),
    SingleR_cells = fct_infreq(SingleR_cells),
    SingleR_clusters = fct_infreq(SingleR_clusters),
  )



# Basic plots -------------------------------------------------------------

nb_data %>% 
  ggplot(aes(tSNE_1, tSNE_2)) +
  geom_point(aes(color = seurat_clusters), size = .1, show.legend = TRUE) +
  facet_wrap(vars(sample), nrow = 3) +
  coord_fixed() +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "t-SNEs after integration via Seurat",
    caption = "Cells colored by clusters calculated via Seurat::FindClusters(resolution = 0.5)"
  ) +
  NULL

ggsave("plots/clusters.pdf",
       units = "mm", width = 297, height = 210)


nb_data %>% 
  ggplot(aes(tSNE_1, tSNE_2)) +
  geom_point(
    aes(color = SingleR_cells),
    size = .1,
    show.legend = TRUE
  ) +
  facet_wrap(vars(sample), nrow = 3) +
  scale_color_discrete(guide = guide_legend(override.aes = list(size = 3))) +
  coord_fixed() +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "t-SNEs after integration via Seurat",
    caption = "Cells colored by cell type determined by SingleR"
  ) +
  NULL

ggsave("plots/celltype_SingleR_cells.pdf",
       units = "mm", width = 297, height = 210)


nb_data %>% 
  ggplot(aes(tSNE_1, tSNE_2)) +
  geom_point(
    aes(color = SingleR_clusters),
    size = .1,
    show.legend = TRUE
  ) +
  facet_wrap(vars(sample), nrow = 3) +
  scale_color_discrete(guide = guide_legend(override.aes = list(size = 3))) +
  coord_fixed() +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "t-SNEs after integration via Seurat",
    caption = "Cells colored by cell type determined by SingleR (using aggregated cluster profiles)"
  ) +
  NULL

ggsave("plots/celltype_SingleR_clusters.pdf",
       units = "mm", width = 297, height = 210)



# SingleR results ---------------------------------------------------------

facet_cell_types <- function(method, nrow = NULL) {
  nb_data %>%
    ggplot(aes(tSNE_1, tSNE_2)) +
    geom_point(
      color = "gray80",
      data = nb_data %>% select(!.data[[method]]),
      size = .1,
    ) +
    geom_point(
      aes(color = .data[[method]]),
      size = .1,
      show.legend = FALSE
    ) +
    facet_wrap(vars(.data[[method]]), nrow = nrow) +
    coord_fixed() +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = "t-SNE after integration via Seurat",
      subtitle = str_glue("Cell types from SingleR are highlighted ({method})")
    ) +
    NULL

  ggsave(str_glue("plots/celltype_{method}_highlight.png"),
         dpi = 300, units = "mm", width = 297, height = 210)
}

facet_cell_types("SingleR_clusters", nrow = 2)
facet_cell_types("SingleR_cells", nrow = 4)



facet_patients <- function(method,
                           cell_type,
                           cell_color = "black",
                           index = 0) {
  message("Plotting ", cell_type, ", method ", method)
  p <- 
    nb_data %>% 
    ggplot(aes(tSNE_1, tSNE_2)) +
    geom_point(
      color = "gray80",
      size = .1,
      show.legend = TRUE
    ) +
    geom_point(
      color = cell_color,
      data = nb_data %>% filter(.data[[method]] == {{cell_type}}),
      size = .1,
    ) +
    facet_wrap(vars(sample), nrow = 3) +
    scale_color_discrete(guide = guide_legend(override.aes = list(size = 3))) +
    coord_fixed() +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = "t-SNEs after integration via Seurat",
      subtitle = str_glue("Cell type '{cell_type}' is highlighted (method {method})")
    ) +
    NULL
  
  filename <- str_glue("plots/celltype_{method}_highlight_{index}_{cell_type}.png")
  ggsave(filename, units = "mm", dpi = 300, width = 297, height = 210)
}

facet_patients("SingleR_clusters", "T_cells")
facet_patients("SingleR_cells", "T_cells")


make_all_highlight_plots <- function(method) {
  cell_types <- levels(nb_data[[method]])
  pal <- scales::hue_pal()(length(cell_types))
  
  list(
    cell_type = cell_types,
    cell_color = pal,
    index = seq_along(cell_types)
  ) %>% 
    pwalk(facet_patients, method = method)  
}

make_all_highlight_plots("SingleR_clusters")
make_all_highlight_plots("SingleR_cells")

