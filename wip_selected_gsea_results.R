plot_gsea_dots <- function(data,
                           comparisons,
                           cell_types,
                           db = "MSigDB_Hallmark_2020",
                           top_n_positive = 5L,
                           top_n_negative = 5L,
                           max_p_adj = 0.05,
                           min_abs_NES = 1,
                           filename = "auto",
                           ...) {
  data <- 
    data %>% 
    mutate(comparison = rename_contrast(comparison)) %>%
    filter(comparison %in% {{comparisons}}, cell_type %in% {{cell_types}})
  
  data_top_terms <-
    data %>%
    filter(
      db == {{db}},
      padj <= max_p_adj, abs(NES) >= min_abs_NES) %>%
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
      cell_type = factor(cell_type, names(CELL_TYPE_ABBREVIATIONS))
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
    ylab(NULL) +
    scale_color_distiller(
      "normalized\nenrichment\nscore",
      palette = "PiYG",
      direction = 1,
      limits = c(-color_limit, color_limit)
    ) +
    scale_size_area(TeX("-log_{10} (p_{adj})")) +
    coord_fixed() +
    facet_wrap(vars(cell_type), nrow = 1) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = .5),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  p
}

plot_gsea_dots(dge$gsea, c("Mc", "Ac", "Sc"), "B")
ggsave_default("dge_mm/gsea_B_TvsC", height = 100)

plot_gsea_dots(dge$gsea, "Mas", "B", top_n_positive = Inf, top_n_negative = Inf)
ggsave_default("dge_mm/gsea_B_MvsN", height = 100)

plot_gsea_dots(dge$gsea, c("Mc", "Ac", "Sc"), "M")
ggsave_default("dge_mm/gsea_M_TvsC", height = 100)

plot_gsea_dots(dge$gsea, "Mas", "M", top_n_positive = Inf, top_n_negative = Inf)
ggsave_default("dge_mm/gsea_M_MvsN", height = 100)
