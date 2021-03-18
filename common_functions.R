# Functions that are used in several analysis scripts



# Visualization -----------------------------------------------------------

#' Save a plot with sensible defaults.
#'
#' @param filename Filename, will be saved in subfolder `plots/`. May contain
#'   additional subfolders, which are possibly created. If `NULL`, exit without
#'   creating a plot.
#' @param type: Type of image file.
#' @param plot If `NULL`, save the `last_plot()` via `ggsave()`. Otherwise, save
#'   the graphics object in `plot` via the `png(); print(); dev.off()` workflow.
#' @param width Width in mm.
#' @param height Height in mm.
#' @param crop If `TRUE`, remove white margins from the saved plot.
#' @param ... Other parameters passed to `ggsave()` or `png()`.
#'
#' @return The filename, invisibly.
ggsave_default <- function(filename,
                           type = "png",
                           plot = NULL,
                           width = 297,
                           height = 210,
                           crop = TRUE,
                           ...) {
  if (is.null(filename))
    return()
  
  filename <- stringr::str_glue("plots/{filename}.{type}")
  filename %>%
    fs::path_dir() %>%
    fs::dir_create()
  
  if (is.null(plot)) {
    ggplot2::ggsave(filename, dpi = 300, units = "mm", limitsize = FALSE,
                    width = width, height = height, ...)  
  } else {
    rlang::exec(type, filename, res = 300, units = "mm",
                width = width, height = height,  ...)
    print(plot)
    dev.off()
  }
  
  if (crop)
    knitr::plot_crop(filename)
  
  invisible(filename)
}

#' Make a Seurat-style dotplot of features vs clusters.
#'
#' @param counts Count matrix whose rownames correspond to features.
#' @param features Vector of features, plotted in the given order.
#' @param groups Factor of groups to which each cell belongs, plotted in the
#'   order of its levels.
#' @param min_exp Lower limit for the scaled average expression.
#' @param max_exp Upper limit for the scaled average expression.
#' @param panel_annotation A dataframe used for drawing rectangles on the panel
#'   background. Must contain three columns "xmin", "xmax" (aesthetics for
#'   `geom_rect()`), and "fill" (color name).
#'
#' @return A ggplot object.
plot_dots <- function(counts, features, groups,
                      min_exp = -2.5, max_exp = 2.5, panel_annotation = NULL) {
  scale_and_limit <- function(x) {
    scale(x)[,1] %>% 
      pmax(min_exp) %>% 
      pmin(max_exp)
  }
  
  vis_data <- 
    counts[features, , drop = FALSE] %>% 
    Matrix::t() %>% 
    as.matrix() %>% 
    as_tibble(rownames = "cell") %>% 
    group_by(id = groups) %>% 
    summarise(
      across(
        where(is.numeric),
        list(
          avg_exp = ~mean(expm1(.)),
          pct_exp = ~length(.[. > 0]) / length(.) * 100
        ),
        .names = "{.col}__{.fn}"
      )
    ) %>% 
    mutate(across(ends_with("avg_exp"), scale_and_limit)) %>%
    pivot_longer(
      !id,
      names_to = c("feature", ".value"),
      names_pattern = "(.+)__(.+)"
    ) %>% 
    mutate(feature = factor(feature, levels = features))
  
  if (!is.null(panel_annotation)) {
    panel_bg <- list(
      geom_point(color = "white"),  # initialize discrete coordinate system
      geom_rect(
        data = panel_annotation,
        aes(xmin = xmin, xmax = xmax,
            ymin = 0.5, ymax = nlevels(vis_data$feature) + 0.5,
            fill = fill),
        show.legend = FALSE,
        inherit.aes = FALSE,
      ),
      scale_fill_identity()
    )
  } else {
    panel_bg <- NULL
  }
  
  ggplot(vis_data, aes(id, feature)) +
    panel_bg +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    scale_x_discrete("cluster", expand = expansion(add = 0.5)) +
    scale_y_discrete("feature", expand = expansion(add = 0.5)) +
    scale_color_scico(
      "scaled\naverage\nexpression",
      palette = "oslo",
      direction = -1,
      aesthetics = "color"
    ) +
    scale_radius("% expressed", range = c(0, 6)) +
    coord_fixed(
      # xlim = c(0.5, nlevels(vis_data$id) + 0.5),
      # ylim = c(0.5, nlevels(vis_data$feature) + 0.5),
      clip = "off"
    ) +
    theme_classic() +
    theme(panel.grid = element_blank())
}




# Logging -----------------------------------------------------------------

default_logger <- log4r::logger(threshold = "DEBUG")

debug <- function(..., .envir = parent.frame()) {
  log4r::debug(default_logger, glue::glue(..., .envir = .envir))
}

info <- function(..., .envir = parent.frame()) {
  log4r::info(default_logger, glue::glue(..., .envir = .envir))
}

warn <- function(..., .envir = parent.frame()) {
  log4r::warn(default_logger, glue::glue(..., .envir = .envir))
}
