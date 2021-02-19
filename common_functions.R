# Functions that are used in several analysis scripts

#' Save a plot with sensible defaults.
#'
#' @param filename Filename, will be saved in subfolder `plots/`. May contain
#'   additional subfolders, which are possibly created. If `NULL`, exit without
#'   creating a plot.
#' @param plot If `NULL`, save the `last_plot()` via `ggsave()`. Otherwise, save
#'   the graphics object in `plot` via the `png(); print(); dev.off()` workflow.
#' @param width Width in mm.
#' @param height Height in mm.
#' @param crop If `TRUE`, remove white margins from the saved plot.
#' @param ... Other parameters passed to `ggsave()` or `png()`.
#'
#' @return The filename, invisibly.
ggsave_default <- function(filename,
                           plot = NULL,
                           width = 297,
                           height = 210,
                           crop = TRUE,
                           ...) {
  if (is.null(filename))
    return()
  
  filename <- stringr::str_glue("plots/{filename}.png")
  filename %>%
    fs::path_dir() %>%
    fs::dir_create()
  
  if (is.null(plot)) {
    ggplot2::ggsave(filename, dpi = 300, units = "mm", limitsize = FALSE,
                    width = width, height = height, ...)  
  } else {
    png(filename, res = 300, width = width, height = height, units = "mm", ...)
    print(plot)
    dev.off()
  }
  
  if (crop)
    knitr::plot_crop(filename)
  
  invisible(filename)
}