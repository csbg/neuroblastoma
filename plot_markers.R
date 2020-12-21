library(Seurat)
library(tidyverse)
library(fs)
library(patchwork)

ggsave_default <- function(filename, width = 297, height = 210,
                           crop = TRUE, ...) {
  if (is.null(filename))
    return()
  
  filename <- str_glue("plots/{filename}.png")
  filename %>% path_dir() %>% dir_create()
  
  ggsave(filename, dpi = 300, units = "mm", limitsize = FALSE,
         width = width, height = height, ...)
  
  if (crop) knitr::plot_crop(filename)
  
  invisible(filename)
}



# Load data ---------------------------------------------------------------

nb <- readRDS("data_generated/all_datasets_current/nb_only_RNA.rds")

nb@meta.data <- 
  nb@meta.data %>% 
  rownames_to_column("cell") %>% 
  left_join(
    read_csv(
      "data_raw/metadata/sample_groups.csv",
      col_types = "cf",
      comment = "#"
    ) %>%
      mutate(group = fct_relevel(group, "I", "II", "III", "IV")),
    by = "sample"
  ) %>%
  left_join(
    read_csv("data_generated/all_datasets_current/nb_singler.csv"),
    by = "cell"
  ) %>% 
  extract(
    cell_type,
    into = "cell_type_broad",
    regex = "([^:]*)",
    remove = FALSE
  ) %>%
  column_to_rownames("cell")

canonical_markers <-
  read_csv("data_raw/metadata/cell_type_marker_combinations.csv")



# Mitochondrial genes -----------------------------------------------------

FeaturePlot(nb, "percent.mt")

ggplot(nb@meta.data, aes(integrated_snn_res.0.5, percent.mt)) +
  geom_violin(aes(fill = integrated_snn_res.0.5), show.legend = FALSE) +
  geom_jitter(alpha = .1) +
  xlab("Cluter") +
  ylab("% mitochondrial genes")



# NB markers --------------------------------------------------------------

FeaturePlot(
  nb,
  c("L1CAM", "MYCN", "PHOX2B", "B4GALNT1", "NCAM1"),
  min.cutoff = "q5",
  max.cutoff = "q95",
  coord.fixed = TRUE,
  ncol = 3,
  order = TRUE
)
ggsave_default("markers/nb_markers")


Idents(nb) <- "integrated_snn_res.0.5"

DotPlot(
  nb,
  features = c("L1CAM", "MYCN", "PHOX2B", "B4GALNT1", "NCAM1"),
  split.by = "group",
  cols = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
)
ggsave_default("markers/nb_dotplot", height = 400, width = 210)



# Canonical cell type markers ---------------------------------------------

plot_markers <- function(cell_type, genes) {
  message("Plotting ", cell_type)
  
  genes <- str_split(genes, ",")[[1]]
  
  # p <-
  #   FeaturePlot(
  #   nb,
  #   genes,
  #   min.cutoff = "q5",
  #   max.cutoff = "q95",
  #   coord.fixed = TRUE,
  #   order = TRUE,
  #   combine = FALSE
  # ) %>% 
  #   wrap_plots() +
  #   plot_annotation(caption = str_glue("Canonical markers for {cell_type}"))
  # 
  # ggsave_default(str_glue("markers/feature_{cell_type}"))
  # p
  
  Idents(nb) <- "cell_type_broad"
  RidgePlot(nb, genes)
  ggsave_default(str_glue("markers/ridge_{cell_type}"))
}

pwalk(
  canonical_markers[1, ],
  plot_markers
)

FindMarkers()
RidgePlot(
  nb,
  c("IL7R", "CCR7")
)
ggsave_default("markers/ridge")


# TODO
# use SCtransformed data for plots?
# make ridgeplot work
# include other plots
# unify data loading (with functions defined in summary_integration.R)