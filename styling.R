library(tidyverse)
library(fs)
library(openxlsx)



# Common definitions ------------------------------------------------------

# use "patients" (instead of "samples")

# UMAP dimensions are called "UMAP1" and "UMAP2"

# use the following abbreviations and long forms for cell types
# factors should be ordered like this vector
CELL_TYPE_ABBREVIATIONS <- c(
  "T" = "T cell",
  NK  = "natural killer cell",
  B   = "B cell",
  M   = "myeloid cell",
  pDC = "plasmacytoid dendritic cell",
  E   = "erythroid lineage cell",
  SC  = "hematopoietic precursor cell",
  NB  = "neuroblastoma cell"
)

# color cell types using ColorBrewer's Set1
CELL_TYPE_COLORS <- c(
  "T"   = "#1f78b4",
  NK    = "#a6cee3",
  B     = "#33a02c",
  M     = "#ff7f00",
  pDC   = "#fdbf6f",
  E     = "#b15928",
  SC    = "#6a3d9a",
  NB    = "#e31a1c",
  other = "black",
  na    = "gray80"
)

GROUP_NAMES_LONG = c(
  C = "control",
  M = "MYCN amplified",
  A = "ATRX deleted",
  S = "sporadic"
)

# color groups using one of the wesanderson color palettes
# (https://wesandersonpalettes.tumblr.com/post/110716093015/ash-should-we-dance)
# factors should be ordered like this vector
GROUP_COLORS <- c(
  C   = "#b9b09f",
  M   = "#a9d8c8",
  A   = "#d16b54",
  S   = "#e8c95d",
  "T" = "#433447"
)

# same as group colors; Ma (M vs A) etc are blends
CONTRAST_COLORS <- c(
  "M vs C" = "#a9d8c8",
  "A vs C" = "#d16b54",
  "S vs C" = "#e8c95d",
  "M vs A" = "#bda28e",
  "M vs S" = "#c9d193",
  "A vs S" = "#dd9a59",
  "M vs A+S" = "#c3b991"
)

PATIENT_ORDER <- c("C1", "C2", "C3", "C4", "C5",
                   "M1", "M2", "M3", "M4",
                   "A1", "A2",
                   "S1", "S2", "S3", "S4", "S5")

# color patients using the cyclic romaO from the Scientific Colour Maps
PATIENT_COLORS <-
  scico::scico(16, palette = "romaO") %>%
  colorspace::lighten(0.3) %>% 
  set_names(c("S1", "A1", "M1", "C1",
              "S2", "C2", "S5", "A2",
              "S3", "M3", "C3", "S4",
              "M4", "C4", "M2", "C5")) %>% 
  magrittr::extract(PATIENT_ORDER)

# tumor subclusters
SUBCLUSTER_COLORS <- c(
  "1" = "#6699ff",
  "2" = "#ff6666",
  "3" = "#669900",
  "4" = "#993366"
)

# colors for adrenal medullary cell types
ADRMED_CELLS_COLORS <- c(
  "late SCPs" = "#921813",
  "SCPs" = "#be202e",
  "cycling SCPs" = "#be6867",
  "Bridge" = "#ef845f",
  "connecting Chromaffin cells" = "#33ac75",
  "Chromaffin cells" = "#006e3b",
  "late Chromaffin cells" = "#686c58",
  "cycling Neuroblasts" = "#aacedc",
  "Neuroblasts" = "#244a94",
  "late Neuroblasts" = "#303d63"
)

# MYCN status
MYCN_STATUS_COLORS <- c(
  "normal" = "gray95",
  "amplified" = "gray50"
)

# consistently differentially expressed genes
CONSISTENT_GENES_COLORS <- c(
  "A, M, and S" = "#205d89",
  "M and S" = "#cf784b",
  "A and S" = "#73a87c",
  "A and M" = "#c1bc78"
)

# Paul Tol's bright
MYELOID_COLORS <- c(
  "classical mono" = "#ccbb44",
  "nonclassical mono" = "#66ccee",
  "mDCs" = "#ee6677",
  other = "#bbbbbb"
)



# ggplot functions --------------------------------------------------------

BASE_TEXT_SIZE_MM = 1.76  # mm, corresponds to 5 pt, use in geom_text()
BASE_TEXT_SIZE_PT = 5 # pt, use in theme()
BASE_LINE_SIZE = 0.25 # pt

#' Common theme for figures in the publication.
#'
#' This theme bases upon `theme_bw()` and ensures
#' - common line widths of `BASE_LINE_SIZE`
#' - common text sizes of `BASE_TEXT_SIZE_PT`
#' - a uniform plot margin of 1 mm
#' - a medium strip text, an empty strip background, and
#'   1 mm padding between strip text and panel
#'
#' @param grid If `TRUE`, show the grid.
#' @param rotate If `TRUE`, rotate x-axis tick labels by 90°.
#' @param ... Other parameters passed to `theme_bw()`.
#'
#' @return A theme object.
theme_nb <- function(grid = TRUE,
                     rotate = FALSE,
                     ...){
  res <-
    theme_bw(...) + 
    theme(
      line = element_line(size = BASE_LINE_SIZE),
      axis.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      axis.title = element_text(color = "black", size = BASE_TEXT_SIZE_PT),
      legend.background = element_blank(),
      legend.text = element_text(color = "black", size = BASE_TEXT_SIZE_PT), 
      legend.title = element_text(size = BASE_TEXT_SIZE_PT),
      panel.border = element_rect(size = BASE_LINE_SIZE * 2),
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      strip.background = element_blank(),
      strip.text = element_text(
        color = "black",
        size = BASE_TEXT_SIZE_PT
      ),
      strip.text.x = element_text(margin = margin(b = 1, unit = "mm")),
      strip.text.y = element_text(margin = margin(l = 1, unit = "mm"))
    )
  
  if (!grid)
    res <-
      res +
      theme(panel.grid = element_blank())
  if (rotate)
    res <-
      res +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  res
}



#' Save a publication-quality plot.
#'
#' @param filename Filename, will be saved in subfolder `plots/final`.
#'                 If `NULL`, exit without creating a plot.
#' @param type: Type of image file.
#' @param plot If `NULL`, save the `last_plot()` via `ggsave()`.
#'   Otherwise, save the graphics object in `plot` via the
#'   `png()/pdf(); print(); dev.off()` workflow.
#' @param legends If `FALSE`, save plot without legends
#'   and append "_noLegends" to the filename.
#' @param dpi Resolution.
#' @param width Width in cm.
#' @param height Height in cm.
#' @param ... Other parameters passed to the plotting function.
ggsave_publication <- function(filename,
                               type = "pdf",
                               plot = NULL,
                               legends = TRUE,
                               dpi = 1200,
                               width = 4,
                               height = 4,
                               ...) {
  if (is.null(filename))
    return()
  
  # construct filename
  if (legends) {
    legendstr <- ""
    make_guides <- NULL
  } else {
    legendstr <- "_noLegends"
    make_guides <- theme(legend.position = "none")
  }
  
  filename <- str_glue("plots/final/{filename}{legendstr}.{type}")
  filename %>%
    path_dir() %>%
    dir_create()
  
  if (is.null(plot)) {
    # if last_plot() is available, use ggsave()
    ggsave(
      filename,
      plot = last_plot() + make_guides,
      dpi = dpi,
      units = "cm",
      limitsize = FALSE,
      width = width,
      height = height,
      ...
    )  
  } else {
    # for non-ggplot objects, use the base R functions directly;
    # only png and pdf backends are supported
    if (type == "png") {
      png(
        filename,
        res = dpi,
        units = "cm",
        width = width,
        height = height,
        ...
      )
    } else if (type == "pdf") {
      pdf(
        filename,
        width = width / 2.54,  # dimensions for pdf() must be inches
        height = height / 2.54,
        ...
      )
    } else {
      stop("Type", type, "cannot be saved.")
    }
    
    print(plot)
    dev.off()
  }
}

# color scale for gene expression dotplots: oslo from Scientific Colour Maps
scale_color_dotplot <- function(...) {
  scico::scale_color_scico(
    palette = "oslo",
    direction = -1,
    ...
  )
}

# color scale for enrichment analyses and gene expression heatmaps:
# ColorBrewer's RdBu
scale_color_gsea <- function(...) {
  scale_color_distiller(
    palette = "RdBu",
    direction = -1,
    oob = scales::oob_squish_any,  # make sure that value capping
    ...                            # yields colored dots
  )
}

# color scale for enrichment analyses (only positive odds ratios):
# ColorBrewer's Reds
scale_color_gsea_2 <- function(...) {
  scale_color_distiller(
    palette = "Reds",
    direction = 1,
    oob = scales::oob_squish_any,  # make sure that value capping
    ...                            # yields colored dots
  )
}

# color scale for CytoTRACE scores:
# ColorBrewer's PiYG with yellow as midpoint
scale_color_cytotrace <- function(...) {
  scale_color_gradient2(
    low = "#1b7837",
    mid = "#fee090",
    high = "#af186f",
    midpoint = 0.5,
    oob = scales::oob_squish_any,
    ...
  )
}



# Renaming functions ------------------------------------------------------

#' Rename a string or factor while retaining level order.
#'
#' @param s String or factor to be renamed.
#' @param nm Tibble with columns "old" and "new", containing old and new names.
#'
#' @return Renamed string or factor.
rename_str_or_fct <- function(s, nm) {
  if (is.factor(s))
    suppressWarnings(fct_recode(s, !!!deframe(nm[c("new", "old")])))
  else
    recode(s, !!!deframe(nm))
}

# use this function to create renaming functions
# for groups, patients, and contrasts

rename_groups <- partial(
  rename_str_or_fct,
  nm = tribble(
    ~old,  ~new,
    "I",   "C",
    "II",  "M",
    "III", "A",
    "IV",  "S"
  )
)

rename_patients <- partial(
  rename_str_or_fct,
  nm = tribble(
    ~old,        ~new,
    "2014_0102", "C1",
    "2016_1853", "C2",
    "2016_2950", "C3",
    "2018_4252", "C4",
    "2020_1288", "C5",
    "2016_4503", "M1",
    "2018_1404", "M2", # R1_NB_2018_1404 is M2a, MF220_NB_BM_Patient1 is M2b 
    "2019_5022", "M3",
    "2019_5754", "M4",
    "2005_1702", "A1",
    "2016_3924", "A2",
    "2006_2684", "S1",
    "2018_1625", "S2",
    "2018_6056", "S3",
    "2019_2495", "S4",
    "2020_1667", "S5"
  )
)

rename_contrast <- partial(
  rename_str_or_fct,
  nm = tribble(
    ~old,       ~new,
    "II_vs_I",  "M vs C",
    "III_vs_I", "A vs C",
    "IV_vs_I",  "S vs C",
    "II_vs_IV", "M vs S",
    "III_vs_IV", "A vs S",
    "II_vs_III", "M vs S",
    "MNA_vs_other", "M vs A+S",
  )
)

rename_contrast_long <- partial(
  rename_str_or_fct,
  nm = tribble(
    ~old,       ~new,
    "II_vs_I",  "MYCN amplified vs control",
    "III_vs_I", "ATRX deleted vs control",
    "IV_vs_I",  "sporadic vs control",
    "II_vs_IV", "MYCN amplified vs sporadic",
    "III_vs_IV", "ATRX deleted vs sporadic",
    "II_vs_III", "MYCN amplified vs ATRX deleted",
    "MNA_vs_other", "MYCN amplified vs MYCN normal"
  )
)

rename_adrenal_cells <- partial(
  rename_str_or_fct,
  nm = tribble(
    ~old,                          ~new,
    "Bridge",                      "Bridge",
    "Chromaffin cells",            "Chromaffin cells",
    "connecting Progenitor cells", "Connecting progenitor cells",
    "cycling Neuroblasts",         "Cycling neuroblasts",
    "cycling SCPs",                "Cycling SCPs",
    "late Chromaffin cells",       "Late chromaffin cells",
    "late Neuroblasts",            "Late neuroblasts",
    "late SCPs",                   "Late SCPs",
    "Neuroblasts",                 "Neuroblasts",
    "SCPs",                        "SCPs"
  )
)

rename_myeloid <- partial(
  rename_str_or_fct,
  nm = tribble(
    ~old,  ~new,
    "classical mono",    "classical monocytes",
    "nonclassical mono", "nonclassical monocytes",
    "mDCs",              "myeloid dendritic cells",
    "other",             "other"
  )
)

rename_myeloid_newline <- partial(
  rename_str_or_fct,
  nm = tribble(
    ~old,  ~new,
    "classical mono",    "classical\nmonocytes",
    "nonclassical mono", "nonclassical\nmonocytes",
    "mDCs",              "myeloid\ndendritic cells",
    "other",             "other"
  )
)



# Tables ------------------------------------------------------------------

#' Save a (list of) data frame as well-formatted XLSX file.
#' 
#' @param tables Data frames to be saved. A list of dataframes is saved into
#'   separate worksheets, whose names equal the list names. If the list is
#'   unnamed, worksheets will be named "Sheet1" etc.
#' @param filename Filename without extension; will be saved to `tables/`.
#' @param sheet_name Sheet name if a single unnamed table should be saved.
#'
#' @return Nothing.
save_table <- function(tables, filename, sheet_name = "Sheet1") {
  filename <- str_glue("tables/{filename}.xlsx")
  wb <- createWorkbook()
  
  # ensure that tables is a named list
  if (inherits(tables, "list")) {
    if (is.null(names(tables)))
      tables <- set_names(tables, paste0("Sheet", seq_along(tables)))
  } else if (inherits(tables, "data.frame")) {
    tables <- list(tables) %>% set_names(sheet_name)
  } else {
    stop("'tables' must be a data frame or list of data frames.")
  }

  # populate Excel file with worksheets  
  iwalk(
    tables,
    function(table, sheet_name) {
      addWorksheet(wb, sheet_name)
      writeData(
        wb,
        sheet_name,
        table,
        headerStyle = createStyle(textDecoration = "bold")
      )
      freezePane(wb, sheet_name, firstRow = TRUE)
      setColWidths(wb, sheet_name, 1:ncol(table), "auto")
    }
  )
  
  saveWorkbook(wb, filename, overwrite = TRUE)
}
