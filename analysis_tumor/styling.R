require(scico)
require(tidyverse)
# styles for publication-quality figures

# Common definitions ------------------------------------------------------

CELL_TYPE_ABBREVIATIONS <- c(
  "T" = "T cell",
  NK  = "natural killer cell",
  B   = "B cell",
  M   = "monocyte",
  pDC = "plasmacytoid dendritic cell",
  E   = "erythroid lineage cell",
  SC  = "hematopoietic precursor cell",
  NB  = "neuron"
)

# ColorBrewer Set1
CELL_TYPE_COLORS <- c(
  "T" = "#1f78b4",
  "NK" = "#a6cee3",
  "B" = "#33a02c",
  "M" = "#ff7f00",
  "pDC" = "#fdbf6f",
  "E" = "#b15928",
  "SC" = "#6a3d9a",
  "NB" = "#e31a1c",
  "other" = "black",
  na = "gray80"
)

# SAS colors Patients NB 
patient_palette <- qualitative_hcl(5)



PATIENT_COLORS <-
  scico::scico(16, palette = "romaO") %>%
  colorspace::lighten(0.3) %>% 
  set_names(c("C1", "C2", "C3", "C4", "C5",
              "M1", "M2", "M3", "M4",
              "A1", "A2",
              "S1", "S2", "S3", "S4", "S5"))

PATIENT_COLORS_NB <-
  scico::scico(16, palette = "romaO") %>%
  colorspace::lighten(0.3) %>% 
  as_tibble() %>%
  cbind(
    ID=
      c( "S1","A1","M1","C1", "S2","C2","S5","A2", "S3", "M3","C3", "S4","M4","C4","M2","C5")
  ) %>%
  arrange(ID)

PATIENT_COLORS_NB <- 
  as.character(PATIENT_COLORS_NB$value) %>%
  set_names(c(PATIENT_COLORS_NB$ID))

# SAS NB-cluster colors

cluster_palette <- colorspace::qualitative_hcl(n=4,c=50,l=60)

NBCLUSTER_COLORS <- c(
  #"1" = cluster_palette[3],
  #"2" = cluster_palette[1],
  #"3" = cluster_palette[2],
  #"4" = cluster_palette[4]
  "1" = "#6699ff",
  "2" = "#ff6666",
  "3" = "#669900",
  "4" = "#993366"
)


# https://wesandersonpalettes.tumblr.com/post/110716093015/ash-should-we-dance
GROUP_COLORS <- c(
  C = "#b9b09f",
  M = "#a9d8c8",
  A = "#d16b54",
  S = "#e8c95d"
)

# SAS Heatmap colors + settings
color_heat <- rev(brewer.pal(8, "RdBu"))
breaks_heat <- seq(0, 1, length.out = 7)

# ggplot functions --------------------------------------------------------

BASE_TEXT_SIZE_MM = 1.76  # corresponds to 5 pt
BASE_TEXT_SIZE_PT = 5
BASE_LINE_SIZE = 0.25

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
    make_guides <- theme(legend.position = "right")
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



scale_color_dotplot <- function(...) {
  scico::scale_color_scico(
    palette = "oslo",
    direction = -1,
    ...
  )
}

scale_color_gsea <- function(...) {
  scale_color_distiller(
    palette = "RdBu",
    direction = -1,
    oob = scales::oob_squish_any,  # make sure that value capping
    ...                            # yields colored dots
  )
}

scale_color_gsea_2 <- function(...) {
  scale_color_distiller(
    palette = "Reds",
    direction = 1,
    oob = scales::oob_squish_any,  # make sure that value capping
    ...                            # yields colored dots
  )
}

scale_color_GrRd <- function(...){
  scale_color_distiller(
    palette = "RdYlGn",
    direction = -1,
    oob = scales::oob_squish_any,  # make sure that value capping
    ...                            # yields colored dots
  )
}

scale_color_PiYG <- function(...){
 scale_color_gradient2(low ="#1b7837",
                       mid = "#fee090",
                       high = "#af186f",
                       midpoint = 0.5,
                       oob = scales::oob_squish_any,
                       ...
                      )
 # scale_color_distiller(
#    palette = "PRGn",
#    direction = -1,
#    oob = scales::oob_squish_any,  # make sure that value capping
#    ...                            # yields colored dots
#  ) 
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
rename_str_or_fct <- function(s, nm) {
  if (is.factor(s))
    suppressWarnings(fct_recode(s, !!!deframe(nm[c("new", "old")])))
  else
    recode(s, !!!deframe(nm))
}

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
    "2018_1404", "M2",
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
    "II_vs_I",  "M",
    "III_vs_I", "A",
    "IV_vs_I",  "S"
  )
)


# SaS: rename literature celltypes in Ḿarker-"Enrichment" Plot
rename_celltype <- partial(
  rename_str_or_fct,
  nm = tribble(
    ~old,       ~new,
    "Bridge",  "Bridge",
    "Chromaffin cells", "Chromaffin cells",
    "connecting Progenitor cells",  "Connecting progenitor cells",
    "cycling Neuroblasts","Cycling neuroblasts",
    "cycling SCPs", "Cycling SCPs",
    "late Chromaffin cells", "Late chromaffin cells",
    "late Neuroblasts", "Late neuroblasts",
    "late SCPs", "Late SCPs",
    "Neuroblasts", "Neuroblasts",
    "SCPs", "SCPs"
  )
)

rename_EnrTerms <- function(str){
  str <- str %>% #(EnrData$Cl4 %>% filter(db == "NCI-Nature_2016"))$Term %>% #NCI-Nature_2016 GO_Biological_Process_2018
    str_replace_all("[:punct:]GO[:punct:]","") %>%
    str_replace_all("[0-9]{7}", "") %>%
    str_replace_all("[:space:][:punct:]$","") %>%
    str_replace_all("[:space:][:punct:]human[:punct:]","") %>%
    str_replace_all("[:space:]Homo sapiens","") %>%
    str_replace_all("[:space:][:graph:]{36}$","")
}
