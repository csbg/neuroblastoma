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

# https://wesandersonpalettes.tumblr.com/post/110716093015/ash-should-we-dance
GROUP_COLORS <- c(
  C = "#b9b09f",
  M = "#a9d8c8",
  A = "#d16b54",
  S = "#e8c95d"
)

PATIENT_COLORS <-
  scico::scico(16, palette = "romaO") %>%
  colorspace::lighten(0.3) %>% 
  set_names(c("C1", "C2", "C3", "C4", "C5",
              "M1", "M2", "M3", "M4",
              "A1", "A2",
              "S1", "S2", "S3", "S4", "S5"))



# ggplot functions --------------------------------------------------------

BASE_TEXT_SIZE_MM = 1.76  # corresponds to 5 pt
BASE_TEXT_SIZE_PT = 5
BASE_LINE_SIZE = 0.25

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
                               guides = TRUE,
                               dpi = 1200,
                               width = 4,
                               height = 4,
                               ...) {
  if (is.null(filename))
    return()
  
  if (guides) {
    guidestr <- ""
    make_guides <- NULL
  } else {
    guidestr <- "_noGuides"
    make_guides <- theme(legend.position = "none")
  }
  filename <- stringr::str_glue("plots/final/{filename}{guidestr}.{type}")
  filename %>%
    fs::path_dir() %>%
    fs::dir_create()
  
  if (is.null(plot)) {
    ggplot2::ggsave(filename, plot = last_plot() + make_guides,
                    dpi = dpi, units = "cm", limitsize = FALSE,
                    width = width, height = height, ...)  
  } else {
    if (type == "png") {
      png(filename, res = dpi, units = "cm",
          width = width, height = height, ...)
    } else if (type == "pdf") {
      pdf(filename, width = width / 2.54, height = height / 2.54, ...)
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
    oob = scales::oob_squish_any,
    ...
  )
}



# Renaming functions ------------------------------------------------------

rename_str_or_fct <- function(s, nm) {
  if (is.factor(s))
    suppressWarnings(fct_recode(s, !!!deframe(nm[c("new", "old")])))
  else
    recode(s, !!!deframe(nm))
}

rename_groups <- function(s) {
  rename_str_or_fct(
    s,
    tribble(
      ~old, ~new,
      "I", "C",
      "II", "M",
      "III", "A",
      "IV", "S"
    )
  )
}


rename_patients <- function(s) {
  rename_str_or_fct(
    s,
    tribble(
      ~old, ~new,
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
}

rename_contrast <- function(s) {
  rename_str_or_fct(
    s,
    tribble(
      ~old, ~new,
      "II_vs_I",  "M",
      "III_vs_I", "A",
      "IV_vs_I",  "S"
    )
  )
}
