# styles for publication-quality figures

# Common definitions ------------------------------------------------------

cell_type_abbreviations <- c(
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
cell_type_colors <- c(
  "T cell" = "#1f78b4",
  "natural killer cell" = "#a6cee3",
  "B cell" = "#33a02c",
  "monocyte" = "#ff7f00",
  "plasmacytoid dendritic cell" = "#fdbf6f",
  "erythroid lineage cell" = "#b15928",
  "hematopoietic precursor cell" = "#6a3d9a",
  "neuron" = "#e31a1c",
  "other" = "black",
  na = "gray80"
)

# https://wesandersonpalettes.tumblr.com/post/110716093015/ash-should-we-dance
group_colors <- c(
  I = "#b9b09f",
  II = "#a9d8c8",
  III = "#d16b54",
  IV = "#e8c95d"
)



# ggplot functions --------------------------------------------------------

theme_nb <- function(base_size = 5,
                     grid = TRUE,
                     rotate = FALSE,
                     ...){
  res <-
    theme_bw(...) + 
    theme(
      axis.text = element_text(colour = "black", size = base_size),
      legend.text = element_text(colour = "black", size = base_size), 
      strip.text = element_text(colour = "black", size = base_size),
      strip.background = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.margin = margin(),
      legend.box.spacing = unit(c(2, 0, 0, 0), "mm")
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
                    dpi = 300, units = "cm", limitsize = FALSE,
                    width = width, height = height, ...)  
  } else {
    rlang::exec(type, filename, res = 300, units = "cm",
                width = width, height = height,  ...)
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
    palette = "PiYG",
    direction = 1,
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
      "IV", "R"
    )
  )
}

rename_groups(c("I", "II", "III", "IV"))
rename_groups(factor(c("I", "II", "III", "IV")))


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
      "2006_2684", "R1",
      "2018_1625", "R2",
      "2018_6056", "R3",
      "2019_2495", "R4",
      "2020_1667", "R5"
    )
    # read_csv("metadata/sample_groups.csv", comment = "#") %>% 
    #   distinct(group, old = sample) %>% 
    #   arrange(group, old) %>% 
    #   mutate(group = rename_groups(group)) %>% 
    #   group_by(group) %>% 
    #   mutate(n = row_number()) %>% 
    #   unite(group, n, col = "new", sep = "") %>% 
    #   select(old, new)
  )
}

rename_patients(c("2014_0102", "2020_1667", "2019_5022", "2005_1702"))
rename_patients(factor(c("2014_0102", "2020_1667", "2019_5022", "2005_1702")))



read_csv("metadata/sample_groups.csv", comment = "#") %>% 
  distinct(group, old = sample) %>% 
  arrange(group, old) %>% 
  mutate(group = rename_groups(group)) %>% 
  group_by(group) %>% 
  mutate(n = row_number()) %>% 
  unite(group, n, col = "new", sep = "") %>% 
  select(old, new)
