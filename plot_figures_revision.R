library(monocle3)
library(scuttle)
library(scran)
library(fgsea)
library(tidyverse)
library(latex2exp)
library(scico)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

nb_cds <-
  readRDS("data_generated/rna_decontaminated.rds") %>%
  logNormCounts(assay.type = "soupx_counts")

nb_metadata <- readRDS("data_generated/metadata.rds")

nb_signatures <-
  read_csv("metadata/nb_signatures.csv", comment = "#") %>% 
  group_by(cell_type) %>% 
  summarise(gene = list(gene)) %>% 
  deframe()



# GSEA of NB signatures ---------------------------------------------------

markers <- findMarkers(logcounts(nb_cds), groups = nb_metadata$cellont_abbr)


gsea <- imap_dfr(
  as.list(markers),
  function(markers, cell_type) {
    info("Calculating markers for {cell_type}")
    ranked_genes <-
      markers %>% 
      as_tibble(rownames = "gene") %>% 
      select(gene, summary.logFC) %>% 
      deframe()
    
    fgseaMultilevel(
      nb_signatures,
      ranked_genes,
      eps = 0,
      nPermSimple = 10000
    ) %>% 
      as_tibble() %>% 
      mutate(cell_type = {{cell_type}}, .before = 1)
  }
)

color_limit <- max(abs(gsea$NES), na.rm = TRUE)

gsea %>% 
  mutate(
    pathway =
      pathway %>%
      fct_recode("NCC-like" = "ncc_like") %>%
      fct_relevel("adrenergic", "noradrenergic", "mesenchymal") %>% 
      fct_rev(),
    cell_type = fct_relevel(cell_type, union("NB", names(CELL_TYPE_ABBREVIATIONS)))
  ) %>% 
  filter(cell_type != "other") %>% 
  ggplot(aes(cell_type, pathway)) +
  geom_point(aes(size = -log10(padj), color = NES)) +
  xlab("cell type") +
  ylab("gene signature") +
  scale_color_gsea(
    "normalized\nenrichment\nscore",
    limits = c(-color_limit, color_limit),
    breaks = c(-color_limit, 0, color_limit),
    labels = function(x) round(x, 2),
    guide = guide_colorbar(
      barheight = unit(10, "mm"),
      barwidth = unit(2, "mm"),
      ticks = FALSE
    )
  ) +
  scale_size_area(
    name = TeX("-log_{10} p_{adj}"),
    max_size = 3,
    limits = c(0, 20),
    breaks = c(0, 10, 20),
    labels = c("0", "10", "20 or higher"),
    oob = scales::oob_squish
  ) +
  coord_fixed() +
  theme_nb(grid = FALSE) +
  theme(
    legend.box = "horizontal",
    legend.box.just = "bottom",
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    legend.margin = margin(),
  )
ggsave_publication("R2_signature_gsea", width = 8, height = 2.2)
