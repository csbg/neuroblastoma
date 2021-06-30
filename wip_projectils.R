library(ProjecTILs)
library(tidyverse)
library(patchwork)
source("common_functions.R")



# Load data ---------------------------------------------------------------

ref <- readRDS("data_wip/ref_TILAtlas_mouse_v1.rds")

state_colors <-
  c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB",
    "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd") %>% 
  set_names(levels(ref$functional.cluster))

nb <- readRDS("data_generated/rna_integrated_monocle.rds")

nb_metadata <- readRDS("data_generated/metadata.rds")

nb_T <-
  CreateSeuratObject(
    counts(nb),
    project = "NB",
    meta.data = nb_metadata %>% column_to_rownames("cell")
  ) %>% 
  subset(subset = cellont_abbr == "T")
nb_T



# Project cells -----------------------------------------------------------

query.projected <- make.projection(
  nb_T,
  ref = ref,
  filter.cells = FALSE,
  ncores = 6
)
query.projected <- cellstate.predict(ref, query.projected)


plot.projection(ref, query.projected)
plot.statepred.composition(ref, query.projected, metric = "Percent")
plot.states.radar(ref, query.projected, min.cells = 20)

saveRDS(query.projected, "data_wip/rna_projectils.rds")



# Compare groups ----------------------------------------------------------

query.list <- SplitObject(query.projected, "group")

abundance_fold_change <- 
  query.projected@meta.data %>% 
  as_tibble(rownames = "cell") %>% 
  count(group, functional.cluster) %>% 
  pivot_wider(names_from = group, values_from = n) %>%
  mutate(
    across(!functional.cluster, ~. / sum(.) * 100),
    II_vs_I = II / I,
    III_vs_I = III / I,
    IV_vs_I = IV / I,
  )

plot_abundance <- function(group) {
  ggplot(abundance_fold_change, aes(functional.cluster, {{group}})) +
    geom_col(aes(fill = functional.cluster), show.legend = FALSE) +
    scale_fill_manual(values = state_colors) +
    xlab(NULL) +
    ylab("abundance (%)") +
    coord_flip()
}

plot_fold_change <- function(group) {
  ggplot(abundance_fold_change, aes(functional.cluster, {{group}})) +
    geom_col(aes(fill = functional.cluster), show.legend = FALSE) +
    scale_fill_manual(values = state_colors) +
    xlab(NULL) +
    ylab("fold change") +
    geom_hline(yintercept = 1) +
    scale_y_continuous(trans = "log2") +
    coord_flip()
}


wrap_plots(
  plot.projection(ref, query.list$I) + ggtitle("I"),
  plot_abundance(I),
  plot.projection(ref, query.list$II) + ggtitle("II"),
  plot_fold_change(II_vs_I),
  plot.projection(ref, query.list$III) + ggtitle("III"),
  plot_fold_change(III_vs_I),
  plot.projection(ref, query.list$IV) + ggtitle("IV"),
  plot_fold_change(IV_vs_I),
  ncol = 2,
  widths = c(6, 4)
)
ggsave_default("wip/projectils_subtypes", width = 400, height = 600)


plot.states.radar(
  ref,
  query.list,
  genes4radar = c("Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7", "Gzmb", "Pdcd1",
                  "Havcr2", "Tox", "Entpd1", "Cxcr5", "Ifng", "Cxcl13",
                  "Xcl1", "Itgae"),
  return.as.list = TRUE
) %>% 
  wrap_plots(guides = "collect")
ggsave_default("wip/projectils_radar")


Hs2Mm.convert.table %>% 
  as_tibble() %>% 
  filter(Gene.MM %in% c("Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7", "Gzmb", "Pdcd1",
                        "Havcr2", "Tox", "Entpd1", "Cxcr5", "Ifng", "Cxcl13",
                        "Xcl1", "Itgae"))
