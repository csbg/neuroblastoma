# Plot QC metrics (% mitochondrial genes, feature and molecule count)

library(tidyverse)
library(Seurat)
library(fs)
library(patchwork)



# Load data ---------------------------------------------------------------

data_dir <- "data_raw"

samples <-
  tibble(
    file = dir_ls(
      data_dir,
      recurse = TRUE,
      regexp = "filtered_feature_bc_matrix.h5"
    ),
    bsf_id = str_match(file, "/(.*)_trans")[, 2]
  ) %>%
  left_join(
    read_csv("data_raw/metadata/sample_groups.csv", comment = "#"),
    by = "bsf_id"
  ) %>%
  drop_na() %>% 
  arrange(group, sample)
  
nb_list <- pmap(
  samples,
  function(file, sample, group, ...) {
    message("Loading ", sample)
    file %>%
      Read10X_h5() %>%
      CreateSeuratObject(str_glue("{sample} ({group})")) %>%
      PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
  }
)



# Features, molecules, %MT genes ------------------------------------------

cutoffs <- tribble(
  ~name, ~value,
  "% mt genes", 10,
  "features",  200,
  "features", 5000
)

p <-
  nb_list %>% 
  map(
    ~.x@meta.data %>% 
      as_tibble(rownames = "cell") %>% 
      pivot_longer(c(nCount_RNA, nFeature_RNA, percent.mt)) %>% 
      mutate(
        name = recode(name,
                      nCount_RNA = "molecules",
                      nFeature_RNA = "features",
                      percent.mt = "% mt genes")
      ) %>% 
      ggplot(aes("1", value)) + 
      geom_violin(aes(fill = name), scale = "width", show.legend = FALSE) +
      geom_jitter(alpha = .02, size = 1) +
      geom_hline(
        data = cutoffs,
        aes(yintercept = value),
        linetype = "dashed",
        color = "red"
      ) +
      facet_wrap(vars(name), scales = "free_y") +
      ylab("") +
      ggtitle(.x@project.name) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank()
      )
  )

wrap_plots(p, nrow = 4)
ggsave("plots/qc_plots.jpg", dpi = 300, units = "mm", width = 400, height = 400)
