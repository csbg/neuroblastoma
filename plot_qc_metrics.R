library(tidyverse)
library(Seurat)
library(fs)
library(patchwork)

data_dir <- "data_raw"

samples <-
  tibble(
    file = dir_ls(
      data_dir,
      recurse = TRUE,
      regexp = "filtered_feature_bc_matrix.h5"
    ),
    sample = str_match(file, "/[^_]*_(.*)_trans")[, 2]
  ) %>%
  filter(sample != "16_4503_Re_DOWN") %>% 
  left_join(
    read_csv("data_raw/metadata/sample_groups.csv", comment = "#"),
    by = "sample"
  ) %>%
  arrange(group, sample)
  
nb_list <- pmap(
  samples,
  function(file, sample, group) {
    message("Loading ", sample)
    file %>%
      Read10X_h5() %>%
      CreateSeuratObject(str_glue("{sample} ({group})")) %>%
      PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
  }
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

wrap_plots(p, nrow = 3)
ggsave("plots/qc_plots.jpg", dpi = 300, units = "mm", width = 550, height = 297)
