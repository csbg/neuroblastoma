library(tidyverse)
library(Seurat)
library(fs)
library(patchwork)

data_dir <- "data_raw"

samples <- tibble(
  sample_file = dir_ls(
    data_dir,
    recurse = TRUE,
    regexp = "filtered_feature_bc_matrix.h5"
  ),
  sample_name = str_match(sample_file, "/[^_]*_(.*)_trans")[, 2]
) %>%
  filter(sample_name != "16_4503_Re_DOWN")
samples

nb_list <- pmap(
  samples,
  function(sample_file, sample_name) {
    message("Loading ", sample_name)
    sample_file %>% 
      Read10X_h5() %>% 
      CreateSeuratObject(sample_name) %>% 
      AddMetaData(sample_name, col.name = "sample") %>%
      PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")
  }
)

p <-
  nb_list %>% 
  map(
    ~.x@meta.data %>% 
      as_tibble(rownames = "cell") %>% 
      pivot_longer(c(nCount_RNA, nFeature_RNA, percent.mt)) %>% 
      ggplot(aes("1", value)) + 
      geom_violin(aes(fill = name), scale = "width", show.legend = FALSE) +
      geom_jitter(alpha = .02, size = 1) +
      facet_wrap(vars(name), scales = "free_y") +
      ggtitle(.x@project.name) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      )
  )

wrap_plots(p, nrow = 3)
ggsave("plots/qc_plots.png", dpi = 300, units = "mm", width = 420, height = 297)
