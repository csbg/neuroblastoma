# cell types in high-res clusters, using marker genes from Baryawno 2020
# run with data from plot_markers.R

cell_type_order <- c(
  "CD4+ Naive",
  "CD8+ Naive",
  "CTL-1",
  "CTL-2",
  "Treg",
  "Th1/17",
  
  "NKT",
  "NK",
  
  "Pro-B",
  "Immature B cells",
  "Mature B",
  "memBcell",  
  
  "Monocytes",
  "mDC",
  "TAM",
  "TIM",
  
  "pDC",
  
  "Erythroid",
  
  "Progenitors",
  
  "NB",
  
  "Endothelial",
  "Osteoblasts",
  "Osteoclasts",    
  "Pericytes"
)

markers2 <-
  read_csv("metadata/microenvironment_markers.csv", comment = "#") %>% 
  pivot_longer(everything(), names_to = "cell_type", values_to = "gene") %>% 
  mutate(cell_type = as_factor(cell_type) %>% fct_relevel(cell_type_order)) %>% 
  arrange(cell_type) %>% 
  filter(!is.na(gene), gene %in% rownames(nb)) %>% 
  distinct(gene, .keep_all = TRUE)
  
y_annotation_data <-
  markers2 %>%
  arrange(desc(row_number())) %>%
  mutate(r = row_number()) %>%
  group_by(label = cell_type) %>%
  summarise(
    yintercept = first(r) - 0.5,
    label_y = mean(r)
  )

cluster_order <- list(
  c("1", "6", "7", "11", "12", "17", "21", "29"), # T
  c("16"), # NK
  c("14", "28"), # pro B
  c("4", "8", "9", "30"), # B
  c("27"), # mem B
  c("2", "10", "13", "24", "3", "18", "20", "22", "25"), # mono
  c("19", "31", "34"), #pDC
  c("23", "26", "33"), # ery
  c("15"), # sc
  c("5", "32") # nb
)

x_annotation_data <-
  tibble(
    n = cluster_order %>% map_int(length) %>% cumsum(),
    xmin = lag(n, default = 0) + 0.5,
    xmax = n + 0.5,
  ) %>% 
  mutate(
    fill = case_when(
      row_number() %% 2 == 0 ~ "white",
      TRUE                   ~ "gray90"
    )
  )

plot_dots(
  logcounts(nb),
  rev(markers2$gene),
  fct_relevel(colData(nb)$cluster_20, flatten_chr(cluster_order)),
  panel_annotation = x_annotation_data
) +
  geom_hline(
    yintercept = y_annotation_data$yintercept,
    linetype = "dashed",
    size = 0.25
  ) +
  geom_text(
    data = y_annotation_data,
    aes(
      x = nlevels(colData(nb)$cluster_20) + 2,
      y = label_y,
      label = label
    ),
    size = 4,
    hjust = 0
  ) +
  theme(legend.position = "left") +
  NULL
ggsave_default("wip/markers_microenvironment", height = 1500)
