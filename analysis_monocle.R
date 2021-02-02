# requires functions and packages of analysis_integration.R

nb_data <-
  nb_data %>% 
  left_join(
    read_csv(
      "data_generated/all_datasets_current/nb_metadata_monocle.csv",
      col_types = "cffdddd"
    ) %>% 
      rename_with(~paste0("monocle_", .x), .cols = !cell) %>% 
      mutate(
        across(where(is.factor), fct_inseq)
      ),
    by = "cell"
  )



# Clusters ----------------------------------------------------------------

plot_clusters_all(nb_data, monocle_UMAP_1, monocle_UMAP_2, monocle_cluster,
                  show_resolution = FALSE, filename = "monocle/clusters_all_UMAP")
plot_clusters_all(nb_data, monocle_tSNE_1, monocle_tSNE_2, monocle_cluster,
                  show_resolution = FALSE, filename = "monocle/clusters_all_tSNE")
plot_clusters_all(nb_data, monocle_UMAP_1, monocle_UMAP_2, monocle_partition,
                  show_resolution = FALSE, filename = "monocle/partitions_all_UMAP")
plot_clusters_all(nb_data, monocle_tSNE_1, monocle_tSNE_2, monocle_partition,
                  show_resolution = FALSE, filename = "monocle/partitions_all_tSNE")


plot_clusters_all(nb_data %>% arrange(monocle_UMAP_1),
                  monocle_UMAP_1, monocle_UMAP_2, sample,
                  show_resolution = FALSE, label_direct = FALSE,
                  filename = "monocle/clusters_all_UMAP_sample")

plot_clusters_all(nb_data %>% arrange(monocle_UMAP_1),
                  monocle_UMAP_1, monocle_UMAP_2, group,
                  show_resolution = FALSE, label_direct = FALSE,
                  filename = "monocle/clusters_all_UMAP_group")


plot_clusters_per_sample(nb_data, monocle_UMAP_1, monocle_UMAP_2,
                         monocle_cluster, sample,
                         nrow = 3, filename = "monocle/clusters_sample_UMAP")

plot_clusters_per_sample(nb_data %>% mutate(group2 = group), monocle_UMAP_1, monocle_UMAP_2,
                         monocle_cluster, group2,
                         nrow = 3, filename = "monocle/clusters_groups_UMAP")



# Cell types --------------------------------------------------------------

plot_clusters_all(nb_data, monocle_UMAP_1, monocle_UMAP_2, cell_type_broad_lumped,
                  label_direct = FALSE, show_resolution = FALSE,
                  filename = "monocle/celltype_all_UMAP")

plot_clusters_per_sample(nb_data, monocle_UMAP_1, monocle_UMAP_2, cell_type_broad_lumped,
                         sample, nrow = 3, show_legend = TRUE,
                         filename = "monocle/celltype_sample_UMAP")

plot_clusters_per_sample(nb_data %>% mutate(group2 = group),
                         monocle_UMAP_1, monocle_UMAP_2, cell_type_broad_lumped,
                         group2, nrow = 3, show_legend = TRUE,
                         filename = "monocle/celltype_group_UMAP")

plot_clusters_selected(nb_data, monocle_UMAP_1, monocle_UMAP_2, cell_type_broad_lumped,
                       folder = "monocle/clusters_highlighted_celltype")


plot_cvt_heatmap(nb_data, cell_type_broad, monocle_cluster,
                 filename = "monocle/cvt_heatmap_broad")
plot_cvt_heatmap(nb_data, cell_type_fine, monocle_cluster,
                 lump_n = 35, filename = "monocle/cvt_heatmap_fine")



plot_cvt_bar(nb_data, cell_type_broad, monocle_cluster,
             lump_n = 5, save_subplots = FALSE,
             width = 420, height = 297,
             filename = "monocle/cvt_bar_broad_cluster")
plot_cvt_bar(nb_data, cell_type_broad, monocle_partition,
             lump_n = 5, save_subplots = FALSE,
             filename = "monocle/cvt_bar_broad_partition")




plot_gene_program(nb_data, adrenergic, mesenchymal, monocle_cluster,
                  ncol = 8,
                  filename = "monocle/gene_programs_am_clusters")
plot_gene_program(nb_data, noradrenergic, ncc_like, monocle_cluster,
                  ncol = 8,
                  filename = "monocle/gene_programs_nn_clusters")



plot_clusters_all(nb_data, monocle_UMAP_1, monocle_UMAP_2, percent.mt,
                  label_direct = FALSE, show_resolution = FALSE,
                  color_scale = scale_color_viridis_c(),
                  filename = "monocle/qc_mtgene_umap")


plot_clusters_all(nb_data %>% arrange(cxds_score),
                  monocle_UMAP_1, monocle_UMAP_2, cxds_score,
                  label_direct = FALSE, show_resolution = FALSE,
                  color_scale = scale_color_viridis_c(),
                  filename = "monocle/qc_doublet_cxds")

plot_clusters_all(nb_data %>% arrange(bcds_score),
                  monocle_UMAP_1, monocle_UMAP_2, bcds_score,
                  label_direct = FALSE, show_resolution = FALSE,
                  color_scale = scale_color_viridis_c(),
                  filename = "monocle/qc_doublet_bcds")

plot_clusters_all(nb_data %>% arrange(hybrid_score),
                  monocle_UMAP_1, monocle_UMAP_2, hybrid_score,
                  label_direct = FALSE, show_resolution = FALSE,
                  color_scale = scale_color_viridis_c(),
                  filename = "monocle/qc_doublet_hybrid")

nb_data %>% 
  mutate(
    monocle_cluster_ordered = fct_reorder(
      monocle_cluster, cxds_score, .desc = TRUE
    )
  ) %>% 
  ggplot(aes(monocle_cluster_ordered, cxds_score)) +
  geom_boxplot(
    aes(fill = monocle_cluster),
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  NULL

nb_data %>% 
  mutate(
    monocle_cluster_ordered = fct_reorder(
      monocle_cluster, bcds_score, .desc = TRUE
    )
  ) %>% 
  ggplot(aes(monocle_cluster_ordered, bcds_score)) +
  geom_boxplot(
    aes(fill = monocle_cluster),
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  NULL

nb_data %>% 
  mutate(
    monocle_cluster_ordered = fct_reorder(
      monocle_cluster, hybrid_score, .desc = TRUE
    )
  ) %>% 
  ggplot(aes(monocle_cluster_ordered, hybrid_score)) +
  geom_boxplot(
    aes(fill = monocle_cluster),
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  NULL
