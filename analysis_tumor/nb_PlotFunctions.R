################################################################################
# UMAP plot function 
################################################################################


plot_UMAP_cluster <- function(cds){
  
  UMAP <- reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
    as_tibble(rownames = "cell") %>%
    inner_join(colData(cds) %>%
                 as_tibble() %>%
                 mutate(cell = colnames(cds)),
               by = "cell") %>%
    select(1:9 | "cluster") %>%
    mutate(group = rename_groups(group),
           patient = rename_patients(sample)
    )
  
  
  ggplot(UMAP,
         aes(umap_1,umap_2)
  ) +
    geom_point( size = .01,
                shape = 16,
                aes(
                  color = cluster
                )) +
    xlab("UMAP1")+
    ylab("UMAP2")+
    annotate(
      geom = "text",
      x = c(-2,2.5,1,-2.7),
      y = c(2,-1.5,1,-2.7),
      #label = c("Cl.-1", "Cl.-2","Cl.-3","Cl.-4"),
      label = c("1", "2","3","4"),
      size = BASE_TEXT_SIZE_PT
      ) +
    theme_nb(grid = F) +
    scale_x_continuous(limits = c(-5,5))+
    scale_y_continuous(limits = c(-4.5,4.5))+
    scale_color_manual(
      guide = "legend",
      values = NBCLUSTER_COLORS
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(shape = 16 ,size = 1)
      )
    ) +
    theme(legend.key.height = unit(0.45,"cm"),
          legend.key.width = unit(0.45,"cm"),
          plot.margin = margin(r = 0.5,t=0.5, unit = "mm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
          )
}


################################################################################
#Group coloring
plot_UMAP_group <- function(cds){
  
UMAP <- reducedDim(cds, "UMAP") %>%
  magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
  as_tibble(rownames = "cell") %>%
  inner_join(colData(cds) %>%
               as_tibble() %>%
               mutate(cell = colnames(cds)),
             by = "cell") %>%
  select(1:9) %>%
  mutate(group = rename_groups(group),
         patient = rename_patients(sample)
  )


ggplot(UMAP,
       aes(umap_1,umap_2)
       ) +
  geom_point( size = BASE_LINE_SIZE,
              shape = 16,
    aes(
      color = group
        )) +
  xlab("UMAP1")+
  ylab("UMAP2")+
  scale_color_manual(
    guide = "legend",
    values = GROUP_COLORS[2:4]
    ) +
  guides(
    colour = guide_legend()
    ) +
  scale_x_continuous(limits = c(-5,5))+
  scale_y_continuous(limits = c(-4.5,4.5))+
  theme_nb(grid = F) +
  theme(
    legend.key.height = unit(0.45,"cm"),
    legend.key.width = unit(0.45,"cm"),
    plot.margin = margin(r = 0.5,t=0.5, unit = "mm"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
    )

}

################################################################################
#Patient coloring
plot_UMAP_patient <- function(cds){
  cds <- unaligned_cds
  UMAP <- reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
    as_tibble(rownames = "cell") %>%
    inner_join(colData(cds) %>%
                 as_tibble() %>%
                 mutate(cell = colnames(cds)),
               by = "cell") %>%
    select(1:9) %>%
    mutate(group = rename_groups(group),
           patient = 
             as.factor(
               rename_patients(sample)
             )
           ) %>%
    droplevels()
  
  
  ggplot(UMAP,
         aes(umap_1,umap_2)
  ) +
    geom_point( size = BASE_LINE_SIZE,
                shape = 16,
                aes(
                  color = patient
                )) +
    theme_nb(grid = F) +
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    scale_color_manual(
      values = PATIENT_COLORS_NB[grepl(names(PATIENT_COLORS_NB), 
                                       pattern=("^[AMS]"))] , 
      name = "Patient:",
      drop = T
      ) +
    guides(
      colour = guide_legend(
        direction = "vertical",
        title.position = "top",
        override.aes = list(shape = 16, 
                            size = BASE_TEXT_SIZE_MM
                            ),
        #nrow = 1,
        keywidth = 0.1,
        keyheight = 0.2,
        label.position = "left"
        )
    ) +
    theme(legend.key.height = unit(0.3,"cm"),
          legend.box = "vertical",
          legend.box.just = "left",
          legend.direction = "vertical",
          legend.justification = "top",
          plot.margin = margin(r = 0.5,t=0.5, unit = "mm")
          #axis.title.x = element_blank(),
          #axis.title.y = element_blank()
    )
}

################################################################################
# Coloring by patients within a specific specified group

plot_UMAP_group_patient <- 
  function(
  cds,
  group_ID
  )
    {
    cds <- cds[,grepl(paste0("(^",group_ID,"$)"), colData(cds)$group)]
    colData(cds) <- droplevels(colData(cds))
    UMAP <- reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(
      c("umap_1", "umap_2")
      ) %>% 
    as_tibble(rownames = "cell") %>%
    inner_join(
      colData(cds) %>%
        as_tibble() %>%
        mutate(cell = colnames(cds)),
      by = "cell"
    ) %>%
    select(1:9) %>%
    mutate(
      group = rename_groups(group),
      patient = rename_patients(sample)
    )
    
    ggplot(UMAP,
       aes(umap_1,umap_2)
    ) +
    geom_point(
      size = BASE_LINE_SIZE,
      shape = 16,
      aes(color = patient
    )) +
    xlab("UMAP1")+
    ylab("UMAP2")+
    geom_text(
      check_overlap = T,
      show.legend = F,
      size=BASE_TEXT_SIZE_MM,
      aes(
        label = paste("Group",group_ID %>% rename_groups()),
        x = -4,
        y = 4.25
        )
    ) +
    theme_nb(grid = F) +
    scale_color_manual(values = PATIENT_COLORS_NB, drop = T) +
    scale_x_continuous(limits = c(-5,5))+
    scale_y_continuous(limits = c(-4.5,4.5))+
    guides(
     colour = guide_legend(
                override.aes = 
                    list(
                      shape = 16, size = 1
                    )
              )
     ) +
    theme(plot.margin = margin(r = 0.5, t=0.5, unit = "mm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
          )
    
}



################################################################################
# UMAP + trajectory & cytotrace coloring

plot_trajectory <- function(cds) {
  # get UMAP coordinates of trajectory nodes
  ica_space_df <-
    cds@principal_graph_aux$UMAP$dp_mst %>%
    t() %>% 
    magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
    as_tibble(rownames = "node_id")
  
  # get edges of trajectory graph and add UMAP coordinates of their nodes
  edge_df <-
    cds@principal_graph$UMAP %>%
    igraph::as_data_frame() %>%
    left_join(
      ica_space_df %>%
        rename(
          source_umap_1 = umap_1,
          source_umap_2 = umap_2
        ),
      by = c(from = "node_id")
    ) %>%
    left_join(
      ica_space_df %>%
        rename(
          target_umap_1 = umap_1,
          target_umap_2 = umap_2
        ),
      by = c(to = "node_id")
    )
  
  colData(cds) %>% 
    as_tibble(rownames = "cell") %>% 
    left_join(
      reducedDim(cds, "UMAP") %>%
        magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
        as_tibble(rownames = "cell"),
      by = "cell"
    ) %>% 
    ggplot(aes(umap_1, umap_2)) +
    geom_point(
      aes(color = cytoTRACE),
      shape = 16,
      size = BASE_LINE_SIZE/2,
      show.legend = TRUE
    ) +
    geom_segment(
      data = edge_df,
      aes(
        x = source_umap_1, 
        y = source_umap_2,
        xend = target_umap_1, 
        yend = target_umap_2
      ),
      size = BASE_LINE_SIZE,
      color = "#525252",
      alpha = 0.75
    ) +
    theme_nb(grid = FALSE)  +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_x_continuous(limits = c(-5,5))+
    scale_y_continuous(limits = c(-4.5,4.5))+
    scale_color_PiYG(name="cytoTrace", 
                     guide = guide_colorbar(
                       barheight = unit(12,"mm"),
                       barwidth = unit(2,"mm")
                       )
                     )  +
    coord_fixed() +
    theme(legend.key.width = unit(0.45,"cm"),
          legend.key.height = unit(0.45,"cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(t=0.5,r = 0.5, unit = "mm")
          )
}

################################################################################
# Enrichment cluster/patients plotting function

plot_enrich_ClusterPatient <- function(data,
                                       circle_significant = FALSE,
                                       sig_p = 0.05) {
  
data <- 
  data %>% 
  mutate(
    is_significant = Adjusted.P.value <= sig_p
    ) 

if (circle_significant) {
  circle_significant <- geom_point(
    data = data %>% filter(is_significant),
    shape = 1,
    stroke = BASE_LINE_SIZE*2,
    alpha = 0.5
  )
} else {
  circle_significant <- NULL
}

color_limit <- log2(10)
#color_limit <- round(log2(max(abs(data$Odds.Ratio))),digits=0)/2
horizontal_grid <- NULL 

#if (nlevels(data$Patient) > 5) {
#  horizontal_grid <-
#    geom_hline(
#      yintercept = seq(5, nlevels(data$Patient), 5),
#      size = BASE_LINE_SIZE,
#      color = "grey92"
#    )
#} else {
#  horizontal_grid <- NULL  
#}

p <- 
  ggplot(data, 
         aes(cluster, 
             Patient, 
             size = -log10(Adjusted.P.value)
             )
         ) +
  scale_y_discrete() +
  horizontal_grid +    
  geom_point(
    aes(color = log2(Odds.Ratio)
        )
    ) +
  circle_significant +
  horizontal_grid+
  xlab("NB subcluster") +
  ylab(NULL) +
  scale_color_gsea(
    name = TeX("log_2 odds ratio"),
    limits = c(-3,3),
    breaks = c(-2,0,2),
    guide = guide_colorbar(
      barheight = unit(1.2, "cm"),
      barwidth = unit(2, "mm"),
      label.position = "left",
      title.vjust = 0.2,
      order = 2
    )
  ) +
  scale_size(
    name = TeX("-log_{10} p_{adj}"),
    range = c(0,2),
    breaks = c(0,1,3,10),
    guide = guide_legend(order = 1,
                         label.position = "left")
  )  +
  coord_fixed() +
  theme_nb(grid = FALSE) +
  theme(
    aspect.ratio = 3,
    legend.box = "vertical",
    legend.box.just = "left",
    legend.direction = "vertical",
    legend.justification = "right",
    legend.title.align = 1,
    legend.title = element_text(angle = 0),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(0.1, "mm"),
    legend.position = "left",
    legend.spacing = unit(0.5, "mm"),
    axis.text.x = element_text(),
    legend.margin = margin(0, 0.1, 0, -3, "mm"),
    plot.margin = margin(t=0.5, r = 0.1, unit= "mm")
  )
}

################################################################################
# Enrichment plot literature markers
plot_enrich_Literature <- function(data,
                                   circle_significant = FALSE,
                                   sig_p = 0.05) {
  data <- 
    data %>% 
    mutate(
      is_significant = Adjusted.P.value <= sig_p
    ) 
  
  if (circle_significant) {
    circle_significant <- geom_point(
      data = data %>% filter(is_significant),
      shape = 1,
      stroke = BASE_LINE_SIZE*2,
      alpha = 0.5
    )
  } else {
    circle_significant <- NULL
  }
  
  color_limit <- log2(10)
  
  if (nlevels(data$Marker) > 5) {
    horizontal_grid <-
      geom_hline(
        yintercept = seq(5, nlevels(data_vis$Term), 5),
        size = BASE_LINE_SIZE,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  
  p <- 
    ggplot(data, 
           aes(cluster, 
               Marker %>% rename_celltype(), 
               size = -log10(Adjusted.P.value)
           )
    ) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(
      aes(color = log2(Odds.Ratio)
      )
    ) +
    circle_significant +
    xlab("NB subcluster") +
    ylab(NULL) +
    scale_color_gsea(
      name = TeX("log_2 odds ratio"),
      limits = c(-3,3),
      breaks = c(-2,0,2),
      guide = guide_colorbar(
        barheight = unit(12, "mm"),
        barwidth = unit(2, "mm"),
        label.position = "left",
        title.vjust = 0.2,
        order = 2
      )
    ) +
    scale_size(
      name = TeX("-log_{10} p_{adj}"),
      range = c(0,2),
      breaks = c(0,1,3,10),
      guide = guide_legend(order = 1,
                           label.position = "left"
                           )
    )  +
    coord_fixed() +
    expand_limits(x = 0.5)+
    scale_x_discrete(limits = factor(1:4)) +
    theme_nb(grid = FALSE) +
    theme(
      aspect.ratio = 3,
      legend.box = "vertical",
      legend.box.just = "left",
      legend.direction = "vertical",
      legend.justification = "center",
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      #legend.position = c(1, 1.07),
      legend.spacing = unit(1, "mm"),
      legend.margin = margin(0, 1, 0, -3, "mm"),
      plot.margin = margin(c(0.5, 0, 0, 0), "mm")
    )
  p  
}

################################################################################
# WES GSEA Enrichments ---------------------------------------------------------

plot_gsea <- function(data,
                      #db,
                      circle_significant = FALSE,
                      max_p_adj = 0.05,
                      min_OR = 5,
                      min_overlap_size = 2) {

  
  #data <- EnrData[[Cl4]]
  #db <- "KEGG_2019_Human"
  
  data <- 
    data %>%
    separate(
      Overlap,
      into = c("overlap_size", "geneset_size"),
      convert = TRUE
    )
  
  top_terms <-
    data %>% 
    filter(
      #Term == Terms,
      #db == {{db}},
      Adjusted.P.value <= max_p_adj,
      Odds.Ratio >= min_OR,
      overlap_size >= min_overlap_size
    ) %>% 
    pull(Term) %>%
    unique()
  
  data_vis <- 
    data %>% 
    filter(
      #db == {{db}},
      Term %in% top_terms
    ) %>% 
    mutate(
      is_significant =
        Adjusted.P.value <= max_p_adj &
        Odds.Ratio >= min_OR &
        overlap_size >= min_overlap_size,
      Term =
        Term %>% 
      #  str_match("(.+) Homo") %>% 
      #  magrittr::extract(,2) %>% 
        as_factor() 
        #reorder(cluster)
      #  relevel(ref = "DNA Repair") %>%
       # relevel(ref = "DNA repair") %>%
        #relevel(ref = "DNA replication") %>%
        #relevel(ref = "Cell cycle")
    ) %>%
    mutate(Term = fct_reorder(Term,log2(Odds.Ratio),  max)) %>%
    mutate(Term = fct_reorder(Term, cluster,  max))
  
  # return(data_vis)
  
  if (nlevels(data_vis$Term) > 5) {
    horizontal_grid <-
      geom_hline(
        yintercept = seq(5, nlevels(data_vis$Term), 5),
        size = BASE_LINE_SIZE,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  
  color_limit <- log2(64)#max(log2(data_vis$Odds.Ratio)) # max(abs(data_vis$NES)) #log2(10)  # 
  
  if (circle_significant) {
    circle_significant <- geom_point(
      data = data_vis %>% filter(is_significant),
      shape = 1,
      stroke = BASE_LINE_SIZE*2,
      alpha = 0.5
      )
  } else {
    circle_significant <- NULL
  }
  
  p <- 
    ggplot(data_vis, aes(cluster, Term, size = -log10(Adjusted.P.value))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = log2(Odds.Ratio))) +
    circle_significant +
    xlab("NB subcluster") +
    ylab(NULL) +
    scale_color_gsea_2(
      name = TeX("log_2 odds ratio"),
      breaks = c(0,2,4,6),
      limit = 
        c(0, color_limit),
      guide = guide_colorbar(
        barheight = unit(1.2, "cm"),
        barwidth = unit(2, "mm"),
        title.vjust = 0.1,
        label.position = "left",
        order = 2
      )
    ) +
    scale_size(
      name = TeX("-log_{10} p_{adj}"),
      range = c(0,2),
      breaks = c(0.1,1.3,3,10),
      guide = guide_legend(order = 1,
                           label.position = "left")
    )  +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      legend.box = "vertical",
      legend.box.just = "left",
      legend.direction = "vertical",
      legend.justification = "left",
      legend.title = element_text(size= BASE_TEXT_SIZE_PT, angle = 0),
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(0.5, "mm"),
      legend.position = "left",
      legend.spacing = unit(2, "mm"),
      legend.margin = margin(0, 0, 0, -3, "mm"),
      panel.spacing = unit(-.1, "pt"),
      plot.margin = margin(t = .5,r=0.1, unit = "mm")
    )
  
  p
}

################################################################################
# Heatmap of differentially expressed genes


plot_DEG_heatmap <- function(data_DEG,
                             data_CDS,
                             legend = F
                             ){
  genes <- 
    data_DEG$gene_short_name %>% 
    unique()
  
  matrix <-
    exprs(data_CDS)
  
  cell_group <- 
    tibble(
      cell = rownames(colData(data_CDS)),
      cell_group = clusters(data_CDS)[colnames(data_CDS)]
      )
  
  gene_group <- 
    tibble(
      id = genes,
      gene_group = genes
      )  
  
  EXPR <- 
    aggregate_gene_expression(
      data_CDS, 
      gene_group,
      cell_group,
      norm_method = "log" ,
      scale_agg_values = F
      ) %>% 
    as.matrix()
    
  EXPR_scale <- t(scale(t(EXPR)))
  
  colnames(EXPR_scale) <- colnames(EXPR)
  

 H =
   Heatmap(EXPR_scale,
           col = color_heat,
           cluster_rows = T,
           show_row_dend = F,
           show_column_dend = T,
           show_column_names = T,
           show_row_names = F,
           #column_order = c("1","2","3","4"),
           column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
           column_names_rot = 0,
           column_dend_height = unit(3,"pt"),
           column_dend_gp = gpar(lwd = BASE_LINE_SIZE),
           clustering_distance_rows = "pearson",
           clustering_distance_columns = "pearson",
           clustering_method_rows = "ward.D2",
           clustering_method_columns = "ward.D2",
           column_title = "NB subcluster",
           column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
           column_title_side = "bottom",
           column_dend_side = "bottom",
           show_heatmap_legend = legend, 
           heatmap_legend_param = list(
             title = latex2exp::TeX("log_{10} expr_{scaled}"),
             direction = "vertical",
             labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
             grid_height = unit(1,"mm"),
             grid_width = unit(1,"mm"),
             legend_height = unit(1.2, "cm"),
             legend_width = unit(2, "mm"),
             title_position = "topcenter"
             )
           ) 
           
    
  
H + ht_opt(
  legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  DENDROGRAM_PADDING = unit(0.5, "mm"),
  TITLE_PADDING = unit(c(0.5,0.5), "mm"),
  DIMNAME_PADDING = unit(0.5, "mm"),
  HEATMAP_LEGEND_PADDING = unit(4,"mm")
  )
  
draw(H, padding = unit(c(0.1, 4, 0.5, 0.5), "mm"), heatmap_legend_side = "left")

}


################################################################################
# Violinplots interesting genes

plot_violin <- function(data_cds,
                        genes
                        ){

  data_cds <- 
    logNormCounts(data_cds, 
                  size.factors = colData(data_cds)$Size_Factor
    )
  
  colData(data_cds)$size_factor <- colData(data_cds)$Size_Factor %>% as.numeric()
  
  plotExpression(data_cds,
                 x = "cluster",
                 features = genes$gene %>% unique(),
                 log2_values = T,
                 colour_by = "cluster",
                 ncol = 4,
                 scales = "free",
                 xlab = "NB Subcluster",
                 one_facet = T,
                 size_by = "Size_Factor",
                 show_median = T
                 ) + 
    theme_nb(grid=F)+
    theme(rect = element_rect(size = unit(0.1,"pt")),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position="none"
          ) +
    scale_colour_manual(name = "Cluster",
                        values = NBCLUSTER_COLORS)+
    scale_size_area(max_size = 0.5) +
    guides(size = "none",
           colour = guide_legend()
           )
    
}

################################################################################
# Violin plots variant 2

plot_violin2 <- function(data_cds,
                         genes,
                         cluster_number){
  
  genes <- genes[[cluster_number]]$gene
  
  data_cds <- scuttle::logNormCounts(
                            data_cds,
                            log = T,
                            transform = "log",
                            size.factors = colData(data_cds)$Size_Factor,
                            pseudo.count = 1,
                            subset.row = genes
                            )

  data <- 
    logcounts(data_cds)[genes, ,drop = FALSE] %>%
    t() %>%
    as.matrix() %>%
    as_tibble(rownames = "cell") %>%
    pivot_longer(!cell, 
                 names_to = "gene", 
                 values_to = "logexp") %>%
    left_join(
      data_cds@colData %>%
        as_tibble(rownames = "cell") %>% 
        select(cell, group, sample, cluster, Size_Factor),
      by = "cell"
    ) %>% 
    group_by(gene) %>% 
    mutate(logexp = logexp / max(logexp)) %>% 
    ungroup()
  
  ggplot(data,
         aes(cluster, logexp)
         ) +
     geom_violin(
      aes(fill = cluster,
          color = NULL),
      scale = "width",
      size = BASE_LINE_SIZE,
      alpha = 0.65,
      trim = T
    ) +
    geom_boxplot(notch = T,
                 aes(fill = NA,
                     size = 1),
                 size = BASE_LINE_SIZE,
                 outlier.shape = NA,
                 width = 1
    ) +
    xlab("NB subcluster") +
    ylab(NULL)+
    scale_y_continuous(
      name = TeX("log_2 expression_{norm.}"),
      limits =  c(0, 1),
      breaks = c(0, 1),
      expand = expansion(add = .025)
    ) +
    scale_fill_manual(values = NBCLUSTER_COLORS) +
    #scale_size_area(max_size = 0.00001)+
    coord_flip() +
    facet_grid(~gene,scales = "free") +
    theme_nb(grid = FALSE) +
    theme(legend.position = "none",
          #axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
          )
}


################################################################################
#UMAP unaligned colored by aligned clusters
plot_UMAP_alignedcluster <- function(cds){
  
  UMAP <- reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
    as_tibble(rownames = "cell") %>%
    inner_join(colData(cds) %>%
                 as_tibble() %>%
                 mutate(cell = colnames(cds)),
               by = "cell") %>%
    select(1:9 | "aligned_cluster") %>%
    mutate(group = rename_groups(group),
           patient = rename_patients(sample)
    )
  
  
  ggplot(UMAP,
         aes(umap_1,umap_2)
  ) +
    geom_point( size = .01,
                shape = 16,
                aes(
                  color = aligned_cluster
                )) +
    xlab("UMAP1")+
    ylab("UMAP2")+
    theme_nb(grid = F) +
    scale_color_manual(
      guide = "legend",
      values = NBCLUSTER_COLORS
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(shape = 16 ,size = 1)
      )
    ) +
    theme(legend.key.height = unit(0.45,"cm"),
          legend.key.width = unit(0.45,"cm"),
          plot.margin = margin(r = 0.5,t=0.5, unit = "mm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    )
}


################################################################################
# Cellcycle bargraph

plot_CC_patient <- function(cds, labs) {
  
  labs <- if (labs == T) {
    xlab("Patients")
  } 
  else {
    xlab(NULL)
  }
    
  
  meta <- colData(cds) %>% 
    as_tibble() %>%
    select("group","phase","sample") %>%
    mutate("patient" = rename_patients(sample),
           "group" = rename_groups(group))
  
  ggplot(meta,
         aes(as.factor(patient))
         ) +
    geom_histogram(stat = "count",aes(fill=patient))+
    facet_grid(rows=vars(phase), scales = "free", space="free_x" )+
    labs+
    ylab(NULL)+
    scale_fill_manual(name = "Patients",values = PATIENT_COLORS_NB)+
    guides(
      fill= guide_legend(
        title.position = "top",
        override.aes = list(shape = 16, 
                            size = BASE_LINE_SIZE
        ),
        keywidth = 0.1,
        keyheight = 0.2,
        label.position = "top"
        )
      )+
    theme_nb(grid = F)+
    theme(axis.text.y = element_blank(),
          legend.position = "none",
          panel.spacing = unit(-.5, "pt"),
          plot.margin = margin(t = .5,r=0.5,l=-0.5, unit = "mm")
          )
  
}

################################################################################
# Violin plot & boxplot cytoTRACE score for each cluster

plot_Violin_cyto <- function(cds ,labs){

  labs_x <- if (labs == T) {
    xlab("NB subcluster")
  } 
  else {
    xlab(NULL)
  }
  
  labs_y <- if (labs == T) {
    ylab("CytoTRACE score") 
  } 
  else {
    ylab(NULL)
  } 
 
  meta <- colData(cds) %>% 
    as_tibble() %>%
    select("group","cluster","cytoTRACE","sample") %>%
    mutate("patient" = rename_patients(sample),
           "group" = rename_groups(group))
  
  ggplot(meta,
         aes(cluster,cytoTRACE)
         )+
    geom_violin(scale = "width",
                trim = T,
                aes(fill=cluster))+
    geom_boxplot(notch = T,
                 alpha=0.1,
                 size=0.5,
                 varwidth= T,
                 aes(fill=cluster))+
    scale_fill_manual(values = NBCLUSTER_COLORS)+
    theme_nb(grid = F)+
    theme(legend.position = "none",
          plot.margin = margin(t = .5,r=0.5,l=0.1, unit = "mm")
          )+
    labs_x +
    labs_y
}

################################################################################
# Violin plot of cell-metadata for each cluster

plot_Violin_cluster <- function(cds ){

  meta <- colData(aligned_cds) %>% 
    as_tibble() %>%
    select("cluster","percent_mito") 

  ggplot(meta, 
         aes(cluster, percent_mito)
  ) +
    geom_violin(scale = "width",
                trim = T,
                aes(fill = cluster))+
    geom_boxplot(notch = T,
                 alpha=0.1,
                 size=0.5,
                 varwidth= T,
                 aes(fill=cluster))+
    xlab("NB subcluster")+
    ylab("% Mitochondrial Reads")+
    scale_fill_manual(values = NBCLUSTER_COLORS)+
    theme_nb(grid = F)+
    theme(legend.position = "none",
          plot.margin = margin(t = .5,r=0.5,l=0.1, unit = "mm")
    )
  
}

################################################################################
# Enrichment dotplot function adapted from W.E.S.
plot_enrichr_dots <- function(data,
                              db,
                              max_p_adj = 0.05,
                              min_odds_ratio = 10,
                              min_overlap_size = 2,
                              log_odds_cap = 2,
                              direction = "up",
                              filename = NULL,
                              guide_show = F,
                              ...) {
  
  data$cluster <- as.factor(data$cluster)
  
  data_filtered <-
    data %>% 
    filter(
      db == {{db}},
      direction == {{direction}},
      overlap_size >= min_overlap_size
    ) %>% 
    group_by(Term) %>%
    filter(
      min(Adjusted.P.value) < max_p_adj,
      max(Odds.Ratio) > min_odds_ratio
    ) %>%
    ungroup() %>%
    select(contrast, cluster, Term, Odds.Ratio, Adjusted.P.value) %>% 
    mutate(
      Term = as_factor(Term) %>% fct_reorder(Odds.Ratio, max),
      is_significant =
        Adjusted.P.value < max_p_adj & Odds.Ratio > min_odds_ratio
    ) %>% 
    mutate(Odds.Ratio = pmin(Odds.Ratio, 10^log_odds_cap))
  
  p <- ggplot(
    data_filtered, 
    aes(cluster, 
        Term, 
        size = -log10(Adjusted.P.value)
    )
  ) +
    geom_point(aes(color = log10(Odds.Ratio))) +
    geom_point(data = data_filtered %>% filter(is_significant), shape = 1) +
    coord_fixed() +
    #labs(
    #  y = "",
    #  title = str_glue( "Enrichr results ({db})"),
    #  caption = str_glue(
    #    "minimum overlap size: {min_overlap_size}; ",
    #    "bordered circles: adjusted p value < {max_p_adj}, ",
    #    "odds ratio > {min_odds_ratio}; ",
    #    "color scale capped at {log_odds_cap}"
    #  )
    #) +
    scale_color_gsea(guide="legend") +
    scale_size_continuous(guide="legend") +
    guides(
      color = guide_legend(title = TeX("log_{10} (odds ratio)")),
      size = guide_legend(title = TeX("log_{10} (padj)")
      )
    ) +
    theme_nb(grid = F)
  
  ggsave_publication(
    filename,
    plot = p,
    type = "pdf",
    legends = guide_show,
    height = 11,
    width = 11
  )
}






