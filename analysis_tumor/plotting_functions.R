# UMAPs -------------------------------------------------------------------

## Group coloring ----

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


## Unaligned colored by aligned clusters ----

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



# Enrichment --------------------------------------------------------------


## GSEA Enrichments ----

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




# Violin plots of selected genes ------------------------------------------

## Variant 1 ----

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



## Variant 2 ----

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



# Cell cycle bar chart ----------------------------------------------------

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



# Violin plot & boxplot cytoTRACE score for each cluster ----

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



# Violin plot of cell metadata for each cluster ----

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



# Enrichment dotplot function adapted from W.E.S. ----

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
