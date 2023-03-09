
# NB UMAP -----------------------------------------------------------------

umap <- 
  DimPlot(
    object = nblast, 
    label = TRUE, 
    label.size = 3, 
    pt.size = 0.2,
    #cols = CELL_TYPE_COLORS,
    #group.by = 'new.ident',
    reduction = "umap_scopen") &
  #scale_colour_brewer(palette = "Set1") &
  #hrbrthemes::theme_ipsum() &
  theme_nb() &
  scale_color_manual(
    name = "cell type",
    values = CELL_TYPE_COLORS,
    labels =
      str_glue(
        "{names(CELL_TYPE_ABBREVIATIONS)} ({CELL_TYPE_ABBREVIATIONS})"
      )) 
#theme(legend.position = "none")
#guides(color = guide_legend(override.aes = list(size=10), ncol=1) ) &
#theme(legend.text=element_text(size=10))

ggsave_publication('umap_wes_theme_nb_annot', 
                   width = 9, height= 6)

DimPlot(
  object = nblast, 
  label = TRUE, 
  group.by = 'seurat_clusters',
  reduction = "umap_scopen") &
  theme_nb() &
  labs(title = "Seurat clusters") &
  umap &
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank())
ggsave_publication('seurat_umap', 
                   width = 8, height= 7)

# scale_colour_brewer(palette = "Set1") &
# hrbrthemes::theme_ipsum()

cowplot::plot_grid(umap,  cov_plot1 / cov_plot2 ) 
(umap |cov_plot1 | cov_plot2 / (p3_tnfa /p3_e2f))

umap + plot_spacer() + 
  (cov_plot1 / cov_plot2) + 
  plot_spacer() + (p3_tnfa /p3_e2f) + plot_spacer()

ggpubr::ggarrange(umap, cov_plot1, theme_void(), theme_void(),
                  theme_void(), cov_plot2, theme_void(),theme_void(),
                  p3_tnfa , theme_void(), theme_void(),theme_void(),
                  p3_e2f,
                  #labels = c("A", "B", "C"),
                  ncol = 4, nrow = 4)
multiplot(umap,
          cov_plot1,cov_plot2,
          p3_tnfa,p3_e2f,ncol=4) 
multiplot(all_promoters, all_others, ncol=2)
ggpubr::ggarrange(all_promoters, all_others)

all_promoters
ggsave_publication('all_promoters_motif', 
                   width = 10, height= 15)


# Annotation with UCell ---------------------------------------------------
# annotate clusters
library(UCell)
load("/media/AGFORTELNY/people/rohit/projects/markers.RData")


nblast <- AddModuleScore_UCell(nblast, features = baryawno_list, 
                               assay = 'RNA', slot = "counts")
baryawno_list.markers <- paste0(names(baryawno_list), "_UCell")
VlnPlot(nblast, features = baryawno_list.markers[11], log = TRUE, pt.size = 0) &
  stat_summary(fun = median, geom='point', size = 2, colour = "black") #&
#scale_y_continuous(limits = c(1,1.015)) &
#theme_nb()
pdf(file = 'plots/violin/ucell_baryawno.pdf',   # The directory you want to save the file in
    width = 8, 
    height = 4)
for (i in 1:length(baryawno_list.markers)){
  print(VlnPlot(nblast, features = baryawno_list.markers[i], 
                log = TRUE, pt.size = 0) &
          stat_summary(fun = median, geom='point', 
                       size = 2, colour = "black"))
  
}
dev.off()
# my_levels <- c(4,3,2,1)

# Relevel object@ident
# Idents(nb) <- factor(Idents(nb), levels = ord)
nblast <- AddModuleScore_UCell(nblast, features = markers, 
                               assay = 'RNA', slot = "counts")
signature.markers <- paste0(names(markers), "_UCell")
# Signature markers

(t_ucell <-
    VlnPlot(nblast, features = signature.markers[1], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    #scale_y_continuous(limits = c(1,1.015)) &
    theme_nb())
(nkt_ucell <-
    VlnPlot(nblast, features = signature.markers[2], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    theme_nb())
(b_ucell <- VlnPlot(nblast, features = signature.markers[3], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    theme_nb())
(monocytes_ucell <- 
    VlnPlot(nblast, features = signature.markers[4], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    theme_nb())
(erythroid_ucell <-
    VlnPlot(nblast, features = signature.markers[5], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    theme_nb())
(nblast_ucell <-
    VlnPlot(nblast, features = signature.markers[6], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    theme_nb())
(stem_ucell<-
    VlnPlot(nblast, features = signature.markers[7], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    theme_nb())
(nk_ucell<-
    VlnPlot(nblast, features = signature.markers[9], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    theme_nb())

(nblast_ucell/ monocytes_ucell/ b_ucell/ erythroid_ucell/ nkt_ucell/ t_ucell ) +
  plot_layout(guides = "collect")
ggsave_publication('ucell_marker_violin', 
                   width = 19, height= 28)

#scale_y_continuous(limits = c(1,1.015)) &
# theme(legend.position = 'none',
#       text = element_text(size=18, face = "bold"),
#       axis.text.x = element_text(angle=0,  hjust = 1, size=12, color="black" ),
#       axis.text.y = element_text( size=12, color="black" ),
#       axis.title    = element_text( size=15, face = "bold" ),
#       legend.text = element_text( size=40, color="black"),
#       legend.title = element_text( size=15, color="black"),
#       axis.line.x = element_line(size=1.5, color="black"),
#       axis.line.y = element_line(size=1.5, color="black"))

Tcell_signature = c("CD2","CD3E","CD3D")
Myeloid_signature = c("SPI1","FCER1G","CSF1R")

Seurat::FeaturePlot(nblast,     
                    reduction = "umap_scopen", 
                    features = signature.markers[1:6], 
                    ncol = 3, 
                    cols = c("grey90", "magenta"),
                    order = T) &
  scale_colour_gradientn(
    colours = rev(brewer.pal(n = 5, 
                             name = "RdYlGn"))) 

Seurat::FeaturePlot(nblast,     
                    reduction = "umap_scopen", 
                    features = signature.markers[1:9], 
                    ncol = 3, 
                    #cols = rev(Blue2Gray8Steps),
                    cols = c("grey90", "#ef8a62", "#2166ac"),
                    order = T) &
  scale_fill_manual(values=ModifiedSpectralScheme11Steps) +
  (sm1 <- Seurat::FeaturePlot(nblast,     
                              reduction = "umap_scopen", 
                              features = signature.markers[1],
                              cols = c("grey90", "#ef8a62", "#2166ac"),
                              order = T) & 
     labs(title = "T cell markers") &
     theme_nb() &
     #theme(legend.position="none") &
     # scale_colour_gradientn(colours = rev(brewer.pal(
     #   n = 5,
     #   name = "RdYlGn"))) &
     theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
           legend.key.height = unit(0.3, 'cm'), #change legend key height
           legend.key.width = unit(0.3, 'cm'), #change legend key width
           legend.title = element_text(size=5), #change legend title font size
           legend.text = element_text(size=5)) )

(sm2 <- Seurat::FeaturePlot(nblast,     
                            reduction = "umap_scopen", 
                            features = signature.markers[2], 
                            cols = c("grey90", "#ef8a62", "#2166ac"),
                            order = T) &
    labs(title = "NKT cell markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) & 
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm3 <- Seurat::FeaturePlot(nblast,     
                            reduction = "umap_scopen", 
                            features = signature.markers[3],
                            cols = c("grey90", "#ef8a62", "#2166ac"),
                            order = T) &
    labs(title = "B cell markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm4 <- Seurat::FeaturePlot(nblast,     
                            reduction = "umap_scopen", 
                            features = signature.markers[4],
                            cols = c("grey90", "#ef8a62", "#2166ac"),
                            order = T) &
    labs(title = "Myeloid markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm5 <- Seurat::FeaturePlot(nblast,     
                            reduction = "umap_scopen", 
                            features = signature.markers[5], 
                            cols = c("grey90", "#ef8a62", "#2166ac"),
                            order = T) &
    labs(title = "Erythroid markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm6 <- Seurat::FeaturePlot(nblast,     
                            reduction = "umap_scopen", 
                            features = signature.markers[6],
                            cols = c("grey90", "#ef8a62", "#2166ac"),
                            order = T) &
    labs(title = "NB markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 

(sm8 <- Seurat::FeaturePlot(nblast,     
                            reduction = "umap_scopen", 
                            features = signature.markers[8],
                            cols = c("grey90", "#ef8a62", "#2166ac"),
                            order = T) &
    labs(title = "NK markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm9 <- Seurat::FeaturePlot(nblast,     
                            reduction = "umap_scopen", 
                            features = signature.markers[9],
                            cols = c("grey90", "#ef8a62", "#2166ac"),
                            order = T) &
    labs(title = "pDC markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 

(sm1 | sm8   | sm3) /(sm4 |sm5 |sm6 ) +
  plot_annotation(tag_levels = 'A')
ggsave_publication('umap_signatures2', 
                   width = 18, height= 9)

ggsave_publication('umap_signatures', 
                   width = 18, height= 9)


# Footprinting ------------------------------------------------------------

fp1 <- PlotFootprint(nblast, features = 'JUNB')
fp2 <- PlotFootprint(nblast, features = 'NFE2')
fp3 <- PlotFootprint(nblast, features = 'ELF5')
fp4 <- PlotFootprint(nblast, features = 'SPIB')
fp5 <- PlotFootprint(nblast, features = 'PHOX2B')
fp6 <- PlotFootprint(nblast, features = 'PHOX2A')
(fp1 | fp2 |fp3 ) /(fp4 |fp5 |fp6) &
  plot_layout(guides = 'collect') &
  theme_nb()


prow <- cowplot::plot_grid(
  fp1 & theme(legend.position="none"),
  fp2 & theme(legend.position="none"),
  fp3 ,
  fp4 & theme(legend.position="none"),
  fp5 & theme(legend.position="none"),
  fp6 ,
  align = 'vh',
  labels = c("A", "B", "C", "D", "E", "F"),
  hjust = -1,
  nrow = 2,
  ncol = 3
)
prow

fp_1 <- PlotFootprint(nblast, features = 'JUNB') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells") 
fp_2 <- PlotFootprint(nblast, features = 'NFE2') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")
fp_3 <- PlotFootprint(nblast, features = 'ELF5') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")
fp_4 <- PlotFootprint(nblast, features = 'SPIB') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")
fp_5 <- PlotFootprint(nblast, features = 'PHOX2B') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")
fp_6 <- PlotFootprint(nblast, features = 'PHOX2A') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")

prow2 <- cowplot::plot_grid(
  fp_1 & theme(legend.position="none") ,
  fp_2 & theme(legend.position="none") ,
  fp_3 ,
  fp_4 & theme(legend.position="none") ,
  fp_5 & theme(legend.position="none") ,
  fp_6 ,
  align = 'vh',
  labels = c("A", "B", "C", "D", "E", "F"),
  hjust = -1,
  nrow = 2,
  ncol = 3
)
prow2
ggsave_publication('cell_type_footprinting', 
                   width = 30, height= 15)

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  fp1 & theme(legend.box.margin = margin(0, 0, 0, 12))
)

cowplot::plot_grid(prow, legend, rel_widths = c(3, .4))

fp1 & 
  guides(fill="none") &
  theme(legend.position="none")
PlotFootprint(nblast, 
              features = c('JUNB', 'NFE2', 'ELF5', 'SPIB', 'PHOX2B','PHOX2A')) &
  scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells") &
  plot_layout(guides = 'collect') &
  theme_nb()


# plot the footprint data for each group of cells
p2 <- PlotFootprint(nblast, features = c("GATA2", "CEBPA", "EBF1"))
p2 + patchwork::plot_layout(ncol = 1)


p1 <- PlotFootprint(nblast, features =  moi[1:3])
p1 + patchwork::plot_layout(ncol = 3)


p3 <- PlotFootprint(nblast, features =  moi[1:3])
p3 + patchwork::plot_layout(ncol = 1)

p4 <- PlotFootprint(nblast, features =  moi[8:10])
p4 + patchwork::plot_layout(ncol = 1)


p5 <- PlotFootprint(nblast, features =  moi[11])
p5 + patchwork::plot_layout(ncol = 1)

names(table(nblast@active.ident))


# Coverage plots ----------------------------------------------------------

patients_col <- c(
  "I-2014_0102" = "navy" ,  "I-2018_4252" = "navy",  "I-2020_1288" = "navy",  
  "I-MF244_GNB" = "navy",   "I-MF244_GNM" = "navy", "II-16_4503"   = "royalblue",
  "II-2018_1404" = "royalblue", "II-2019_5022" = "royalblue", "II-2019_5754" = "royalblue", 
  "III-2005_1702" = "red",  "III-2016_3924" = "red", "IV-2006_2684"  = "seagreen",
  "IV-2018_1625" = "seagreen", "IV-2018_6056" = "seagreen", "IV-2019_2495" = "seagreen", 
  "IV-2020_1667" = "seagreen")

# add patient group information
groups <- list(c('2', '3', '6', '15', '16'),
               c('1', '4', '5', '11'),
               c('7','12'),
               c('8', '9', '10','13', '14'))

nblast@meta.data$group <- nblast@meta.data %>% 
  dplyr::mutate(
    group = dplyr::case_when(orig.ident %in% groups[[1]] ~"C",
                             orig.ident %in% groups[[2]] ~"M", 
                             orig.ident %in% groups[[3]] ~"A",
                             orig.ident %in% groups[[4]] ~"S")) %>% 
  pull(group)
nblast$group

group_col <- c('C' = 'navy',
               'M' = 'royalblue',
               'A' = 'red',
               'S' = 'seagreen')

fragment.path <- '/media/AGFORTELNY/PROJECTS/Neuroblastoma/data_raw/bsf/COUNT/AGGR_ALL_ATAC/fragments.tsv.gz'
nblast <- Seurat::SetFragments(object = nblast,file = fragment.path)

no_gene <-  c("chr12-92603948-92604523", 
              "chr12-92551861-92552230",
              "chr2-152200475-152200993",
              "chr1-66464406-66465145",
              "chr12-92559040-92559218",
              "chr17-34203317-34203629")
nfkb1 <- 'chr4-102521897-102522522'
manba <- "chr4-102627693-102628315"
itpr2 <- "chr12-26334745-26335048"
hipk2 <- "chr7-139748889-139749340"
cdc37 <- "chr19-10386069-10386429"
oxsr1 <- "chr3-38185669-38186114"

# tile_plot <- TilePlot(
#   object = nblast,
#   group.by = 'group',
#   region = 'chr4-102521897-102522522',
# )
# tile_plot
# 
# expr_plot <- ExpressionPlot(
#   object = nblast,
#   features = "NFKB1",
#   group.by = 'group',
#   assay = "RNA"
# )
# expr_plot

peak <-'chr4-102521897-102522522'
ranges.show <- StringToGRanges(peak)
ranges.show$color <- "coral"

cov_plot_group <- CoveragePlot(
  object = nblast,
  region.highlight = ranges.show,
  region = peak,
  #features = "NFKB1",
  annotation = TRUE,
  extend.upstream = 1000,
  extend.downstream  = 1000,
  group.by = 'group', 
  peaks = FALSE
) &
  scale_fill_manual( values = group_col ) &
  theme(text = element_text(size=16))

peak1 <-'chr4-102521897-102522522'
ranges.show1 <- StringToGRanges(peak1)
ranges.show1$color <- "coral"

cov_plot1 <- CoveragePlot(
  object = nblast,
  region.highlight = ranges.show1,
  region = peak1,
  #features = "NFKB1",
  annotation = TRUE,
  extend.upstream = 1000,
  extend.downstream  = 1000,
  group.by = 'group', 
  peaks = FALSE
) &
  scale_fill_manual( values = group_col ) &
  theme(text = element_text(size=10)) &
  theme_nb()
cov_plot1

peak2 <-'chr5-132603930-132604445' 
ranges.show2 <- StringToGRanges(peak2)
ranges.show2$color <- "royalblue"
cov_plot2 <- CoveragePlot(
  object = nblast,
  region.highlight = ranges.show2,
  region = peak2,
  #features = "RAD50",
  annotation = TRUE,
  extend.upstream = 1000,
  extend.downstream  = 1000,
  group.by = 'group', 
  peaks = FALSE
) &
  scale_fill_manual( values = group_col ) &
  theme(text = element_text(size=10)) &
  theme_nb()
cov_plot2

chr12-66212986-66213443
chr16-79379258-79379415
il15 <- chr4-141645655-141645941
GACAT3 <- chr2-16148619-16230261
kdm2b <- chr12-121570550-121570768
cflar <- chr2-201135337-201135742
ano6 <- chr12-45232871-45233355

peak3 <-'chr12-121570550-121570768' 
ranges.show3 <- StringToGRanges(peak3)
ranges.show3$color <- "coral"

cov_plot3 <- CoveragePlot(
  object = nblast,
  region.highlight = ranges.show3,
  region = peak3,
  #features = "RAD50",
  annotation = TRUE,
  extend.upstream = 1000,
  extend.downstream  = 1000,
  group.by = 'group', 
  peaks = FALSE
) &
  scale_fill_manual( values = group_col ) &
  theme(text = element_text(size=10)) &
  theme_nb()
cov_plot3

cowplot::plot_grid(cov_plot1 , 
                   cov_plot2, 
                   ncol = 1, 
                   nrow = 2)
ggsave_publication('coverage_plot_nb2', 
                   width = 12, height= 15)
(cov_plot1 & theme_nb() )|
  (cov_plot2 & theme_nb() ) 
ggsave_publication('coverage_plot_nb2', 
                   width = 20, height= 8)

(cov_plot1 & theme(legend.position = "none")) |(cov_plot3  & theme(legend.position = "none"))
ggsave_publication('coverage_plot_nb_kdm2b', 
                   width = 10, height= 7)


cov_plot_sample <- CoveragePlot(
  object = nblast,
  region.highlight = ranges.show,
  region = peak,
  #features = "NFKB1",
  annotation = TRUE,
  extend.upstream = 5000,
  extend.downstream  = 5000,
  group.by = 'sample_id', 
  peaks = FALSE
) &
  scale_fill_manual( values = patients_col ) &
  theme(text = element_text(size=40))


.draw_coverage_plot(nblast, "chr17-60609440-60609846", 'group', group_col)



# Bulk Coverage Plots -----------------------------------------------------

save_coverage_plots <- function(dataframe, path ){
  .draw_coverage_plot <- function(object, peak, group_by, color_scheme) {
    peak <- peak
    ranges.show <- StringToGRanges(peak)
    ranges.show$color <- "coral"
    cov_plt <- 
      CoveragePlot(
        object = object,
        region.highlight = ranges.show,
        region = peak,
        annotation = TRUE,
        extend.upstream = 10000,
        extend.downstream  = 10000,
        group.by = group_by , #'group' 
        peaks = FALSE
      ) &  
      scale_fill_manual( values = color_scheme ) &
      theme(text = element_text(size=28)) 
    return(cov_plt)
  }
  genes_up  <- 
    dataframe %>% 
    filter(!is.na(pathway)) %>% 
    group_by(group, pathway) %>% 
    slice_max(order_by = logFC_ATAC, n = 10) 
  genes_up_list <- genes_up$peak
  genes_down <- 
    dataframe %>% 
    filter(!is.na(pathway)) %>%
    group_by(group, pathway) %>% 
    slice_min(order_by = logFC_ATAC, n = 10) 
  genes_down_list <- genes_down$peak
  
  pat = path
  
  for (i in 1:length(genes_down_list)){
    peak_dn <- genes_down_list[i]
    print(peak_dn)
    filename <- 
      paste0(genes_down[i, ]$pathway, '_' , genes_down[i, ]$gene ,"_dn.pdf")
    peak_dn_aggr <- .draw_coverage_plot(nblast, peak_dn, 'group', group_col)
    peak_dn_sample <- .draw_coverage_plot(nblast, peak_dn, 'sample_id', patients_col)
    
    pdf(file.path(pat, filename), width = 40, height = 20)
    print(peak_dn_aggr)
    print(peak_dn_sample)
    dev.off()
  }
}
save_coverage_plots(monocytes_atac_df,
                    "/media/AGFORTELNY/people/rohit/projects/neuroblastoma/browser_tracks")
# filename_aggr <- 
#   paste0(genes_up[i, ]$pathway, '_' , genes_up[i, ]$gene ,"_aggr.pdf")
# filename_sample <- 
#   paste0(genes_up[i, ]$pathway, '_' ,genes_up[i, ]$gene ,"_sample.pdf")



# PatchWork ---------------------------------------------------------------

tmp_cov_plot_1 <- 
  cov_plot1 & 
  theme_nb() & 
  theme(legend.position = "none")
tmp_cov_plot_2 <- 
  cov_plot2 & 
  theme_nb()
cov_plot1 | cov_plot2  
tmp_umap <- umap + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank())
ggsave_publication('tmp_umap', 
                   width = 9, height= 6)


layout <- "
AA####
AA####
BC####
DE####
"
umap  + scatter + 
  IRF1 +theme(legend.position = "none")  + 
  NFKB1 + theme(legend.position = "none") + 
  fos +theme(legend.position = "none") + 
  bach1 + theme(legend.position = "none") + 
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')
# umap +  promoter_dotplot  + distal_dotplot + 

ggsave_publication('FigureLayout_patchwork_trial', 
                   width = 10, height= 20)

scatter <- 
  p3_e2f / p3_tnfa + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom') 
tmp_umap + p3_e2f + p3_igr + 
  tmp_cov_plot_1 + tmp_cov_plot_2 +  
  IRF1 + IRF2 + 
  fos + bach1 + 
  distal_dotplot + promoter_dotplot +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')

ggsave_publication('FigureLayout_patchwork', 
                   width = 20, height= 29)

(tmp_umap |(p3_e2f/p3_igr))/
  (tmp_cov_plot_1 | tmp_cov_plot_2 | plot_spacer() |plot_spacer() ) / 
  ((IRF1 & theme(legend.position = "none") | IRF2 |plot_spacer() |plot_spacer()) /
     (fos & theme(legend.position = "none")| bach1 | plot_spacer()|plot_spacer())) +
  #plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')

(umap /(p3_e2f / p3_tnfa/(IRF1 |NFKB1)/(fos|bach1))) | ((promoter_dotplot | distal_dotplot))

ggsave_publication('FigureLayout_patchwork_trial', 
                   width = 20, height= 26)

