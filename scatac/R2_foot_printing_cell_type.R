pacman::p_load(Signac, Seurat, 
               scater, tidyverse, ggplot2, 
               patchwork, hrbrthemes, ggbeeswarm, 
               enrichR, UCell, RColorBrewer)
source('/media/AGFORTELNY/people/rohit/projects/neuroblastoma/styling_atac.R')
nb  <- 
  readRDS(
    #file = '/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized.rds')
    file = '/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized_motifs_added.rds')

groups <- list(c('2', '3', '6', '15', '16'),
               c('1', '4', '5', '11'),
               c('7','12'),
               c('8', '9', '10','13', '14'))

nb$group <- nb@meta.data %>% 
  dplyr::mutate(
    group = dplyr::case_when(
      orig.ident %in% groups[[1]] ~"C",
      orig.ident %in% groups[[2]] ~"M", 
      orig.ident %in% groups[[3]] ~"A",
      orig.ident %in% groups[[4]] ~"S")) %>% 
  pull(group)

nb$activ.ident <- nb@meta.data$group

# gather the footprinting information for sets of motifs
nb <- Footprint(
  object = nb,
  motif.name = c("JUNB", "NFE2", "ELF5",
                 "SPIB", "PHOX2B","PHOX2A" ),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(nb, 
                    group.by = 'activ.ident',
                    features = c("GATA2")) &
  theme_nb() &
  scale_color_manual(values = c('navy', 'red', 'royalblue', 'green')) 
  
p2 + patchwork::plot_layout(ncol = 1)

nb@meta.data %>% View()
