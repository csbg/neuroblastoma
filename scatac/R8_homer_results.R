pacman::p_load(tidyverse, data.table, janitor)

homer_df <- function( file, cell_type ) {
  homer <- data.table::fread(file) %>% clean_names() 
  homer %>% 
    data.frame() %>% 
    separate(motif_name,c('c1','source', 'c3'),'/',extra ='merge') %>% 
    separate('c1',  c('motif_name', 'motif_class'), '\\(') %>% 
    mutate(motif_class = str_replace(motif_class,"\\)",""),
           cell_type = cell_type ) %>% 
    mutate(
      percent_of_target_sequences_with_motif = as.numeric(str_replace(
        percent_of_target_sequences_with_motif, '%', '' ))) %>% 
    mutate(
      percent_of_background_sequences_with_motif = as.numeric(str_replace(
        percent_of_background_sequences_with_motif, '%', '' ))) %>% 
    dplyr::select(motif_name, motif_class, source, 
           consensus, p_value, log_p_value,
           percent_of_target_sequences_with_motif,
           percent_of_background_sequences_with_motif,
           q_value_benjamini, cell_type)
}


# Monocytes fold ----------------------------------------------------------


path = '/Volumes/GFS_BIOMEDBIO_AGFORTELNY/people/rohit/projects/clusterObj'

homer_list <- list(
  'nb' = homer_df(file.path(path, 'nb/knownResults.txt'), "NB"),
  'tcell' = homer_df(file.path(path, 'tcell/knownResults.txt'), "T-cell"),
  'bcell' = homer_df(file.path(path, 'bcell/knownResults.txt'), "B-cell"),
  'monocytes' = homer_df(file.path(path, 'monocytes/knownResults.txt'), "Monocytes"),
  'basophils' = homer_df(file.path(path, 'basophils/knownResults.txt'), "Basophils"),
  'cluster11' = homer_df(file.path(path, 'cluster11/knownResults.txt'), "Cluster11"),
  'cluster22' = homer_df(file.path(path, 'cluster22/knownResults.txt'), "Cluster22"),
  'cluster32' = homer_df(file.path(path, 'cluster32/knownResults.txt'), "Cluster32"),
  'cluster33' = homer_df(file.path(path, 'cluster33/knownResults.txt'), "Cluster33"),
  'cluster37' = homer_df(file.path(path, 'cluster37/knownResults.txt'), "Cluster37")
 )


my_path = '/Volumes/ROHIT32/group_motif'
my_path = '/Volumes/GFS_BIOMEDBIO_AGFORTELNY/people/rohit/projects/clusterObj/group_microenv'
homer_fold <- list(
  'mono_atrx_dn' = homer_df(file.path(my_path, 
                  'monocytes_down_atrx_30_12/knownResults.txt'), "Atrx Close"),
  'mono_atrx_up' = homer_df(file.path(my_path, 
                 'monocytes_up_atrx_30_12/knownResults.txt'), "Atrx Open"),
  'mono_mycn_dn' = homer_df(file.path(my_path, 
                  'monocytes_down_mycn_30_12/knownResults.txt'), "MYCN Close"),
  'mono_mycn_up' = homer_df(file.path(my_path, 
                  'monocytes_up_mycn_30_12/knownResults.txt'), "MYCN Open"),
  'mono_sporadic_dn' = homer_df(file.path(my_path, 
                      'monocytes_down_sporadic_30_12/knownResults.txt'), "Sporadic Close"),
  'mono_sporadic_up' = homer_df(file.path(my_path, 
                      'monocytes_up_sporadic_30_12/knownResults.txt'), "Sporadic Open")
)

my_path = '/Volumes/GFS_BIOMEDBIO_AGFORTELNY/people/rohit/projects/clusterObj/monocytes/distal'
homer_fold <- list(
  'mono_atrx_dn' = homer_df(file.path(my_path, 
                                      'atrx_promoter300_down/knownResults.txt'), "Atrx Close"),
  'mono_atrx_up' = homer_df(file.path(my_path, 
                                      'atrx_promoter300_up/knownResults.txt'), "Atrx Open"),
  'mono_mycn_dn' = homer_df(file.path(my_path, 
                                      'mycn_promoter300_down/knownResults.txt'), "MYCN Close"),
  'mono_mycn_up' = homer_df(file.path(my_path, 
                                      'mycn_promoter300_up/knownResults.txt'), "MYCN Open"),
  'mono_sporadic_dn' = homer_df(file.path(my_path, 
                                          'sporadic_promoter300_down/knownResults.txt'), "Sporadic Close"),
  'mono_sporadic_up' = homer_df(file.path(my_path, 
                                          'sporadic_promoter300_up/knownResults.txt'), "Sporadic Open")
)
my_path <- '/media/AGFORTELNY/people/rohit/projects/clusterObj/monocytes/distal'
homer_fold_distal <- list(
  'mono_atrx_dn' = homer_df(file.path(my_path, 
                                      'atrx_distal_down/knownResults.txt'), "Atrx Close"),
  'mono_atrx_up' = homer_df(file.path(my_path, 
                                      'atrx_distal_up/knownResults.txt'), "Atrx Open"),
  'mono_mycn_dn' = homer_df(file.path(my_path, 
                                      'mycn_distal_down/knownResults.txt'), "MYCN Close"),
  'mono_mycn_up' = homer_df(file.path(my_path, 
                                      'mycn_distal_up/knownResults.txt'), "MYCN Open"),
  'mono_sporadic_dn' = homer_df(file.path(my_path, 
                                          'sporadic_distal_down/knownResults.txt'), "Sporadic Close"),
  'mono_sporadic_up' = homer_df(file.path(my_path, 
                                          'sporadic_distal_up/knownResults.txt'), "Sporadic Open")
)



test_homer_dn <- do.call( rbind, list(homer_fold$mono_atrx_dn, 
                              homer_fold$mono_mycn_dn, 
                              homer_fold$mono_sporadic_dn))

test_homer_dn %>% 
  dplyr::filter( q_value_benjamini < 0.05 ) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log2(odds_ratio)) %>% 
  dplyr::select(motif_name, p_value, log_p_value, cell_type, log_odds_ratio) %>% 
  dplyr::filter(abs(log_odds_ratio) >1) %>%
  #dplyr::mutate(log_p_value = -log_p_value) %>% 
  pivot_wider(!c( p_value ),
              names_from = cell_type, values_from =log_odds_ratio,  
              values_fn = max, values_fill = 0) %>%
#dplyr::filter(motif_name %in% detect_moi) %>% 
tibble::column_to_rownames(., var = 'motif_name')  %>% 
  as.matrix() %>% 
  ComplexHeatmap::Heatmap()
pheatmap::pheatmap()

homer_fold_db <- do.call(rbind, homer_fold)
homer_fold_distal_db <- do.call(rbind, homer_fold_distal)

hierarch.ordering <- function(dt, toOrder, orderBy, value.var, aggregate = FALSE){
  
  if(!aggregate){orMT <- t(as.matrix(dcast.data.table(dt, get(orderBy) ~ get(toOrder), value.var=value.var)[,-"orderBy",with=F]))}
  else{orMT <- t(as.matrix(dcast.data.table(dt, get(orderBy) ~ get(toOrder), value.var=value.var, fun.aggregate=mean)[,-"orderBy",with=F]))}
  orMT[is.na(orMT)] <- 1
  orMT[orMT == Inf] <- max(orMT[orMT != Inf])
  hclustObj <- hclust(dist(orMT))
  dt[[toOrder]] <- factor(dt[[toOrder]], levels=hclustObj$labels[hclustObj$order])
  return(dt)
  
}


htmap <- 
  homer_fold_db %>% 
  dplyr::filter( q_value_benjamini < 0.05 ) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log(odds_ratio)) %>% 
  dplyr::select(motif_name, p_value, log_p_value, cell_type, log_odds_ratio) %>% 
  dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >=0.75) %>% 
  dplyr::mutate(log_p_value = -log_p_value)  %>%  #test_homer_db
  pivot_wider(!c( p_value, log_p_value ),
              names_from = cell_type, values_from =log_odds_ratio,  
              values_fn = max, values_fill = 0) %>% 
  #dplyr::filter(motif_name %in% detect_moi) %>% 
  tibble::column_to_rownames(., var = 'motif_name')  %>% 
  as.matrix() %>%
  ComplexHeatmap::Heatmap(show_column_dend = FALSE, 
                          show_row_dend = FALSE,
                          col = c('black', 'orange'),
                          name = 'log(odds.ratio)',
                          heatmap_legend_param = list(
                            legend_direction = "horizontal", 
                            legend_width = unit(5, "cm")
                          ))
ComplexHeatmap::draw(htmap, heatmap_legend_side = "top")

pheatmap::pheatmap()
ord <- c("MYCN Close", "Atrx Close", "Sporadic Close", 
         "MYCN Open" , "Atrx Open" , "Sporadic Open")
test_homer_db  %>% 
  ggplot(.) +
  aes(x= factor(cell_type, levels = ord), y = motif_name) +
  geom_point(aes(size = log_odds_ratio,  color= as.numeric(p_value))) +
  hrbrthemes::theme_ipsum() + ylab('') + xlab('') 
  
  
  lemon::facet_rep_grid(Name ~ ., scales = "free", space="free", 
                        repeat.tick.labels = 'bottom')
homer_fold_distal_db %>% 
  dplyr::filter( q_value_benjamini < 0.05 ) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log(odds_ratio)) %>% 
  dplyr::select(motif_name, p_value, log_p_value, cell_type, log_odds_ratio) %>% 
  dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >=0.75) %>% 
  dplyr::mutate(log_p_value = -log_p_value)  %>%  #test_homer_db
  pivot_wider(!c( p_value, log_p_value ),
              names_from = cell_type, values_from =log_odds_ratio,  
              values_fn = max, values_fill = 0) %>% 
  #dplyr::filter(motif_name %in% detect_moi) %>% 
  tibble::column_to_rownames(., var = 'motif_name')  %>% 
  as.matrix() %>%
  ComplexHeatmap::Heatmap(show_column_dend = FALSE, 
                          show_row_dend = FALSE,
                          col = c('black', 'orange'),
                          name = 'log(odds.ratio)',
                          heatmap_legend_param = list(
                            legend_direction = "horizontal", 
                            legend_width = unit(5, "cm")
                          ))
# T-cell fold -------------------------------------------------------------

my_path = '/Volumes/ROHIT32/group_motif'
t_homer_fold <- list(
  't_atrx_dn' = homer_df(file.path(my_path, 
                                      't-cells_atrx_down/knownResults.txt'), "T Atrx Down"),
  't_atrx_up' = homer_df(file.path(my_path, 
                                      't-cells_atrx_up/knownResults.txt'), "T Atrx Up"),
  't_mycn_dn' = homer_df(file.path(my_path, 
                                      't-cells_mycn_down/knownResults.txt'), "T MYCN Down"),
  't_mycn_up' = homer_df(file.path(my_path, 
                                      't-cells_mycn_up/knownResults.txt'), "T MYCN Up"),
  't_sporadic_dn' = homer_df(file.path(my_path, 
                                          't-cells_sporadic_down/knownResults.txt'), "T Sporadic Down"),
  't_sporadic_up' = homer_df(file.path(my_path, 
                                          't-cells_sporadic_up/knownResults.txt'), "T Sporadic Up")
)


t_homer_fold_db <- do.call(rbind, t_homer_fold)

t_homer_fold_db %>% 
  dplyr::filter( q_value_benjamini < 0.05 ) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log(odds_ratio)) %>% 
  dplyr::select(motif_name, p_value, log_p_value, cell_type, log_odds_ratio) %>% 
  dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >=0.9) %>% 
  dplyr::mutate(log_p_value = -log_p_value) %>% 
  pivot_wider(!c( p_value, log_p_value ),
              names_from = cell_type, values_from =log_odds_ratio,  
              values_fn = max, values_fill = 0) %>%
  #dplyr::filter(motif_name %in% detect_moi) %>% 
  tibble::column_to_rownames(., var = 'motif_name')  %>% 
  as.matrix() %>% 
  ComplexHeatmap::Heatmap()


# B-cell fold -------------------------------------------------------------

b_homer_fold <- list(
  'b_atrx_dn' = homer_df(file.path(my_path, 
                                   'b-cells_atrx_down/knownResults.txt'), "B Atrx Down"),
  'b_atrx_up' = homer_df(file.path(my_path, 
                                   'b-cells_atrx_up/knownResults.txt'), "B Atrx Up"),
  'b_mycn_dn' = homer_df(file.path(my_path, 
                                   'b-cells_mycn_down/knownResults.txt'), "B MYCN Down"),
  'b_mycn_up' = homer_df(file.path(my_path, 
                                   'b-cells_mycn_up/knownResults.txt'), "B MYCN Up"),
  'b_sporadic_dn' = homer_df(file.path(my_path, 
                                       'b-cells_sporadic_down/knownResults.txt'), "B Sporadic Down"),
  'b_sporadic_up' = homer_df(file.path(my_path, 
                                       'b-cells_sporadic_up/knownResults.txt'), "B Sporadic Up")
)


b_homer_fold_db <- do.call(rbind, b_homer_fold)

b_homer_fold_db %>% 
  dplyr::filter( q_value_benjamini < 0.05 ) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log(odds_ratio)) %>% 
  dplyr::select(motif_name, p_value, log_p_value, cell_type, log_odds_ratio) %>% 
  dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >=0.9) %>% 
  dplyr::mutate(log_p_value = -log_p_value) %>% 
  pivot_wider(!c( p_value, log_p_value ),
              names_from = cell_type, values_from =log_odds_ratio,  
              values_fn = max, values_fill = 0) %>%
  #dplyr::filter(motif_name %in% detect_moi) %>% 
  tibble::column_to_rownames(., var = 'motif_name')  %>% 
  as.matrix() %>% 
  ComplexHeatmap::Heatmap()

# Monocytes combined ------------------------------------------------------

homerdb %>% 
  dplyr::filter(q_value_benjamini < 0.05 & percent_of_target_sequences_with_motif >2 ) %>% 
  dplyr::select(motif_name, p_value, log_p_value, cell_type) %>% 
  dplyr::mutate(log_p_value = -log_p_value) %>% 
  pivot_wider(!c( p_value ),
              names_from = cell_type, values_from =log_p_value,  
              values_fn = max, values_fill = 0) %>% 
  #dplyr::filter(motif_name %in% detect_moi) %>% 
  tibble::column_to_rownames(., var = 'motif_name')  %>% 
  as.matrix() %>% 
  ComplexHeatmap::Heatmap()
  pheatmap::pheatmap()

  motifdb[motifdb > 0] = 1

  motifdb %>% 
    as.matrix() %>% 
    #pheatmap::pheatmap()
    ComplexHeatmap::Heatmap(row_km = 5)

  moi <- c('AP', 'NFAT', 'AP', 'YY-1', 'STAT', 'GATA', 
           'NFKB', 'IL2Ra', 'ATF-2', 'CREB', 'CREB', 'Hlf', 
           'Ik', 'POU2F1', 'POU2F1', 'POU2F1b', 'POU2F1c', 
           'EGR1',  'BRG', 'LITAF', 'ETS', 'SP1', 
           'c-fos', 'c-jun', 'c-Fos', 'c-Jun', 'p53', 'TBP', 
           'CTCF', 'MAZ1', 'hnRNPK', 'CNBP',  'hnRNPA1', 'FBP', 
           'FIR', 'TFIIH', 'GCNF', 'GCNF', 'GCNF', 'AML1a', 
           'Myb', 'Myc', 'Max1', 'Sox9', 'Tal1beta', 'USF-1', 
           'USF1', 'E2F',  'REST','NRSF',
           'IRF', 'MyoD',  'HOXA9', 'HOXA9B', 'Meis', 
           'Pax-4a', 'POU3F2', 'Arnt','HNF', 'Max', 
           'MZF-1', 'SRY', 'TFIID', 'E47', 'Evi', 'FOXD3', 'Ik2', 'Pax')
#  delta, CBP/p300, C/EBP, nm23/puf60, NRSF form 1, NRSF form 2,
  
  motif_names <- rownames(motif_mtx)
  detect_moi <- motif_names[str_detect(motif_names, paste(moi, collapse = "|"))]



pacman::p_load(motifStack, universalmotif)

pcm <- read.table(file.path(find.package("motifStack"), 
                              "extdata", "bin_SOLEXA.pcm"))


homer_known1_monocytes <- universalmotif::read_homer(
  '/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/monocytes/knownResults/known1.motif')
homer_known1_monocytes_Atf3 <- universalmotif::convert_motifs(homer_known1_monocytes, 'motifStack-pcm')


motif <- new("pcm", 
             mat=homer_known1_monocytes_Atf3@mat, 
             name=homer_known1_monocytes_Atf3@name)
##pfm object
#motif <- pcm2pfm(pcm)
#motif <- new("pfm", mat=motif, name="bin_SOLEXA")
plot(motif, ic.scale=TRUE, ylab="probability")

motifStack(motif)

mpath = '/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/monocytes/knownResults'
mnames = list.files(mpath, pattern=".motif")



for (i in file.path(mpath, mnames)) {
  outfile <- file.path(mpath, 
                      paste0(strsplit(basename(i), '\\.')[[1]][1], '.pcm'))
  homer_known <- universalmotif::read_homer(i)
  #universalmotif::write_jaspar(homer_known, file= outfile, overwrite = TRUE)
  universalmotif::write_matrix(homer_known, file=outfile, type = 'PCM', overwrite = TRUE)
}


pcms<-importMatrix(list.files(mpath, pattern = "*.pwm"),
                   format="jaspar", to="pfm")
motifs<-AAmotifAlignment(pcms)

pwms <- importMatrix(dir('/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/monocytes/knownResults', 
                ".pwm", full.names = TRUE), format = 'jaspar')

motifs <- DNAmotifAlignment(pwms)

motifStack(motifs[1:10], layout = "radialPhylog")

str(motifs[1:10])

len <- length(motifs[1:10])
df <- data.frame(x=.5, y=(seq.int(len)-.5)/len, 
                 width=.75, height=1/(len+1))
df$motif <- motifs[1:10]
library(ggplot2)
ggplot(df, aes(x=x, y=y, width=width, height=height, motif=motif)) +
  geom_motif(use.xy = TRUE) + theme_bw() + xlim(0, 1) + ylim(0, 1) + 
  coord_flip()


monofile <- '/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/monocytes/knownResults.txt'
tmp_mono <- data.table::fread(monofile) %>% 
  clean_names() 

tmp_mono %>% 
  dplyr::filter(q_value_benjamini < 0.05) %>% 
  mutate(percent_of_target_sequences_with_motif = as.numeric(str_replace(percent_of_target_sequences_with_motif, 
                                                              '%', '' ))) %>% 
  mutate(percent_of_background_sequences_with_motif = as.numeric(str_replace(percent_of_background_sequences_with_motif, 
                                                                         '%', '' ))) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log2(odds_ratio)) %>% 
  dplyr::select(motif_name, consensus, p_value, log_p_value, q_value_benjamini,odds_ratio, log_odds_ratio) %>% 
  View()


dplyr::filter(q_value_benjamini < 0.05 & percent_of_target_sequences_with_motif >5) %>% View()
 
col_ord <- c('NB', 'Monocytes', 'Basophils', 'T-cell', 'B-cell', 
             'Cluster11', 'Cluster22', 'Cluster37') 
# PLot odds ratio
homerdb %>% 
  dplyr::filter(q_value_benjamini < 0.05) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log2(odds_ratio)) %>% 
  dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >1) %>% 
  dplyr::select(motif_name, p_value, log_odds_ratio, cell_type) %>% 
  pivot_wider(!c( p_value ),
                names_from = cell_type, values_from =log_odds_ratio,  
                values_fn = max, values_fill = 0) %>% 
  tibble::column_to_rownames(., var = 'motif_name')  %>% 
  dplyr::select(col_ord) %>% 
  as.matrix() %>% 
  #pheatmap::pheatmap()
  ComplexHeatmap::Heatmap(name = 'log(odds.ratio)',
                            cluster_columns = FALSE)
 
# plot log(p-value) 
homerdb %>% 
  dplyr::filter(q_value_benjamini < 0.05) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log2(odds_ratio)) %>% 
  dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >1) %>%
  dplyr::select(motif_name, p_value, log_p_value, cell_type) %>% 
  dplyr::mutate(p_value  = ifelse(p_value < 1e-10, 1e-10, p_value )) %>%
  dplyr::mutate(log_p_value = -log2(p_value)) %>% 
  pivot_wider(!c( p_value ),
              names_from = cell_type, values_from =log_p_value,  
              values_fn = max, values_fill = 0) %>% 
  tibble::column_to_rownames(., var = 'motif_name')  %>% 
  as.matrix() %>%
  #pheatmap::pheatmap() %>% 
  ComplexHeatmap::Heatmap(name = 'log(p-value)',
                          cluster_columns = FALSE,
                          )
# Patient group analysis --------------------------------------------------

library(circlize)
col_fun = colorRamp2(c(0, 2, 4), c("royalblue", "white", "red"))
col_fun(seq(-3, 3))
  
path_group = '/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/group_microenv'

homer_group_list <- list(
  'atrx' = homer_df(file.path(path_group, 'monocytes_atrx/knownResults.txt'), "ATRX"),
  'mycn' = homer_df(file.path(path_group, 'monocytes_mycn/knownResults.txt'), "MYCN"),
  'sporadic' = homer_df(file.path(path_group, 'monocytes_sporadic/knownResults.txt'), "Sporadic")
)

homer_group_db <- do.call(rbind, homer_group_list)

homer_group_db %>% 
  dplyr::filter(q_value_benjamini < 0.05) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log2(odds_ratio)) %>% 
  dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >1) %>% 
  dplyr::select(motif_name, p_value, log_odds_ratio, cell_type) %>% 
  pivot_wider(!c( p_value ),
              names_from = cell_type, values_from =log_odds_ratio,  
              values_fn = max, values_fill = 0) %>% 
  tibble::column_to_rownames(., var = 'motif_name')  %>% 
  as.matrix() %>%
  ComplexHeatmap::Heatmap(name = 'log(odds.ratio)', 
                          cluster_columns = FALSE)
  



homerdb %>% 
  dplyr::filter(q_value_benjamini < 0.05 & cell_type =="Monocytes") %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log2(odds_ratio)) %>% 
  dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >1) %>% 
  arrange(p_value) %>% 
  head(20) %>% 
  pull(motif_name) ->moi_monocytes

filtered_list <- list()

for (i in motifs) {
  new_name <- strsplit(i$name,"\\(")[[1]][1]
  if (new_name %in% moi_monocytes){
    print(new_name)
    filtered_list[[new_name]] <- i
  }
}

motifStack(filtered_list, layout = 'radialPhylog')


monocytes_enr <- read.table('/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/Monocytes.MASTenriched_peaks.txt',
           sep = ',', header = TRUE)
monocytes_enr2 <- read.table('/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/Monocytes.enriched_peaks.txt',
                            sep = ',', header = TRUE)
tcell_enr <- read.table('/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/T-cells.MASTenriched_peaks.txt',
                            sep = ',', header = TRUE)
tcell_enr2 <- read.table('/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/T-cells.enriched_peaks.txt',
                             sep = ',', header = TRUE)




# homer and homer goi -----------------------------------------------------

hierarch.ordering <- 
  function(dt, toOrder, orderBy, value.var, aggregate = FALSE){
  if(!aggregate){orMT <- t(as.matrix(dcast.data.table(dt, get(orderBy) ~ get(toOrder), value.var=value.var)[,-"orderBy",with=F]))}
  else{orMT <- t(as.matrix(dcast.data.table(dt, get(orderBy) ~ get(toOrder), value.var=value.var, fun.aggregate=mean)[,-"orderBy",with=F]))}
  orMT[is.na(orMT)] <- 1
  orMT[orMT == Inf] <- max(orMT[orMT != Inf])
  hclustObj <- hclust(dist(orMT))
  dt[[toOrder]] <- factor(dt[[toOrder]], levels=hclustObj$labels[hclustObj$order])
  return(dt)
  
}

hierarch.ordering_rev <- 
  function(dt, toOrder, orderBy, value.var, aggregate = FALSE){
    if(!aggregate){orMT <- t(as.matrix(dcast.data.table(dt, get(orderBy) ~ get(toOrder), value.var=value.var)[,-"orderBy",with=F]))}
    else{orMT <- t(as.matrix(dcast.data.table(dt, get(orderBy) ~ get(toOrder), value.var=value.var, fun.aggregate=mean)[,-"orderBy",with=F]))}
    orMT[is.na(orMT)] <- 1
    orMT[orMT == Inf] <- max(orMT[orMT != Inf])
    hclustObj <- hclust(dist(orMT))
    dt[[toOrder]] <- factor(dt[[toOrder]], levels=rev(hclustObj$labels[hclustObj$order]))
    return(dt)
  }

pacman::p_load(tidyverse, data.table, janitor)

homer_df <- function( file, cell_type ) {
  homer <- data.table::fread(file) %>% clean_names() 
  homer %>% 
    data.frame() %>% 
    separate(motif_name,c('c1','source', 'c3'),'/',extra ='merge') %>% 
    separate('c1',  c('motif_name', 'motif_class'), '\\(') %>% 
    dplyr::mutate(motif_class = str_replace(motif_class,"\\)",""),
           cell_type = cell_type ) %>% 
    dplyr::mutate(
      percent_of_target_sequences_with_motif = as.numeric(str_replace(
        percent_of_target_sequences_with_motif, '%', '' ))) %>% 
    dplyr::mutate(
      percent_of_background_sequences_with_motif = as.numeric(str_replace(
        percent_of_background_sequences_with_motif, '%', '' ))) %>% 
    dplyr::select(motif_name, motif_class, source, 
                  consensus, p_value, log_p_value,
                  percent_of_target_sequences_with_motif,
                  percent_of_background_sequences_with_motif,
                  q_value_benjamini, cell_type)
}

homer_df_id <- function( file, cell_type, location ) {
  homer <- data.table::fread(file) %>% clean_names() 
  homer <- homer %>% 
    data.frame() %>% 
    dplyr::mutate(location = location) %>% 
    separate(motif_name,c('c1','source', 'c3'),'/',extra ='merge') %>% 
    separate('c1',  c('motif_name', 'motif_class'), '\\(') %>% 
    dplyr::mutate(motif_class = str_replace(motif_class,"\\)",""),
                  cell_type = cell_type) %>% 
    dplyr::mutate(
      percent_of_target_sequences_with_motif = as.numeric(str_replace(
        percent_of_target_sequences_with_motif, '%', '' ))) %>% 
    dplyr::mutate(
      percent_of_background_sequences_with_motif = as.numeric(str_replace(
        percent_of_background_sequences_with_motif, '%', '' ))) %>% 
    dplyr::select(
      location, motif_name, motif_class, source,
                  consensus, p_value, log_p_value,
                  percent_of_target_sequences_with_motif,
                  percent_of_background_sequences_with_motif,
                  q_value_benjamini, cell_type)
}

draw_heatmap <- function(data_list, title) {
  require(ComplexHeatmap)
  require(tidyverse)
  
  dataframe <- do.call(rbind, data_list)
  heatmap <- 
    dataframe %>% 
    dplyr::filter( q_value_benjamini < 0.05 ) %>% 
    mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
    mutate(log_odds_ratio = log(odds_ratio)) %>% 
    dplyr::select(motif_name, p_value, log_p_value, cell_type, log_odds_ratio) %>% 
    dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >=1) %>% 
    dplyr::mutate(log_p_value = -log_p_value)  %>%
    pivot_wider(!c( p_value, log_p_value ),
                names_from = cell_type, values_from =log_odds_ratio,  
                values_fn = max, values_fill = 0) %>% 
    tibble::column_to_rownames(., var = 'motif_name')  %>% 
    as.matrix() %>%
    ComplexHeatmap::Heatmap(
      show_column_dend = FALSE, 
      show_row_dend = FALSE, 
      col = c('black', 'orange'),
      name = 'log(odds.ratio)', 
      column_title = title,
      heatmap_legend_param = list( legend_direction = "horizontal", legend_width = unit(5, "cm")))
  heatmap
}

my_path = '/media/AGFORTELNY/people/rohit/projects/clusterObj/monocytes'

#my_path = '/Volumes/GFS_BIOMEDBIO_AGFORTELNY/people/rohit/projects/clusterObj/monocytes'

all1 <- draw_heatmap(homer_fold_promoter, 'promoter')
all2 <- draw_heatmap(homer_fold_intermediate,'intermediate')
all3 <- draw_heatmap(homer_fold_distal, 'distal')
all1 + all2 + all3

all_promoters <- draw_heatmap(homer_fold_promoter, 'promoter')
all_others <- draw_heatmap(homer_fold_others, "others")

# DotPlot order -----------------------------------------------------------

penguins %>% 
  group_by(sex) %>%
  summarize(ave_bill_length_mm=mean(bill_length_mm))
#function(dt, toOrder, orderBy, value.var, aggregate = FALSE)
do.call(rbind, homer_fold_promoter)  %>% 
  dplyr::filter( q_value_benjamini < 0.05 ) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log(odds_ratio)) %>% 
  setDT(.)  ->all_promoters_dt 

do.call(rbind, homer_fold_others)  %>% 
  dplyr::filter( q_value_benjamini < 0.05 ) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log(odds_ratio)) %>% 
  filter(log_odds_ratio >1) %>% 
  setDT(.)  ->all_others_dt 

all_others_df <- 
  do.call(rbind, homer_fold_others)  %>% 
  dplyr::filter( q_value_benjamini < 0.05 ) %>% 
  mutate(odds_ratio = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) %>% 
  mutate(log_odds_ratio = log(odds_ratio)) 
  
#a-zA-Z]{4} |\d{4}
str_match(s, )
stringr::str_match(all_others_df$motif_name, "[a-zA-Z]{3}[\\d]")
all_others_df$group_mot <- 
  str_extract(all_others_df$motif_name, 
              substr(all_others_df$motif_name, 1, 3))

all_others_dt <-  
  all_others_df %>% setDT(.)
ordered_distal <-  hierarch.ordering(all_others_dt, 
                                     "motif_name", 
                                     "cell_type",
                                     "log_odds_ratio",
                                     aggregate = FALSE) 
ordered_distal_rev <-  hierarch.ordering_rev(all_others_dt, 
                                     "motif_name", 
                                     "cell_type",
                                     "log_odds_ratio",
                                     aggregate = FALSE) 
hierarch.ordering_rev

ordered_distal %>% 
  group_by(group_mot, cell_type, q_value_benjamini) %>% 
  summarize(log_odds_ratio_grp=sum(log_odds_ratio)) %>% 
  ggplot(., aes(factor(cell_type, levels = ord_l), y = group_mot)) +
  geom_point(aes(size = log_odds_ratio_grp, color = q_value_benjamini)) +
  hrbrthemes::theme_ipsum()
  
  

# Motifs for footprinting -------------------------------------------------

file_dir <- '/media/AGFORTELNY/people/rohit/projects/bigwig/monocytes/hint_patients/hint_results/Lineplots'
foi <- list()
for (i in ordered_promoters$motif_name){
  search <- paste0('*',toupper(i),'*.txt') #substring(i, 1, 4)
  found <- Sys.glob(file.path(file_dir, search))
  if (!identical(found, character(0))) {
    foi[[length(foi) + 1]] <- found
    print(found)
  } 
}

for (i in unlist(foi)) {
  print(i)
}

Sys.glob(file.path(file_dir, "*ATF*.txt")) ## file_dir = file containing directory
files <- 
  list.files(file_dir)
print(grep(i, files, fixed=T))


# Motifs DotPlot ----------------------------------------------------------
temp_other_dt <- 
  all_others_dt %>% 
  data.frame() %>% View()
  filter(log_odds_ratio >= 0.5) %>% 
  data.table::setDT()

ordered_temp_dt <- 
  hierarch.ordering(temp_other_dt, 
                  "motif_name", 
                  "cell_type",
                  "log_odds_ratio",
                  aggregate = FALSE) 
ordered_temp_dt %>% 
ggplot(., aes(x= factor(cell_type, levels = ord_l), y = motif_name)) +
  geom_point(aes(color = log_odds_ratio, size = -log(q_value_benjamini))) +
  labs(title = "distal",x = 'Patient group wrt to Ctrl', y = 'Motif', 
       size = "Adj_pvalue", color = "Log(Odds_ratio)") +
  scale_colour_gradient(high = "red", low = "blue", na.value = NA) +
  theme_nb() +
  theme(axis.text.x = element_text(angle = 90)) 


ordered_distal <-  hierarch.ordering(all_others_dt, 
                    "motif_name", 
                    "cell_type",
                    "log_odds_ratio",
                    aggregate = FALSE) 
  dplyr::select(motif_name, p_value, log_p_value, cell_type, log_odds_ratio) %>% 
  dplyr::filter(is.finite(log_odds_ratio) & abs(log_odds_ratio) >=0.75) %>% 
  dplyr::mutate(log_p_value = -log_p_value)  
  ord_l <- c("MYCN Close", "Atrx Close", "Sporadic Close", 
             "MYCN Open", "Atrx Open" ,"Sporadic Open" )

ordered_promoters <- 
  hierarch.ordering(all_promoters_dt, 
                    "motif_name", 
                    "cell_type",
                    "log_odds_ratio",
                    aggregate = FALSE) 
ordered_promoters_rev <- 
  hierarch.ordering_rev (all_promoters_dt, 
                    "motif_name", 
                    "cell_type",
                    "log_odds_ratio",
                    aggregate = FALSE) 

(promoter_dotplot_rev <- 
    ordered_promoters_rev  %>% 
    mutate(q_value_benjamini = case_when(q_value_benjamini == 0 ~ 0.00001,
                                         TRUE ~ q_value_benjamini )) %>% 
    ggplot(., aes(x= factor(cell_type, levels = rev(ord_l)), y = motif_name)) +
    geom_point(aes(color = log_odds_ratio, size = -log(q_value_benjamini))) +
    labs(title = "promoters",x = 'Patient group wrt to Ctrl', y = 'Motif', 
         size = "-log10(Adj_pvalue)", color = "Log(Odds_ratio)") +
    scale_colour_gradient(high = "red", low = "grey80", na.value = NA) +
    theme_nb() +
    #scale_y_discrete(position = "right") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=6), #change legend title font size
          legend.text = element_text(size=6)) + #change legend text font size
    # theme(strip.text.x = element_blank(),
    #       #legend.background = element_rect(fill = "white", color = "grey80"),
    #       strip.background = element_rect(colour="white", fill="white"),
    #       legend.position=c(.1,.85))
    theme(legend.position="left")
)
(distal_dotplot_rev <- 
    ordered_distal_rev  %>% 
    mutate(q_value_benjamini = case_when(q_value_benjamini == 0 ~ 0.00001,
                                         TRUE ~ q_value_benjamini )) %>% 
    ggplot(., aes(x= factor(cell_type, levels = rev(ord_l)), y = motif_name)) +
    geom_point(aes(color = log_odds_ratio, size = -log10(q_value_benjamini))) +
    labs(title = "distal",x = 'Patient group wrt to Ctrl', y = 'Motif', 
         size = "-log10(Adj_pvalue)", color = "Log(Odds_ratio)") +
    scale_colour_gradient(high = "red", low = "grey80", na.value = NA) +
    theme_nb() + #coord_flip() +
    scale_y_discrete(position = "right") +
    # theme(axis.text.x = element_text(angle = 90))
    #
    theme(axis.text.x = element_text(angle = 90))  +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=6), #change legend title font size
          legend.text = element_text(size=6)) + #change legend text font size
    # theme(strip.text.x = element_blank(),
    #       #legend.background = element_rect(fill = "white", color = "grey80"),
    #       strip.background = element_rect(colour="white", fill="white"),
    #       legend.position=c(0.9,0.85))
    theme(legend.position="right") 
)

promoter_dotplot <- 
promoter_dotplot_rev +
  theme(legend.position="none") +
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=10)) 
distal_dotplot <- 
  distal_dotplot_rev +
  theme(legend.position="none") +
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=10)) 

ggsave_publication('motif_dotplot_rev', 
                   width = 20, height= 25)

promoter_dotplot_rev + distal_dotplot_rev 
ggsave_publication('motif_dotplot_rev', 
                   width = 15, height= 18)



ggsave_publication('motif_dotplot2', 
                   width = 16, height= 28)

ggsave_publication('promoter_dotplot', 
                   width = 10, height= 20)
ggsave_publication('promoter_dotplot2', 
                   width = 10, height= 20)
ggsave_publication('promoter_dotplot3', 
                   width = 10, height= 20)

((example_plot/plot_spacer()/plot_spacer()) | example_plot2) + 
  plot_layout(guides = 'collect')
ggsave_publication('promoter_dotplot_combined', 
                   width = 20, height= 60)
# Pathway specific motif analysis -----------------------------------------

goi1 <- draw_heatmap(homer_fold_promoter_goi, 'promoter')
goi2 <- draw_heatmap(homer_fold_intermediate_goi, 'intermediate')
goi3 <- draw_heatmap(homer_fold_distal_goi, 'distal')

goi1 + goi2 + goi3

goi1 +
  labs(title = "Promoter IGR")

homer_fold_promoter <- list(
  'mono_atrx_dn' = homer_df_id(
    file.path(my_path, '/promoter/atrx_promoter_down/knownResults.txt'), "Atrx Close", 'promoter'),
  'mono_atrx_up' = homer_df_id(
    file.path(my_path, '/promoter/atrx_promoter_up/knownResults.txt'), "Atrx Open",'promoter'),
  'mono_mycn_dn' = homer_df_id(
    file.path(my_path, '/promoter/mycn_promoter_down/knownResults.txt'), "MYCN Close",'promoter'),
  'mono_mycn_up' = homer_df_id(
    file.path(my_path, '/promoter/mycn_promoter_up/knownResults.txt'), "MYCN Open",'promoter'),
  'mono_sporadic_dn' = homer_df_id(
    file.path(my_path,  '/promoter/sporadic_promoter_down/knownResults.txt'), "Sporadic Close",'promoter'),
  'mono_sporadic_up' = homer_df_id(
    file.path(my_path, '/promoter/sporadic_promoter_up/knownResults.txt'), "Sporadic Open",'promoter')
)

homer_fold_intermediate <- list(
  'mono_atrx_dn' = homer_df_id(
    file.path(my_path, '/intermediate/atrx_intermediate_down/knownResults.txt'), "Atrx Close", 'intermediate'),
  'mono_atrx_up' = homer_df_id(
    file.path(my_path, '/intermediate/atrx_intermediate_up/knownResults.txt'), "Atrx Open", 'intermediate'),
  'mono_mycn_dn' = homer_df_id(
    file.path(my_path, '/intermediate/mycn_intermediate_down/knownResults.txt'), "MYCN Close", 'intermediate'),
  'mono_mycn_up' = homer_df_id(
    file.path(my_path, '/intermediate/mycn_intermediate_up/knownResults.txt'), "MYCN Open", 'intermediate'),
  'mono_sporadic_dn' = homer_df_id(
    file.path(my_path,  '/intermediate/sporadic_intermediate_down/knownResults.txt'), "Sporadic Close", 'intermediate'),
  'mono_sporadic_up' = homer_df_id(
    file.path(my_path, '/intermediate/sporadic_intermediate_up/knownResults.txt'), "Sporadic Open", 'intermediate')
)

homer_fold_distal <- list(
  'mono_atrx_dn' = homer_df_id(
    file.path(my_path, '/distal/atrx_distal_down/knownResults.txt'), "Atrx Close", 'distal'),
  'mono_atrx_up' = homer_df_id(
    file.path(my_path, '/distal/atrx_distal_up/knownResults.txt'), "Atrx Open", 'distal'),
  'mono_mycn_dn' = homer_df_id(
    file.path(my_path, '/distal/mycn_distal_down/knownResults.txt'), "MYCN Close", 'distal'),
  'mono_mycn_up' = homer_df_id(
    file.path(my_path, '/distal/mycn_distal_up/knownResults.txt'), "MYCN Open", 'distal'),
  'mono_sporadic_dn' = homer_df_id(
    file.path(my_path, '/distal/sporadic_distal_down/knownResults.txt'), "Sporadic Close", 'distal'),
  'mono_sporadic_up' = homer_df_id(
    file.path(my_path, '/distal/sporadic_distal_up/knownResults.txt'), "Sporadic Open", 'distal')
)
#/Volumes/GFS_BIOMEDBIO_AGFORTELNY/people/rohit/projects/clusterObj/monocytes/
homer_fold_others <- list(
  'mono_atrx_dn' = homer_df_id(
    file.path(my_path, '/others/atrx_down_others/knownResults.txt'), "Atrx Close", 'promoter'),
  'mono_atrx_up' = homer_df_id(
    file.path(my_path, '/others/atrx_up_others/knownResults.txt'), "Atrx Open",'promoter'),
  'mono_mycn_dn' = homer_df_id(
    file.path(my_path, '/others/mycn_down_others/knownResults.txt'), "MYCN Close",'promoter'),
  'mono_mycn_up' = homer_df_id(
    file.path(my_path, '/others/mycn_up_others/knownResults.txt'), "MYCN Open",'promoter'),
  'mono_sporadic_dn' = homer_df_id(
    file.path(my_path,  '/others/sporadic_down_others/knownResults.txt'), "Sporadic Close",'promoter'),
  'mono_sporadic_up' = homer_df_id(
    file.path(my_path, '/others/sporadic_up_others/knownResults.txt'), "Sporadic Open",'promoter')
)

homerdb <- list('homer_promoter_db' = do.call(rbind, homer_fold_promoter),
               'homer_intermediate_db' = do.call(rbind, homer_fold_intermediate), 
               'homer_distal_db' = do.call(rbind, homer_fold_distal))

motifdb <- do.call(rbind, homerdb)

motifdb$motif_class <- gsub(',', ':', motifdb$motif_class)
write.csv(motifdb, 
          file = file.path(my_path, 'homer_motif_analyses.txt'), 
          row.names=FALSE, quote = FALSE, col.names = TRUE, sep = "\t" )


# scratch -----------------------------------------------------------------



homer_fold_promoter_goi <- list(
  'mono_atrx_dn' = homer_df(
    file.path(my_path, '/promoter/atrx_promoter_down_goi/knownResults.txt'), "Atrx Close"),
  'mono_atrx_up' = homer_df(
    file.path(my_path, '/promoter/atrx_promoter_up_goi/knownResults.txt'), "Atrx Open"),
  'mono_mycn_dn' = homer_df(
    file.path(my_path, '/promoter/mycn_promoter_down_goi/knownResults.txt'), "MYCN Close"),
  'mono_mycn_up' = homer_df(
    file.path(my_path, '/promoter/mycn_promoter_up_goi/knownResults.txt'), "MYCN Open"),
  'mono_sporadic_dn' = homer_df(
    file.path(my_path,  '/promoter/sporadic_promoter_down_goi/knownResults.txt'), "Sporadic Close"),
  'mono_sporadic_up' = homer_df(
    file.path(my_path, '/promoter/sporadic_promoter_up_goi/knownResults.txt'), "Sporadic Open")
)

homer_fold_intermediate_goi <- list(
  'mono_atrx_dn' = homer_df(
    file.path(my_path, '/intermediate/atrx_intermediate_down_goi/knownResults.txt'), "Atrx Close"),
  'mono_atrx_up' = homer_df(
    file.path(my_path, '/intermediate/atrx_intermediate_up_goi/knownResults.txt'), "Atrx Open"),
  'mono_mycn_dn' = homer_df(
    file.path(my_path, '/intermediate/mycn_intermediate_down_goi/knownResults.txt'), "MYCN Close"),
  'mono_mycn_up' = homer_df(
    file.path(my_path, '/intermediate/mycn_intermediate_up_goi/knownResults.txt'), "MYCN Open"),
  'mono_sporadic_dn' = homer_df(
    file.path(my_path,  '/intermediate/sporadic_intermediate_down_goi/knownResults.txt'), "Sporadic Close"),
  'mono_sporadic_up' = homer_df(
    file.path(my_path, '/intermediate/sporadic_intermediate_up_goi/knownResults.txt'), "Sporadic Open")
)
homer_fold_distal_goi <- list(
  'mono_atrx_dn' = homer_df(
    file.path(my_path, '/distal/atrx_distal_down_goi/knownResults.txt'), "Atrx Close"),
  'mono_atrx_up' = homer_df(
    file.path(my_path, '/distal/atrx_distal_up_goi/knownResults.txt'), "Atrx Open"),
  'mono_mycn_dn' = homer_df(
    file.path(my_path, '/distal/mycn_distal_down_goi/knownResults.txt'), "MYCN Close"),
  'mono_mycn_up' = homer_df(
    file.path(my_path, '/distal/mycn_distal_up_goi/knownResults.txt'), "MYCN Open"),
  'mono_sporadic_dn' = homer_df(
    file.path(my_path, '/distal/sporadic_distal_down_goi/knownResults.txt'), "Sporadic Close"),
  'mono_sporadic_up' = homer_df(
    file.path(my_path, '/distal/sporadic_distal_up_goi/knownResults.txt'), "Sporadic Open")
)


network_nodes_promoters <- 
  ordered_promoters_rev[,c('motif_name', 
                           'p_value', 
                           'cell_type',
                           'log_p_value', 
                           'odds_ratio', 
                           'log_odds_ratio' )]
network_nodes_promoters$tss <- "promoter"
fwrite(network_nodes_promoters, "network_nodes_promoters.csv")


network_nodes_distal <- 
  ordered_distal_rev[,c('motif_name', 
                        'p_value', 
                        'cell_type',
                        'log_p_value', 
                        'odds_ratio', 
                        'log_odds_ratio' )]
network_nodes_distal$tss <- "distal"

network_nodes <- rbind(
  network_nodes_promoters, 
  network_nodes_distal)
fwrite(network_nodes, "network_nodes.csv")









