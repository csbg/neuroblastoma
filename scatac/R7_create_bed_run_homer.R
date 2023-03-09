# Motif analysis of significant peaks -------------------------------------

read_mult_csv <- function(pattern, cluster) {
  pacman::p_load(plyr,readr)
  mydir = "/media/AGFORTELNY/people/rohit/projects/clusterObj/"
  myfiles = list.files(path=mydir, pattern=pattern, full.names=TRUE)
  dat_csv = ldply(myfiles, read_csv)
  union_peaks <- StringToGRanges(unique(dat_csv$X1))
  enriched.motifs$cluster <- cluster
  return(enriched.motifs)
}

get_enriched_motifs <- function() {
  pacman::p_load(plyr,readr)
  mydir = "/media/AGFORTELNY/people/rohit/projects/clusterObj/"
  myfiles = list.files(path=mydir, pattern=pattern, full.names=TRUE)
  dat_csv = ldply(myfiles, read_csv)
  #union_peaks <- StringToGRanges(unique(dat_csv$X1))
  enriched.motifs <- FindMotifs(
    object = nb,
    features = unique(dat_csv$X1)
  )
  enriched.motifs$cluster <- cluster
  
  return(enriched.motifs)
}


get_peak_bed <- function(pattern, outfile){
  pacman::p_load(plyr,readr)
  mydir = "/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/"
  myfiles = list.files(path=mydir, pattern=pattern, full.names=TRUE)
  dat_csv = ldply(myfiles, read_csv)
  names(dat_csv) <- c("peak", "p_val", "avg_log2FC", 
                      "pct.1", "pct.2", "p_val_adj", "cell_type" )
  union_peaks <- StringToGRanges(unique(dat_csv$peak))
  peak_bed <- data.frame(
    seqnames = seqnames(union_peaks),
    start = start(union_peaks),
    end = end(union_peaks),
    strand = strand(union_peaks))
  peak_bed$strand <- ifelse(peak_bed$start < peak_bed$end, '+','-')
  peak_bed$id <- paste0(as.character(peak_bed$seqnames), '_',
                        as.character(peak_bed$start),'_',
                        as.character(peak_bed$end) )
  peak_bed <- peak_bed %>% 
    dplyr::select(id, seqnames, start, end, strand )
  # write.table(peak_bed, file.path(mydir, outfile), 
  #             quote = FALSE, sep="\t", row.names = FALSE)
  return(peak_bed)
}

pk11_peak <- 
  get_peak_bed("11.*enriched_peaks.txt",  "11.enriched_peaks.bed" )
pk22_peak <- 
  get_peak_bed("22.*enriched_peaks.txt", "22.enriched_peaks.bed")
pk26_peak <- 
  get_peak_bed("26.*enriched_peaks.txt", "26.enriched_peaks.bed")
pk32_peak <- 
  get_peak_bed("32.*enriched_peaks.txt", "32.enriched_peaks.bed")
pk33_peak <- 
  get_peak_bed("33.*enriched_peaks.txt", "33.enriched_peaks.bed")
pk37_peak <- 
  get_peak_bed("37.*enriched_peaks.txt", "37.enriched_peaks.bed")
pk_tcell_peak <- 
  get_peak_bed("T-cells.*enriched_peaks.txt", "tcell.enriched_peaks.bed")
pk_bcell_peak <- 
  get_peak_bed("B-cells.*enriched_peaks.txt", "bcell.enriched_peaks.bed")
pk_bcell_mem_peak <- 
  get_peak_bed("B_cells_memory_34.*enriched_peaks.txt", "bcell_mem.enriched_peaks.bed")
pk_baso_peak <- 
  get_peak_bed("Basophils.*enriched_peaks.txt", "baso.enriched_peaks.bed")
pk_nb_peak <- 
  get_peak_bed("NB.*enriched_peaks.txt", "nb.enriched_peaks.bed")
pk_mono_peak <- 
  get_peak_bed("Monocytes.*enriched_peaks.txt", "Monocytes.enriched_peaks.bed")


cell_type_bed_files <- function(cell_file) {
  pacman::p_load(plyr,readr, tools)
  mydir = "/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj"
  .create_save_bed <- function(dataframe, outfile) {
    df_gr <- 
      Signac::StringToGRanges(unique(dataframe$gene))
    
    peak_bed <- data.frame(
      seqnames = seqnames(df_gr),
      start = start(df_gr),
      end = end(df_gr),
      strand = strand(df_gr))
    
    peak_bed$strand <- 
      ifelse(peak_bed$start < peak_bed$end, '+','-')
    
    peak_bed$id <- 
      paste0(as.character(peak_bed$seqnames), '_',
             as.character(peak_bed$start),'_',
             as.character(peak_bed$end) )
    peak_bed <- 
      peak_bed %>% 
      dplyr::select(id, seqnames, start, end, strand )
    write.table(peak_bed, file.path(mydir, outfile), 
                quote = FALSE, sep="\t", row.names = FALSE)
  }
  
  dat_csv1 <-
    read.csv(file.path(mydir, cell_file ), sep = ' ', header = TRUE)
  outstr <- 
    paste0('/group_microenv/',
           tolower(file_path_sans_ext(cell_file)))
  dat_csv1 %>% 
    filter(logFC_groupII < 0) %>% 
    .create_save_bed(., 
                     outfile = paste0(outstr, '_mycn_down.bed'))
  dat_csv1 %>% 
    filter(logFC_groupIII < 0) %>% 
    .create_save_bed(., 
                     outfile = paste0(outstr, '_atrx_down.bed'))
  dat_csv1 %>% 
    filter(logFC_groupIV < 0) %>% 
    .create_save_bed(., 
                     outfile= paste0(outstr, '_sporadic_down.bed'))
  dat_csv1 %>% 
    filter(logFC_groupII > 0) %>% 
    .create_save_bed(., 
                     outfile= paste0(outstr, '_mycn_up.bed'))
  dat_csv1 %>% 
    filter(logFC_groupIII > 0) %>% 
    .create_save_bed(., 
                     outfile = paste0(outstr, '_atrx_up.bed'))
  dat_csv1 %>% 
    filter(logFC_groupIV > 0) %>% 
    .create_save_bed(., 
                     outfile = paste0(outstr,'_sporadic_up.bed'))
}
for (i in c('Monocytes.csv', 'Basophils.csv', 'T-cells.csv',
            'B_cells_memory_34.csv', 'B-cells.csv', '37.csv', 
            '36.csv', '33.csv', '32.csv', 
            '26.csv', '22.csv', '11.csv')) {
  cell_type_bed_files(i) 
}


# Run Homer ---------------------------------------------------------------

# Homer run from terminal 
#!/bin/bash
# for i in *_down.bed
# do
# findMotifsGenome.pl $i hg38 $(basename $i .bed) -size 100 -mask
# done
# 
# for i in *_up.bed
# do
# findMotifsGenome.pl $i hg38 $(basename $i .bed) -size 100 -mask
# done


pk11 <- read_mult_csv("11.*enriched_peaks.txt", 'Cluster_11')
pk22 <- read_mult_csv("22.*enriched_peaks.txt", 'Cluster_22')
pk26 <- read_mult_csv("26.*enriched_peaks.txt", 'Cluster_26')
pk32 <- read_mult_csv("32.*enriched_peaks.txt", 'Cluster_32')
pk33 <- read_mult_csv("33.*enriched_peaks.txt", 'Cluster_33')
pk37 <- read_mult_csv("37.*enriched_peaks.txt", 'Cluster_37')
pk_tcell <- read_mult_csv("T-cells.*enriched_peaks.txt", 'T-cell')
pk_bcell <- read_mult_csv("B-cells.*enriched_peaks.txt", 'B-cell')
pk_bcell_mem <- read_mult_csv("B_cells_memory_34.*enriched_peaks.txt", 'B-cell-memory')
pk_baso <- read_mult_csv("Basophils.*enriched_peaks.txt", 'Basophils')
pk_mono <- read_mult_csv("Monocytes.*enriched_peaks.txt", 'Monocytes')
pk_nb <- read_mult_csv("NB.*enriched_peaks.txt", 'NB')

big_list <- list(pk11, pk22,pk26, pk32, pk33,pk37,
                 pk_tcell, pk_bcell, pk_bcell_mem, 
                 pk_baso, pk_mono, pk_nb
)
big_df <- do.call(rbind, big_list )  

write.csv(big_df,"/media/AGFORTELNY/people/rohit/projects/cell_type_motif_db.csv", 
          row.names = FALSE)

big_df <- read.csv("/media/AGFORTELNY/people/rohit/projects/cell_type_motif_db.csv",
                   header = TRUE) 

# moi <- c('AP', 'NFAT', 'AP', 'YY-1', 'STAT', 'GATA', 
#          'NFKB', 'IL2Ra', 'ATF-2', 'CREB', 'CREB', 'Hlf', 
#          'Ik', 'POU2F1', 'POU2F1', 'POU2F1b', 'POU2F1c', 
#          'EGR1',  'BRG', 'LITAF', 'ETS', 'SP1', 
#          'c-fos', 'c-jun', 'c-Fos', 'c-Jun', 'p53', 'TBP', 
#          'CTCF', 'MAZ1', 'hnRNPK', 'CNBP',  'hnRNPA1', 'FBP', 
#          'FIR', 'TFIIH', 'GCNF', 'GCNF', 'GCNF', 'AML1a', 
#          'Myb', 'Myc', 'Max1', 'Sox9', 'Tal1beta', 'USF-1', 
#          'USF1', 'E2F',  'REST','NRSF',
#          'IRF', 'MyoD',  'HOXA9', 'HOXA9B', 'Meis', 
#          'Pax-4a', 'POU3F2', 'Arnt','HNF', 'Max', 
#          'MZF-1', 'SRY', 'TFIID', 'E47', 'Evi', 'FOXD3', 'Ik2', 'Pax')
# delta, CBP/p300, C/EBP, nm23/puf60, NRSF form 1, NRSF form 2,

# detect_moi <- motif_names[str_detect(motif_names, paste(moi, collapse = "|"))]
# # Heatmap with -log10(pvalue)
# big_df %>% 
#   dplyr::select(motif, fold.enrichment, pvalue, motif.name, cluster ) %>% 
#   mutate(pvalue=replace(pvalue, pvalue==0, '1e-321')) %>% View()
# mutate(logp = -log10(as.numeric(pvalue))) %>%
#   pivot_wider(c('motif.name', 'cluster', 'logp' ),
#               names_from = cluster, values_from =logp) %>% 
#   dplyr::filter(motif.name %in% detect_moi ) %>% 
#   tibble::column_to_rownames(., var = 'motif.name') %>% 
#   as.matrix() %>% 
#   pheatmap::pheatmap()
# 
# # Heatmap with logFC
# big_df %>% 
#   dplyr::select(motif, fold.enrichment, pvalue, motif.name, cluster ) %>% 
#   mutate(pvalue=replace(pvalue, pvalue==0, '1e-321')) %>% 
#   mutate(logp = -log10(as.numeric(pvalue))) %>%
#   pivot_wider(c('motif.name', 'cluster', 'fold.enrichment' ),
#               names_from = cluster, values_from = fold.enrichment) %>% 
#   dplyr::filter(motif.name %in% detect_moi ) %>% 
#   dplyr::select(-Cluster_32) %>% 
#   tibble::column_to_rownames(., var = 'motif.name') %>% 
#   as.matrix() %>% 
#   pheatmap::pheatmap()


