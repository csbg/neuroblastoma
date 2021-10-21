######################################################
# Step 1: 
####################################################

db <- dbs[1:22]

EnrData_Cl8 <- list(
  "1" = DEG_2_1,
  "2" = DEG_2_2,
  "3" = DEG_2_3,
  "4" = DEG_2_4,
  "5" = DEG_2_5,
  "6" = DEG_2_6,
  "7" = DEG_2_7,
  "8" = DEG_2_8
) %>% 
  map_dfr(~enrich_genes(.$gene_short_name, databases = db), .id = "cluster")

EnrData_Cl8_altMeth <-list(
  "1" = cluster_altmeth_expression %>% filter(cell_group==1 ),
  "2" = cluster_altmeth_expression %>% filter(cell_group==2 ),
  "3" = cluster_altmeth_expression %>% filter(cell_group==3 ),
  "4" = cluster_altmeth_expression %>% filter(cell_group==4 ),
  "5" = cluster_altmeth_expression %>% filter(cell_group==5 ),
  "6" = cluster_altmeth_expression %>% filter(cell_group==6 ),
  "7" = cluster_altmeth_expression %>% filter(cell_group==7 ),
  "8" = cluster_altmeth_expression %>% filter(cell_group==8 ) 
) %>% 
  map_dfr(~enrich_genes(.$gene_short_name, databases = db), .id = "cluster")


EnrData_Cl4 <- list(
  "1" = DEG_1,
  "2" = DEG_2,
  "3" = DEG_3,
  "4" = DEG_4 
  )%>% 
  map_dfr(~enrich_genes(.$gene_short_name, databases = db), .id = "cluster")


EnrData_Cl4_altMeth <-list(
  "1" = cluster_altmeth_expression_2 %>% filter(cell_group==1 ),
  "2" = cluster_altmeth_expression_2 %>% filter(cell_group==2 ),
  "3" = cluster_altmeth_expression_2 %>% filter(cell_group==3 ),
  "4" = cluster_altmeth_expression_2 %>% filter(cell_group==4 ) 
)%>% 
  map_dfr(~enrich_genes(.$gene_short_name, databases = db), .id = "cluster")


#Helper function enrichClusterPlot top3 / top5 (enrichment data, db number )
#Top3 by EST Enrichment Dotplots

## Enrichment plots using plotenrich_dotplot() function from W.E.S.
#####################################################################
# aggregate enrichment data & transform to output format

EnrData <- list(Cl8 = EnrData_Cl8 %>% # Cluster 8 standard
                  mutate(
                    "contrast" = 1,
                    "direction" = "up",
                    "Overlap" = EnrData_Cl8$Overlap 
                    ), 
                Cl8a = EnrData_Cl8_altMeth %>% # Cluster 8 alternative Method
                  mutate("contrast"=1,
                         "direction"="up",
                         "Overlap"=EnrData_Cl8_altMeth$Overlap
                         ),   
                Cl4 = EnrData_Cl4 %>% # Cluster 4 standard
                  mutate(
                    "contrast" = 1,
                    "direction" = "up",
                    "Overlap" = EnrData_Cl4$Overlap
                    ),            
                Cl4a = EnrData_Cl4_altMeth %>% # Cluster 4 alternative Method
                  mutate("contrast" = 1,"direction"="up",
                         "Overlap" = EnrData_Cl4_altMeth$Overlap
                         ),
                Selected= 
                  sel <- EnrData_Cl4 %>% 
                    mutate(
                      "contrast" = 1,
                      "direction" = "up",
                      "Overlap" = EnrData_Cl4$Overlap,
                      "Term" = rename_EnrTerms(Term),
                      ) %>%
                    filter(db %in% c(
                         "NCI-Nature_2016",
                         "CORUM",
                         "GO_Biological_Process_2018",
                         "KEGG_2019_Human",
                         "MSigDB_Hallmark_2020"
                         )
                    ) %>%
                  filter(
                    Term %in% 
                      c("TNF-alpha signalling via NF-kB",
                        "dopaminergic neuron differentiation",
                        "IL2 effects mediated by PI3K",
                        "IL3-mediated signaling events",
                        "IL6-mediated signaling events",
                        "Trk receptor signaling mediated by PI3K and PLC-gamma",
                        "PDGFR-beta signaling pathway",
                        "ERG-JUN-FOS DNA-protein complex",
                        "JUND-FOSB-SMAD3-SMAD4 complex",
                        "CALM1-FKBP38-BCL2 complex",
                          
                        "Cell Cycle",
                        "DNA replication",
                        "DNA repair",
                        "PCNA-DNA polymerase delta complex",
                        "nucleoside metabolic process",
                        "DNA synthesome complex",
                        "6S methyltransferase and RG-containing Sm proteins complex",
                        "BRD4-RFC complex",
                        "histone mRNA metabolic process",
                        "DNA synthesome complex (17 subunits)",
                        
                        
                        "intermediate filament bundle assembly",
                        "negative regulation of glial cell proliferation",
                        "peripheral nervous system neuron development",
                        "noradrenergic neuron differentiation",
                        "myelin maintenance",
                        "BRAF-RAF1-13-3-3 complex",
                        "kinase maturation complex 1", 
                        "IKBKB-CDC37-KIAA1967-HSP90 complex", 
                        "ITGA3-ITB1-BSG complex",
                        
                        "Regulation of Telomerase",
                        "p53-SP1 complex",
                        "CAND1-CUL4A-RBX1 complex", 
                        "Er-alpha-p53-hdm2 compex", 
                        "NUMB-TP53-MDM2 complex"
                        )
                  )
                )


Method <- c("_Cl8",
            "_Cl8_AltMethod",
            "_Cl4",
            "_Cl4_AltMethod")

folder <- c("/pngEnrichments/CL8/stdMethod/",
            "/pngEnrichments/CL8/altMethod/",
            "/pngEnrichments/CL4/stdMethod/",
            "/pngEnrichments/CL4/altMethod/")

# Plotting Loop:
#######################################################
for(i in 1:length(EnrData)){  
data <- EnrData[i] %>% as_tibble()
colnames(data) <- "y"
data <- data$y  

clusterTemp <- Method[i]
folderTemp <- folder[i]
 
data$cluster <- as.factor(data$cluster)


dbs[1] #GO Molecular Function
plot_enrichr_dots(data,
                  db = dbs[1],
                  min_odds_ratio = 20,
                  filename =  path(paste0(folderTemp,dbs[2],clusterTemp)))

dbs[2] #GO Biological Process
plot_enrichr_dots(data,
                  db = dbs[2],
                  min_odds_ratio = 50,
                  filename = path(paste0(folderTemp,dbs[2],clusterTemp)),
                  height=500)

dbs[3] # GO Cellular Component
plot_enrichr_dots(data,
                  db = dbs[3],
                  min_odds_ratio = 20,
                  filename = path(paste0(folderTemp,dbs[3],clusterTemp)))

dbs[4] #WikiPathways_2019_Human
plot_enrichr_dots(data,
                  db = dbs[4],
                  min_odds_ratio = 10,
                  filename = path(paste0(folderTemp,dbs[4],clusterTemp)),
                  height=280)

dbs[5] #Panther_2016
plot_enrichr_dots(data,
                  db = dbs[5],
                  min_odds_ratio = 1,
                  filename = path(paste0(folderTemp,dbs[5],clusterTemp)))

dbs[6] #KEGG_2019_Human
plot_enrichr_dots(data,
                  db = dbs[6],
                  min_odds_ratio = 10,
                  filename = path(paste0(folderTemp,dbs[6],clusterTemp)))

dbs[7] #NCI-Nature_2016
plot_enrichr_dots(data,
                  db = dbs[7],
                  min_odds_ratio = 10,
                  filename = path(paste0(folderTemp,dbs[7],clusterTemp)))

dbs[8] #TRRUST_Transcription_Factors_2019
plot_enrichr_dots(data,
                  db = dbs[8],
                  min_odds_ratio = 10,
                  filename = path(paste0(folderTemp,dbs[8],clusterTemp)))

dbs[9]#ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X
plot_enrichr_dots(data,
                  db = dbs[9],
                  min_odds_ratio = 5,
                  filename = path(paste0(folderTemp,dbs[9],clusterTemp)))

dbs[10]#TRANSFAC_and_JASPAR_PWMs
plot_enrichr_dots(data,
                  db = dbs[10],
                  min_odds_ratio = 1,
                  filename = path(paste0(folderTemp,dbs[10],clusterTemp)))

dbs[11]#MSigDB_Hallmark_2020
plot_enrichr_dots(data,
                  db = dbs[11],
                  min_odds_ratio = 1,
                  filename = path(paste0(folderTemp,dbs[11],clusterTemp)))

dbs[12]#Corum
plot_enrichr_dots(data,
                  db = dbs[12],
                  min_odds_ratio = 50,
                  filename = path(paste0(folderTemp,dbs[12],clusterTemp)),
                  height=340)

dbs[13] #NCI-60_Cancer_Cell_Lines
plot_enrichr_dots(data,
                  db = dbs[13],
                  min_odds_ratio = 5,
                  filename = path(paste0(folderTemp,dbs[13],clusterTemp)))

dbs[14] #Chromosome_Location
plot_enrichr_dots(data,
                  db = dbs[14],
                  min_odds_ratio = 0,
                  filename = path(paste0(folderTemp,dbs[14],clusterTemp)))

dbs[15]#ClinVar_2019
plot_enrichr_dots(data,
                  db = dbs[15],
                  min_odds_ratio = 5,
                  filename = path(paste0(folderTemp,dbs[15],clusterTemp)))

dbs[16]#Human_Gene_Atlas
plot_enrichr_dots(data,
                  db = dbs[16],
                  min_odds_ratio = 5,
                  filename = path(paste0(folderTemp,dbs[16],clusterTemp)))

dbs[17]#Jensen_TISSUES
plot_enrichr_dots(data,
                  db = dbs[17],
                  min_odds_ratio = 20,
                  filename = path(paste0(folderTemp,dbs[17],clusterTemp)))

dbs[18] #Jensen_COMPARTMENTS
plot_enrichr_dots(data,
                  db = dbs[18],
                  min_odds_ratio = 50,
                  filename = path(paste0(folderTemp,dbs[18],clusterTemp)))

dbs[19] #Jensen_DISEASES
plot_enrichr_dots(data,
                  db = dbs[19],
                  min_odds_ratio = 20,
                  filename = path(paste0(folderTemp,dbs[19],clusterTemp)))

dbs[20]#Reactome_2016
plot_enrichr_dots(data,
                  db = dbs[20],
                  min_odds_ratio = 50,
                  filename = path(paste0(folderTemp,dbs[20],clusterTemp)))

dbs[21] #ChEA_2016
plot_enrichr_dots(data,
                  db = dbs[21],
                  min_odds_ratio = 3,
                  filename = path(paste0(folderTemp,dbs[21],clusterTemp)))

dbs[22] #ENCODE_TF_ChIP-seq_2015
plot_enrichr_dots(data,
                  db = dbs[22],
                  min_odds_ratio = 5,
                  filename = path(paste0(folderTemp,dbs[22],clusterTemp)))
}

source("~/code/resources/RFunctions/Basics.R")


##############################################################################
# Enrichment of patients in clusters = Relative Abundance enrichment
metadata <- readRDS(file = paste0(path,"/NB_CellMetadata_Cl4.rds"))
metadata <- data.table(as.data.frame(metadata))

metadata <- data.table(as.data.frame(vCDS1@colData)) # 8 clusters
metadata <- data.table(as.data.frame(colData(aligned_cds))) #4 clusters

pDT <- setNames(metadata[,c("sample", "group", "recluster"), with=F], c("Patient", "Group", "Cluster"))

for(colx in colnames(pDT)){
  pDT[[colx]] <- as.character(pDT[[colx]])
}
pDT <-pDT[,.N, by=c("Patient", "Group", "Cluster")]
str(pDT)
res <- data.table()
cx <- pDT$Cluster[1]
patx <- pDT$Patient[1]
for(cx in unique(pDT$Cluster)){
  for(patx in unique(pDT$Patient)){
    cont.matrix <- matrix(nrow=2, ncol=2, c(
      sum(pDT[Patient == patx & Cluster == cx]$N),
      sum(pDT[Patient != patx & Cluster == cx]$N),
      sum(pDT[Patient == patx & Cluster != cx]$N),
      sum(pDT[Patient != patx & Cluster != cx]$N)
    ))
    fish <- fisher.test(cont.matrix)
    res <- rbind(res, data.table(Cluster=cx, Patient=patx, p=fish$p.value, OR=fish$estimate))
  }
}
res[,padj := p.adjust(p, method="BH")]
res <- merge(res, unique(pDT[,c("Patient", "Group"),with=F]), by="Patient")

Enr_Patients <- res %>% 
  transmute(
  Patient = rename_patients(Patient),
  group = rename_groups(Group),
  cluster = Cluster,
  Adjusted.P.value = padj,
  Odds.Ratio = OR
  )





#####

Enr_PatientCluster <- ggplot(res, aes(
  x=factor(Cluster),
  y=Patient, 
  size=pmin(-log10(padj), 5))
) + 
  geom_point(aes(color =log2(OR + min(res[OR != 0]$OR))),shape=16) +
  geom_point(data = res %>% filter(sig == 1), shape = 1) +
  scale_color_gsea(name= TeX("log_{2} odds ratio"))+
  xlab("Cluster")+
  #scale_color_gradient2(name="log2(Odds.Ratio)", low="blue", high="red",mid="white",midpoint=0) +
  scale_size_continuous(name=TeX("-log_{10} p_{adj}")) + 
  facet_grid(Group ~ ., scales = "free", space = "free") + 
  theme_nb(grid=F)+
  theme(legend.key.height = unit(0.5,"cm"))


ggsave(plot=Enr_PatientCluster,filename=paste0(path_plot,"/MarkerEnrichment/nb_Patients_Clusters_Enrich_Cl4.pdf"),
       height=8, width=8, 
       units = "cm")




# Enrichment of cell-type markers in cluster
###############################################################################

# Enrichment of external markers in DE ------------------------------------------
#Select background
bg <- list(Cl4=unique((DEG_1%>%rbind(DEG_2,DEG_3,DEG_4))$gene_short_name), #Cl4
           Cl4a=unique(cluster_altmeth_expression_2$gene_short_name), # Cl4 Alt.Meth.
          Cl8=unique((DEG_2_1%>%rbind(DEG_2_2,DEG_2_3,DEG_2_4,DEG_2_5,DEG_2_6,DEG_2_7,DEG_2_8))$gene_short_name), #Cl8
          Cl8a=unique(cluster_altmeth_expression$gene_short_name)) # Cl8 Alt.Meth.

enrDT <- list(
  Cl4 = DEG_1 %>% #EnrichmentData Method 1-Cl4
    rbind(DEG_2,DEG_3,DEG_4) %>% 
    ungroup(), 
  Cl4a = cluster_altmeth_expression_2 %>% #EnrichmentData Method 2-Cl4
    rename("cluster"="cell_group") %>% 
    as_tibble(), 
  Cl8 = DEG_2_1 %>%  #EnrichmentData Method 1-Cl8
    rbind(DEG_2_2,DEG_2_3,DEG_2_4,DEG_2_5,DEG_2_6,DEG_2_7,DEG_2_8) %>% 
    ungroup(), 
  Cl8a = cluster_altmeth_expression %>%  #EnrichmentData Method 2-Cl8
    rename("cluster"="cell_group") %>% 
    as_tibble() 
  ) 
###############################################################################################################
 #    Markers Kildisiute
###

name <- c("Cl4","Cl4a","Cl8","Cl8a")

table(markers_1$dataset)

for(i in 1:length(bg)){
  temp_bg <- bg[i] %>% unlist()
  temp_enrDT <- enrDT[i] %>% as_tibble() 
  colnames(temp_enrDT) <- "y"
  temp_enrDT <- temp_enrDT$y

  enrRes <- data.table() #empty Results table
  
  #Markers Kildisiute et al.
  marker.kild <- unique(markers_1%>%filter(dataset=="adrenal_gland"|dataset=="stringent_adrenal"))%>%
    filter(gene %in% temp_bg) %>% 
    rename("ID" = "cluster")%>%
    filter(ID!="Podocytes")  # remove....just 1 gene from Podos
  
  cx <- marker.kild[1]
  mx <- marker.kild$ID[1]
  
  for(cx in unique(temp_enrDT$cluster)){
    for(mx in unique(marker.kild$ID)){
      g.cx <- (temp_enrDT %>% filter(cluster == cx))$gene_short_name
      g.mx <- (marker.kild%>% filter(ID == mx))$gene
      
      A <- length(intersect(g.cx, g.mx))
      B <- length(intersect(g.cx, temp_bg)) - A
      C <- length(intersect(g.mx, temp_bg)) - A
      D <- length(temp_bg) - A - B - C
      stopifnot(A + B + C + D == length(temp_bg))
      
      cont.matrix <- matrix(nrow=2, ncol=2, c(A,B,C,D))
      fish <- fisher.test(cont.matrix)
      enrRes <- rbind(enrRes, data.table(Cluster=cx, MarkerList=mx, p=fish$p.value, OR=fish$estimate, Genes=paste(intersect(g.cx, g.mx), collapse=",")))
    }
  }
  
  enrRes[,padj := p.adjust(p, method="BH")]
  enrRes[padj < 0.05]
  
  enrRes <- enrRes%>%mutate(sig=if_else(padj<0.05,1,0))
  
  p1 <- ggplot(enrRes, aes(
    x=factor(Cluster),
    y=MarkerList, 
    #OR + min(enrRes[OR != 0]$OR)
    size=pmin(-log10(padj), 5)
  )) + 
    geom_point(aes(color =log2(OR + min(enrRes[OR != 0]$OR)))) +
    geom_point(data = enrRes %>% filter(sig == 1), shape = 1) +
    scale_color_gradient2(name="log2OR", low="blue", high="red") +
    scale_size_continuous(name="padj") + 
    theme_bw(12)+
    xlab("Cluster") 
#    scale_shape_manual(values=c(22,21),guide = NULL)
  
ggsave(paste0(path_plot,"/MarkerEnrichment/nb_Enrich_DEG_in_markers_KildisiuteAdrenal_",name[i],".pdf"),plot=p1, w=5, h=4)

}


##############################
# Markers Jansky et al.

temp_bg <- bg[1] %>% unlist() # Cl4-DE background

temp_enrDT <- enrDT[[1]] %>% as_tibble() # Cl4- DE results

table(temp_enrDT$cluster)
  
enrRes <- data.table() #empty Results table
  
marker.jansky <- 
  markers_2 %>% 
  filter(gene %in% temp_bg) %>% 
  rename("ID" = "cluster")
  
cx <- marker.jansky[1]
mx <- marker.jansky$ID[1]
  
for(cx in unique(temp_enrDT$cluster)){
    for(mx in unique(marker.jansky$ID)){
      g.cx <- (temp_enrDT %>% filter(cluster == cx))$gene_short_name
      g.mx <- (marker.jansky%>% filter(ID == mx))$gene
      
      A <- length(intersect(g.cx, g.mx))
      B <- length(intersect(g.cx, temp_bg)) - A
      C <- length(intersect(g.mx, temp_bg)) - A
      D <- length(temp_bg) - A - B - C
      stopifnot(A + B + C + D == length(temp_bg))
      
      cont.matrix <- matrix(nrow=2, ncol=2, c(A,B,C,D))
      fish <- fisher.test(cont.matrix)
      enrRes <- rbind(enrRes, data.table(Cluster=cx, MarkerList=mx, p=fish$p.value, OR=fish$estimate, Genes=paste(intersect(g.cx, g.mx), collapse=",")))
    }
  }
  
enrRes[,padj := p.adjust(p, method="BH")]
enrRes[padj < 0.05]
  
Enr_Jansky_Cl4 <- enrRes %>%
    transmute(
      cluster = Cluster,
      Marker = MarkerList,
      Odds.Ratio = OR,
      Adjusted.P.value = padj,
      genes = Genes
    )
  
table(Enr_Jansky_Cl4$Marker)  

#### 
  p1 <- ggplot(enrRes, aes(
    x=factor(Cluster),
    y=MarkerList, 
    size=pmin(-log10(padj), 5)
    ))+
    geom_point(aes(color =log2(OR + min(enrRes[OR != 0]$OR)))) +
    geom_point(data = enrRes %>% filter(sig == 1), shape = 1) +
    scale_color_gsea(name= TeX("log_{2} odds ratio"))+
    scale_size_continuous(name=TeX("-log_{10} p_{adj}")) + 
    theme_nb(grid=F)+
    theme(legend.key.height = unit(0.5,"cm"))+
    xlab("Cluster")+
    ylab("")

  
  ggsave(paste0(path_plot,"/MarkerEnrichment/nb_Enrich_DEG_in_markers_Jansky_",name[i],".pdf"),plot=p1, w=4, h=4)


  
  
  
  
  
  
  
  
  
  
##########################
# unused
# Plotting of DE results --------------------------------------------------

ff <- list.files("results/", pattern="NB_Regression", full.names = TRUE)
ff <- setNames(lapply(ff, fread), gsub(".*_(CL\\d+)\\.csv", "\\1", basename(ff)))
sapply(ff, nrow)/7

str(ff[[1]])

ff[[1]][,.N, by=c("term", "gene_short_name")][N > 1]
ff[[2]][,.N, by=c("term", "gene_short_name")][N > 1]
ff[[1]][gene_short_name == "HSPA14"]

cx <- names(ff)[1]
res <- data.table()
for(cx in names(ff)){
  fx <- copy(ff[[cx]])
  fx[, q_value2 := as.numeric(gsub(",", ".", q_value))]
  fx[is.na(q_value2)]
  fx[, normalized_effect2 := as.numeric(gsub(",", ".", estimate))]
  fx[is.na(normalized_effect2)]
  fx <- fx[,.(sig=sum(q_value2 < 0.05 & normalized_effect2 > 0, na.rm = TRUE), effect=mean(normalized_effect2, na.rm=TRUE)), by="gene_short_name"]
  res <- rbind(res, data.table(fx, cluster=cx))
}


res[sig >= 3][,.N, by="cluster"]
res[sig >= 4][,.N, by="cluster"]

with(res[sig >= 3], split(gene_short_name, cluster))

x <- jaccard(with(res[sig >= 3], split(gene_short_name, cluster)))
diag(x) <- NA
pheatmap(x)

res[sig >= 3][,.N, by="cluster"]
write.tsv(res[sig >= 3], "nf/markers.tsv")



enrRes_jansky <- data.table()
cx <- marker.jansky[1]
mx <- marker.jansky$ID[1]
for(cx in unique(enrDT$cluster)){
  for(mx in unique(marker.jansky$ID)){
    g.cx <- (enrDT %>% filter(cluster == cx))$gene_short_name
    g.mx <- (marker.jansky%>% filter(ID == mx))$gene
    
    A <- length(intersect(g.cx, g.mx))
    B <- length(intersect(g.cx, bg)) - A
    C <- length(intersect(g.mx, bg)) - A
    D <- length(bg) - A - B - C
    stopifnot(A + B + C + D == length(bg))
    
    cont.matrix <- matrix(nrow=2, ncol=2, c(A,B,C,D))
    fish <- fisher.test(cont.matrix)
    enrRes_jansky <- rbind(enrRes_jansky, data.table(Cluster=cx, MarkerList=mx, p=fish$p.value, OR=fish$estimate, Genes=paste(intersect(g.cx, g.mx), collapse=",")))
  }
}

enrRes_jansky[,padj := p.adjust(p, method="BH")]
enrRes_jansky[padj < 0.05]

ggplot(enrRes_jansky, aes(
  x=factor(Cluster),
  y=MarkerList, 
  color=log2(OR + min(enrRes[OR != 0]$OR)), 
  size=pmin(-log10(padj), 5))
) + 
  geom_point(shape=16) +
  scale_color_gsea(name= TeX("log_{2} odds ratio"))+
  scale_size_continuous(name=TeX("-log_{10} p_{adj}")) +
  theme_nb(grid=F)+
  xlab("Cluster")

ggsave(paste0(path_plot,"/MarkerEnrichment/nb_Enrich_DEG_in_markers_Jansky_Cl8a.pdf"), w=5, h=4)



##############################################
##############################################
#Module Enrichment

#group genes into modules (UMAP run on genes + knn) &
# enrich gene modules
#################################################################
#Find gene modules
Reg_Mod <- find_gene_modules(vCDS[cluster_aggr_expression$gene_short_name,],
                             resolution = 1e-2)
#join with cluster_aggr_expr
Reg_mod_join <- Reg_Mod %>%rename("gene_short_name"=id) %>%
  left_join(cluster_aggr_expression,by="gene_short_name")

#Set DB to enrich
(dbs)
(db <- dbs[4])

#enrich found gene modules (1. data, 2. number of modules)
Reg_Mod_Enr <- enrichWrite(Reg_mod_join %>% rename(id=gene_short_name),14)

#     Summarize all duplicated genes-terms (same genes from different clusters) 
#     term is chosen according to max(Estimate)

Reg_ModGenes_Enr <- Reg_Mod_Enr 

Reg_Mod_Enr_Top3 <- tibble(id=Reg_ModGenes_Enr$id, 
                           enr=Reg_ModGenes_Enr$enr.1, 
                           Odds.Ratio=Reg_ModGenes_Enr$Odds.Ratio.1, 
                           Adj.P.val=Reg_ModGenes_Enr$Adjusted.P.value.1) %>% 
  add_row(Reg_ModGenes_Enr %>% select(id, enr.2, Odds.Ratio.2, Adjusted.P.value.2) %>% 
            rename(enr = enr.2, Odds.Ratio=Odds.Ratio.2, Adj.P.val=Adjusted.P.value.2))  %>% 
  add_row(Reg_ModGenes_Enr %>% select(id, enr.3, Odds.Ratio.3, Adjusted.P.value.3)%>% 
            rename(enr = enr.3, Odds.Ratio=Odds.Ratio.3, Adj.P.val=Adjusted.P.value.3)) %>% 
  left_join(Reg_ModGenes_Enr%>% select(id,module,cluster,norm_effect,q_val,est),copy = TRUE)

table(duplicated(Reg_Mod_Enr_Top3$enr))
table(duplicated(Reg_Mod_Enr_Top3$id))

