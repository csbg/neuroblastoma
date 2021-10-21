# Import Data
source("nb_Packages.R")
source("styling.R")
source("nb_common_functions.R")

cds <- readRDS(file = paste0(path,"/nb_CDS.rds"))
metadata <- readRDS(file = paste0(path,"/nb_Metadata_tumor.rds"))
DEG <- readRDS(file = paste0(path,"/nb_DEG.rds"))

sig.jansky <- 
  readxl::read_excel(
  paste0(path, "/Jansky_et_al_SuppMat_Tables.xlsx"),
  sheet = "Suppl. Table 5",
  skip = 1
  )
#      Jansky et al. supplementary tables
# https://www.nature.com/articles/s41588-021-00806-1

##############################################################################
# Enrichment of patients in clusters ~ "Relative Abundance testing"

metadata <- data.table(as.data.frame(metadata))

pDT <- 
  setNames(
  metadata[ , c("sample", "group", "cluster"), with = F], 
  c("Patient", "Group", "Cluster")
  )

for (colx in colnames(pDT)) {
  pDT[[colx]] <- as.character(pDT[[colx]])
}

pDT <- pDT[ ,.N,
           by = c("Patient", "Group", "Cluster")
           ]

res <- data.table()
cx <- pDT$Cluster[1]
patx <- pDT$Patient[1]

# create contigency matrix
for (cx in unique(pDT$Cluster)) {
  for (patx in unique(pDT$Patient)) {
    cont.matrix <- 
      matrix(
        nrow = 2,
        ncol = 2, 
        c(
          sum(pDT[Patient == patx & Cluster == cx]$N),
          sum(pDT[Patient != patx & Cluster == cx]$N),
          sum(pDT[Patient == patx & Cluster != cx]$N),
          sum(pDT[Patient != patx & Cluster != cx]$N)
          )
        )
    # fisher's exact test
    fish <- fisher.test(cont.matrix)
    
    res <- 
      rbind(res,
            data.table(
              Cluster = cx,
              Patient = patx,
              p = fish$p.value,
              OR = fish$estimate)
            )
                 
  }
}

# Adjust p-value with Benjamini-Hochberg method
res[,padj := p.adjust(p, method="BH")]

res <- 
  merge(res,
        unique(pDT[ , c("Patient", "Group"), with = F]),
        by = "Patient"
        )

Patient_Cluster_Enrichment <- res %>% 
  transmute(
    Patient = rename_patients(Patient),
    Group = rename_groups(Group),
    Cluster = Cluster,
    Adjusted.P.value = padj,
    Odds.Ratio = OR
  ) %>%
  as_tibble()

saveRDS(Patient_Cluster_Enrichment,
        file = paste0(path, "/nb_Patient_Enrichment_Cluster.rds")
        )


################################################################################
# cluster specific DEG Enrichment in Literature Gene Signature

#Select background
bg <- unique(DEG$gene_short_name)
          
enrDT <- DEG %>% 
  as.data.table()

sig.jansky <- 
  sig.jansky %>% 
  filter(gene %in% bg) 

cx <- sig.jansky[1]
mx <- sig.jansky$ID[1]
enrRes <- data.table()

for(cx in unique(enrDT$cluster)) {
  for(mx in unique(sig.jansky$ID)) {
    g.cx <- (enrDT %>% filter(cluster == cx))$gene_short_name
    g.mx <- (sig.jansky %>% filter(ID == mx))$gene
    
    A <- length(intersect(g.cx, g.mx))
    B <- length(intersect(g.cx, bg)) - A
    C <- length(intersect(g.mx, bg)) - A
    D <- length(bg) - A - B - C
    
    stopifnot(A + B + C + D == length(bg))
    
    cont.matrix <- 
      matrix(
        nrow = 2, 
        ncol = 2, 
        c(A,B,C,D)
        )
    fish <- fisher.test(cont.matrix)
    enrRes <- rbind(enrRes, 
                    data.table(Cluster = cx, 
                               MarkerList = mx, 
                               p = fish$p.value, 
                               OR = fish$estimate, 
                               Genes = paste(
                                        intersect(g.cx, g.mx), 
                                        collapse = ","
                                        )
                               )
                    )
  }
}
# Adjust p-value by B.-H. method
enrRes[,padj := p.adjust(p, method = "BH")]

nb_DEG_Enrichment_SigJansky <- 
  enrRes %>%
  transmute(
    Cluster = Cluster,
    Signature = MarkerList,
    Odds.Ratio = OR,
    Adjusted.P.value = padj,
    Genes = Genes
  ) %>% 
  as_tibble()

saveRDS(nb_DEG_Enrichment_SigJansky,
        file = paste0(path, "/nb_DEG_Enrichment_Jansky.rds")
)

################################################################################
# Gene Set Enrichment Analysis (GSEA) on cluster specific DEG

db <- dbs[1:26]

GSEA <- 
  list(
  "1" = DEG %>% filter(cluster == 1),
  "2" = DEG %>% filter(cluster == 2),
  "3" = DEG %>% filter(cluster == 3),
  "4" = DEG %>% filter(cluster == 4)
  ) %>% 
  map_dfr(
    ~enrich_genes(.$gene_short_name, databases = db),
    .id = "cluster"
    )

saveRDS(GSEA,
        file = paste0(path, "/nb_GSEA.rds")
)