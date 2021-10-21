source("analysis_tumor/packages.R") 
source("analysis_tumor/common_functions.R") 
source("styling.R") 



# Import data -------------------------------------------------------------

cds <- readRDS("analysis_tumor/data_generated/cds.rds")
metadata <- readRDS("analysis_tumor/data_generated/metadata_tumor.rds")
DEG <- readRDS("analysis_tumor/data_generated/deg.rds")

sig_jansky <- read_excel(
  "analysis_tumor/metadata/Jansky_SuppMat_Tables.xlsx",
  sheet = "Suppl. Table 5",
  skip = 1
)



# Enrichment of patients in clusters --------------------------------------

metadata <- data.table(as.data.frame(metadata))

pDT <- setNames(
  metadata[, c("sample", "group", "cluster"), with = FALSE], 
  c("Patient", "Group", "Cluster")
)

for (colx in colnames(pDT)) {
  pDT[[colx]] <- as.character(pDT[[colx]])
}

pDT <- pDT[ , .N, by = c("Patient", "Group", "Cluster")]

res <- data.table()
cx <- pDT$Cluster[1]
patx <- pDT$Patient[1]

# create contigency matrix
for (cx in unique(pDT$Cluster)) {
  for (patx in unique(pDT$Patient)) {
    cont.matrix <- matrix(
      nrow = 2,
      ncol = 2, 
      c(
        sum(pDT[Patient == patx & Cluster == cx]$N),
        sum(pDT[Patient != patx & Cluster == cx]$N),
        sum(pDT[Patient == patx & Cluster != cx]$N),
        sum(pDT[Patient != patx & Cluster != cx]$N)
        )
      )
    
    # Fisher's exact test
    fish <- fisher.test(cont.matrix)
    
    res <- rbind(
      res,
      data.table(
        Cluster = cx,
        Patient = patx,
        p = fish$p.value,
        OR = fish$estimate
      )
    )
  }
}

# adjust p-value with Benjamini-Hochberg method
res[, padj := p.adjust(p, method = "BH")]

res <- merge(
  res,
  unique(pDT[ , c("Patient", "Group"), with = FALSE]),
  by = "Patient"
)

patient_cluster_enrichment <-
  res %>% 
  transmute(
    Patient = rename_patients(Patient),
    Group = rename_groups(Group),
    Cluster = Cluster,
    Adjusted.P.value = padj,
    Odds.Ratio = OR
  ) %>%
  as_tibble()

saveRDS(
  patient_cluster_enrichment,
  "analysis_tumor/data_generated/enrichment_patient_cluster.rds"
)



# Cluster specific DEG enrichment in literature gene signature ------------

# select background
bg <- unique(DEG$gene_short_name)
          
enrDT <- as.data.table(DEG)

sig_jansky <- 
  sig_jansky %>% 
  filter(gene %in% bg) 

cx <- sig_jansky[1]
mx <- sig_jansky$ID[1]
enrRes <- data.table()

for (cx in unique(enrDT$cluster)) {
  for (mx in unique(sig_jansky$cluster)) {
    g.cx <- (enrDT %>% filter(cluster == cx))$gene_short_name
    g.mx <- (sig_jansky %>% filter(cluster == mx))$gene
    
    A <- length(intersect(g.cx, g.mx))
    B <- length(intersect(g.cx, bg)) - A
    C <- length(intersect(g.mx, bg)) - A
    D <- length(bg) - A - B - C
    
    stopifnot(A + B + C + D == length(bg))
    
    cont.matrix <- matrix(
      nrow = 2, 
      ncol = 2, 
      c(A, B, C, D)
    )
    fish <- fisher.test(cont.matrix)
    enrRes <- rbind(
      enrRes, 
      data.table(
        Cluster = cx, 
        MarkerList = mx, 
        p = fish$p.value, 
        OR = fish$estimate, 
        Genes = paste(intersect(g.cx, g.mx), collapse = ",")
      )
    )
  }
}

# adjust p-value by B.-H. method
enrRes[, padj := p.adjust(p, method = "BH")]

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

saveRDS(
  nb_DEG_Enrichment_SigJansky,
  "analysis_tumor/data_generated/enrichment_jansky.rds"
)



# GSEA on cluster specific DEG --------------------------------------------

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

saveRDS(GSEA, "analysis_tumor/data_generated/gsea.rds")
