visium_c2 <- readRDS("~/Misc/sindhi/sindhiAll/FINALANALYSIS/visiumc2/integrated_C2_subset_sct (1).RDS")
source("./1.library.R")

visium_c2 <- readRDS("./visium_c2/integrated_C2_subset_sct.RDS")

Idents(visium_c2) ="label_scRNA_integ"
visiumc2_marker = FindAllMarkers(visium_c2,logfc.threshold = .25,only.pos = TRUE )%>%dplyr::filter(p_val<0.05)
write.csv(visiumc2_marker,"./DEGList/Markers_visium_c2_LFC0.25_pval0.05.csv")

Idents(visium_c2)="Outcome"
visium_c2_rejector = subset(visium_c2, idents = c( "Rejector"))

visium_c2_rejector <- SCTransform(visium_c2_rejector,
                                  vars.to.regress = "percent.mt", verbose = FALSE, assay = "Spatial")%>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:10, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:10, verbose = FALSE)
# %>% FindClusters(resolution = 0.5, verbose = FALSE)

visium_c2_rejector = PrepSCTFindMarkers(visium_c2_rejector)

Idents(visium_c2_rejector) = "label_scRNA_integ"
visiumc2_rejector_marker = FindAllMarkers(visium_c2_rejector,logfc.threshold = .25,only.pos = TRUE )%>%dplyr::filter(p_val<0.05)
write.csv(visiumc2_rejector_marker,"./DEGList/Markers_visium_c2_rejection_LFC0.25_pval0.05.csv")

visium_c2_rejector$anno_dsa = paste0(visium_c2_rejector$label_scRNA_integ,"_",visium_c2_rejector$DSA)
Idents(visium_c2_rejector)="anno_dsa"
DEGs_visium_c2_rejector_Glul_DSAvNoDSA = FindMarkers(visium_c2_rejector, ident.1 ="GLUL_yes", ident.2 = "GLUL_no" )%>%
  dplyr::filter(abs(avg_log2FC)>0.25)%>%
  filter(p_val<0.05)%>%rownames_to_column()%>%
  mutate(condition = "visium_c2_rejector_Glul_DSAvNoDSA")

DEGs_visium_c2_rejector_Midzone_DSAvNoDSA = FindMarkers(visium_c2_rejector, ident.1 ="Midzone_yes", ident.2 = "Midzone_no" )%>%
  # dplyr::filter(abs(avg_log2FC)>0.25)%>%
  filter(p_val<0.05)%>%rownames_to_column()%>%
  mutate(condition = "visium_c2_rejector_Midzone_DSAvNoDSA")

DEGs_visium_c2_rejector_Cholangiocytes_DSAvNoDSA = FindMarkers(visium_c2_rejector, ident.1 ="Cholangiocytes_yes", ident.2 = "Cholangiocytes_no" )%>%
  dplyr::filter(abs(avg_log2FC)>0.25)%>%
  filter(p_val<0.05)%>%rownames_to_column()%>%
  mutate(condition = "visium_c2_rejector_Cholangiocytes_DSAvNoDSA")

DEGs_visium_c2_rejector_DSAvNoDSA =  rbind(DEGs_visium_c2_rejector_Glul_DSAvNoDSA,rbind(DEGs_visium_c2_rejector_Midzone_DSAvNoDSA,DEGs_visium_c2_rejector_Cholangiocytes_DSAvNoDSA))

source("./GeneList.R")
Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <- DEGs_visium_c2_rejector_DSAvNoDSA %>%
  mutate(
    TFH = ifelse(rowname %in% TFH, "TFH", NA),
    GC = ifelse(rowname %in% GCMarker, "GC", NA),
    HumoralImmunity = ifelse(rowname %in% humoralgenes, "HumoralImmunity", NA),
    GC_TFs = ifelse(rowname %in% GC_TFs, "GC_TFs", NA),
    TF_TFs = ifelse(rowname %in% TF_TFs, "TF_TFs", NA),
    GCZones = ifelse(rowname %in% GCZones, "GCZones", NA),
    cd40Signaling = ifelse(rowname %in% cd40Signaling, "cd40Signaling", NA),
    nfkbactivation = ifelse(rowname %in% nfkbactivation, "nfkbactivation", NA),
    ProliferationdarkZoneBcells = ifelse(rowname %in% ProliferationdarkZoneBcells, "ProliferationdarkZoneBcells", NA),
    otherMarkers = ifelse(rowname %in% otherMarkers, "otherMarkers", NA),
    calciumsignaling = ifelse(rowname %in% calciumsignaling, "calciumsignaling", NA),
    BCRgenes = ifelse(rowname %in% BCRgenes, "BCRgenes", NA),
    memoryBcells = ifelse(rowname %in% memoryBcells, "memoryBcells", NA)
  )

# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-   filter(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,!if_all(c("memoryBcells","BCRgenes","calciumsignaling","otherMarkers","ProliferationdarkZoneBcells","cd40Signaling","GCZones","TF_TFs","GC_TFs","HumoralImmunity","GC","TFH","nfkbactivation"), is.na))
write.csv(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,"./DEGList/Restricted_DEGs_visium_c2_rejector_DSAvNoDSA.csv")


Idents(visium_c2_rejector)="SeuAnno"
Markers_seuAnno = FindAllMarkers(visium_c2_rejector, only.pos = TRUE,logfc.threshold= .25 )%>%dplyr::filter(p_val<0.05)
write.csv(Markers_seuAnno,"./DEGList/Markers_visium_c2_seuAnno_LFC0.25_pval0.05.csv")


# Create a new identifier based on seu_label and DSA
visium_c2_rejector$seu_label_dsa <- paste0(visium_c2_rejector$SeuAnno, "_", visium_c2_rejector$DSA)
Idents(visium_c2_rejector) <- "seu_label_dsa"

# Define the list of cell types for comparison
cell_types <-unique(visium_c2_rejector$SeuAnno)

# Initialize an empty list to store DEGs
DEG_list <- list()

# Loop through cell types and perform differential expression analysis
for (cell_type in cell_types) {
  deg <- FindMarkers(visium_c2_rejector,
                     ident.1 = paste0(cell_type, "_yes"),
                     ident.2 = paste0(cell_type, "_no")) %>%
    filter(p_val < 0.05, abs(avg_log2FC) > 0.25) %>%
    rownames_to_column() %>%
    add_column(cluster = paste0("visium_c2_rejector_", cell_type, "_DSAvNoDSA"))

  DEG_list[[cell_type]] <- deg
}

# Combine all DEGs into one dataframe
DEG_all_seu_label <- do.call(rbind, DEG_list)


source("./GeneList.R")
Restricted_DEGs <- DEG_all_seu_label %>%
  mutate(
    TFH = ifelse(rowname %in% TFH, "TFH", NA),
    GC = ifelse(rowname %in% GCMarker, "GC", NA),
    HumoralImmunity = ifelse(rowname %in% humoralgenes, "HumoralImmunity", NA),
    GC_TFs = ifelse(rowname %in% GC_TFs, "GC_TFs", NA),
    TF_TFs = ifelse(rowname %in% TF_TFs, "TF_TFs", NA),
    GCZones = ifelse(rowname %in% GCZones, "GCZones", NA),
    cd40Signaling = ifelse(rowname %in% cd40Signaling, "cd40Signaling", NA),
    nfkbactivation = ifelse(rowname %in% nfkbactivation, "nfkbactivation", NA),
    ProliferationdarkZoneBcells = ifelse(rowname %in% ProliferationdarkZoneBcells, "ProliferationdarkZoneBcells", NA),
    otherMarkers = ifelse(rowname %in% otherMarkers, "otherMarkers", NA),
    calciumsignaling = ifelse(rowname %in% calciumsignaling, "calciumsignaling", NA),
    BCRgenes = ifelse(rowname %in% BCRgenes, "BCRgenes", NA),
    memoryBcells = ifelse(rowname %in% memoryBcells, "memoryBcells", NA)
  )

# Save the results to a CSV file
write.csv(Restricted_DEGs, "./DEGList/DEGs_visium_c2_rejector_seu_label_DSAvNoDSA_LFC.25_pval.05.csv")
