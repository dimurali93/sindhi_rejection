# integrated_C2_subset_sct <- readRDS("~/Misc/sindhi/sindhiAll/FINALANALYSIS/visiumc2/integrated_C2_subset_sct (1).RDS")
source("./1.library.R")

visium_c2 <- readRDS("./visium_c2/integrated_C2_subset_sct.RDS")

Idents(visium_c2)="label_scRNA_integ"
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
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-mutate(DEGs_visium_c2_rejector_DSAvNoDSA,TFH = case_when(
#   rowname %in% TFH ~ "TFH"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,GC = case_when(
#   rowname %in% GCMarker ~ "GC"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,HumoralImmunity = case_when(
#   rowname %in% humoralgenes ~ "HumoralImmunity"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,GC_TFs = case_when(
#   rowname %in% GC_TFs ~ "GC_TFs"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,TF_TFs = case_when(
#   rowname %in% TF_TFs ~ "TF_TFs"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,GCZones = case_when(
#   rowname %in% GCZones ~ "GCZones"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,cd40Signaling = case_when(
#   rowname %in% cd40Signaling ~ "cd40Signaling"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,nfkbactivation = case_when(
#   rowname %in% nfkbactivation ~ "nfkbactivation"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,ProliferationdarkZoneBcells = case_when(
#   rowname %in% ProliferationdarkZoneBcells ~ "ProliferationdarkZoneBcells"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,otherMarkers= case_when(
#   rowname %in% otherMarkers ~ "otherMarkers"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,calciumsignaling= case_when(
#   rowname %in% calciumsignaling ~ "calciumsignaling"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,BCRgenes= case_when(
#   rowname %in% BCRgenes ~ "BCRgenes"))
#
# Restricted_DEGs_visium_c2_rejector_DSAvNoDSA <-  mutate(Restricted_DEGs_visium_c2_rejector_DSAvNoDSA,memoryBcells = case_when(
#   rowname %in% memoryBcells ~ "memoryBcells"))


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




