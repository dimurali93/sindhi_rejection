integrated_C2_subset_sct <- readRDS("~/Misc/sindhi/sindhiAll/FINALANALYSIS/visiumc2/integrated_C2_subset_sct (1).RDS")

Idents(integrated_C2_subset_sct)="Outcome"
integrated_C2_subset_sct_rejector = subset(integrated_C2_subset_sct, idents = c( "Rejector"))
remove(integrated_C2_subset_sct)

integrated_C2_subset_sct_rejector <- SCTransform(integrated_C2_subset_sct_rejector, vars.to.regress = "percent.mt", verbose = FALSE, assay = "Spatial")%>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:10, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:10, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

integrated_C2_subset_sct_rejector = PrepSCTFindMarkers(integrated_C2_subset_sct_rejector)
Idents(integrated_C2_subset_sct_rejector) = "label_scRNA_integ"

visiumc2_marker = FindAllMarkers(integrated_C2_subset_sct_rejector,logfc.threshold = .25 )%>%dplyr::filter(p_val<0.05)

integrated_C2_subset_sct_rejector$anno_dsa = paste0(integrated_C2_subset_sct_rejector$label_scRNA_integ,"_",integrated_C2_subset_sct_rejector$DSA)

Idents(integrated_C2_subset_sct_rejector)="anno_dsa"

DEGs_Visiumc2_RejectionGlul = FindMarkers(integrated_C2_subset_sct_rejector, ident.1 ="GLUL_yes", ident.2 = "GLUL_no" )%>%
  dplyr::filter(abs(avg_log2FC)>0.25)%>%
  filter(p_val<0.05)%>%rownames_to_column()%>%
  mutate(condition = "Glul_RejectorDSA")

DEGs_Visiumc2_RejectionMidzone = FindMarkers(integrated_C2_subset_sct_rejector, ident.1 ="Midzone_yes", ident.2 = "Midzone_no" )%>%
  # dplyr::filter(abs(avg_log2FC)>0.25)%>%
  filter(p_val<0.05)%>%rownames_to_column()%>%
  mutate(condition = "Midzone_RejectorDSA")

DEGs_Visiumc2_RejectionCholangiocytes = FindMarkers(integrated_C2_subset_sct_rejector, ident.1 ="Cholangiocytes_yes", ident.2 = "Cholangiocytes_no" )%>%
  dplyr::filter(abs(avg_log2FC)>0.25)%>%
  filter(p_val<0.05)%>%rownames_to_column()%>%
  mutate(condition = "Cholangiocytes_RejectorDSA")


DEGs_Visiumc2_Rejection =  rbind(DEGs_Visiumc2_RejectionGlul,rbind(DEGs_Visiumc2_RejectionMidzone,DEGs_Visiumc2_RejectionCholangiocytes))


Restricted_DEGs_Visiumc2_Rejection <-mutate(DEGs_Visiumc2_Rejection,TFH = case_when(
  rowname %in% TFH ~ "TFH"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,GC = case_when(
  rowname %in% GCMarker ~ "GC"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,HumoralImmunity = case_when(
  rowname %in% humoralgenes ~ "HumoralImmunity"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,GC_TFs = case_when(
  rowname %in% GC_TFs ~ "GC_TFs"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,TF_TFs = case_when(
  rowname %in% TF_TFs ~ "TF_TFs"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,GCZones = case_when(
  rowname %in% GCZones ~ "GCZones"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,cd40Signaling = case_when(
  rowname %in% cd40Signaling ~ "cd40Signaling"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,nfkbactivation = case_when(
  rowname %in% nfkbactivation ~ "nfkbactivation"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,ProliferationdarkZoneBcells = case_when(
  rowname %in% ProliferationdarkZoneBcells ~ "ProliferationdarkZoneBcells"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,otherMarkers= case_when(
  rowname %in% otherMarkers ~ "otherMarkers"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,calciumsignaling= case_when(
  rowname %in% calciumsignaling ~ "calciumsignaling"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,BCRgenes= case_when(
  rowname %in% BCRgenes ~ "BCRgenes"))

Restricted_DEGs_Visiumc2_Rejection <-  mutate(Restricted_DEGs_Visiumc2_Rejection,memoryBcells = case_when(
  rowname %in% memoryBcells ~ "memoryBcells"))


Restricted_DEGs_Visiumc2_Rejection <-   filter(Restricted_DEGs_Visiumc2_Rejection,!if_all(c("memoryBcells","BCRgenes","calciumsignaling","otherMarkers","ProliferationdarkZoneBcells","cd40Signaling","GCZones","TF_TFs","GC_TFs","HumoralImmunity","GC","TFH","nfkbactivation"), is.na))

write.csv(RestrictedListvisiumc2,"/home/divya/Misc/sindhi/DEGs_RestrictedListvisiumc2_expanded.csv")



