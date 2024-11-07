visiumc1 <- readRDS("~/Misc/sindhi/sindhiAll/FINALANALYSIS/Sindhi_Visium_FFPE_5R&2NR_C1/visiumc1.rds")
visiumc2 <- readRDS("~/Misc/sindhi/sindhiAll/FINALANALYSIS/visiumc2/integrated_C2_subset_sct (1).RDS")
DefaultAssay(visiumc1)="Spatial"
DefaultAssay(visiumc2)="Spatial"
visiumc1@assays$SCT<-NULL
visiumc2@assays$SCT<-NULL
visiumc1$batch ="visiumc1"
visiumc2$batch ="visiumc2"
library(Seurat)
DimPlot(visiumc1)
DimPlot(visiumc2)
visiumc1.list <- SplitObject(visiumc1, split.by = "orig.ident")
visiumc2.list <- SplitObject(visiumc2, split.by = "orig.ident")
visiumall <- c(visiumc1.list, visiumc2.list)
merged <- merge(visiumall[[1]],
                visiumall[2:length(visiumall)])
merged <- RenameCells(merged,
                      new.names =paste("cell", 1:ncol(merged),
                                       sep = "_"))
library(dplyr)
merged=  NormalizeData(merged)
merged =   FindVariableFeatures(merged,nfeatures = 3000)
merged = ScaleData(merged)
Idents(merged) <- "Outcome"
merged_rejector = subset(merged, idents ="Rejector")
Idents(merged_rejector) <- "orig.ident"
merged_rejector <- JoinLayers(merged_rejector)
markers <- FindAllMarkers(merged_rejector,
                          verbose = TRUE,
                          only.pos = TRUE)
markers <- dplyr::filter(markers,
                         p_val_adj < 0.01)

visiumall <- visiumall[names(visiumall) %in% unique(merged_rejector$orig.ident)]


visiumall_2 <- lapply(X = visiumall, FUN = function(data) {
  data <- data %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000)
})
features <- SelectIntegrationFeatures(visiumall_2,
                                      nfeatures = 2000)
batch_markers=markers
batch_markers <-batch_markers %>%
  group_by(cluster) %>%
  top_n(n = 150, wt = avg_log2FC)
table(batch_markers$gene %in% features)
features <- setdiff(features, batch_markers$gene)
print(length(features))
library(tidyverse)
visiumall_3 <- lapply(X = visiumall_2, FUN = function(data) {
  print(unique(data@meta.data$source))
  data <- data %>%
    ScaleData(features = features) %>%
    RunPCA(features = features)
})
# load selected integration features
anchors <- FindIntegrationAnchors(object.list = visiumall_3,
                                  normalization.method = "LogNormalize",
                                  anchor.features = features,
                                  reduction = "cca")
integrated <- IntegrateData(anchorset = anchors,k.weight = 70)
DefaultAssay(integrated) <- "integrated"
# # run the standard workflow on the integrated assay
# pct <- integrated[["pca"]]@stdev / sum(integrated[["pca"]]@stdev) * 100
# cumu <- cumsum(pct)
# co1 <- which(cumu > 90 & pct < 5)[1]
# co2 <-sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
#            decreasing = T)[1] + 1
# ndims <- min(co1, co2)
integrated <- integrated %>%
  ScaleData() %>%
  RunPCA(npcs = 10) %>%
  RunUMAP(reduction = "pca", dims = 1:10) %>%
  FindNeighbors(reduction = "pca", dims = 1:10)

Idents(integrated)="seurat_clusters"

DimPlot(integrated,split.by = "seurat_clusters")
FeaturePlot(integrated, features = c("GLUL","ASS1"))
SpatialDimPlot(integrated,crop = FALSE,ncol = 3)
DefaultAssay(integrated)<-"Spatial"
library(ggplot2)
library(ggpubr)
Idents(integrated)="seurat_clusters"
Idents(integrated)="label"
DotPlot(integrated, features =c( "GLUL", "ASS1","KRT7","SOX9","SCTR","AQP1","CFTR"))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")
SpatialDimPlot(integrated,crop = FALSE,ncol = 3)
integrated$seurat_clusters <- as.factor(as.numeric(as.character(integrated$seurat_clusters)))
Idents(integrated) <- "seurat_clusters"
levels(integrated)
new_ids <- c("GLUL",
             "mid",
             "mid",
             "mid",
             "Cholangiocyte",
             "GLUL",
             "mid",
             "mid")
names(new_ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new_ids)
integrated$label <- Idents(integrated)
DimPlot(integrated,label = T)
Idents(integrated)="DSA"
integrated = JoinLayers(integrated)
integrated$dsa_label = paste0(integrated$DSA,"_",integrated$label)

integrated$seu_label = paste0(integrated$seurat_clusters,"_",integrated$label)
Idents(integrated)="seu_label"
allmarkers =FindAllMarkers(integrated,logfc.threshold = .25)%>%filter(p_val<0.05)
write.csv(allmarkers,"./Misc/sindhi/AllMarkers_Combinedvisiumc1c2_LFC.25_pval0.05.csv")


Idents(integrated)="label"
allmarkers =FindAllMarkers(integrated,logfc.threshold = .25)%>%filter(p_val<0.05)
write.csv(allmarkers,"./Misc/sindhi/AllMarkers_AnnotationCombinedvisiumc1c2_LFC.25_pval0.05.csv")





Idents(integrated)="dsa_label"
DEG_Glul=FindMarkers(integrated,ident.1 = "yes_GLUL",ident.2 = "no_GLUL")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="GLUL_RDSAversusNoDSA")
DEG_Cholangiocyte=FindMarkers(integrated,ident.1 = "yes_Cholangiocyte",ident.2 = "no_Cholangiocyte")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="cholangiocyte_RDSAversusNoDSA")
DEG_mid=FindMarkers(integrated,ident.1 = "yes_mid",ident.2 = "no_mid")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="mid_RDSAversusNoDSA")
DEG_all=rbind(DEG_Glul,rbind(DEG_Cholangiocyte,DEG_mid))


GCMarker = c( "CD38", "CD83","CXCR4","CXCR5","AICDA","IRF8","CD40","IGHG1")
TFH = c("BCL6", "CXCR5", "ICOS", "IL21", "SH2D1A", "IRF4","PDCD1")
library(dplyr)
humoralgenes<- read_csv("./Misc/sindhi/GO_term_summary_20241031_140920.csv")
humoralgenes=data.frame(humoralgenes)
humoralgenes = humoralgenes[,"...12"]


humoralgenes=(unlist(humoralgenes))

RestrictedList <- DEG_all %>%
  mutate(group = case_when(
    rowname %in% humoralgenes ~ "HumoralImmunity"))

RestrictedList <- RestrictedList %>%
  mutate(group2 = case_when(
    rowname %in% TFH ~ "TFH"))

RestrictedList <- RestrictedList %>%
  mutate(group3 = case_when(
    rowname %in% GCMarker ~ "GC"))

RestrictedList <- RestrictedList %>%
  filter(!if_all(all_of(c("group3","group2","group")), is.na))

colnames(RestrictedList)=c("Gene","p_val","avg_log2FC","pct.1", "pct.2","p_val_adj","cluster", "HumoralImmunity","TFH","GC")

# write.csv(RestrictedList,"./Misc/sindhi/DEG_all_visumc1c2_rdsaVsNoDSA_restrictedList.csv")
write.csv(RestrictedList,"./Misc/sindhi/DEG_all_visumc1c2_rdsaVsNoDSA_restrictedList_LFC0.25_pval0.05.csv")

# write.csv(DEG_all,"./Misc/sindhi/DEG_all_visumc1c2_rdsaVsNoDSA.csv")
write.csv(DEG_all,"./Misc/sindhi/DEG_all_visumc1c2_rdsaVsNoDSA_LFC0.25_pval0.05.csv")


Idents(integrated)="label"
DotPlot(integrated, features =c( "GLUL", "ASS1","KRT7","SOX9","SCTR","AQP1","CFTR"))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")
Idents(integrated)="dsa_label"
DotPlot(integrated, features =c( "CD38", "CD83","CXCR4","CXCR5","AICDA","IRF8","CD40","IGHG1"))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

DoHeatmap(integrated, features = c("CD38", "CD83","CXCR4","CXCR5","AICDA","IRF8","CD40","IGHG1"),assay = "integrated")

GCMarker = c( "CD38", "CD83","CXCR4","CXCR5","AICDA","IRF8","CD40","IGHG1")
TFH = c("BCL6", "CXCR5", "ICOS", "IL21", "SH2D1A", "IRF4","PDCD1")
library(dplyr)
humoralgenes<- read_csv("Misc/sindhi/GO_term_summary_20241031_140920.csv")
humoralgenes=data.frame(humoralgenes)
humoralgenes = humoralgenes%>%
  # filter(Annotated.Term =="humoral immune response mediated by circulating immunoglobulin")%>%
  select("...12")



Obj <- ScaleData(integrated, features = rownames(integrated))




DotPlot(integrated, features =unique(humoralgenes))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")
DoHeatmap(Obj, features = unlist(unique(humoralgenes)))
DoHeatmap(Obj, features = GCMarker)
DoHeatmap(Obj, features = TFH)





##################################
DEG_all_visumc1c2_rdsaVsNoDSA_LFC0_25_pval0_05 <- read.csv("Misc/sindhi/DEG_all_visumc1c2_rdsaVsNoDSA_LFC0.25_pval0.05.csv")
DEG_all = DEG_all_visumc1c2_rdsaVsNoDSA_LFC0_25_pval0_05
GCMarker = c( "CD38", "CD83","CXCR4","CXCR5","AICDA","IRF8","CD40","IGHG1")
TFH = c("BCL6", "CXCR5", "ICOS", "IL21", "SH2D1A", "IRF4","PDCD1")
library(dplyr)
humoralgenes<- read_csv("./Misc/sindhi/GO_term_summary_20241031_140920.csv")
humoralgenes=data.frame(humoralgenes)
humoralgenes = humoralgenes[,"...12"]


humoralgenes=(unlist(humoralgenes))

TF_TFs = as.set(c("CD4", "IL21", "PDCD1", "CXCR3", "CCR6", "ICOS", "IRF4", "BCL6", "CXCR5", "CXCL13", "IL4","ASCL2", "TCF", "MAF", "BATF"))
GC_TFs = as.set(c("BCL6", "FOXO1", "FOXP1" ))
GCZones= as.set(c("CD83" , "CXCR4"))
cd40Signaling = as.set(c("CD40", "TRAF1", "ICAM1"))
nfkbactivation = as.set(c("NFKB2", "RELB", "NFKB1", "REL"))
ProliferationdarkZoneBcells =as.set(c("PCNA", "MKI67", "CDK1", "CDC20", "MYC"))
memoryBcells = as.set(c("PRDM1", "CCR6", "CD27", "XBP1", "IRF4", "MZB1", "TNFRSF17"))
otherMarkers = as.set(c("AICDA" , "BCL2A1", "EZH2", "E2F1"))
calciumsignaling =as.set(c("PTPN6", "CD22"))
BCRgenes=as.set(c("BTK", "BLK", "BLNK"))


RestrictedListvisiumc1c2 <-mutate(DEG_all,TFH = case_when(
  rowname %in% TFH ~ "TFH"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,GC = case_when(
  rowname %in% GCMarker ~ "GC"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,HumoralImmunity = case_when(
  rowname %in% humoralgenes ~ "HumoralImmunity"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,GC_TFs = case_when(
  rowname %in% GC_TFs ~ "GC_TFs"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,TF_TFs = case_when(
  rowname %in% TF_TFs ~ "TF_TFs"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,GCZones = case_when(
  rowname %in% GCZones ~ "GCZones"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,cd40Signaling = case_when(
  rowname %in% cd40Signaling ~ "cd40Signaling"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,nfkbactivation = case_when(
  rowname %in% nfkbactivation ~ "nfkbactivation"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,ProliferationdarkZoneBcells = case_when(
  rowname %in% ProliferationdarkZoneBcells ~ "ProliferationdarkZoneBcells"))


RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,otherMarkers= case_when(
  rowname %in% otherMarkers ~ "otherMarkers"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,calciumsignaling= case_when(
  rowname %in% calciumsignaling ~ "calciumsignaling"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,BCRgenes= case_when(
  rowname %in% BCRgenes ~ "BCRgenes"))

RestrictedListvisiumc1c2 <-  mutate(RestrictedListvisiumc1c2,memoryBcells = case_when(
  rowname %in% memoryBcells ~ "memoryBcells"))


RestrictedListvisiumc1c2 <-   filter(RestrictedListvisiumc1c2,!if_all(c("memoryBcells","BCRgenes","calciumsignaling","otherMarkers","ProliferationdarkZoneBcells","cd40Signaling","GCZones","TF_TFs","GC_TFs","HumoralImmunity","GC","TFH","nfkbactivation"), is.na))

write.csv(RestrictedListvisiumc1c2,"~/Misc/sindhi/DEGs_RestrictedListvisiumc1c2_expanded.csv")
