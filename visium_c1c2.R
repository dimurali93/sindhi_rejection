source("./1.library.R", echo=TRUE)

visium_c1 <- readRDS("~/Misc/sindhi/sindhiAll/FINALANALYSIS/Sindhi_Visium_FFPE_5R&2NR_C1/visiumc1.rds")
visium_c2 <- readRDS("~/Misc/sindhi/sindhiAll/FINALANALYSIS/visiumc2/integrated_C2_subset_sct (1).RDS")

visium_c1 <- readRDS("/ix1/rsindhi/dim95/Rprojects/sindhi_rejection/visium_c1/visiumc1.rds")
visium_c2 <- readRDS("/ix1/rsindhi/dim95/Rprojects/sindhi_rejection/visium_c2/integrated_C2_subset_sct.RDS")


DefaultAssay(visium_c1)="Spatial"
DefaultAssay(visium_c2)="Spatial"

visium_c1@assays$SCT<-NULL
visium_c2@assays$SCT<-NULL

visium_c1$batch ="visium_c1"
visium_c2$batch ="visium_c2"

DimPlot(visium_c1)
DimPlot(visium_c2)

visium_c1.list <- SplitObject(visium_c1, split.by = "orig.ident")
visium_c2.list <- SplitObject(visium_c2, split.by = "orig.ident")

visiumc1c2_rejection <- c(visium_c1.list, visium_c2.list)

merged <- merge(visiumc1c2_rejection[[1]],
                visiumc1c2_rejection[2:length(visiumc1c2_rejection)])
merged <- RenameCells(merged,
                      new.names =paste("cell", 1:ncol(merged),
                                       sep = "_"))

merged=  NormalizeData(merged)
merged =   FindVariableFeatures(merged,nfeatures = 3000)
merged = ScaleData(merged)

merged@meta.data$Outcome <- gsub(
  pattern = "^Non-Rejector$",
  replacement = "NonRejector",
  merged@meta.data$Outcome
)
spe_cat <-as.SingleCellExperiment(merged)
spe_cat <- merged
spe_cat <- JoinLayers(spe_cat, assay = "Spatial")
spe_cat$sample_id <- factor(merged$orig.ident)
spe_cat$condition <- factor(merged$Outcome)
spe_cat$cluster_id <- factor(merged$batch)
spe_cat$DSA <- factor(merged$DSA)
spe_cat1=as.SingleCellExperiment(spe_cat)
# # Add celltype information to metadata
metadata(spe_cat1)$cluster_codes <- data.frame(celltype = factor(spe_cat1$orig.ident))
C2_MDS_1 = pbMDS(spe_cat1,
                 by = "sample_id",assay = "counts" ,
                 shape_by = "DSA", size_by = TRUE,
                 label_by = "sample_id",  dims = c(1,2),
                 k = "celltype") +
  scale_color_manual(values = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  theme(legend.position="right", aspect.ratio=1)
C2_MDS_1
# C2_MDS_1 = pbMDS(spe_cat1,
#                  by = "sample_id",assay = "counts" ,
#                  label_by = "sample_id",
#                  k = "celltype") +
#   scale_color_manual(values = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) +
#   theme(legend.position="right")
Idents(merged) <- "Outcome"
merged_rejector = subset(merged, idents ="Rejector")
Idents(merged_rejector) <- "orig.ident"
merged_rejector <- JoinLayers(merged_rejector)
markers <- FindAllMarkers(merged_rejector,
                          verbose = TRUE,
                          only.pos = TRUE)
markers <- dplyr::filter(markers,
                         p_val_adj < 0.01)


visiumc1c2_rejection <- visiumc1c2_rejection[names(visiumc1c2_rejection) %in% unique(merged_rejector$orig.ident)]


visiumc1c2_rejection_2 <- lapply(X = visiumc1c2_rejection, FUN = function(data) {
  data <- data %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000)
})
features <- SelectIntegrationFeatures(visiumc1c2_rejection_2,
                                      nfeatures = 2000)
batch_markers=markers
batch_markers <-batch_markers %>%
  group_by(cluster) %>%
  top_n(n = 150, wt = avg_log2FC)
table(batch_markers$gene %in% features)
features <- setdiff(features, batch_markers$gene)
print(length(features))

visiumc1c2_rejection_3 <- lapply(X = visiumc1c2_rejection_2, FUN = function(data) {
  print(unique(data@meta.data$source))
  data <- data %>%
    ScaleData(features = features) %>%
    RunPCA(features = features)
})
# load selected integration features
anchors <- FindIntegrationAnchors(object.list = visiumc1c2_rejection_3,
                                  normalization.method = "LogNormalize",
                                  anchor.features = features,
                                  reduction = "cca")

integrated <- IntegrateData(anchorset = anchors)
# DefaultAssay(integrated) <- "integrated"
integrated <- integrated %>%
  ScaleData() %>%
  RunPCA(npcs = 10) %>%
  RunUMAP(reduction = "pca", dims = 1:10) %>%
  FindNeighbors(reduction = "pca", dims = 1:10)%>%
  FindClusters(resolution = 0.5)

DotPlot(integrated, features =c( "GLUL", "ASS1","KRT7","SOX9","SCTR","AQP1","CFTR"),group.by = "seurat_clusters", assay ="Spatial")+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")


DimPlot(integrated,group.by  = "seurat_clusters")
FeaturePlot(integrated, features = c("GLUL","ASS1"))

Idents(integrated)="seurat_clusters"
SpatialDimPlot(integrated,crop = FALSE,ncol = 3)

RidgePlot(integrated, features = "GLUL", ncol = 2,group.by = "seurat_clusters")+ geom_vline(aes(xintercept = 2.5),colour = "red")
VlnPlot(integrated, features = "GLUL", ncol = 2,group.by = "seurat_clusters")+geom_boxplot()+ geom_hline(aes(yintercept = 2.5),colour = "red")

RidgePlot(integrated, features = c("KRT7","SOX9","SCTR","AQP1","CFTR"), ncol = 2,group.by = "seurat_clusters")
VlnPlot(integrated, features = c("KRT7","SOX9","SCTR","AQP1","CFTR"), ncol = 2,group.by = "seurat_clusters",  combine = TRUE)+geom_boxplot()

RidgePlot(integrated, features = "ASS1", ncol = 2,group.by = "seurat_clusters")
VlnPlot(integrated, features = "ASS1", ncol = 2,group.by = "seurat_clusters")+geom_boxplot()

integrated$seurat_clusters <- as.factor(as.numeric(as.character(integrated$seurat_clusters)))
Idents(integrated) <- "seurat_clusters"
levels(integrated)
new_ids <- c("IntraLobularRegion",
             "IntraLobularRegion",
             "IntraLobularRegion",
             "GLUL",
             "Cholangiocyte","Cholangiocyte","IntraLobularRegion",
             "Cholangiocyte")
names(new_ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new_ids)
integrated$label <- Idents(integrated)
DimPlot(integrated,label = T, group.by = "seurat_clusters",split.by = "label")

Idents(integrated)="DSA"
integrated = JoinLayers(integrated)


integrated$seu_label = paste0(integrated$seurat_clusters,"_",integrated$label)
Idents(integrated)="seu_label"

markers_visiumc1c2_rejection_seu_label =FindAllMarkers(integrated,logfc.threshold = .25, only.pos = TRUE)%>%filter(p_val<0.05)
write.csv(markers_visiumc1c2_rejection_seu_label,"./DEGList/Markers_Visiumc1c2_rejection_seu_label_LFC.25_pval0.05.csv")


integrated$seu_label_dsa = paste0(integrated$seu_label,"_",integrated$DSA)
Idents(integrated)="seu_label_dsa"
DEG_5_Glul=FindMarkers(integrated,ident.1 = "5_GLUL_yes",ident.2 = "5_GLUL_no")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_5_glul_DSAvNoDSA")

DEG_3_mid=FindMarkers(integrated,ident.1 = "3_mid_yes",ident.2 = "3_mid_no")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_3_mid_DSAvNoDSA")


DEG_2_mid=FindMarkers(integrated,ident.1 = "2_mid_yes",ident.2 = "2_mid_no")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_2_mid_DSAvNoDSA")

DEG_4_Cholangiocyte=FindMarkers(integrated,ident.1 = "4_Cholangiocyte_yes",ident.2 = "4_Cholangiocyte_no")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_4_Cholangiocyte_DSAvNoDSA")


DEG_0_GLUL=FindMarkers(integrated,ident.1 = "0_GLUL_yes",ident.2 = "0_GLUL_no")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_0_GLUL_DSAvNoDSA")

DEG_1_mid=FindMarkers(integrated,ident.1 = "1_mid_yes",ident.2 = "1_mid_no")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_1_mid_DSAvNoDSA")


DEG_6_mid=FindMarkers(integrated,ident.1 = "6_mid_yes",ident.2 = "6_mid_no")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_6_mid_DSAvNoDSA")

DEG_7_mid=FindMarkers(integrated,ident.1 = "7_mid_yes",ident.2 = "7_mid_no")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_7_mid_DSAvNoDSA")





DEG_all_seu_label <- rbind(DEG_0_GLUL, DEG_1_mid, DEG_2_mid, DEG_3_mid, DEG_4_Cholangiocyte,DEG_5_Glul, DEG_6_mid, DEG_7_mid)
source("./GeneList.R")
Restricted_DEGs_visium_c1c2_rejector_DSAvNoDSA_seu_label <- DEG_all_seu_label %>%
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



# write.csv(RestrictedList,"./Misc/sindhi/DEG_all_visumc1c2_rdsaVsNoDSA_restrictedList.csv")
write.csv(Restricted_DEGs_visium_c1c2_rejector_DSAvNoDSA_seu_label,"./DEGList/DEGs_visium_c1c2_rejector_seu_label_DSAvNoDSA_LFC.25_pval.05.csv")



Idents(integrated)="label"
markers_visiumc1c2_rejection_label =FindAllMarkers(integrated,logfc.threshold = .25, only.pos = TRUE)%>%filter(p_val<0.05)
write.csv(markers_visiumc1c2_rejection_label,"./DEGList/Markers_Visiumc1c2_rejection_label_LFC.25_pval0.05.csv")

integrated$dsa_label = paste0(integrated$DSA,"_",integrated$label)
Idents(integrated)="dsa_label"
DEG_Glul=FindMarkers(integrated,ident.1 = "yes_GLUL",ident.2 = "no_GLUL")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_glul_DSAvNoDSA")
DEG_Cholangiocyte=FindMarkers(integrated,ident.1 = "yes_Cholangiocyte",ident.2 = "no_Cholangiocyte")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_cholangiocyte_DSAvNoDSA")
DEG_mid=FindMarkers(integrated,ident.1 = "yes_mid",ident.2 = "no_mid")%>%
  # filter(p_val_adj<0.05)%>%
  filter(p_val<0.05)%>%
  filter(abs(avg_log2FC)>0.25)%>%rownames_to_column()%>%add_column(cluster ="visium_c1c2_rejector_mid_DSAvNoDSA")
DEG_all_label=rbind(DEG_Glul,rbind(DEG_Cholangiocyte,DEG_mid))

source("./GeneList.R")
Restricted_DEGs_visium_c1c2_rejector_DSAvNoDSA_label <- DEG_all_label %>%
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


# RestrictedList <- RestrictedList %>%
#   filter(!if_all(all_of(c("group3","group2","group")), is.na))

# colnames(RestrictedList)=c("Gene","p_val","avg_log2FC","pct.1", "pct.2","p_val_adj","cluster", "HumoralImmunity","TFH","GC")

# write.csv(RestrictedList,"./Misc/sindhi/DEG_all_visumc1c2_rdsaVsNoDSA_restrictedList.csv")
write.csv(Restricted_DEGs_visium_c1c2_rejector_DSAvNoDSA_label,"./DEGList/DEGs_visium_c1c2_rejector_lable_DSAvNoDSA_LFC.25_pval.05.csv")

# write.csv(DEG_all,"./Misc/sindhi/DEG_all_visumc1c2_rdsaVsNoDSA.csv")
# write.csv(DEG_all,"./Misc/sindhi/DEG_all_visumc1c2_rdsaVsNoDSA_LFC0.25_pval0.05.csv")


Idents(integrated)="label"
DotPlot(integrated, features =c( "GLUL", "ASS1","KRT7","SOX9","SCTR","AQP1","CFTR"))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")
Idents(integrated)="dsa_label"
DotPlot(integrated, features =c( "CD38", "CD83","CXCR4","CXCR5","AICDA","IRF8","CD40","IGHG1"))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

DoHeatmap(integrated, features = c("CD38", "CD83","CXCR4","CXCR5","AICDA","IRF8","CD40","IGHG1"),assay = "integrated")
