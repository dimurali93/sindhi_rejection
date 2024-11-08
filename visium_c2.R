source("/ix1/rsindhi/dim95/Rprojects/sindhi_rejection/1.library.R")
SP_RD_Path= "/ix1/rsindhi/dim95/RawData/Visium_Transplant_c2/SpaceRanger_OUT/"

sampleMeta = read.csv("./manifest_c2_simplified.csv", fill =T)
SampleNames <-list.dirs(SP_RD_Path, recursive = F, full.names = F)


#Read all the files for N and NR conditions and add names to the list with their image names or filenames
Read10x_spatial_mod = function(SP_RD_Path,numSamples){
  # get entire file loaction for the files presesnt in the spaceranger output folder
  fileLocation = list.files(path=SP_RD_Path, full.names = TRUE)
  # get sample names from the folder where all raw data is present
  SampleNames <-list.dirs(SP_RD_Path, recursive = F, full.names = F)
  images <- listData <- list()
  for (i in 1:length(SampleNames)) {
    images[[i]] <- Read10X_Image(image.dir =paste0(fileLocation[[i]],"/outs/spatial"))
    listData[[i]] = Load10X_Spatial(data.dir= paste0(fileLocation[[i]],"/outs"), filename = "filtered_feature_bc_matrix.h5",
                                    assay = "Spatial", slice = SampleNames[[i]], image = images[[i]])
    listData[[i]]$orig.ident <- SampleNames[[i]]
    # save inital seurat object
    saveRDS(listData, "seu_list_visium_c2.RDS")
  }
  names(listData)<- SampleNames
  return(listData)
}
listData = Read10x_spatial_mod(SP_RD_Path,numSamples)



# Filter Dead cells
filterDeadCells <- function(object) {
  seu = object
  seu[["percent.mt"]] <- PercentageFeatureSet(seu,pattern = "^MT-*", assay = "Spatial")
  seu[["percent_hb"]] <- PercentageFeatureSet(seu, pattern ="^HB.-*", assay = "Spatial")
  seu <- seu[, seu$nFeature_Spatial > 500 &
               seu$nCount_Spatial > 200 &
               seu$percent.mt < 8 &
               seu$percent_hb < 10]
  seu <- seu[!grepl("^ALB", rownames(seu)), ]
  seu <- seu[!grepl("^HB.*-", rownames(seu)), ]
  return(seu)
}
listData_Filtered<- lapply(listData, filterDeadCells)

# process seurat and SCT
# processNewSeurat <- function(parent.object,  res = 0.5) {
#   seu = parent.object
#   seu <- SCTransform(seu,
#                      assay = "Spatial",
#                      vst.flavor = "v2",
#                      vars.to.regress = "percent.mt",
#                      verbose = TRUE) %>%
#     RunPCA(npcs = 20,verbose = TRUE)
#   pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100
#   cumu <- cumsum(pct)
#   co1 <- which(cumu > 90 & pct < 5)[1]
#   co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
#   ndims <- min(co1, co2)
#   print(ndims)
#   seu = RunUMAP(seu, reduction = "pca",
#                 dims = 1:ndims,
#                 verbose = TRUE) %>%
#     FindNeighbors(reduction = "pca",
#                   dims = 1:ndims,
#                   verbose = TRUE) %>%
#     FindClusters(resolution = 0.5,
#                  verbose = TRUE)
#   return(seu)
# }
listData_Filtered$`CHS19-3245_02`<-NULL
listData_Filtered$`CHS20-5330_02`<-NULL
listData_Filtered$`CHS21-1905_02`<-NULL
listData_Filtered$`CHS21-5518_01`<-NULL
listData_Filtered$`CHS22-1363`<-NULL

merged <- merge(listData_Filtered[[1]],
                listData_Filtered[2:length(listData_Filtered)])
merged <- RenameCells(merged,
                      new.names =paste("cell", 1:ncol(merged),
                                       sep = "_"))

DefaultAssay(merged)="Spatial"
merged =  NormalizeData(merged)
merged =   FindVariableFeatures(merged,nfeatures = 3000)
merged = ScaleData(merged)
Idents(merged) <- "orig.ident"
merged<- JoinLayers(merged)
markers <- FindAllMarkers(merged,
                          verbose = TRUE,
                          only.pos = TRUE)
markers <- dplyr::filter(markers,
                         p_val_adj < 0.01)



listData_Filtered_sct <- lapply(X = listData_Filtered_sct, FUN = function(data) {
  data <- data %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 5000)
})

features <- SelectIntegrationFeatures(listData_Filtered_sct,
                                      nfeatures = 5000)

batch_markers=markers
batch_markers <-batch_markers %>%
  group_by(cluster) %>%
  top_n(n = 150, wt = avg_log2FC)
table(batch_markers$gene %in% features)
features <- setdiff(features, batch_markers$gene)
print(length(features))
listData_Filtered_sct <- lapply(X = listData_Filtered_sct, FUN = function(data) {
  print(unique(data@meta.data$source))
  data <- data %>%
    ScaleData(features = features) %>%
    RunPCA(features = features)
})

# load selected integration features
anchors <- FindIntegrationAnchors(object.list = listData_Filtered_sct,
                                  normalization.method = "LogNormalize",
                                  anchor.features = features,
                                  reduction = "cca")
visiumc2_integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(visiumc2_integrated) <- "integrated"

visiumc2_integrated <- visiumc2_integrated %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:15) %>%
  FindNeighbors(reduction = "pca", dims = 1:15)

Idents(visiumc2_integrated)="seurat_clusters"
DefaultAssay(visiumc2_integrated)="Spatial"
DotPlot(visiumc2_integrated,
                          features = c("GLUL", "ASS1","KRT7","SOX9","SCTR","AQP1","CFTR","EPCAM"))+
               scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
               theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))



visiumc2_integrated$seurat_clusters <- as.factor(as.numeric(as.character(visiumc2_integrated$seurat_clusters)))
Idents(visiumc2_integrated) <- "seurat_clusters"
levels(visiumc2_integrated)
new_ids <- c("GLUL",
             "mid",
             "mid",
             "Cholangiocyte",
             "Cholangiocyte")
names(new_ids) <- levels(visiumc2_integrated)
visiumc2_integrated <- RenameIdents(visiumc2_integrated, new_ids)
visiumc2_integrated$label <- Idents(visiumc2_integrated)

DotPlot(visiumc2_integrated,
        features = c("GLUL", "ASS1","KRT7","SOX9","SCTR","AQP1","CFTR","EPCAM"))+
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))



# DimPlot(visiumc2_integrated,reduction = "umap",group.by = "seurat_clusters")
SpatialDimPlot(visiumc2_integrated,ncol = 3,crop = FALSE)


integrated_C2_subset_sct <- readRDS("~/Misc/sindhi/sindhiAll/FINALANALYSIS/visiumc2/integrated_C2_subset_sct (1).RDS")

Idents(integrated_C2_subset_sct)="Outcome"
integrated_C2_subset_sct_rejector = subset(integrated_C2_subset_sct, idents = c( "Rejector"))
remove(integrated_C2_subset_sct)

integrated_C2_subset_sct_rejector <- SCTransform(integrated_C2_subset_sct_rejector, vars.to.regress = "percent.mt", verbose = FALSE, assay = "Spatial")%>%
  RunPCA( verbose = FALSE) %>%
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



