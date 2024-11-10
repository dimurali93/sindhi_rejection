rm(list = ls())
# source library
source("/ix1/rsindhi/dim95/Rprojects/sindhi_rejection/1.library.R")

# set directory for raw files
SP_RD_Path = "/ix1/rsindhi/dim95/RawData/Visium_Transplant_c2/SpaceRanger_OUT/"

# read meta data
sampleMeta = read.csv(
  "/ix1/rsindhi/dim95/Rprojects/sindhi_rejection/metadata/manifest_c2_simplified.csv",
  fill = T)
SampleNames <- list.dirs(SP_RD_Path, recursive = F, full.names = F)

# read raw 10x spaceranger fiels
Read10x_spatial_mod = function(SP_RD_Path, numSamples) {
  # get entire file loaction for the files presesnt in the spaceranger output folder
  fileLocation = list.files(path = SP_RD_Path, full.names = TRUE)
  # get sample names from the folder where all raw data is present
  SampleNames <- list.dirs(SP_RD_Path, recursive = F, full.names = F)
  images <- listData <- list()
  for (i in 1:length(SampleNames)) {
    images[[i]] <-
      Read10X_Image(image.dir = paste0(fileLocation[[i]], "/outs/spatial"))
    listData[[i]] = Load10X_Spatial(
      data.dir = paste0(fileLocation[[i]], "/outs"),
      filename = "filtered_feature_bc_matrix.h5",
      assay = "Spatial",
      slice = SampleNames[[i]],
      image = images[[i]]
    )
    listData[[i]]$orig.ident <- SampleNames[[i]]
    # save inital seurat object
    # saveRDS(listData, "seu_list1.RDS")
  }
  names(listData) <- SampleNames
  return(listData)
}
listData = Read10x_spatial_mod(SP_RD_Path, numSamples)


# get spots per samples
spotsPerSamplelist = function(listdata) {
  data_df <- data.frame(Sample = character(length(listdata)),
                        SpotsWithData = numeric(length(listdata)))
  for (i in 1:length(listdata)) {
    seurat_obj <- listdata[[i]]
    counts_matrix <- seurat_obj@assays$Spatial@cells
    num_spots_with_data <- sum(rowSums(counts_matrix) > 0)
    data_df$Sample[i] <-
      paste("Sample", i)  # Replace with your sample names
    data_df$SpotsWithData[i] <- num_spots_with_data
  }
  data_df$Sample = names(listdata)
  spotsPerSample = ggplot(data_df, aes(x = Sample, y = SpotsWithData, fill = "grey")) +
    geom_bar(stat = "identity") + # Create a bar plot
    xlab("Sample") +
    ylab("Number of Spots per sample") +
    ggtitle("Number of Spots  per Sample") +
    theme_minimal()
  return(spotsPerSample)
}
spotsPerSample = spotsPerSamplelist(listData)


# merge all samples and remove outlier discussed earlier
merged_C2 =  merge(
  listData[[1]],
  y = c(
    listData[[2]],
    listData[[3]],
    listData[[4]],
    listData[[5]],
    listData[[6]],
    listData[[7]],
    listData[[8]],
    listData[[9]],
    listData[[10]],
    listData[[11]],
    listData[[12]],
    listData[[13]]
  ),
  add.cell.ids = SampleNames[1:13]
)

# add metadata
seurat_metadata <- merged_C2@meta.data
metadata_combined <-
  merge(seurat_metadata, sampleMeta, by = "orig.ident")
rownames(metadata_combined) <- rownames(merged_C2@meta.data)
merged_C2@meta.data <- metadata_combined

# adding mito and hb percent
merged_C2[["percent.mt"]] <- PercentageFeatureSet(merged_C2,
                                                  pattern = "^MT-*", assay = "Spatial")

merged_C2[["percent_hb"]] <- PercentageFeatureSet(merged_C2,
                                                  pattern = "^HB.-*", assay = "Spatial")

# filtering dead cells
merged_C2 <-
  merged_C2[, merged_C2$nFeature_Spatial > 300 &
              merged_C2$nCount_Spatial > 150 &
              merged_C2$percent.mt < 8 &
              merged_C2$percent_hb < 10]

# removing hb and albumin
merged_C2 <- merged_C2[!grepl("^ALB", rownames(merged_C2)),]
merged_C2 <- merged_C2[!grepl("^Hb.*-", rownames(merged_C2)),]

# normalization using SCT
merged_C2_sct <- SCTransform(
  merged_C2,
  assay = "Spatial",
  vst.flavor = "v2",
  vars.to.regress = "percent.mt",
  verbose = FALSE
) %>%
  RunPCA(npcs = 30,
         verbose = FALSE) %>%
  RunUMAP(reduction = "pca",
          dims = 1:30,
          verbose = FALSE) %>%
  FindNeighbors(reduction = "pca",
                dims = 1:30,
                verbose = FALSE) %>%
  FindClusters(resolution = 0.5,
               verbose = FALSE)

###########################################################################
# subset samples and renormalize
###########################################################################
# outliers from MDS plots
# 8       1039     Rejector    CHS21-4785 NO DSA >>> dont remove
# 6       1029     Rejector CHS21-1905_02 NO DSA >>> dont remove
# 12      1073 Non-Rejector    CHS22-1363 >>>>>>>>>>>>>>>>>>>>>>>>REMOVE
# 4       1011     Rejector CHS20-5330_02 matched removal>>> dont remove
# 1        977     Rejector CHS19-3245_02>>> dont remove
# 5       1029     Rejector CHS21-1905_01 NO DSA >>> dont remove

# samples removed from scRNA
# 991 - NA spatial
# 1011 - CHS20-5330_01
# 1015 - NA spatial
# 1025 -CHS21-5518_01,CHS21-5518_02
# 1027 - NA spatial
# 1031 - NA spatial
# 1046 - NA spatial

#
# # remove samples as outliers
Idents(merged_C2_sct) = "orig.ident"
merged_C2_subset = subset(
  merged_C2_sct ,
  idents = c(
    "CHS19-3245_02",
    "CHS20-5330_02",
    'CHS21-1905_02',
    "CHS21-2579",
    "CHS21-5518_01",
    "CHS22-1363"
  ),
  invert = TRUE
)

# Normalize after subset
merged_C2_subset_sct <- SCTransform(
  merged_C2_subset,
  assay = "Spatial",
  vst.flavor = "v2",
  verbose = FALSE) %>%
  RunPCA(npcs = 30,
         verbose = FALSE)

pct <-
  merged_C2_subset_sct[["pca"]]@stdev / sum(merged_C2_subset_sct[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <-
  sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)


merged_C2_subset_sct = RunUMAP(
  merged_C2_subset_sct,
  reduction = "pca",
  dims = 1:pcs,
  verbose = FALSE) %>%
  FindNeighbors(reduction = "pca",
                dims = 1:pcs,
                verbose = FALSE) %>%
  FindClusters(resolution = 0.5,
               verbose = FALSE)

###########################################################################
# Identify Markers
###########################################################################

SCTResults(object = merged_C2_subset_sct, slot = "umi.assay")
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[2]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[3]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[4]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[5]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[6]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[7]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[8]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[9]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[10]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[11]],
     name = "umi.assay") <- "Spatial"
slot(object = merged_C2_subset_sct@assays$SCT@SCTModel.list[[12]],
     name = "umi.assay") <- "Spatial"
SCTResults(object = merged_C2_subset_sct, slot = "umi.assay")

integrated_C2_subset_sct <- IntegrateLayers(
  object = merged_C2_subset_sct,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  verbose = F
)

pct <-
  merged_C2_subset_sct[["pca"]]@stdev / sum(merged_C2_subset_sct[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <-
  sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

integrated_C2_subset_sct <-
  FindNeighbors(integrated_C2_subset_sct, dims = 1:pcs)
integrated_C2_subset_sct <-
  FindClusters(integrated_C2_subset_sct, resolution = .5)
integrated_C2_subset_sct <-
  RunUMAP(integrated_C2_subset_sct, dims = 1:pcs)

DimPlot(integrated_C2_subset_sct,
        group.by = "seurat_clusters")
DotPlot(
  integrated_C2_subset_sct,
  features = c(
    "SOX9",
    "KRT7",
    "EPCAM",
    "BICC1",
    "DCDC2" ,
    "GPNMB1",
    "C1QB1",
    "SPP1",
    "GLUL",
    "CD74",
    "HLA-DRA",
    "LTB"
  )
) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red") +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5,
    hjust = 1))

Idents(integrated_C2_subset_sct) = "seurat_clusters"
integrated_C2_subset_sct_clusters <-
  Idents(integrated_C2_subset_sct)
integrated_C2_subset_sct$label_scRNA_integ <-
  integrated_C2_subset_sct$seurat_clusters

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <-
  c(
    "GLUL",
    "Midzone",
    "Midzone",
    "Midzone",
    "Cholangiocytes",
    "GLUL",
    "Midzone",
    "Midzone"
  )
integrated_C2_subset_sct$label_scRNA_integ <-
  plyr::mapvalues(x = integrated_C2_subset_sct$label_scRNA_integ,
                  from = current.cluster.ids,
                  to = new.cluster.ids)

umap_seurat = DimPlot(integrated_C2_subset_sct, group.by = "seurat_clusters")
umap_labels = DimPlot(integrated_C2_subset_sct, group.by = "label_scRNA_integ")
umap_outcome = DimPlot(integrated_C2_subset_sct, group.by = "Outcome")
umap_DSA = DimPlot(integrated_C2_subset_sct, group.by = "DSA")
umap_orig = DimPlot(integrated_C2_subset_sct, group.by = "orig.ident")
umap_OutcomeDSA = DimPlot(integrated_C2_subset_sct, group.by = "OutcomeDSA")
umap_orig | umap_seurat | umap_labels
umap_DSA | umap_outcome | umap_OutcomeDSA



Idents(integrated_C2_subset_sct) = "label_scRNA_integ"

SCTResults(object = integrated_C2_subset_sct, slot = "umi.assay")
integrated_C2_subset_sct <-
  JoinLayers(integrated_C2_subset_sct, assay = "Spatial")
integrated_C2_subset_sct <-
  PrepSCTFindMarkers(object = integrated_C2_subset_sct)

integrated_C2_subset_sct$SeuAnno = paste0(
  integrated_C2_subset_sct$seurat_clusters,
  "_",
  integrated_C2_subset_sct$label_scRNA_integ
)
integrated_C2_subset_sct$SeuAnno_Outcome = paste0(integrated_C2_subset_sct$SeuAnno,
                                                  "_",
                                                  integrated_C2_subset_sct$Outcome)
Idents(integrated_C2_subset_sct) = "SeuAnno"

visiumc2Markers = FindAllMarkers(integrated_C2_subset_sct, only.pos = TRUE)
visiumc2Markers = visiumc2Markers %>% filter(p_val_adj < 0.05)
# write.csv(visiumc2Markers, "Markersvisiumc2.csv")
Idents(integrated_C2_subset_sct) = "SeuAnno_Outcome"
control = c("Non-Rejector")
conditions = c("Rejector")
#
allcombs = c(unique(integrated_C2_subset_sct$SeuAnno_Outcome))
annotations = c(unique(integrated_C2_subset_sct$SeuAnno))
DEconditionsversusNormal = data.frame()
temp = data.frame()
for (i in annotations) {
  for (j in conditions) {
    ident1 = paste0(i, "_", j)
    idnet2 = paste0(i, "_", control)
    print(ident1)
    print(idnet2)
    if (idnet2 %in% allcombs) {
      temp = FindMarkers(
        integrated_C2_subset_sct,
        ident.1 = ident1,
        ident.2 = idnet2,
        # logfc.threshold = 0.236,
        # min.pct = 0.5,
        only.pos = F
      ) %>% rownames_to_column() %>% add_column(grp = paste0(ident1, "_", idnet2))
      DEconditionsversusNormal = rbind(DEconditionsversusNormal, temp)
    }
  }
}
DEconditionsversusNormal = DEconditionsversusNormal %>% filter(p_val_adj <
                                                                 0.05)
# write.csv(DEconditionsversusNormal,
#           "DEGs_VisiumCohort2_LFC0.5_min.pct0.01.csv")

commonwithCohot1c5 = DEconditionsversusNormal %>% filter(rowname %in% unique(c(c4_c5Genes$C5))) %>%
  filter(p_val < 0.05) %>% filter(abs(avg_log2FC) > .25)
# write.csv(commonwithCohot1c5, "commonwithCohot1c5.csv")

commonwithCohot1c4 = DEconditionsversusNormal %>% filter(rowname %in% unique(c(c4_c5Genes$C4))) %>%
  filter(p_val < 0.05) %>% filter(abs(avg_log2FC) > .25)
# write.csv(commonwithCohot1c4, "commonwithCohot1c4.csv")

gitcreds::gitcreds_set()
library("usethis")

library(gitcreds)
git config --global user.email "dimurali@ucsd.edu"
git config --global user.name "dimurali93"
