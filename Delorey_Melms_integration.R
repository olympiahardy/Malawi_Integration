library(Seurat)
library(tidyverse)
library(harmony)
library("data.table")
library("Matrix")
library("Seurat")
delorey <- readRDS("/datastore/Olympia/Delorey_Lung_Atlas/lung.rds")
lung_all <- readRDS("/home1/2505621h/lung_atlas_cosmic.rds")
lung.cosmic <- readRDS("/home1/2505621h/lung_atlas_cosmic.rds")
melms <- Read10X("/datastore/Olympia/Melms_Lung_Atlas/")
melms.counts <- readMM("/datastore/Olympia/Melms_Lung_Atlas/matrix.mtx.gz")
melms.counts <- as.matrix(melms.counts)
rownames(melms.counts) <- melms.genes$V1
colnames(melms.counts) <- melms.cells$V1

melms.counts <- as(as.matrix(melms.counts), "sparseMatrix")
melms.meta <- read.csv("/datastore/Olympia/Melms_Lung_Atlas/lung_metaData.txt", header = T,
                       sep = "\t", row.names = 1)
melms.meta <- melms.meta[-1,]
melms <- CreateSeuratObject(melms.counts, meta.data = melms.meta)

saveRDS(melms, "/datastore/Olympia/Melms_Lung_Atlas/Melms.rds")
melms <- readRDS("/datastore/Olympia/Melms_Lung_Atlas/Melms.rds")

lung_covid <- subset(lung_all, subset = Group == "COVID-19")
delorey <- subset(delorey, subset = doublet == "False")
delorey[["percent.ribo"]] <- PercentageFeatureSet(delorey, pattern = "^RP[SL]")

melms_covid <- subset(melms, subset = disease__ontology_label == "COVID-19")
melms_covid[["percent.ribo"]] <- PercentageFeatureSet(melms_covid, pattern = "^RP[SL]")
melms_covid[["percent.mt"]] <- PercentageFeatureSet(melms_covid, pattern = "^MT-")
names(melms_covid@meta.data)[names(melms_covid@meta.data) == 'donor_id'] <- 'Sample'
names(melms_covid@meta.data)[names(melms_covid@meta.data) == 'disease__ontology_label'] <- 'Group'
melms_covid$Origin <- "Melms"

lung.list <- list(Malawi = lung_covid, Delorey = delorey, Melms = melms_covid )
lung.list$Malawi$Origin <- "Malawi"
lung.list$Delorey$Origin <- "Delorey"
names(lung.list$Delorey@meta.data)[names(lung.list$Delorey@meta.data) == 'donor'] <- 'Sample'
names(lung.list$Delorey@meta.data)[names(lung.list$Delorey@meta.data) == 'percent_mito'] <- 'percent.mt'
names(lung.list$Delorey@meta.data)[names(lung.list$Delorey@meta.data) == 'disease'] <- 'Group'

lung.list$Delorey <- SCTransform(lung.list$Delorey, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
lung.list$Malawi <- SCTransform(lung.list$Malawi, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
lung.list$Melms <- SCTransform(lung.list$Melms, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)



features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 3000)

lung.list$Malawi@assays$SCT@var.features <- features
lung.list$Delorey@assays$SCT@var.features <- features
lung.list$Melms@assays$SCT@var.features <- features

lung_all <- merge(lung.list$Malawi, y = c(lung.list$Delorey, lung.list$Melms), add.cell.ids = c("Malawi", "Delorey", "Melms"), project = "COVID_lung")

lung_all <- lung_all %>% 
  RunPCA(npcs = 40, assay = "SCT", features = features)
ElbowPlot(lung_all, ndims = 40)
lung_all <- RunUMAP(lung_all, reduction = "pca", dims = 1:35)
lung_all <- FindNeighbors(lung_all, reduction = "pca", dims = 1:35, verbose=FALSE)
lung_all <- FindClusters(lung_all, verbose = FALSE, res=0.5)
DimPlot(lung_all, pt.size = 1, label = T) + NoLegend()
DimPlot(lung_all, pt.size = 1, label = T, group.by = "Origin") + NoLegend()
DimPlot(lung_all, pt.size = 1, group.by = "Sample")

lung_all <- RunHarmony(lung_all, 
                       group.by.vars = c("Sample", "Origin"), 
                       reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
ElbowPlot(lung_all, ndims = 40)

lung_all <- RunUMAP(lung_all, reduction = "harmony", dims = 1:38)
lung_all <- FindNeighbors(lung_all, reduction = "harmony", dims = 1:38, verbose=FALSE)
lung_all <- FindClusters(lung_all, verbose = FALSE, res=0.5)

p1 <- DimPlot(lung_all, pt.size = 1, label = T, group.by = "seurat_clusters") + NoLegend()
p2 <- DimPlot(lung_all, pt.size = 1, group.by = "Origin")
p3 <- DimPlot(lung_all, pt.size = 1, group.by = "celltype_annotation")
p1 + p2 + p3

saveRDS(lung_all, "lung_malawi_delorey_melms.rds")
lung_all <- readRDS("lung_malawi_delorey_melms.rds")
library(future)
availableCores()
plan()
plan("multisession", workers = 20)
lung_all.markers <- FindAllMarkers(lung_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
top.lung_lung_all.markers <- lung_all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

# Removing junk and RBC - 15, 21, 22
lung_all <- subset(lung_all, idents = c(0:14, 16:20, 23:28))

# Macrophage and neutrophil subsets
lung_all <- SCTransform(lung_all, vars.to.regress = c("percent.mt", "percent.ribo"), 
                        verbose = T, variable.features.n = 3000, vst.flavor="v2")
lung_all <- lung_all %>% 
  RunPCA(npcs = 40, assay = "SCT", features = features)

lung.list <- SplitObject(lung_all, split.by = "Origin")
features <- SelectIntegrationFeatures(object.list = lung_all, nfeatures = 3000)

lung_all <- RunHarmony(lung_all, 
                                group.by.vars = c("Sample", "Origin"), 
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
ElbowPlot(lung_all, ndims = 40)

lung_all <- RunUMAP(lung_all, reduction = "harmony", dims = 1:25)
lung_all <- FindNeighbors(lung_all, reduction = "harmony", dims = 1:25, verbose=FALSE)
lung_all <- FindClusters(lung_all, verbose = FALSE, res=0.3)

DimPlot(lung_all, label = T, pt.size = 1)

#26,31 
lung_all <- subset(lung_all, idents = c(0:25, 27:30))
saveRDS(lung_all, file = "final_integrated_malawi_melms_delorey.rds")
DimPlot(lung_all, label = T, pt.size = 1, split.by = "Origin")

write_csv(lung_all.markers, "lung_atlas_markers_dmm.csv")
lung_all.markers <- read.csv("lung_atlas_markers_dmm.csv", header = T)

lung_all <- readRDS(file = "final_integrated_malawi_melms_delorey.rds")
library(pochi)
AbundancePlot(lung_all, group.by = "final.celltype", split.by = "Origin", replicate.by = "Sample")

# Macrophage and neutrophil subsets
lung_macs <- subset(lung_all, idents = c(0,4,6,21,22,20))

lung.macs.list <- SplitObject(lung_macs, split.by = "Origin")

lung.macs.list$Delorey <- SCTransform(lung.macs.list$Delorey, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
lung.macs.list$Malawi <- SCTransform(lung.macs.list$Malawi, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
lung.macs.list$Melms <- SCTransform(lung.macs.list$Melms, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)

features <- SelectIntegrationFeatures(object.list = lung.macs.list, nfeatures = 2000)

lung_macs <- SCTransform(lung_macs, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)

lung_macs <- lung_macs %>% 
  RunPCA(npcs = 40, assay = "SCT")

lung_macs <- RunHarmony(lung_macs, 
                       group.by.vars = c("Sample", "Origin"), 
                       reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
ElbowPlot(lung_macs, ndims = 40)

lung_macs <- RunUMAP(lung_macs, reduction = "harmony", dims = 1:25)
lung_macs <- FindNeighbors(lung_macs, reduction = "harmony", dims = 1:25, verbose=FALSE)
lung_macs <- FindClusters(lung_macs, verbose = FALSE, res=0.1)

DimPlot(lung_macs, label = T)
DimPlot(lung_macs, label = T, split.by = "Origin")

Idents(lung_macs) <- "seurat_clusters"
lung_macs.markers.clus <- FindAllMarkers(lung_macs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
top.lung_macs.markers.clus <- lung_macs.markers.clus %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

DoHeatmap(subset(lung_macs, downsample = 1500), features = top.lung_macs.markers.clus$gene)


lung_macs <- RenameIdents(lung_macs, c("0" = "Alveolar macrophages",
                                     "1" ="Monocyte-derived macrophages",
                                     "2" ="Neutrophils",
                                     "3" ="Alveolar macrophages",
                                     "4" ="Monocytes",
                                     "5" ="Alveolar macrophages",
                                     "6" ="Monocyte-derived macrophages",
                                     "7" ="Alveolar macrophages"))

VlnPlot(lung_macs, features = c("LYVE1", "SPP1", "MERTK", "APOE", "PPARG", "NAMPT", "FCGR3A", 
                                "FCGR3B", "CD14", "CD74", "CD163", "CD68"), assay = "RNA", pt.size = 0)
# MDM
VlnPlot(lung_macs, features = c("LYZ", "ACP5", "TYROBP", "LGALS1", "CD68", "AIF1", "CTSL", 
                                "EMP3", "FCER1G", "LAPTM5"), assay = "RNA", pt.size = 0)

# IntM
VlnPlot(lung_macs, features = c("LYZ", "CD163", "TYROBP", "LGALS1", "CD68", "CTSL", "EMP3", 
"AIF1", "FCER1G", "LAPTM5"), assay = "RNA", pt.size = 0)

# Alv.M
VlnPlot(lung_macs, features = c("HLA-DPB1", "HLA-DPA1", "CTSZ", "ACP5", "COTL1", "FCER1G", "C1QC", 
"LAPTM5", "CTSS", "HLA-DQA1", "APOE"), assay = "RNA", pt.size = 0)

lung_macs <- NormalizeData(lung_macs)
lung_macs <- ScaleData(lung_macs)

mono.macs_barcodes <- WhichCells(lung_macs, idents = "Monocyte-derived macrophages")
alv.macs_barcodes <- WhichCells(lung_macs, idents = "Alveolar macrophages")
neut_barcodes <- WhichCells(lung_macs, idents = "Neutrophils")
mono_barcodes <- WhichCells(lung_macs, idents = "Monocytes")


saveRDS(lung_macs, "lung_macs_mmd.rds")
lung_macs <- readRDS("lung_macs_mmd.rds")


lung_macs2 <- subset(lung_macs, idents = c("Alveolar macrophages", "Interstital macrophages"))

lung.macs.list2 <- SplitObject(lung_macs2, split.by = "Origin")

lung.macs.list2$Delorey <- SCTransform(lung.macs.list2$Delorey, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
lung.macs.list2$Malawi <- SCTransform(lung.macs.list2$Malawi, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
lung.macs.list2$Melms <- SCTransform(lung.macs.list2$Melms, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)

features <- SelectIntegrationFeatures(object.list = lung.macs.list2, nfeatures = 1500)

lung_macs2 <- SCTransform(lung_macs2, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)

lung_macs2 <- lung_macs2 %>% 
  RunPCA(npcs = 40, assay = "SCT", features = features)

lung_macs2 <- RunHarmony(lung_macs2, 
                        group.by.vars = c("Sample", "Origin"), 
                        reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
ElbowPlot(lung_macs2, ndims = 40)

lung_macs2 <- RunUMAP(lung_macs2, reduction = "harmony", dims = 1:25)
lung_macs2 <- FindNeighbors(lung_macs2, reduction = "harmony", dims = 1:25, verbose=FALSE)
lung_macs2 <- FindClusters(lung_macs2, verbose = FALSE, res=0.1)

DimPlot(lung_macs2, label = T)
DimPlot(lung_macs2, label = T, split.by = "Origin")

lung_macs2.markers <- FindAllMarkers(lung_macs2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
top.lung_macs2.markers <- lung_macs2.markers %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC)

# T-cells
lung_tcell <- subset(lung_all, idents = c(3,10))
lung_tcell.list <- SplitObject(lung_tcell, split.by = "Origin")

lung_tcell.list$Delorey <- SCTransform(lung_tcell.list$Delorey, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
lung_tcell.list$Malawi <- SCTransform(lung_tcell.list$Malawi, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
lung_tcell.list$Melms <- SCTransform(lung_tcell.list$Melms, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)

features <- SelectIntegrationFeatures(object.list = lung_tcell.list, nfeatures = 1500)

lung_tcell <- lung_tcell %>% 
  RunPCA(npcs = 40, assay = "SCT", features = features)

saveRDS(lung_tcell, "lung_tcells_mmd.rds")
lung_tcell <- readRDS("lung_tcells_mmd.rds")
lung_tcell <- RunHarmony(lung_tcell, 
                        group.by.vars = c("Sample", "Origin"), 
                        reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
ElbowPlot(lung_tcell, ndims = 40)

lung_tcell <- RunUMAP(lung_tcell, reduction = "harmony", dims = 1:30)
lung_tcell <- FindNeighbors(lung_tcell, reduction = "harmony", dims = 1:30, verbose=FALSE)
lung_tcell <- FindClusters(lung_tcell, verbose = FALSE, res=0.1)

DimPlot(lung_tcell, label = T)
DimPlot(lung_tcell, label = T, split.by = "Origin")

lung_tcell.markers <- FindAllMarkers(lung_tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
top.lung_tcell.markers <- lung_tcell.markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

lung_tcell <- RenameIdents(lung_tcell, c("0" = "CD4+ T cell",
                                       "1" ="NK cell",
                                       "2" ="CD4+ T cell",
                                       "3" ="CD8+ T cell",
                                       "4" ="CD8+ T cell",
                                       "5" ="CD8+ T cell",
                                       "6" ="CD4+ T cell",
                                       "7" ="NK cell",
                                       "8" ="CD4+ T cell",
                                       "9" ="CD4+ T cell",
                                       "10" ="CD4+ T cell",
                                       "11" ="CD4+ T cell"))

cd4_barcodes <- WhichCells(lung_tcell, idents = "CD4+ T cell")
cd8_barcodes <- WhichCells(lung_tcell, idents = "CD8+ T cell")
nk_barcodes <- WhichCells(lung_tcell, idents = "NK cell")

saveRDS(lung_tcell, "lung_tcells_mmd.rds")
lung_tcell <- readRDS("lung_tcells_mmd.rds")

# Overall annotation
lung_all <- RenameIdents(lung_all, c("0" = "Myeloid",
                                                       "1" ="Endothelium",
                                                       "2" ="Fibroblasts",
                                                       "3" ="T cells",
                                                       "4" ="Myeloid",
                                                       "5" ="AT2",
                                                       "6" ="Myeloid",
                                                       "7" ="AT1",
                                             "8" = "Plasma cells",
                                             "9" = "AT2",
                                             "10" = "T cells",
                                             "11" = "Fibroblasts",
                                             "12" = "Fibroblasts",
                                             "13" = "Fibroblasts",
                                             "14" = "AT1",
                                             "15" = "Ciliated cells",
                                             "16" = "Smooth muscle cells",
                                             "17" = "Secretory cells",
                                             "18" = "Lymphatic endothelium",
                                             "19" = "Myofibroblasts",
                                             "20" = "Myeloid",
                                             "21" = "Myeloid",
                                             "22" = "Myeloid",
                                             "23" = "B cells",
                                             "24" = "Mast cells",
                                             "25" = "Smooth muscle cells",
                                             "27" = "Fibroblasts",
                                             "28" = "Mesothelial",
                                     "29" = "Lipofibroblasts",
                                     "30" = "AT2"))

DimPlot(lung_all, label = T, pt.size = 1)
DimPlot(lung_all, split.by = "Origin", pt.size = 1)
lung_all$broad.celltype <- Idents(lung_all)

# # Just have a look at immune
# lung_immune <- subset(lung_all, idents = c("T cells", "Myeloid", "Plasma cells", "B cells"))
# 
# lung_immune.list <- SplitObject(lung_immune, split.by = "Origin")
# 
# lung_immune.list$Delorey <- SCTransform(lung_immune.list$Delorey, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
# lung_immune.list$Malawi <- SCTransform(lung_immune.list$Malawi, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
# lung_immune.list$Melms <- SCTransform(lung_immune.list$Melms, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = T)
# 
# features <- SelectIntegrationFeatures(object.list = lung_immune.list, nfeatures = 2000)
# 
# lung_immune <- lung_immune %>% 
#   RunPCA(npcs = 40, assay = "SCT", features = features)
# 
# saveRDS(lung_immune, "lung_immunes_mmd.rds")
# lung_immune <- readRDS("lung_immunes_mmd.rds")
# lung_immune <- RunHarmony(lung_immune, 
#                          group.by.vars = c("Sample", "Origin"), 
#                          reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
# ElbowPlot(lung_immune, ndims = 40)
# 
# lung_immune <- RunUMAP(lung_immune, reduction = "harmony", dims = 1:36)
# lung_immune <- FindNeighbors(lung_immune, reduction = "harmony", dims = 1:36, verbose=FALSE)
# lung_immune <- FindClusters(lung_immune, verbose = FALSE, res=0.1)
# 
# DimPlot(lung_immune, label = T)
# DimPlot(lung_immune, label = T, split.by = "Origin")
# 
# lung_immune.markers <- FindAllMarkers(lung_immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
# top.lung_immune.markers <- lung_immune.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 30, order_by = avg_log2FC)

write.csv(cosmic.markers, "cosmic_all_markers.csv")
cosmic.markers <- read.csv("cosmic_all_markers.csv", header = T, row.names = 1)
top.cosmic.markers <- cosmic.markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

lung_all.cells <- Cells(lung_all)
Idents(object = lung_all, cells = lung_all.cells[lung_all.cells %in% cd4_barcodes]) <- "CD4+ T cells"
Idents(object = lung_all, cells = lung_all.cells[lung_all.cells %in% cd8_barcodes]) <- "CD8+ T cells"
Idents(object = lung_all, cells = lung_all.cells[lung_all.cells %in% nk_barcodes]) <- "NK cells"

Idents(object = lung_all, cells = lung_all.cells[lung_all.cells %in% mono.macs_barcodes]) <- "Monocyte-derived macrophages"
Idents(object = lung_all, cells = lung_all.cells[lung_all.cells %in% alv.macs_barcodes]) <- "Alveolar macrophages"
Idents(object = lung_all, cells = lung_all.cells[lung_all.cells %in% neut_barcodes]) <- "Neutrophils"
Idents(object = lung_all, cells = lung_all.cells[lung_all.cells %in% mono_barcodes]) <- "Monocytes"

lung_all$final.celltype <- Idents(lung_all)


DimPlot(lung_all, group.by = "final.celltype", label = T, pt.size = 1)

library(pochi)
AbundancePlot(lung_all, group.by = "final.celltype", split.by = "Origin", replicate.by = "Sample")

delorey_UMAP_colours <- c("#aec7e8", "#2ca02c", "#1f77b4", "#9467bd", "#ff7f0e", "#e07b91", "#90A37F", "#d6bcc0", "#864E87", "#DBBA27", "#75AB6A",
                         "#B959BA", "#11AD82", "#bcbd22", "#D4674A", "#F7867C", "#8c564b",
                         "#f7b6d2", "#e377c2", "#CC998D", "#CCD6EB")
names(delorey_UMAP_colours) <- c("NK cells","CD4+ T cells","CD8+ T cells","Neutrophils", "Alveolar macrophages", 
                                "Endothelium","AT1","AT2","Fibroblasts", "Monocyte-derived macrophages",
                                "Myofibroblasts", "Ciliated cells", "Smooth muscle cells", "Secretory cells",
                                "Lymphatic endothelium", "Mast cells", "Plasma cells", "B cells", "Lipofibroblasts", 
                                "Mesothelial", "Monocytes")

DimPlot(lung_all, label = T, cols = delorey_UMAP_colours, pt.size = 1, label.size = 6) + NoLegend()

DimPlot(lung_all, cols = delorey_UMAP_colours, split.by = "Origin", pt.size = 1)

saveRDS(lung_all, file = "lung_malawi_delorey_melms_annotated.rds")

lung_all.markers <- FindAllMarkers(lung_all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
top.lung_all.markers <- lung_all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC)

DoHeatmap(subset(lung_all, downsample = 800), features = top.lung_all.markers$gene)

lung_all$hemisphere <- ifelse(lung_all$Origin == "Malawi", "Malawi", "Northern hemisphere")
## DE ##
########### Lung Malawi vs Lung Delorey ############
lung_all$celltype.origin <- paste(lung_all$final.celltype, lung_all$hemisphere, sep = "_")
# get the amount of clusters
num_cluster=length(unique(lung_all$final.celltype))
Idents(lung_all) <- lung_all$final.celltype

cluster.ids <- unique(Idents(lung_all))
Idents(lung_all) <- "celltype.origin"
DefaultAssay(lung_all) <- "RNA"

origin.covid.list <- list()
origin.covid.list.up <- list()
origin.covid.list.down <- list()
library(future)
availableCores()
plan()
plan("multisession", workers = 40)
for (x in 1:num_cluster) {
  n=cluster.ids[x]
  print(n)
  check_sanity = FALSE
  ## diff expression \
  COVID.response <- FindMarkers(lung_all,ident.1 = paste(n,"_Malawi", sep=""), assay = "RNA", slot = "data", ident.2 = paste(n,"_Northern hemisphere", sep=""), verbose = T, logfc.threshold = 0, min.pct = 0.1)
  COVID.response$gene <- rownames(COVID.response)
  
  COVID.response <- COVID.response[-c(grep("^RPS", COVID.response$gene),
                                      grep("^RPL", COVID.response$gene),
                                      grep("^RP", COVID.response$gene),
                                      grep("^HB", COVID.response$gene),
                                      grep("^MT-", COVID.response$gene),
                                      grep("MALAT1", COVID.response$gene),
                                      grep("^MTR", COVID.response$gene),
                                      grep("^RNA18S5", COVID.response$gene),
                                      grep("^RNA28S5", COVID.response$gene),
                                      grep("^RNA1\\dS5", COVID.response$gene)),]
  origin.covid.list[[n]] <- COVID.response
  COVID.response1 <- subset(COVID.response, COVID.response$avg_log2FC > 0.5 & COVID.response$p_val_adj < 0.05)
  COVID.response1 <- COVID.response1[order(COVID.response1$avg_log2FC, decreasing = T),]
  COVID.response2 <- subset(COVID.response, COVID.response$avg_log2FC < -0.5 & COVID.response$p_val_adj < 0.05)
  COVID.response2 <- COVID.response2[order(COVID.response2$avg_log2FC),]
  
  origin.covid.list.up[[n]] <- COVID.response1
  origin.covid.list.down[[n]] <- COVID.response2
}
names(origin.covid.list) <- cluster.ids
names(origin.covid.list.up) <- cluster.ids
names(origin.covid.list.down) <- cluster.ids




