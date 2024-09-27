# Integrate single-cell data sets of healthy samples from three studies

library(Seurat)
# Load data
Sys.time()
raw <- read.table("rawdata/01.UMI.txt", header=T, row.names = 1)
Sys.time()
dim(raw)
raw[1:5,1:5]

# The following is the symbol annotation of the gene matrix
library(org.Hs.eg.db)
library(clusterProfiler)
ids <- bitr(rownames(raw), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db )
dim(ids)
head(ids)

ids=ids[ids$SYMBOL != '',]
ids=ids[ids$ENSEMBL %in%  rownames(raw),]
dat=raw[ids$ENSEMBL,]
dim(dat)

#The following is the operation of taking the mean of the repeated symbols
ids$median=apply(dat,1,median) 
ids=ids[order(ids$SYMBOL,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$SYMBOL),]
dat=dat[ids$ENSEMBL,] 
rownames(dat)=ids$SYMBOL
head(dat[1:4,1:4])  
dim(dat)

metadata <- read.csv("rawdata/03.Cell.Barcodes.csv",  header=T)
dim(metadata)
head(metadata)

pbmc <- CreateSeuratObject(counts = dat, min.cells = 3, min.features = 200)
pbmc
row.names(metadata) <- metadata$Index
# add metadata to Seurat object
pbmc <- AddMetaData(object = pbmc, metadata = metadata)
head(pbmc@meta.data)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes)

# Quality control standards were set
minGene=500
maxGene=8000
maxUMI=22000
pctMT=20
pctHB=1

pbmc <- subset(pbmc, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                 nCount_RNA < maxUMI & percent.mt < pctMT & percent.HB < pctHB)
pbmc
median(pbmc$nFeature_RNA)


pbmc <- NormalizeData(pbmc)  
DefaultAssay(pbmc) <- "RNA"
pbmc <- CellCycleScoring(object = pbmc, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

colnames(pbmc@meta.data)
table(pbmc$Phase)

all.genes <- rownames(pbmc)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 3000) %>% ScaleData(features = all.genes,vars.to.regress = c("percent.mt","percent.rb"))
pbmc <- SCTransform(pbmc, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"))

pbmc <- RunPCA(pbmc, verbose = F)
ElbowPlot(pbmc, ndims = 50)
pc.num=1:30
pbmc <- pbmc %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
pbmc <- pbmc %>% FindNeighbors(dims = pc.num) %>% FindClusters(resolution=0.8) %>% FindClusters(resolution=0.5)

table(pbmc$SCT_snn_res.0.5)
dir.create("Data")
saveRDS(pbmc, file = "Data/subN.rds")


dir.create("CellType/MarkerPlots")
source("~/project/Resource/sc_function.R")
markers <- ReadMarker("~/project/Resource/hgPBMC_Marker1.txt")

#DefaultAssay(pbmc) <- "RNA"
for(i in seq_along(markers)){
  #i=1
  marker_i <- markers[[i]]
  marker_i <- CaseMatch(marker_i, rownames(pbmc))
  if(length(marker_i)>=1){
    tmp <- length(marker_i)
    tmp <- ceiling(tmp/3)
    p <- FeaturePlot(pbmc, reduction = "umap", features = marker_i, ncol = 3)
    ggsave(paste0("CellType/MarkerPlots/featureplot_", names(markers)[i], ".png"), p,
           width = 15, height = 4.5*tmp, limitsize = F)
  }
}

### ViolinPlot
Idents(pbmc) <- "seurat_clusters"
for(i in seq_along(markers)){
  #i=1
  marker_i <- markers[[i]]
  marker_i <- CaseMatch(marker_i, rownames(pbmc))
  w <- 0.5*length(levels(pbmc))
  h <- 1.5*length(marker_i)
  h <- ifelse(h<3, 3, h)
  if(length(marker_i)<2){
    try({
      p <- VlnPlot(pbmc, features = marker_i, pt.size = 0.01) + NoLegend()
      ggsave(paste0("CellType/MarkerPlots/vlnplot_", names(markers)[i], ".png"), p, width = w, height = h)
    }, silent = T)
  } else{
    p <- VlnPlot(pbmc, features = marker_i, stack = T, flip = T)
    ggsave(paste0("CellType/MarkerPlots/vlnplot_", names(markers)[i], ".png"), p, width = w, height = h)
  }
}


###### cell annotation
table(pbmc$SCT_snn_res.0.5)
pbmc$celltype.main <- recode(pbmc$seurat_clusters,
                             "0" = "CD8_T",
                             "1" = "NK",
                             "2" = "Monocytes",
                             "3" = "CD4_T",
                             "4" = "CD4_T",
                             "5" = "CD8_T",
                             "6" = "CD4_T",
                             "7" = "NK",
                             "8" = "B_cells",
                             "9" = "Monocytes",
                             "10" = "NKT",
                             "11" = "Monocytes",
                             "12" = "NK",
                             "13" = "B_cells",
                             "14" = "Platelets",
                             "15" = "DCs",
                             "16" = "Monocytes",
                             "17" = "B_cells")
pbmc$batch <- recode(pbmc$orig.ident,
                     "A1" = "Cohort2",
                     "A2" = "Cohort2",
                     "A3" = "Cohort2",
                     "HC1" = "Cohort1",
                     "HC2" = "Cohort1",
                     "HC3" = "Cohort1",
                     "CT1" = "Cohort3",
                     "CT2" = "Cohort3",
                     "CT3" = "Cohort3",
                     "CT4" = "Cohort3",
                     "CT5" = "Cohort3",
                     "SC4" = "Cohort3",
                     "SC5" = "Cohort3",
                     "SC6" = "Cohort3",
                     "SC7" = "Cohort3",
                     "SC1" = "Cohort3",
                     "SC2" = "Cohort3",
                     "SC3" = "Cohort3",
                     "P1" = "Cohort2",
                     "P2" = "Cohort2",
                     "P3" = "Cohort2",
                     "P4" = "Cohort2",
                     "P5" = "Cohort2",
                     "P6" = "Cohort2")
table(pbmc$orig.ident)
Idents(pbmc) <- "celltype.main"


rm(list=ls())
#subT.rds subB.rds subB.rds subM.rds subN.rds 
pbmc <- readRDS("Data/panage_pbmc.rds") 
pbmc
#saveRDS(pbmc, "subRPCA_1103.rds")
head(pbmc@meta.data)
table(pbmc$tissue)
table(pbmc$tissue,pbmc$main)
min(pbmc$nFeature_RNA)
max(pbmc$percent.mt)

pbmc1 <- CreateSeuratObject(pbmc@assays$RNA@counts, meta.data = pbmc@meta.data)

meta <- as.data.frame(pbmc@meta.data)
meta[1:5,1:5]; dim(meta)
write.csv(meta, "tables/Ncell_meta.csv")

VlnPlot(pbmc, features = c("ND1","IL1B"), pt.size = 0, group.by = "celltype")+NoLegend() #, stack = T
DimPlot(pbmc, reduction = "tsne", group.by = "celltype.main", label = T, label.size = 5)
DimPlot(pbmc, reduction = "umap", group.by = "batch",  label = T, label.size = 0) #split.by="Group",
p3 = DimPlot(pbmc, reduction = "umap", split.by = "Class", label = T, label.size = 5,group.by = "celltype.main")+NoLegend()
p3
pc = p1 + p2 + plot_layout(guides = "collect")
pc
ggsave("CellType/celltype_main_manual.png", pc, width = 12, height = 5)

Idents(pbmc) <- "celltype.main"
levels(pbmc)
range(pbmc@assays$RNA@counts)

markers <- c("CD14", "S100A8", "LYZ",
             'CD4',"IL7R",
             'CD8A',"GZMK","CD8B",
             'MS4A1','CD79A',
             'GNLY','NKG7',"KLRD1",
             'IL3RA',
             'FOXP3',
             'PPBP',
             'FCER1A','CD1C',"HLA-DQA1",
             'CD34')
markers <- CaseMatch(markers, rownames(pbmc1))
markers <- as.character(markers)


p <- FeaturePlot(pbmc1, reduction = "umap", features = markers, ncol = 4)
p
ggsave("CellType/Markers_featureplot_umap.png", p, width = 14, height = 15)
p <- FeaturePlot(pbmc, reduction = "tsne", features = markers, ncol = 4)
p
ggsave("CellType/Markers_featureplot_tsne.png", p, width = 14, height = 15)


DefaultAssay(pbmc) <- "RNA"
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc, features = union(VariableFeatures(pbmc), markers))
p <- DoHeatmap(pbmc, features = markers, size = 2.5)
p
ggsave("CellType/Markers_heatmap.png", width = 15, height = 8)


markers <- c("IGF2", "H1FX", "GUK1", "GLUL", "SPTB")
p <- DotPlot(pbmc, features = markers, group.by = "celltype1") + RotatedAxis()
p
ggsave("CellType/Markers_dotplot.pdf", width = 15, height = 9)


p <- VlnPlot(pbmc, features = markers, pt.size = 0, stack = T)+NoLegend()
p
ggsave("CellType/Markers_vlnplot.pdf", width = 15, height = 9)
table(pbmc$Patient)
pbmc@meta.data$celltype <- paste(pbmc@meta.data$Class, pbmc@meta.data$celltype.main, sep = "_")
table(pbmc$celltype)

library(Seurat)
library(SingleR)
library(tidyverse)
library(NMF)
library(patchwork)
library(dplyr)


head(pbmc@meta.data)
table(pbmc$celltype.main)
pbmc <- subset(pbmc, celltype.main %in% c("CD4_T", "CD8_T", "NKT"))
#sco.tmp <- subset(pbmc, seurat_clusters %in% c(1,2,4,8,10,12))
colnames(sco@meta.data)

metadata.sub <- c("orig.ident", "percent.mt", "percent.rb","percent.HB","S.Score","G2M.Score","Phase","Patient","Class",
                  "predicted.celltype.l2","celltype.main")#"doublet_scores","predicted_doublets","scrublet_pred_adj",
pbmc <- CreateSeuratObject(pbmc@assays$RNA@counts, meta.data = pbmc@meta.data[,metadata.sub])
pbmc
#saveRDS(sco, "sco0.rds")

pbmc <- NormalizeData(pbmc)
pbmc <- SCTransform(pbmc, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
pbmc <- RunPCA(pbmc)
ElbowPlot(pbmc, ndims = 50)
pc.num=1:30
pbmc <- pbmc %>% RunUMAP(dims=pc.num) %>% RunTSNE(dims=pc.num)
pbmc <- pbmc %>% FindNeighbors(dims = pc.num) %>% FindClusters(resolution=1) %>% FindClusters(resolution=0.5)
pbmc <- pbmc %>% FindClusters(resolution=0.8) 


## 查看Markers
markers <- c("CD3D","CD3E","CD3G","KLRD1","CD4","FOXP3","IL2RA","CD40LG","CD8A","CD8B","CCR7","SELL","LEF1","GZMK","IFNG","GZMB","NKG7","IL7R","KLRG1","CD28")
p1 <- DotPlot(sco, features = markers)
p1
p <- FeaturePlot(sco, features = markers, ncol = 3,reduction = "umap")
p
ggsave("CellType/markers_featureplot.pdf", p, width = 12, height = 20)
table(pbmc$SCT_snn_res.1)
table(Idents(pbmc)) 
sco <-subset(pbmc, idents = c('0',"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))
sco <- sco %>% RunUMAP(dims=pc.num) %>% RunTSNE(dims=pc.num)
sco <- sco %>% FindNeighbors(dims = pc.num) %>% FindClusters(resolution=1) %>% FindClusters(resolution=0.5)
table(Idents(sco))

DimPlot(pbmc, label = T,group.by = "Class", label.size = 4) + NoLegend()

p2 <- DimPlot(sco, group.by = "celltype.fine", label = T, repel = T, label.size = 4,reduction = "umap")
p2
#p3 <- DimPlot(sco, label = T,group.by = "SCT_snn_res.0.5",reduction = "umap", label.size = 4) + NoLegend()
p <- p1|p2
p
table(sco$predicted.celltype.l2)
ggsave("CellType/Celltype_predicted.celltype.dot.pdf", p, width = 10, height = 5)

##  T cell annotation
sco$celltype.fine <- recode(sco$seurat_clusters,
                            '0' = 'CD4TEM',
                            '1' = 'CD8CTL',
                            '2' = 'CD4Naive',
                            '3' = 'CD8TEM',
                            '4' = 'CD4Naive',
                            '5' = 'CD8TEM',
                            '6' = 'CD4TCM',
                            '7' = 'NKT',
                            '8' = 'CD4CTL',
                            '9' = 'CD8CTL',
                            '10' = 'CD8CTL',
                            '11' = 'Treg',
                            '12' = 'gdT',
                            '13' = 'CD4TEM',
                            '14' = 'CD8TEM',
                            '15' = 'CD4CTL',
                            '16' = 'CD4CTL',
                            '17' = 'NKT',
                            '18' = 'gdT',
                            '19' = 'CD8TEM',
                            '20' = 'CD4Naive',
                            '21' = 'CD8Naive',
                            '22' = 'NKT')
saveRDS(pbmc, "Data/subB.rds")

setwd("/data1/project/super_CD4/")
rm(list=ls())
pbmc <- readRDS("Data/subT.rds")

Idents(pbmc) <- pbmc$celltype.fine
table(Idents(pbmc))
pbmc <- subset(pbmc, idents = c("CD4TEM", 'CD4TCM', "CD4Naive", 'CD4CTL', "Treg"))

DimPlot(pbmc, group.by = "celltype", label = T) + NoLegend()
DimPlot(pbmc, group.by = "celltype.fine", label = T) + NoLegend()

ggsave("CellType/Celltype_Cluster_fine.pdf", width = 12, height = 6)

sco@meta.data$celltype <- paste(sco@meta.data$Class, sco@meta.data$celltype.fine, sep = "_")
table(sco$celltype)
markers <- c("CD3D","CD3E","CD3G","CD4","FOXP3","IL2RA","CD40LG","CD8A","CD8B","CCR7","SELL","LEF1","TCF7","PRF1","GZMA","GZMB","GZMK","IFNG","IL7R","KLRG1","CD28")
p1 <- DotPlot(sco, features = markers,group.by = "celltype.fine")
p1
ggsave("CellType/markers_finedot.pdf", p1, width = 16, height = 6)

table(sco$celltype)
marker<-FindMarkers(sco, group.by="celltype",ident.1 = "SC_Treg", ident.2 = "CT_Treg", logfc.threshold = 0.1, min.pct = 0.1)
dim(marker)
write.xlsx(marker,file="SC-CT_Treg.xlsx",quote=F,row.names = T)

Idents(sco) <- "celltype.fine"
ClusterMarker <- FindAllMarkers(sco, assay = "SCT", slot = "data", only.pos = T)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.xlsx(ClusterMarker,'ClusterMarkers_fine.xlsx', quote=F,row.names=F)


ClusterMarker_filter <- ClusterMarker[!grepl("^RP[SL]", ClusterMarker$gene, ignore.case = F),]
ClusterMarker_filter <- ClusterMarker_filter[!grepl("^MT-", ClusterMarker_filter$gene, ignore.case = F),]
top = 30   
TopMarkers1 <- ClusterMarker_filter %>% filter(p_val_adj == 0) %>%
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_filter %>% filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarker_filter <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarker_filter,'Top30Marker_noRibo_noMito.csv', row.names=F)


TopMarkers2Lines(ClusterMarker_filter, output = "CellType")
table(sco$Class)
table(sco$celltype.fine)

sco_clusters <- FetchData(pbmc, vars = c("celltype.main", "Class"))
head(sco_clusters)
count_SC <- sco_clusters %>% filter(Class == "SC") %>% dplyr::count() %>% as.numeric()
count_CT <- sco_clusters %>% filter(Class == "CT") %>% dplyr::count() %>% as.numeric()

sco_counts <- sco_clusters %>% group_by(Class) %>% dplyr::count(celltype.main)

proportion_SC <- sco_counts %>% filter(Class == "SC") %>% mutate(proportion = n/count_SC)
proportion_CT <- sco_counts %>% filter(Class == "CT") %>% mutate(proportion = n/count_CT)

proportion_sco <- full_join(proportion_CT, proportion_SC, by = "celltype.main") %>%
  mutate(proportion.x = ifelse(is.na(proportion.x), 0,  proportion.x)) %>%
  mutate(proportion.y = ifelse(is.na(proportion.y), 0,  proportion.y)) %>%
  mutate(Class.x = "CT") %>%
  mutate(Class.y = "SC") %>%
  mutate(cluster_type = ifelse(proportion.x > proportion.y, "CT", "SC"))

cluster_type_data <- left_join(x = sco_clusters, y = proportion_sco, by = "celltype.main")
rownames(cluster_type_data) <- rownames(sco_clusters)

pbmc <- AddMetaData(pbmc, dplyr::select(cluster_type_data, cluster_type))
dim(cluster_type_data)
head(sco@meta.data)
### Bar plot for figure 2
n1 <-  dplyr::select(proportion_sco, c(Class.x, celltype.main, proportion.x)) %>%
  mutate(Class = Class.x) %>%
  mutate(proportion = proportion.x) %>%
  mutate(Class.x = NULL) %>%
  mutate(proportion.x = NULL)
t1 <-  dplyr::select(proportion_sco, c(Class.y, celltype.main, proportion.y)) %>%
  mutate(Class = Class.y) %>%
  mutate(proportion = proportion.y) %>%
  mutate(Class.y = NULL) %>%
  mutate(proportion.y = NULL)
proportion_sco2 <- rbind(n1, t1)
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")

ggplot(proportion_sco2, aes(fill = Class, y = proportion, x = celltype.main)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("#46ACC8","#F2300F"))

ggsave2("Tproportion.pdf",  width = 10, height = 6)
table(sco$Patient)

p1 <- DimPlot(sco, label = T,group.by = "Class", label.size = 4)
p1
p2 <- DimPlot(sco, group.by = "Patient", label = T, repel = T, label.size = 4,reduction = "umap")
p2
#p3 <- DimPlot(sco, label = T,group.by = "SCT_snn_res.0.5",reduction = "umap", label.size = 4) + NoLegend()
p <- p1|p2
p
ggsave("CellType/Celltype_class.celltype.pdf", p, width = 10, height = 5)
#
