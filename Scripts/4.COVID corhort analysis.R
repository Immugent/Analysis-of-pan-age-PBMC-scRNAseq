#MISC-COVID2019
#Immunopathological signatures in multisystem inflammatory syndrome in children and pediatric COVID-19
rm(list=ls())
pbmc <- readRDS("Data/misc_final.rds")

head(before@meta.data)
table(pbmc$Class, pbmc$Age)
table(pbmc$age, pbmc$predicted.Group)
table(before$celltype2); table(after$celltype1)

DimPlot(pbmc, group.by =  "age",label = T, pt.size = 1,label.size = 5)
p1 <- DimPlot(pbmc, group.by =  "Class",label = T, pt.size = 1,label.size = 5) + NoLegend()
p2 <- DimPlot(pbmc, group.by = "mergedcelltype", pt.size = 1,label = T, label.size = 5)
p1|p2

metadata.sub <- c("orig.ident", "nCount_RNA", "nFeature_RNA","nCount_CITE","nFeature_CITE","nCount_HTO","Batch","Age",  "Gender","sample",
                  "Class", "mergedcelltype","predicted.id","subject_code","sample_code")
pbmc@meta.data = pbmc@meta.data[,metadata.sub]
head(pbmc@meta.data)
table(pbmc$mergedcelltype)

pbmc$Group <- recode(pbmc$Age,
                     "9" = "Tg",
                     "15" = "Tg",
                     "0.6" = "Cd",
                     "7" = "Tg",
                     "6" = "Cd",
                     "10.4" = "Tg",
                     "4" = "Cd",
                     "0.6" = "Cd",
                     "5" = "Tg",
                     "11" = "Tg",
                     "8" = "Tg",
                     "12" = "Tg",
                     "14" = "Tg",
                     "0.3" = "Cd",
                     "10" = "Tg",
                     "2" = "Cd")
table(pbmc$Group)
DimPlot(pbmc, reduction = "umap",group.by =  "Group",pt.size = 1)

pbmc$main <- recode(pbmc$mergedcelltype,
                    'B_Mem' = 'B cell',
                    'B_Naive' = 'B cell',
                    'CD4_isobinding' = 'T cell',
                    'CD4_Mem' = 'T cell',
                    'CD4_Naive' = 'T cell',
                    'CD8_Mem' = 'T cell',
                    'CD8_Naive' = 'T cell',
                    'cDC' = 'Myeloid',
                    'cKit+CD3-activated' = 'others',
                    'dim' = 'others',
                    'DNT' = 'T cell',
                    'HSC' = 'Hsc',
                    'MAIT' = 'T cell',
                    'Mono_Classical' = 'Myeloid',
                    'Mono_Intermediate' = 'Myeloid',
                    'Mono_NonClassical' = 'Myeloid',
                    'NK_CD16hi' = 'NK',
                    'NK_CD56hiCD16lo' = 'NK',
                    'NKT' = 'T cell',
                    'pDC' = 'Myeloid',
                    'Plasmablast' = 'B cell',
                    'Platelet' = 'Platelet',
                    'RBC' = 'others',
                    'T_Vd2' = 'T cell')
table(pbmc$main)
DimPlot(pbmc, group.by = "main", pt.size = 1,label = T, label.size = 5)
saveRDS(pbmc, "Data/misc_final.rds") 

subT <-subset(pbmc, subset = main =='T cell')
subN <-subset(pbmc, subset = main =='NK')
subM <-subset(pbmc, subset = main =='Myeloid')
subB <-subset(pbmc, subset = main =='B cell')

saveRDS(subB, "subB//misc_subB.rds")
saveRDS(subT, "subT//misc_subT.rds")
saveRDS(subM, "subM//misc_subM.rds")
saveRDS(subN, "subN//misc_subN.rds") 

rm(list=ls())
raw <- read.delim("GSE158055_covid19_counts.mtx.gz")
dim(raw)

#mapping
rm(list=ls())
setwd("/data1/project/super_CD4/")
reference <- readRDS("Data/panage_pbmc.rds")
head(reference@meta.data)
pbmc <- readRDS("mapping/chuanqi_predicted.rds")
head(pbmc@meta.data); pbmc
pbmc<- readRDS("mapping/chuanqi_merged_1112.rds")
after <- readRDS("mapping/chuanqi_after_filter.rds")
head(before@meta.data); before
head(after@meta.data)
before$group <- "before"
after$group <- "after"
table(before$orig.ident,before$age)
table(after$age, after$orig.ident)
pbmc <- merge(before, after)
saveRDS(pbmc, "mapping/chuanqi_merged_1112.rds")

VlnPlot(reference,features = c("percent.mt","percent.rb"),group.by = "Group",pt.size = 0)
DimPlot(reference, reduction = "umap", group.by = "celltype", label = TRUE,  pt.size = 0.5)

pbmc <- subset(pbmc, subset = Class %in% c("COVID","MIS-C"))
pbmc1 <- subset(pbmc, subset = Class %in% c("HC"))
head(pbmc@meta.data)
table(pbmc$celltype2)
table(pbmc$Donor, pbmc$Age)

[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# ??[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")
# ??enes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes)

####NormalizeData(pbmc)  
DefaultAssay(pbmc) <- "RNA"
pbmc <- CellCycleScoring(object = pbmc, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000) %>% ScaleData(vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"))
pbmc <- RunPCA(pbmc, verbose = F)
DimPlot(pbmc, split.by = "orig.ident", reduction = "pca", pt.size = 0.2, label = T, label.size = 5)

pc.num=1:30
pbmc <- RunUMAP(pbmc, dims = 1:30)
pbmc <- pbmc %>% RunUMAP(dims=pc.num) %>% RunTSNE(dims=pc.num) 
pbmc <- FindNeighbors(pbmc, dims = 1:30) %>% FindClusters(resolution=0.1)

DimPlot(pbmc, group.by = "celltype", pt.size = 0.2,label = T, label.size = 5)

head(pbmc@meta.data)
use_colors <- c(Cd ="deepskyblue2",Tg="steelblue", Ad="#0B775E")
p1 <- DimPlot(pbmc, group.by =  "predicted.Group",label = T, pt.size = 0.2,label.size = 0,cols = use_colors)
p1
use_colors <- c( Naive_B = "#FD6467", Memory_B = "darkgoldenrod2",Transitional_B = "steelblue", CD83_B = "#E2D200",
                 Plasma_cell = "#46ACC8", Activated_memoryB = "#E58601")
use_colors <- c( CD56b = "#9C964A", Inflamed_CD56d = "#F2300F",HLA_CD56d = "#5B1A18",Classical_CD56d = "#7294D4")
p2 <- DimPlot(pbmc, group.by = "predicted.celltype", pt.size = 0.2,label = T, label.size = 4,cols = use_colors)+NoLegend()
p2
use_colors <- c(COVID ="brown2", HC="#0B775E", MISC = "#B40F20")
p3 <- DimPlot(pbmc, group.by = "Class", pt.size = 0.2,label = T, label.size = 0,cols = use_colors)
p3
p1|p2|p3

anchors <- FindTransferAnchors(reference = reference, query = pbmc,
                               normalization.method = "LogNormalize",   
                               reduction = "pcaproject",
                               reference.reduction = "pca",  
                               dims = 1:50)


head(reference@meta.data)
pbmc <- MapQuery(anchorset = anchors, query = pbmc,  reference = reference,
                 refdata = list(celltype = "celltype1",                   Group = "Group" ),  reference.reduction = "pca")

head(pbmc@meta.data)
table(pbmc$Age, pbmc$predicted.celltype)
markers <- c('PTPRC',"MS4A1",'CD19','CD79A',"CD79B","CD3D","CD3E","CD3G",'NCR1',"KLRD1","KLRB1",'KLRF1',
             "CD14", "LYZ","VCAN","FCGR3A", "PPBP", "PF4","KIT","CD34")
markers <- c('PTPRC',"CD3D","CD3E","MS4A1",'CD79A',"KLRD1","KLRB1",
             "LYZ","VCAN", "PPBP", "PF4",'CD34')

FeaturePlot(pbmc, reduction = "umap", features = markers, ncol = 4,cols = c("#ccccca", "#B40F20"))
FeaturePlot(pbmc, reduction = "umap", features = c("CD34","KIT"), ncol = 2,cols = c("#ccccca", "#e61a2e"))

use_colors <- c(T_cell = "brown2",NK = "#B40F20", Myeloid = "steelblue",B_cell = "#E58601",Platelet = "deepskyblue2",
                Hpc = "seagreen", Cd = "deepskyblue2", Tg = "#0B775E", Ad = "#E2D200", Ed = "brown2", Sc = "#B40F20")
use_colors <- c(Shift = "#B40F20", Preserve = "steelblue")
DimPlot(pbmc, group.by = "celltype",label = T, label.size = 0,cols = use_colors)
DimPlot(pbmc, group.by = "celltype",label = T, label.size = 0)
p1 <- DimPlot(pbmc, group.by = "celltype1",label = T, label.size = 5,cols = use_colors) + NoLegend()
p2 <- DimPlot(pbmc, group.by = "predicted.Group", label = F, label.size = 4,cols = use_colors)
p <- p1|p2
p
DimPlot(pbmc, group.by = "predicted.Group",split.by = "Class", label = F, label.size = 4,cols = use_colors)
ggsave("CellType/Celltype_SingleR.png", p, width = 12, height = 5)

table(pbmc$Time)
pbmc1 <- subset(pbmc, subset = Time %in% c("Post"))
table(pbmc$celltype,pbmc$predicted.Group)
table(pbmc1$Response,pbmc1$predicted.Group)
table(reference$Group)
head(pbmc@meta.data)
saveRDS(pbmc, "mapping/Feldman_predicted.rds")

table(pbmc$RNA_snn_res.0.1)
pbmc$age <- recode(pbmc$orig.ident,'P1' = '1.6','P2' = '5.4','P3' = '3.3','P4' = '2.1','P5' = '1.9','P6' = '4.7')
pbmc$group <- recode(pbmc$predicted.Group,'Cd' = 'Preserve','Tg' = 'Shift','Ad' = 'Shift','Ed' = 'Shift','Sc' = 'Shift')
pbmc$group <- recode(pbmc$predicted.Group,'Cd' = 'Preserve','Tg' = 'Preserve','Ad' = 'Shift','Ed' = 'Shift','Sc' = 'Shift')


meta <- as.data.frame(pbmc@meta.data)
head(meta)
#write.table(meta,"Kawasaki_meta.txt")

meta$celltype1 <- recode(meta$main, 'B cell' = 'B_cells', 'T cell' = 'T_cells', 'Myeloid' = 'Myeloids', 
                         'NK' = 'NKs', 'Platelet' = 'Platelets')
table(meta$celltype1,meta$predicted.Group)
meta$celltype1 <- factor(meta$celltype1, levels = c('B_cells',"Hsc","Myeloids",'others',"NKs","T_cells","Platelets"))
meta$age <- factor(meta$age, levels = c('1.6',"1.9","2.1",'3.3',"4.7","5.4"))
meta$predicted.Group <- factor(meta$predicted.Group, levels = c('Cd',"Tg","Ad",'Ed',"Sc"))

pbmc$Age <- factor(pbmc$Age, levels = c('0.3','0.6',"2",'4',"5","6","7" ,"8",'9',"10","10.4","11","12",'14',"15"))
meta$predicted.Group <- factor(meta$predicted.Group, levels = c('Cd',"Tg","Ad","Sc"))

use"deepskyblue2","steelblue","#0B775E", "#E2D200", "brown2","#B40F20")
ggplot(meta, aes(y = nCount_RNA, axis1 = age, axis2 = celltype1, axis3 = predicted.Group)) +
  geom_alluvium(aes(fill = celltype1),  width = 1/8, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = c(B_cells = "#70493D", T_cells = "#E2AC76", Myeloids = "#3F752B", NKs = "brown2", 
                               Platelets = "#81B0E4", Hpc = "brown2")) +  guides(fill = FALSE) +
  geom_stratum(alpha = .25, width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("age", "predicted.Group", "celltype")) +
  ggtitle("PBMCs of Kawasaki disease mapping to TIAs")+ theme_classic()  #+ coord_flip() 

ggp(y = nCount_RNA, axis1 = age, axis2 = predicted.Group, axis3 = celltype)) +
  geom_alluvium(aes(fill = predicted.Group),  width = 1/8, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = c(Cd = "#81B0E4", Tg = "#70493D", Ad = "#FF4500",
                               Ed = "brown2", Sc = "#B40F20")) +  guides(fill = FALSE) +
  geom_stratum(alpha = .25, width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("age", "predicted.Group", "celltype")) +
  ggtitle("PBMCs of Kawasaki disease mapping to TIAs")+ theme_classic()  #+ coord_flip() 

# COVID
table(meta$celltype1)
ggplot(meta, aes(y = nCount_RNA, axis1 = Age, axis2 = celltype1, axis3 = predicted.Group)) +
  geom_alluvium(aes(fill = celltype1),  width = 1/8, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = c(B_cells = "#70493D", T_cells = "#E2AC76", Myeloids = "#3F752B",NKs = "brown2", 
                               Platelets = "#81B0E4", others = "#0B775E", Hpc = "steelblue")) +  guides(fill = FALSE) +
  geom_stratum(alpha = .25, width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("age", "predicted.Group", "celltype")) +
  ggtitle("PBMCs of Kawasaki disease mapping to TIAs")+ theme_classic()  #+ coord_flip() 

ggplot(meta, aes(y = nCount_RNA, axis1 = Age, axis2 = predicted.Group, axis3 = celltype1)) +
  geom_alluvium(aes(fill = predicted.Group),  width = 1/8, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = c(Cd = "#70493D", Tg = "#81B0E4", Ad = "#FF4500",Sc = "#B40F20")) +  guides(fill = FALSE) +
  geom_stratum(alpha = .25, width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("age", "predicted.Group", "celltype")) +
  ggtitle("PBMCs of Kawasaki disease mapping to TIAs")+ theme_classic()  #+ coord_flip() 


