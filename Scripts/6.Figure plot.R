# figure plot

pbmc <- readRDS("Data/panage_pbmc.rds") 

head(pbmc@meta.data)
table(pbmc$celltype2)
object.markers <- FindMarkers(pbmc, ident.1 = 'Shift_Myeloid',ident.2 = 'Preserve_Myeloid',
                              group.by = 'celltype2',logfc.threshold = 0, min.pct = 0, pseudocount.use = 0.01)
object.markers$names <- rownames(object.markers)
head(object.markers)
#sig_dge.all <- subset(object.markers, p_val_adj<0.05&abs(avg_log2FC)>0.15) #???в???????
#View(sig_dge.all)
library(dplyr)
object.markers <- object.markers %>%
  mutate(Difference = pct.1 - pct.2)
head(object.markers)
library(ggplot2)
library(ggrepel)

ggplot(object.markers, aes(x=Difference, y=avg_log2FC)) +
  geom_point(size=0.5, color="#999999") +
  geom_label_repel(data=subset(object.markers, avg_log2FC >= 1 & Difference >= 0.2 & pct.2 >= 0.05),
                   aes(label=names), label.padding = 0.1, fill="tomato2", segment.size = 0.25, size=2.5)+
  theme_classic()
ggsave("TopMarkerVol1.pdf", height=8, width=8)

object.markers$group=0
for (i in 1:nrow(object.markers)){
  if (object.markers$avg_log2FC[i] >= 1 & object.markers$Difference[i] >= 0.2 & object.markers$pct.2[i] >= 0.05){
    object.markers$group[i]='up'
  }
  else if(object.markers$avg_log2FC[i] <= -1 & object.markers$Difference[i] <= -0.02 & object.markers$pct.1[i] >= 0.05){
    object.markers$group[i]='down'
  }
  else {
    object.markers$group[i]='no'
  }
}

ggplot(object.markers, aes(x=Difference, y=avg_log2FC)) +
  geom_point(size=1,aes(color=group)) +
  scale_color_manual(values=c('blue','grey','red'))+
  geom_label_repel(data=subset(object.markers, group !='no'), aes(label=names), segment.size = 0.25, size=3)+
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_classic()
ggsave("TopMarkerVol2.pdf", height=8, width=8)


DotPlot_2 <- function(object,
                      features,
                      assay = NULL, 
                      scale = T,
                      dot.range.min = 0,
                      dot.range.max = 3.5,
                      Combine=T,
                      legend.position = "right",
                      label.size = 4,
                      label_widths = 0.1,
                      x.lab = NULL,
                      y.lab = NULL,
                      title = NULL,
                      text.size = 8,
                      text.angle = 90,
                      text.vjust = 0.5,
                      text.hjust = 1,
                      group.by = NULL,
                      color.use = NULL,
                      cols = c("lightgrey", 
                               "blue"),
                      legend.key.size = 0.5,
                      col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                      idents = NULL, split.by = NULL, cluster.idents = FALSE, 
                      scale.by = "radius", scale.min = NA, scale.max = NA,
                      ...
){
  library(dplyr)
  library(ggplot2)
  library(aplot)
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust, 
                                              vjust = text.vjust), #,vjust = 0.5
                   panel.grid=element_blank(), # ȥ??????
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )
  
  if(!is.null(group.by)){
    Idents(object) = group.by
  }
  
  if(Combine){
    if(is.null(color.use)){
      color.use <- alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 1)
    }
    df <- data.frame(x = 0, y = levels(object), stringsAsFactors = F )
    df$y <- factor(df$y, levels = df$y )
    p1 <- ggplot(df, aes(x, y, color = factor(y))) +
      geom_point(size = label.size, show.legend = F) +
      scale_color_manual(values = color.use) +
      theme_classic() +
      scale_x_continuous(expand = c(0,0)) + mytheme + 
      theme(
        plot.margin = margin(r=0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 0,color="white"),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  }
  
  p2 = DotPlot(object = object,cols = cols,
               assay = assay,
               col.min = col.min, col.max = col.max, dot.min = dot.min, dot.scale = dot.scale, 
               idents = idents,split.by = split.by, cluster.idents = cluster.idents, 
               scale.by = scale.by, scale.min = scale.min, scale.max = scale.max,
               features = features,scale = scale)+theme_bw()+
    mytheme + labs(x = x.lab,y = y.lab,title = title)  +
    scale_size(range = c(dot.range.min,dot.range.max))+
    theme(legend.key.size = unit(legend.key.size, "cm"))+ 
    guides(size = guide_legend(title = "Per.Exp"),
           color = guide_colorbar(title = paste("Aver.Exp.","Scaled",sep = "\n")))
  if(Combine){
    p3 = p2 %>% insert_left(p1, width=label_widths)
    return(p3)
  }else{
    return(p2) 
  }
}

check_genes = c("IFI16", "RBPJ", "TOX",  "ZBED2", "VDR", "BATF", "IKZF4","TOX2","PRDM1", "STAT3", 
                "TIGIT","CTLA4", "ENTPD1", "PDCD1", "LAYN", "HAVCR2", "LAG3", "BTLA","TNFRSF9","CD27", 
                "IFNG","GZMK", "GZMA","NKG7", "PRF1", "TNF", "TBX21", "IL2")  
#B cell
markers <- c('PTPRC',"MS4A1",'CD19','CD1C','CD79A','CD79B',"IL4R","TCL1A","XBP1",'IGHG1','MZB1','SDC1',"MME","CD24","PAX5","CD34","CD38","TOP2A","MKI67"
             ,"FAS","CD80","CD86","ITGAX","CD27","CD28",'CR2',"CD22","FCER2","CXCR5","CD40","FCRL3","TNFRSF13B","TNFRSF17","PRDM1","CXCR4","CD83","CD93","CD5")
#T cell
markers <- c('PTPRC',"CD3D","CD3E","CD3G","CD4","CD8A","CD8B","FOXP3","IL2RA","CD40LG","CCR7","SELL","LEF1","GZMA","GZMB","GZMK","IFNG","PRF1","CCL5",
             "NKG7","IL7R","CD44","CD69", "CD27","CD28","TOP2A", "MKI67","TARP","TRGV9","CD160","FEZ1","TOX","PDCD1", "TIGIT","HAVCR2","CTLA4","LAG3")
#NK
markers <- c("IL32",'KLRF1',"NCAM1",'NCR1',"PDCD1","CX3CR1", 'GNLY','NKG7',"FCGR3A","XCL1","XCL2","KLRD1","KLRB1","CD160","TOP2A","MKI67",
             "CD52","IFI6","IFI44L","ISG15","MX1","CST7","GZMA","GZMB","GZMK","GZMH","IFNG","PRF1","CCL3","CCL4","CCL5","FGFBP2","LIF","SELL","HLA-DPB1","HLA-DRB1","HLA-A","HLA-C","HLA-E")
#Myeloid
markers <- c("LYZ","FCN1","MS4A6A", "CD1E", "LAMP3", "IL3RA","VCAN","MNDA","MX1","CD68","CSF1R","CD163","CD14","FCGR3A","AIF1","C1QA","C1QB","C1QC","MARCO","CD1C",
             "ITGAX","CD83","CD80","CD86","CR1","CR2","FCER2","FCER1A","CLEC9A","CLEC4C","LILRA4","KLRD1","NRP1","JCHAIN","CXCR3",  "MZB1","ITM2C","SSR4",
             "TLR2","LTF","FUT4", "PDIA4","NLRP3","PECAM1","CSF3R","FPR1","MPO","CEACAM8","S100A8","S100A9","CXCR2","FCGR3B","SELL","C5AR1")

genes <- c("CD4","CD8A","CD8B","CCR7","SELL","LEF1","TCF7","FOXP3","IL2RA","CTLA4","GZMA","GZMB","GZMK","IFNG","PRF1","CCL5",
           "NKG7","IL7R","CD44","CD69","TOP2A", "MKI67","LTB","TRAT1","TLE5")#T
genes <- c("CD19","MS4A1","IL4R","TCL1A","IGHD","CD83","CXCR4","CXCR5",'NR4A1',"NR4A2","TNFRSF13B","CLECL1",'S100A8','S100A9',"XBP1","PRDM1",'MZB1',"TNFRSF17","CD27", 
           "CD5","IL32")#B
genes <- c("IL32","CCL5","GZMH","CD52","KLRD1","GZMA","CCL3","KLRB1",'KLRF1','NKG7',"FCGR3A","PRF1","GZMB",
           "CST7",'GNLY',"HLA-E","HLA-A","HLA-C","XCL1","XCL2","NCAM1","SELL","GZMK")#NK
genes <- c("NLRP3","C5AR1","FPR1","CSF3R","CD14","LYZ", "VCAN","MNDA","FCGR3A","CD68","CSF1R","PECAM1","ITGAX","CD86","CD1C",
           "CD83","FCER1A","IL3RA","CLEC4C","LILRA4","JCHAIN","CXCR3","MZB1")#myeloid , "S100A4","S100A8","S100A9"

p0 = DotPlot_2(object = pbmc,
               group.by = "celltype",
               Combine = F,
               dot.range.max = 3,
               legend.key.size = 0.4,
               dot.range.min = 0.5,label.size = 2.8,
               features = genes)&
  scale_color_distiller(palette = 'RdYlBu')&coord_flip()&
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1));p0

ggsave("B_dotplot.pdf", path = "./",  width = 3.5, height = 4) #2x3=5,7
#

library(ape)
library(dplyr)
library(psych)
library(ComplexHeatmap)
library(edgeR)
library(limma)
library(grid)
library(circlize)
library("GetoptLong")
library(corrplot)

df1<-as.data.frame(AverageExpression(pbmc,features = rownames(pbmc),group.by = "predicted.celltype")[["RNA"]])
colnames(df1) <- paste(colnames(df1),"Query",sep = "_")

#secondary_type #primary_type
df2<-as.data.frame(AverageExpression(reference,features = rownames(pbmc),group.by = "celltype")[["RNA"]])
colnames(df2) <- paste(colnames(df2),"Ref",sep = "_")
df1$gene <- rownames(df1)
df2$gene <- rownames(df2)
head(df1)
head(df2)
df <- merge(df1,df2,by = "gene")#[,c(-3,-4)]
head(df)

rownames(df) <- df$gene
df <- df[,-1]
dim(df); head(df)

library(factoextra)
df.s <- scale(df,center = T,scale = T)
# 03. ???ξ???
df.t <- as.matrix(t(df.s))
#d <- dist(df.t,method = "euclidean")
#?get_dist
d <- get_dist(df.t,method = "pearson") 
hc <- hclust(d,method = "ward.D2")
py <- as.phylo(hc)
row_dend = as.dendrogram(hc)
plot(py,type="phylogram",use.edge.length = 0,
     cex = 0.8,edge.width = 2,font = 3,label.offset=0.01,
     adj = 0.1)

df.pearson <- cor(df.s,method="pearson")
df.pearson
n <- dim(df)[2]-1
colnames(df.pearson) <- gsub("X","",colnames(df.pearson))
colnames(df.pearson) <- gsub("_Mean","",colnames(df.pearson))
rownames(df.pearson) <- gsub("X","",rownames(df.pearson))
rownames(df.pearson) <- gsub("_Mean","",rownames(df.pearson))

f1 <- colorRamp2(breaks = c(0.7,0.8,0.9,1),
                 c("#114A85", "#FFFFFF","#B61F2E","#8D1114"))

p <- Heatmap(df.pearson, name ="Correlation",
             column_title = " ",row_title = " ",
             # rect_gp = gpar(col = "white"),
             col = f1,
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(df.pearson[i, j] > -1)
                 grid.text(sprintf("%.2f", df.pearson[i, j]),
                           x, y, gp = gpar(fontsize = 10))
             }, 
             cluster_rows = row_dend,cluster_columns = row_dend,
             show_column_names = TRUE,show_row_names = TRUE,row_names_side = "right",
             gap = unit(1,"mm"),row_dend_side = "left",row_title_side = "left",
             column_dend_side = "top",column_title_rot = 0,column_title_side = "top")
p

#OR analysis
library(plyr)

pbmc <- readRDS("misc_final.rds")

out.prefix <- "figs"
#meta.tb <- read.table("tables/meta.txt", header = T)
meta.tb <- as.data.frame(pbmc@meta.data)
meta.tb <- as.data.frame(all_counts_rel)
meta.tb[1:10,1:10]

head(meta.tb)
colnames(meta.tb)

table(pbmc$Class);table(pbmc$sample)
table(meta.tb$type);table(meta.tb$celltype)
#meta.tb$Group <- factor(meta.tb$Group,levels = c("Cd","Tg","Ad","Ed","Sc")) #自定义顺序

do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$celltype,
                          colname.patient = "subject_code",
                          loc = cellInfo.tb$group, 
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]}
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = pdf.width, pdf.height = pdf.height)
  if(verbose==1){ return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                              "p.dist.tb"=p.dist.tb,
                              "OR.dist.tb"=OR.dist.tb,
                              "OR.dist.mtx"=OR.dist.mtx))
  }else{    return(OR.dist.mtx)}}

test.dist.table <- function(count.dist,min.rowSum=0)
{count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
sum.col <- colSums(count.dist)
sum.row <- rowSums(count.dist)
count.dist.tb <- as.data.frame(count.dist)
setDT(count.dist.tb,keep.rownames=T)
count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
colnames(count.dist.melt.tb) <- c("rid","cid","count")
count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
  this.row <- count.dist.melt.tb$rid[i]
  this.col <- count.dist.melt.tb$cid[i]
  this.c <- count.dist.melt.tb$count[i]
  other.col.c <- sum.col[this.col]-this.c
  this.m <- matrix(c(this.c,
                     sum.row[this.row]-this.c,
                     other.col.c,
                     sum(sum.col)-sum.row[this.row]-other.col.c),
                   ncol=2)
  res.test <- fisher.test(this.m)
  data.frame(rid=this.row,
             cid=this.col,
             p.value=res.test$p.value,
             OR=res.test$estimate)}))
count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                by=c("rid","cid"))
count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value, "BH")]
return(count.dist.melt.ext.tb)}

OR.CD8.list <- do.tissueDist(cellInfo.tb=meta.tb,
                             out.prefix=sprintf("BC", out.prefix),
                             pdf.width=8, pdf.height=10, verbose=1)
OR

# progeny
library(progeny)
library(gplots)
library(tidyr)
#BiocManager::install("progeny")

pbmc <- progeny(pbmc, scale = F, organism= "Mouse", top=500, perm=1, return_assay=T) #Human Mouse
pbmc <- ScaleData(pbmc, assay = "progeny")

progeny_scores <- as.data.frame(t(GetAssayData(pbmc, assay = "progeny", slot = "scale.data")))
progeny_scores$cell_id <- rownames(progeny_scores)
progeny_scores <- gather(progeny_scores, Pathway, Activity, -cell_id)

head(pbmc@meta.data)
cells_clusters <- FetchData(pbmc, c("Group", "celltype"))
cells_clusters$cell_id <- rownames(cells_clusters)

progeny_scores <- inner_join(progeny_scores, cells_clusters)
head(progeny_scores)

summarized_progeny_scores <- progeny_scores %>%
  group_by(Pathway, celltype) %>%
  summarise(avg = mean(Activity), std = sd(Activity)) %>%
  pivot_wider(id_cols = Pathway, names_from = celltype, values_from = avg) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

pdf("subN/progeny.pdf", width = 6, height = 8)
heatmap.2(summarized_progeny_scores, trace = "none", density.info = "none", col = bluered(100))
dev.off()
