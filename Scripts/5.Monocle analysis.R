#monocle2 analysis

rm(list=ls())
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(openxlsx)

pbmc <- readRDS("new_CD8T.rds")
head(pbmc@meta.data)
Idents(pbmc) <- pbmc$main1
table(Idents(pbmc))
pbmc <- subset(pbmc, idents = c("Effector_CD8T",'Naive_CD8T',"Exhausted_CD8T","IFN_CD8T")) #
pbmc <- subset(pbmc, subset = main1 != "Dims")

data <- GetAssayData(pbmc, assay = "RNA", slot = "counts")
#data <- as(data, 'sparseMatrix')   
pd <- new('AnnotatedDataFrame', data = pbmc@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())


mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=8)

plotdf = pData(mycds)
head(plotdf); dim(plotdf)
table(plotdf$main1)


disp_table <- dispersionTable(mycds)
r.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
length(order.genes)
mycds <- setOrderingFilter(mycds, order.genes)
p3 <- plot_ordering_genes(mycds)
p3

##??ds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')

##??ds <- orderCells(mycds) #??�ds <- orderCells(mycds, root_state = 3)  #???adlot_cell_trajectory(mycds, color_by = "main1", show_branch_points = F, cell_size = 1) +
  scale_color_manual(values = c('#f8766d','#cd9600','#7cae00','#00be67',
                                '#00bfc4','#00a9ff','#c77cff', '#ff61cc'))
table(mycds$main1)
mycds$main1 <- factor(mycds$main1, levels = c("Naive_CD8T","IFN_CD8T","Effector_CD8T","Exhausted_CD8T"))
plot_cell_trajectory(mycds, color_by = "main1")


#Cel2 <- plot_cell_trajectory(mycds, color_by = "State")
plot2

#Pse3 <-  plot_cell_trajectory(mycds, color_by = "Pseudotime", show_branch_points = F, cell_size = 1) + scale_color_viridis_c()
plot3

#?ϲc <- plot1|plot2|plot3

ggsave("CHY/Trajectory_Combination.pdf", plot = plotc, width = 15, height = 5)

#???s$main1 <- factor(mycds$main1, levels = c("CD8_TPEX1","Cycling_T","CD8_TPEX2","CD8_effector",
                                              "CD8_TEX","CD8_TEM-like","CD8_Naive","CD8_IFN-induced"))

p <- plot_cell_trajectory(mycds, color_by = "main1", show_branch_points = F) + facet_wrap(~main1, nrow = 2)
p
ggsave("CHY/Trajectory_Facet.pdf", plot = p, width = 10, height = 6)

##??F){
  mycds <- orderCells(mycds, root_state = 3)
  plot1 <- plot_cell_trajectory(mycds, color_by = "State")
  ggsave("Monocle/Trajectory_State.pdf", plot = plot1, width = 10, height = 6.5)
  
  plot2 <- plot_cell_trajectory(mycds, color_by = "celltype.fine")
  ggsave("Monocle/Trajectory_Celltype.pdf", plot = plot2, width = 10, height = 6.5)
  
  plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
  ggsave("Monocle/Trajectory_Pseudotime.pdf", plot = plot3, width = 10, height = 6.5)
  
  plotc <- plot1|plot2|plot3
  ggsave("Monocle/Trajectory_Combination.pdf", plot = plotc, width = 10, height = 3.5)
  
  p <- plot_cell_trajectory(mycds, color_by = "celltype.fine") + facet_wrap(~celltype.fine, nrow = 3)
  ggsave("Monocle/Trajectory_Facet.pdf", plot = p, width = 10, height = 10)
}

##??ta <- pData(mycds)
s.cells <- subset(pdata, State=="3") %>% rownames() #??�e(s.cells, file = "Monocle/state3.rda")

##??te.csv(pData(mycds), "Monocle/pseudotime.csv")
save(mycds, file = "xijiCD4CTLds_mo2")

##ָenes <- c("LEF1","TCF7","SELL","CCR7","IL7R","IL2RA", "FOXP3","CCL5","GZMB","GZMK")
s.genes <- c("Pdcd1","Atf4")

p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State",ncol = 2)
p1
p2 <- monocle::plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State",ncol = 2)
p2
p3 <- monocle::plot_genes_in_pseudotime(mycds[s.genes,], color_by = "celltype",ncol = 2)
p3
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 20, height = 8)

genes <- c("TCF7","SELL")#"TCF7","SELL","CCR7","IL7R","CCL5","GZMB","GZMK", "NKG7")
genes <- c("Ccr7","Atf4","Gzmb","Ccl5","Pdcd1","Lag3")
p11 <- monocle::plot_genes_in_pseudotime(mycds[genes,],  color_by = "celltype",ncol = 2) #???
ggsave("GENEs_Jitterplot.pdf", plot = p11, width = 20, height = 8)

#???<- plot_cell_trajectory(mycds, color_by = "group") + facet_wrap(~group, nrow = 1)
p4
ggsave("celltype_plot.pdf", plot = p4, width = 20, height = 5)

#???subset <- mycds[c("Ccr7","Sell","Gzmb","Ccl5","Pdcd1","Lag3"),] #"CCR7","SELL","GZMB","CCL5"
monocle::plot_genes_in_pseudotime(cds_subset, ncol = 3, color_by = 'Pseudotime') + scale_color_viridis_c() #???�?l_trajectory(mycds, markers=c("CCR7","LEF1","IL2RA","FOXP3","PRF1","GZMB"), use_color_gradient=T, show_backbone = F,backbone_color = "#F2F2F2",cell_size=0.5, show_branch_points=F) #?????ʺϿ??߱????Ļ???

#Ѱ????� differentialGeneTest(mycds[order.genes,], cores = 1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")#??ʱ
head(Tiff); dim(Time_diff)
Time_diff <- Time_diff[,c(5,2,3,4,1,6)]
head(Time_diff)
Time_diff_filter <- Time_diff[!grepl("^Rp[sl]", Time_diff$gene_short_name, ignore.case = F),]
dim(Time_diff_filter)
Time_diff_filter <- Time_diff_filter[!grepl("^mt-", Time_diff_filter$gene_short_name, ignore.case = F),]
dim(Time_diff_filter)
write.xlsx(Time_diff_filter,"Lixin/tables/Time_diff_all.xlsx",quote=F,row.names = T)

#???????????<- Time_diff_filter[order(Time_diff_filter$qval), "gene_short_name" ][1:100] #??ȡqval??<- c("Ccr7","Sell","Gzmb","Ccl5","Pdcd1","Lag3")
Time_genes <- c("Ccr7","Sell","Gzmb","Ccl5","Pdcd1","Lag3")
Time_genes <- c("Lef1","Sell","Foxp1","Tcf7","Il7r","Bcl2","Xcl1","Irf8","Irf1","Il2rb","Klrb1c","Mki67",
                "Ctla4","Casp3","Ifng","Gzmb","Ccl5","Cxcr6","Lag3","Havcr2","Prf1","Pdcd1","Gzmk","Tox","Gzma") # "Xcl1",
gene <- read.table("Lixin/genes.txt", header = F)
Time_genes <- gene$V1
head(pData(mycds))
p = plot_pseudotime_heatmap(mycds[Time_genes,], num_clusters=3, show_rownames=T, 
                            return_heatmap=T)
p

ggsave("CHY/Time_heatmap.pdf", p, width = 5, height = 15)
dev.off()

#???????????p$tree_row$labels[p$tree_row$order]
head(hp.genes)
Time_diff_sig <- Time_diff_filter[hp.genes, c("gene_short_name", "pval", "qval")]
write.xlsx(Time_diff_sig,"CHY/200-Time_diff_sig.xlsx",quote=F,row.names = T)

##BEAM????
bBEAM(mycds[order.genes,], branch_point = 1, cores = 1)#??ʱ
head(beam_res)
beam_res  <- beam_res[,c(5,2,3,4,1,6)]
dim(beam_res)
beam_res_filter <- beam_res[!grepl("^RP[SL]", beam_res$gene_short_name, ignore.case = F),]
beam_res_filter <- beam_res_filter[!grepl("^MT-", beam_res_filter$gene_short_name, ignore.case = F),]
dim(beam_res_filter)
write.xlsx(beam_res_filter,"CHY/BEAM_all.xlsx",quote=F,row.names = T)

BEAM_genes <- beam_res_filter[order(beam_res_filter$qval), "gene_short_name"][1:200] #??ȡqval??�enes_branched_heatmap(mycds[BEAM_genes,],  branch_point = 1,
                                 num_clusters = 5, show_rownames = T, return_heatmap = T)
ggsave("CHY/BEAM_heatmap.pdf", p$ph_res, width = 5, height = 15)

#???????????p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- beam_res[hp.genes, c("gene_short_name", "pval", "qval")]
write.xlsx(BEAM_sig,"CHY/200-BEAM_sig.xlsx",quote=F,row.names = T)

save(mycds,file="superold/monocle2.Robj")


#չʾĳ??ͨ