#Cytotrace analysis

rm(list=ls())
library(CytoTRACE)

pbmc <- readRDS("subN.rds")
head(pbmc@meta.data)

mat <- as.matrix(pbmc@assays$RNA@counts)
mat[1:4,1:4];dim(mat)

results <- CytoTRACE(mat = mat)  # enableFast = FALSE 

dir.create("CytoTRACE/")
plotCytoGenes(results, numOfGenes = 10, outputDir = "CytoTRACE/")
table(pbmc$age)
phe <- pbmc$main1 
phe = as.character(phe)
names(phe) <- rownames(pbmc@meta.data)

plotCytoTRACE(results, phenotype = phe, gene = c("Mki67"), outputDir = "CytoTRACE/")



