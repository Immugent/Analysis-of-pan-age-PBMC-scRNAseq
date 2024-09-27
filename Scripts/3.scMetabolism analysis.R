#scMetabolism analysis

rm(list=ls())
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(Seurat)

pbmc <- readRDS("Data/panage_data.rds")
head(pbmc@meta.data)
table(pbmc$celltype)

DefaultAssay(pbmc) <- "RNA"
#method supports VISION, AUCell, ssgsea, and gsva, which VISION is the default method.
pbmc <- sc.metabolism.Seurat(obj = pbmc, method = "VISION", imputation = F, ncores = 1, metabolism.type = "KEGG")
head(pbmc@meta.data)
input.pathway <-rownames(pbmc@assays$METABOLISM$score)
input.pathway <- c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
DotPlot.metabolism(obj = pbmc, pathway = input.pathway,  phenotype = "celltype2", norm = "y")

BoxPlot.metabolism(obj = pbmc, pathway = "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate", phenotype = "celltype1", ncol = 1)


# All metabolic pathways are displayed in seven categories 
Amino_acid <- c("Tyrosine metabolism", "Valine, leucine and isoleucine degradation", "Alanine, aspartate and glutamate metabolism",
                "Cysteine and methionine metabolism", "Phosphonate and phosphinate metabolism", "Histidine metabolism", "Taurine and hypotaurine metabolism"
                , "Linoleic acid metabolism", "Lysine degradation", "Arginine biosynthesis", "Arginine and proline metabolism","beta-Alanine metabolism"
                ,"Glutathione metabolism","Glycine, serine and threonine metabolism","Selenocompound metabolism","Tryptophan metabolism")

Lipid <- c("Fatty acid biosynthesis", "Fatty acid degradation", "Biosynthesis of unsaturated fatty acids", "Arachidonic acid metabolism",
           "Steroid biosynthesis", "Steroid hormone biosynthesis", "alpha-Linolenic acid metabolism","Ether lipid metabolism"
           ,"Fatty acid elongation","Glycerolipid metabolism","Glycerophospholipid metabolism","Primary bile acid biosynthesis","Sphingolipid metabolism")

Carbohydrate <- c("Citrate cycle (TCA cycle)", "Fructose and mannose metabolism", "Galactose metabolism","Pyruvate metabolism","Ascorbate and aldarate metabolism",
                  "Glycolysis / Gluconeogenesis", "Glyoxylate and dicarboxylate metabolism", "Pentose phosphate pathway","Amino sugar and nucleotide sugar metabolism",
                  "Butanoate metabolism","Inositol phosphate metabolism","Pentose and glucuronate interconversions","Propanoate metabolism","Pyruvate metabolism",
                  "Starch and sucrose metabolism")

Glycan <- c("Glycosphingolipid biosynthesis - lacto and neolacto series", "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis", "Mannose type O-glycan biosynthesis",
            "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate", "Glycosphingolipid biosynthesis - ganglio series", "Glycosphingolipid biosynthesis - globo and isoglobo series"
            ,"Glycosaminoglycan biosynthesis - heparan sulfate / heparin","Glycosaminoglycan biosynthesis - keratan sulfate","Glycosaminoglycan degradation",
            "Mucin type O-glycan biosynthesis","N-Glycan biosynthesis","Other glycan degradation","Other types of O-glycan biosynthesis")

Cofactor.vitaman <- c("Porphyrin and chlorophyll metabolism", "Retinol metabolism", "Riboflavin metabolism","Folate biosynthesis","Thiamine metabolism",
                      "One carbon pool by folate", "Vitamin B6 metabolism", "Nicotinate and nicotinamide metabolism","Pantothenate and CoA biosynthesis",
                      "Ubiquinone and other terpenoid-quinone biosynthesis")

Energy <- c("Sulfur metabolism", "Oxidative phosphorylation", "Nitrogen metabolism")
#Energy <- c("Glycolysis / Gluconeogenesis","Oxidative phosphorylation")

Nucleotide <- c("Purine metabolism", "Pyrimidine metabolism")

DotPlot.metabolism.2 = function (obj, pathway, phenotype, norm = "y") 
{ input.norm = norm
input.pathway <- as.character(pathway)
input.parameter <- phenotype
metadata <- obj@meta.data
metabolism.matrix <- obj@assays$METABOLISM$score
metadata[, input.parameter] <- as.character(metadata[, input.parameter])
metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, 
])

gg_table <- c()
for (i in 1:length(input.pathway)) {
  gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], 
                                    input.pathway[i], metabolism.matrix_sub[, i]))
}
gg_table <- data.frame(gg_table)
gg_table_median <- c()
input.group.x <- unique(as.character(gg_table[, 1]))
input.group.y <- unique(as.character(gg_table[, 2]))
for (x in 1:length(input.group.x)) {
  for (y in 1:length(input.group.y)) {
    gg_table_sub <- subset(gg_table, gg_table[, 1] == 
                             input.group.x[x] & gg_table[, 2] == input.group.y[y])
    gg_table_median <- rbind(gg_table_median, cbind(input.group.x[x], 
                                                    input.group.y[y], median(as.numeric(as.character(gg_table_sub[, 
                                                                                                                  3])))))
  }
}
gg_table_median <- data.frame(gg_table_median)
gg_table_median[, 3] <- as.numeric(as.character(gg_table_median[, 
                                                                3]))
gg_table_median_norm <- c()
input.group.x <- unique(as.character(gg_table[, 1]))
input.group.y <- unique(as.character(gg_table[, 2]))
range01 <- function(x) {
  (x - min(x))/(max(x) - min(x))
}

if (input.norm == "y") 
  for (y in 1:length(input.group.y)) {
    gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 2] == input.group.y[y])
    norm_value <- range01(as.numeric(as.character(gg_table_median_sub[,3])))
    gg_table_median_sub[, 3] <- norm_value
    gg_table_median_norm <- rbind(gg_table_median_norm, 
                                  gg_table_median_sub)
  }

if (input.norm == "x") 
  for (x in 1:length(input.group.x)) {
    gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 1] == input.group.x[x])
    norm_value <- range01(as.numeric(as.character(gg_table_median_sub[,3])))
    gg_table_median_sub[, 3] <- norm_value
    gg_table_median_norm <- rbind(gg_table_median_norm, 
                                  gg_table_median_sub)}

if (input.norm == "na") 
  gg_table_median_norm <- gg_table_median
gg_table_median_norm <- data.frame(gg_table_median_norm)
gg_table_median_norm[, 3] <- as.numeric(as.character(gg_table_median_norm[,3]))

library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")
if(is.factor(pathway)){
  gg_table_median_norm$X2 = factor(gg_table_median_norm$X2 ,levels = levels(pathway))
}

if(is.factor(countexp.Seurat@meta.data[,phenotype])){
  gg_table_median_norm$X1 = factor(gg_table_median_norm$X1 ,
                                   levels = levels(countexp.Seurat@meta.data[,phenotype]))
}

ggplot(data = gg_table_median_norm, aes(x = gg_table_median_norm[, 1], y = gg_table_median_norm[, 2], color = gg_table_median_norm[,  3])) + 
  geom_point(data = gg_table_median_norm, aes(size = gg_table_median_norm[, 
                                                                          3])) + ylab("Metabolic Pathway") + xlab(input.parameter) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                hjust = 1), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  scale_color_gradientn(colours = pal) + labs(color = "Value", 
                                              size = "Value") + NULL
}

#input.pathway <- factor(c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)"),
#levels = c("Oxidative phosphorylation", "Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)")) #调整纵轴通路

pbmc$Group <- factor(pbmc$Group,levels = c("Cd","Tg","Ad","Ed","Sc")) #调整横轴分组
head(pbmc@meta.data)
table(pbmc$celltype2)
countexp.Seurat <- pbmc
DotPlot.metabolism.2(obj = pbmc, pathway = Amino_acid,  phenotype = "celltype1", norm = "y")
ggsave("figs/Mo.Amino_acid.pdf",width = 18,height = 5)

DotPlot.metabolism.2(obj = pbmc, pathway = Lipid,  phenotype = "celltype1", norm = "y")
ggsave("figs/Mo.Lipid.pdf",width = 18,height = 4)

DotPlot.metabolism.2(obj = pbmc, pathway = Carbohydrate,  phenotype = "celltype1", norm = "y")
ggsave("figs/Mo.Carbohydrate.pdf",width = 18,height = 4)

DotPlot.metabolism.2(obj = pbmc, pathway = Glycan,  phenotype = "celltype1", norm = "y")
ggsave("figs/Mo.Glycan.pdf",width = 18,height = 4)

DotPlot.metabolism.2(obj = pbmc, pathway = Cofactor.vitaman,  phenotype = "celltype1", norm = "y")
ggsave("figs/Mo.vitaman.pdf",width = 18, height = 4)

DotPlot.metabolism.2(obj = pbmc, pathway = Energy,  phenotype = "celltype1", norm = "y")
ggsave("figs/Mo.Energy.pdf",width = 18,height = 3)

DotPlot.metabolism.2(obj = pbmc, pathway = Nucleotide,  phenotype = "celltype1", norm = "y")
ggsave("figs/Mo.Nucleotide.pdf",width = 18,height = 2)
