############## Venn diagram for DEGs of each group ##############
library (VennDiagram)
library(gridSVG)

query1 <- read.csv("28750res/patient-control_diff.xls", sep="\t", header = T)
query2 <- read.csv("57065res/patient-control_diff.xls", sep="\t", header = T)
query3 <- read.csv("95233res/patient-control_diff.xls", sep="\t", header = T)

V1 = query1$gene
V2 = query2$gene
V3 = query3$gene

venn.plot <- venn.diagram(x= list("GSE28750" = V1,
                                  "GSE57065" = V2,
                                  "GSE95233" = V3),
                          filename = "3_venn.tiff", lwd = 0.5,lty = 1, imagetype = "tiff",
                          height = 600, width = 600,resolution =500, 
                          col=c("#FF7582","#A163F7","#45E3FF"),
                          fill=c("#FF7582","#A163F7","#45E3FF"),
                          alpha = 0.50, cex=0.3, cat.cex=0.3) 