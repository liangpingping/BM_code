########################################
################################
################################
#######  GO  ################
################################
################################
########################################
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05        
qvalueFilter=0.05           
pair <- "patient-control"

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


sig.gene<-read.table(file=paste0(pair,"_diff.xls"), sep='\t',check.names=F,header = T)
head(sig.gene)
symbol<-sig.gene[,1] #ID所在的列
entrezid <- bitr(symbol,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
head(entrezid)

kk=enrichGO(gene = entrezid$ENTREZID,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
write.table(GO,file=paste0(pair,"_GO_all.txt"),sep="\t",quote=F,row.names = F)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file=paste0(pair,"_GO.txt"),sep="\t",quote=F,row.names = F)

showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file=paste0(pair,"_go-bubble.pdf"))
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", 
            color = colorSel,label_format = 100) + 
  facet_grid(ONTOLOGY~., scale='free')+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank())
print(bub)
dev.off()

library(DOSE)
rt=read.table(paste0(pair,"_GO.txt"),sep="\t",check.names=F,header=T)    
test1 <- rt[1:10,]
p <- ggplot(test1,aes(x = parse_ratio(GeneRatio),y = forcats::fct_reorder(Description,p.adjust,.desc=TRUE)))+
  geom_point(aes(color = p.adjust,size = Count))+
  labs(size="Count", x="GeneRatio",y="")+
  scale_color_gradient(low = "red", high = "blue")+ 
  theme_bw()+guides(color = guide_colorbar(reverse = TRUE))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank())
p
pdf(paste0(pair,"_go_bp.pdf"))
p
dev.off()

mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14, hjust = 0.5,face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

p <- ggplot(test1, aes(x = parse_ratio(GeneRatio), y = reorder(Description, GeneRatio))) + #将Description按照RichFactor进行排序
  geom_segment(aes(xend=0, yend = Description),size = 1.1,color = "grey") +
  geom_point(aes(color=-log10(p.adjust), size = Count)) +
  scale_color_distiller(palette = "Spectral",direction = -1)+
  scale_size_continuous(range=c(2, 8)) + #气泡大小范围
  theme_minimal() +
  labs(x = "GeneRatio",
       y = NULL,
       title = "GO terms Enrichment") + mytheme +theme_classic()
pdf(paste0(pair,"_go_bp_lollipopplot.pdf"), width = 7, height = 4)
p
dev.off()

########################################
################################
################################
#######  KEGG  ################
################################
################################
########################################
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05        
qvalueFilter=0.05           

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

sig.gene<-read.table(file=paste0(pair,"_diff.xls"), sep='\t',check.names=F,header = T)
head(sig.gene)
symbol<-sig.gene[,1] #ID所在的列
entrezid <- bitr(symbol,
                 fromType = "SYMBOL",#现有的ID类型
                 toType = "ENTREZID",#需转换的ID类型
                 OrgDb = "org.Hs.eg.db")
head(entrezid)

kk <- enrichKEGG(gene = entrezid$ENTREZID, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(entrezid$SYMBOL[match(strsplit(x,"/")[[1]],as.character(entrezid$ENTREZID))],collapse="/")))

write.table(KEGG,file=paste0(pair,"_KEGG_all.txt"),sep="\t",quote=F,row.names = F)
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file=paste0(pair,"_KEGG.txt"),sep="\t",quote=F,row.names = F)

showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file=paste0(pair,"_ko-bubble.pdf"),width = 8,height = 6)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",label_format = 100, color = colorSel)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank())
dev.off()

library(ggplot2)
library(DOSE)
test <- read.table(paste0(pair,"_KEGG.txt"), sep="\t", header = T)
test1 <- test[1:20,]
p <- ggplot(test1,aes(x = parse_ratio(GeneRatio),y = forcats::fct_reorder(Description,p.adjust,.desc=TRUE)))+
  geom_point(aes(color = p.adjust,size = Count))+
  labs(size="Count", x="GeneRatio",y="")+
  scale_color_gradient(low = "red", high = "blue")+ 
  theme_bw()+guides(color = guide_colorbar(reverse = TRUE))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank())
p
ggsave(paste0(pair,"_KEGG_20_dot.pdf"),width = 6, height=5)

p <- ggplot(test1,aes(x = Count,y = forcats::fct_reorder(Description,p.adjust,.desc=TRUE),fill=p.adjust))+
  geom_bar(stat="identity")+
  labs(x="",y="")+
  scale_fill_gradient(low = "red", high = "blue")+ 
  theme_bw()+guides(color = guide_colorbar(reverse = TRUE))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank())
p
ggsave(paste0(pair,"_KEGG_20_bar.pdf"),height=5)