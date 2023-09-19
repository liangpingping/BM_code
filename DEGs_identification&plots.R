inputFile="D:/geoData/sepsis/28750/input_expr_mean.csv"   
typeFile="D:/geoData/sepsis/28750/pro.txt"   

fdrFilter=0.05                
logFCfilter=1 

rt=read.table(inputFile,sep=",",header=T,check.names=F, row.names = 1)
range(rt, na.rm = T)
qx <- as.numeric(quantile(rt, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  rt[which(rt <= 0)] <- NaN
  rt <- log2(rt)
  write.table(rt, "input_expr_log.txt", row.names = T, sep="\t", quote = F)}

########################################
# box-and-whisker plot
boxplot(rt,outline = F) 
boxplot(rt, boxwex=0.7, notch=T, outline=FALSE, las=2)
plotDensities(rt, legend=F)
# mean-variance trend
ex <- na.omit(rt) # eliminate rows with NAs
plotSA(lmFit(ex))
########################################
express.norm<- normalizeBetweenArrays(as.matrix(rt)) # normalize data
#express.norm <- normalize.quantiles(as.matrix(rt))
boxplot(express.norm, outline=F, col="darkgreen") #均一化后
rt <- express.norm
type=read.table(typeFile,sep="\t",header=T,check.names=F, row.names = 1)
sharedID <- intersect(rownames(type),colnames(rt))

data <- rt[,sharedID]
norm_data <- cbind(as.data.frame(rownames(data)),data)
write.table(norm_data,file="normExp.txt",sep="\t",quote=F,row.names=F)

modType <- type[sharedID,] 
design <- model.matrix(~0+factor(modType))
colnames(design)
colnames(design) <- c("control","patient") #  colnames(design) <- levels(factor(modType))

fit <- lmFit(data,design)
pair <- "patient-control"
cont.matrix<-makeContrasts(pair,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
rownames <- as.data.frame(rownames(allDiff))
colnames(rownames) <-'gene'
all <- cbind(rownames,allDiff)
write.table(all,file=paste0(pair,"_allDEGs.xls"),sep="\t",quote=F,row.names=F)

#write table
diffSig <- all[with(all, (abs(logFC)>logFCfilter & adj.P.Val < fdrFilter )), ]
write.table(diffSig,file=paste0(pair,"_diff.xls"),sep="\t",quote=F,row.names=F)
diffUp <- all[with(allDiff, (logFC>logFCfilter & adj.P.Val < fdrFilter )), ]
write.table(diffUp,file=paste0(pair,"_up.xls"),sep="\t",quote=F,row.names=F)
diffDown <- all[with(allDiff, (logFC<(-logFCfilter) & adj.P.Val < fdrFilter )), ]
write.table(diffDown,file=paste0(pair,"_down.xls"),sep="\t",quote=F,row.names=F)

#write expression level of diff gene
hmExp=data[as.vector(rownames(diffSig)),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file=paste0(pair,"_diffExp.txt"),sep="\t",quote=F,col.names=F)

#volcano
allDiff <- as.data.frame(read.table(paste0(pair,"_allDEGs.xls"),sep="\t",header=T))

pdf(file=paste0(pair,"_vol.pdf"))
xMax=max(-log10(allDiff$adj.P.Val))
yMax=max(abs(allDiff$logFC))
plot(-log10(allDiff$adj.P.Val), allDiff$logFC, xlab="-log10(adj.P.Val)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(allDiff, adj.P.Val<fdrFilter & logFC>logFCfilter)
points(-log10(diffSub$adj.P.Val), diffSub$logFC, pch=20, col="red",cex=0.8)
diffSub=subset(allDiff, adj.P.Val<fdrFilter & logFC<(-logFCfilter))
points(-log10(diffSub$adj.P.Val), diffSub$logFC, pch=20, col="green",cex=0.8)
abline(h=0,lty=2,lwd=3)
dev.off()

#pheatmap
rt=read.table(paste0(pair,"_diffExp.txt"),sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#rt=log2(rt+1)
#rt[rt>15]=15

library(pheatmap)
pheno=read.table(typeFile,sep="\t",header=T,check.names=F)
a=strsplit(pair,"-")
a <- a[[1]]
pheno <- pheno[which(pheno$type %in% a),]
rownames(pheno) <- pheno$id
ann <- as.data.frame(pheno[,2])
rownames(ann) <- rownames(pheno)
colnames(ann) <- 'type'
rt <- as.data.frame(rt[,rownames(ann)])

pdf(paste0(pair,"_pheatmap.pdf"))
pheatmap(rt, annotation=ann, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         #         scale="row",
         fontsize_row=3,
         fontsize_col=5)
dev.off()