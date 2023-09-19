################################
################################
#######  ROC    ################
################################
library(pROC) 

getGenesExp <- function(inFile,typeFile,geneFile){
  rt=read.table(inFile,sep="\t",header=T,check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>0.5,]
  data <- t(data)
  
  type=read.table(typeFile,sep="\t",header=T,check.names=F, row.names = 1)
  type=as.data.frame(type[which(rownames(type)==rownames(data)),]$type)
  colnames(type) <- 'Type'
  data <- cbind(data,type)
  
  genes=read.table(geneFile,header=T,sep="\t",check.names=F)
  genename=t(genes[,1])
  Exp  <- data[,c(genename,"Type")]
  rownames <- rownames(data)
  
  prefix=unlist(strsplit(geneFile,"\\.|\\_"))[1]
  write.table(cbind(rownames,Exp),file=paste0(prefix,"_Exp.txt"),sep="\t",row.names=F,quote=F)}

getGenesExp("normExp.txt","sepsis/65682/pro.txt", "target")


run_proc_for_specified_query <- function(query,db,output){
  qt=read.table(query,header=T,sep="\t",check.names=F)
  genename <- t(qt[,1])
  rt=read.table(db,header=T,sep="\t",check.names=F,row.names=1)
  y=colnames(rt)[ncol(rt)]
  for (i in genename){
    outFile=paste0(i,'_roc.pdf',sep='')
    rocobj1=roc(rt[,y], as.numeric(rt[,i]), ci=TRUE, levels=c("control","patient"))
    auc=rocobj1$auc
    info=data.frame(i,auc)
    write.table(info,file=output,append=TRUE,row.names=FALSE,col.names=FALSE)
    pdf(file=outFile,width=5,height=5)
    plot(rocobj1, print.auc=TRUE, col="red")
    dev.off()}}

run_proc_for_specified_query("target_genes.txt","target_Exp.txt","auc-test")


roc1 <- roc(Type ~ CD3E, smooth = TRUE, percent = TRUE, ci=TRUE, levels=c("control","patient"),rt)
colnames(rt)
colnames(rt)[2] <- "HLADRA"
colnames(rt)
roc2 <- roc(Type ~ HLADRA, smooth = TRUE, percent = TRUE, ci=TRUE, levels=c("control","patient"),rt)
roc3 <- roc(Type ~ IL2RB, smooth = TRUE, percent = TRUE, ci=TRUE, levels=c("control","patient"),rt)
roc4 <- roc(Type ~ ITK, smooth = TRUE, percent = TRUE, ci=TRUE, levels=c("control","patient"),rt)
roc5 <- roc(Type ~ LAT, smooth = TRUE, percent = TRUE, ci=TRUE, levels=c("control","patient"),rt)

auc1=signif(as.numeric(roc1$auc),digits=4)
auc2=signif(as.numeric(roc2$auc),digits=4)
auc3=signif(as.numeric(roc3$auc),digits=4)
auc4=signif(as.numeric(roc4$auc),digits=4)
auc5=signif(as.numeric(roc5$auc),digits=4)

ci1=signif(as.numeric(roc1$ci),digits=4)
ci2=signif(as.numeric(roc2$ci),digits=4)
ci3=signif(as.numeric(roc3$ci),digits=4)
ci4=signif(as.numeric(roc4$ci),digits=4)
ci5=signif(as.numeric(roc5$ci),digits=4)

name1=paste('CD3E_',auc1,"(",ci1[1],"-",ci1[3],")")
name2=paste('HLA-DRA_',auc2, "(",ci2[1],"-",ci2[3],")")
name3=paste('IL2RB_',auc3,"(",ci3[1],"-",ci3[3],")")
name4=paste('ITK_',auc4, "(",ci4[1],"-",ci4[3],")")
name5=paste('LAT_',auc5, "(",ci5[1],"-",ci5[3],")")

plot(roc1)                            
lines.roc(roc2, col= "red")           
lines.roc(roc3, col= "steelblue")        
lines.roc(roc4, col= "goldenrod1")     
lines.roc(roc5, col= "darkgreen")     

legend(x = 30, y = 50,
       fill = c("black", "red", "steelblue", "goldenrod1","darkgreen"),
       #legend = c("CD3E","HLA-DRA","IL2RB","ITK","LAT"),
       legend=c(name1,name2,name3,name4,name5),
       cex = 1)