################################
#######  UpSetR  ################
################################
#install.packages("UpSetR")

library(UpSetR)
outFile="upset_intersectGenes.txt"        #输出交集基因文件
outPic="upset.pdf"                  #输出图片

files=dir()                          #获取目录下所有文件
files=grep("_top30.csv$",files,value=T)     #提取.txt结尾的文件
geneList=list()

for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile==outFile){next}
  rt=read.table(inputFile,header=F, sep=",",skip = 2)         #读取输入文件
  geneNames=as.vector(rt[,2])               #提取基因名称
  geneNames=gsub("^ | $","",geneNames)      #去掉基因首尾的空格
  uniqGene=unique(geneNames)                #基因取unique，唯一基因列表
  header=unlist(strsplit(inputFile,"\\.|\\_"))
  geneList[[header[1]]]=uniqGene
  uniqLength=length(uniqGene)
  print(paste(header[1],uniqLength,sep=" "))
}

#绘制UpSet
upsetData=fromList(geneList)
pdf(file=outPic,onefile = FALSE,width=7,height=4)
upset(upsetData,
      nsets = length(geneList),               #展示多少个数据
      nintersects = 50,                       #展示基因集数目
      order.by = "freq",                      #按照数目排序
      show.numbers = "yes",                   #柱状图上方是否显示数值
      number.angles = 0,                     #字体角度
      point.size = 1.2,                         #点的大小
      mb.ratio = c(0.5,0.5),
      matrix.color="black",                     #交集点颜色
      line.size = 0.6,                        #线条粗细
      mainbar.y.label = "Gene Intersections", 
      sets.x.label = "Set Size")
dev.off()

#保存交集文件
intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)