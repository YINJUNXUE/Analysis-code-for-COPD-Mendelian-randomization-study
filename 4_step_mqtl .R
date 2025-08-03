
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#安装包

#BiocManager::install("GEOquery")

#BiocManager::install("minfi")

#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

###操作：改工作目录  改基因名

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

setwd("/media/desk16/iyun3811/R/jiajihua/step0")

mygene="IRF1"#改为自己的基因

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


methylation <- subset(ann450k, UCSC_RefGene_Name == mygene)


gc=as.data.frame(methylation$Name)
write.table(gc,file = paste0(mygene,"cg.txt"),sep ="\t",quote = F,col.names = F,row.names = F)


