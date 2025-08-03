####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com

library(data.table)
library(foreach)

setwd("/media/desk16/iyun3811/R/jiajihua/step1")

all_dfd=fread("GoDMCmQTL.txt")
head(all_dfd)
all_dfd$cpg


cglist=read.table("SDK1cg.txt",header = F,sep = "\t")

cglist1=as.vector(cglist[,1])

foreach(i=cglist1, .errorhandling = "pass") %do%{
  
  mycg=all_dfd[all_dfd$cpg==i,]
  
  fwrite(mycg,file = paste0(i,".txt"),sep = "\t",quote = F)
}