####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com

##操作：把step0 step1 的文件全部移动过来，更改setwd目录，以及cglist基因文件名。

setwd("/media/desk16/iyun3811/R/jiajihua/step2")

library(data.table)
library(foreach)
source("ld_clump.R")
source("ld_matrix.R")
source("afl2.r")
source("api.R")
source("backwards.R")
source("query.R")
source("utils-pipe.R")
source("variants.R")
source("zzz.R")


cglist=read.table("SDK1cg.txt",header = F,sep = "\t")

cglist1=as.vector(cglist[,1])

foreach(i=cglist1, .errorhandling = "pass") %do%{
  expo_rt=fread(file = paste0(i,".txt"), header = T)
  
  
  #数据与代码声明
  #如果没有购买SCI狂人团队或者生信狂人团队的正版会员
  #没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
  #如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿
  
  ####关注微信公众号生信狂人团队
  ###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
  ###作者邮箱：sxkrteam@shengxinkuangren.com
  
  expo_rt=expo_rt[expo_rt$pval<5e-8,]#5e-8
  
  expo_rt2=expo_rt[,c("rsids2","pval")]
  colnames(expo_rt2)=c("rsid", "pval")
  expo_rt2=expo_rt2[!expo_rt2$rsid=="",]
  
  #数据与代码声明
  #如果没有购买SCI狂人团队或者生信狂人团队的正版会员
  #没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
  #如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿
  
  ####关注微信公众号生信狂人团队
  ###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
  ###作者邮箱：sxkrteam@shengxinkuangren.com
  
  clumdf <- ld_clump_local(dat = expo_rt2, clump_kb = 10000, clump_r2 = 0.1,clump_p=1,
                           bfile ="/media/desk16/iyun3811/R/jiajihua/step2/data_maf0/data_maf0.01_rs_ref", 
                           plink_bin = "/media/desk16/iyun3811/R/jiajihua/step2/plink_linux_x86_64_20240804/plink")
  
  expo_rt3=expo_rt[which(expo_rt$rsids2%in%clumdf$rsid),]
  
  write.table(expo_rt3,file = paste0("clump/",i,"_clump.txt"),row.names = F,sep = "\t",quote = F)
}