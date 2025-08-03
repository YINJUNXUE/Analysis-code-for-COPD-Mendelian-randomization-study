####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑

library(dplyr)
library(coloc)
library(data.table)
library(foreach)
setwd("/media/desk16/iyun3811/R/gongdingwei/gongdignwei_eqtl")


genelist=read.table("Venn.txt",header = F,sep = "\t")

genelist1=as.vector(genelist[,1])

result=data.frame()


foreach(i=genelist1, .errorhandling = "pass") %do%{
  
  
  df1=fread(file = paste0("eqtlclump/",i,".txt"))
  
  head(df1)
  
  df1=data.frame(SNP=df1$SNP,chrom=df1$chr,
                 pos=df1$pos,A1=df1$A1,A2=df1$A2,beta=df1$beta,se=df1$se,MAF=df1$eaf,N=df1$n)
  
  df2=fread("GCST004744_LUAD.tsv",header = T)
  head(df2)
  df2 <- data.frame(
    SNP = df2$variant_id,                    
    chrom = df2$chromosome,              
    pos = df2$base_pair_location,                     
    A1 = df2$effect_allele,                     
    A2 = df2$other_allele,                      
    beta = df2$beta,                   
    se = df2$standard_error                    
  )
  
  
  
  
  dfall=merge(df1,df2,by="SNP")
  dfall = dfall[!duplicated(dfall$SNP),]
  
  
  dfall = dfall %>% filter((A1.x==A1.y&A2.x==A2.y)|(A1.x==A2.y&A2.x==A1.y)) 
  dfall = dfall %>% mutate(beta.y = ifelse(A1.x==A1.y,beta.y,-beta.y))
  
  
  dfall$VAR1 = dfall$se.x^2
  dfall$VAR2 =dfall$se.y^2
  dfall = dfall[dfall$VAR1!=0 & dfall$VAR2!=0 ,]

  cdf1 =data.frame(beta=dfall$beta.x,varbeta=dfall$VAR1,snp=dfall$SNP,MAF=dfall$MAF,N=dfall$N)
  cdf2= data.frame(beta=dfall$beta.y,varbeta=dfall$VAR2,snp=dfall$SNP)
  
  cdf1=na.omit(cdf1)
  cdf2=na.omit(cdf2)
  
  cdf1 = as.list(cdf1)
  cdf2 = as.list(cdf2)
 
  cdf1$type = "quant"#quant
  cdf2$type = "cc"
  
  
  res = coloc.abf(cdf1,cdf2,p1=1e-4,p2=1e-4,p12=1e-5)
  color_result1=res$summary
  
  result=rbind(result,cbind(id=i,nsnps=color_result1[[1]],
                            PP.H0.abf=color_result1[[2]],
                            PP.H1.abf=color_result1[[3]],
                            PP.H2.abf=color_result1[[4]],
                            PP.H3.abf=color_result1[[5]],
                            PP.H4.abf=color_result1[[6]]))
  
}

write.table(result,"cloc_res_eqtl.txt",sep = "\t",row.names = F,quote = F)
