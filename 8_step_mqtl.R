

#install.packages("ggplot2")
rm(list = ls())
#library
library(TwoSampleMR)
library(ggplot2)
library(foreach)
setwd("/media/desk16/iyun3811/R/jiajihua/step4")

genelist=read.table("mqtl_result.txt",header = T,sep = "\t")

genelist1=as.vector(genelist$id)

result=data.frame()
foreach(i=genelist1, .errorhandling = "pass") %do%{
  expo_rt<- read_exposure_data(
    filename = paste0("clump/",i,"_clump.txt"),
    sep = "\t",
    snp_col = "rsids2",
    beta_col = "beta_a1",
    se_col = "se",
    effect_allele_col = "allele1",
    other_allele_col = "allele2",
    eaf_col = "freq_a1",
    pval_col = "pval",
    samplesize_col = "samplesize")
  
  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = "CST7cg.txt.gz",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col ="A2",
    eaf_col = "eaf",
    pval_col = "p")
  
  harm_rt <- harmonise_data(
    exposure_dat =  expo_rt, 
    outcome_dat = outc_rt,action=2)
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure*harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf<- mean(harm_rt$f)
  harm_rt<-harm_rt[harm_rt$f>10,]
  
  mr_result<- mr(harm_rt)
  result_or=generate_odds_ratios(mr_result) 
  filename=paste0("mqtl_gene/",i)
  dir.create(filename)
  write.table(harm_rt, file =paste0(filename,"/harmonise.txt"),row.names = F,sep = "\t",quote = F)
  write.table(result_or[,5:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
  pleiotropy=mr_pleiotropy_test(harm_rt)
  write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
  heterogeneity=mr_heterogeneity(harm_rt)
  write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)
  
  result=rbind(result,cbind(id=i,nsnp=result_or$nsnp[3],beta=result_or$b[3],
                            se=result_or$se[3],
                            pvalue=result_or$pval[3],
                            or=result_or$or[3],
                            or_lci95=result_or$or_lci95[3],
                            or_uci95=result_or$or_uci95[3],
                            FDR=result_or$FDR[3],
                            Heter.Stat=heterogeneity$Q[2],Heter.p=heterogeneity$Q_pval[2],
                            egger_intercept=pleiotropy$egger_intercept,
                            egger_intercept_se=pleiotropy$se,
                            egger_intercept_pval=pleiotropy$pval
  ))

  p1 <- mr_scatter_plot(mr_result, harm_rt)
  ggsave(p1[[1]], file=paste0(filename,"/scatter.pdf"), width=8, height=8)
  
  
  singlesnp_res<- mr_singlesnp(harm_rt)
  singlesnpOR=generate_odds_ratios(singlesnp_res)
  write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
  p2 <- mr_forest_plot(singlesnp_res)
  ggsave(p2[[1]], file=paste0(filename,"/forest.pdf"), width=8, height=8)
  sen_res<- mr_leaveoneout(harm_rt)
  p3 <- mr_leaveoneout_plot(sen_res)
  ggsave(p3[[1]], file=paste0(filename,"/sensitivity-analysis.pdf"), width=8, height=8)
  res_single <- mr_singlesnp(harm_rt)
  p4 <- mr_funnel_plot(singlesnp_res)
  ggsave(p4[[1]], file=paste0(filename,"/funnelplot.pdf"), width=8, height=8)
  presso=run_mr_presso(harm_rt,NbDistribution = 1000)
  capture.output(presso,file = paste0(filename,"/presso.txt"))
}

write.table(result,"mqtl_gene.txt",sep = "\t",quote = F,row.names = F)

