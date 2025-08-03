
setwd("/media/desk16/iyun3811/R/jiajihua/step5")

data1=read.table("cg_outcome.txt",header = T,sep = "\t")
data2=read.table("cg_gene.txt",header = T,sep = "\t")
data3=read.table("gene_outcome.txt",header = T,sep = "\t")
md2=merge(data1,data2,by="id")
#####mediation effect

md2$beta_all=md2$beta.x


md2$beta1=md2$beta.y



md2$beta2=data3$beta

#中介效应（间接效应）
md2$beta12=md2$beta1*md2$beta2

#直接效应
md2$beta_dir=md2$beta_all-md2$beta12

md2$se=sqrt(md2$beta.y^2*md2$se.y^2+data3$beta^2*data3$se^2)


#Z
md2$Z=md2$beta12/md2$se

#p
md2$P=2*pnorm(q=abs(md2$Z), lower.tail=FALSE)

#####95%可信区间
md2$lci=md2$beta12-1.96*md2$se

md2$uci=md2$beta12+1.96*md2$se

#中介效应所占的比例
md2$beta12_p=md2$beta12/md2$beta_all
md2$lci_p=md2$lci/md2$beta_all
md2$uci_p=md2$uci/md2$beta_all

write.table(md2,"mediation effect.txt",sep = "\t",quote = F,row.names = F)

