rm(list=ls())
## this is a function to test the mediation relationship between g, x, y ##
## g is genotype, x is methylation, y is gene expression ##


score_reg <- function(g,x,y) {
  xx <- x- mean(x)
  yy <- y- mean(y)
  gg<- g-mean(g)

  out <- rep(0,6)
  fit1 <- glm(xx~gg-1,x=T,y=T)
  error1 <- fit1$x*(fit1$y-fit1$fitted)

  fit2 <- glm(yy~xx+gg-1,x=T,y=T)
  error2 <- xx*(fit2$y-fit2$fitted)

  fit3 <- glm(yy~gg-1,x=T,y=T)
  error3 <- fit3$x*(fit3$y-fit3$fitted)

  fit4 <- glm(xx~yy-1,x=T,y=T) 
  error4 <- fit4$x*(fit4$y-fit4$fitted)

  out[1] <- drop(sum(error3*error2))/drop(sum(error3*error3))
  out[2] <- drop(sum(error1*error3))/drop(sum(error1*error1))
  out[3] <- drop(sum(error4*error3))/drop(sum(error3*error3))
  
  n=200
  bb1 <- rep(0,n)
  bb2 <- rep(0,n)
  bb3 <- rep(0,n)

  for (i in 1:n) {
    sid <- sample(1:n,replace=T)
    xx1 <- x[sid]
    yy1 <- y[sid]
    gg1 <- g[sid]
    xx1 <- xx1- mean(xx1)
    yy1 <- yy1- mean(yy1)
    gg1 <- gg1- mean(gg1)    
    fit1 <- glm(xx1~gg1-1,x=T,y=T)
    err1 <- fit1$x*(fit1$y-fit1$fitted)
    fit2 <- glm(yy1~xx1+gg1-1,x=T,y=T)
    err2 <- xx1*(fit2$y-fit2$fitted)
    fit3 <- glm(yy1~gg1-1,x=T,y=T)
    err3 <- fit3$x*(fit3$y-fit3$fitted)
    fit4 <- glm(xx1~yy1-1,x=T,y=T)
    err4 <- fit4$x*(fit4$y-fit4$fitted)
  
    bb1[i] <- drop(sum(err3*err2))/drop(sum(err3*err3))
    bb2[i] <- drop(sum(err1*err3))/drop(sum(err1*err1))
    bb3[i] <- drop(sum(err4*err3))/drop(sum(err3*err3))
  }


  out[4] <- var(bb1)
  out[5] <- var(bb2)
  out[6] <- var(bb3)

  p1 <- 2*(1-pnorm(abs(out[1]/sqrt(out[4]))))
  p2 <- 2*(1-pnorm(abs(out[2]/sqrt(out[5]))))
  p3 <- 2*(1-pnorm(abs(out[3]/sqrt(out[6]))))
  
  outp <- max(p1,p2,p3)
  return(outp)
}



### below is the simulation example and compare to CIT test ###

library(cit)
n <- 200
nsim <- 100
out <- matrix(0,nsim,2)
for (l in 1:nsim) {
  set.seed(47384+l)
  g1 <- rbinom(n,1,0.2) + rbinom(n,1,0.2)
  u1 <- rnorm(n,0,1)
  u2 <- rnorm(n,0,1)
  
  ## scenario 1 ##
  #x <- 0.5*g1 + 1*u1
  #y <- 0.3*g1 + 1*u2
  
  ## scenario 2 
  ##x <- 0.7*g1 + 1*u1
  ##y <- 0.3*x +  1*u2 
  
  ## scenario 3 
  x <- 1*g1 + 1*u1
  y <- 0.4*x +  1*u2 - 0.2*g1
  
  ## scenario 4 
  #y <- 0.3*g1 + 1*u1
  #x <- y + 0.5*u2
  
  ## scenario 5 
  #y <- 0.3*g1 + 1*u1
  #x <- y + 0.5*u2 +0.2*g1
  
  # scenario 6 
  #x <- 1*u1
  #y <- 0.3*g1 + 1*u2 +0.5*x
  
  # scenario 7
  #y <- 1*u1
  #x <- y + 0.5*u2 +0.2*g1
  
  # scenario 8
  #y <- 1*u1
  #x <- 0.5*u2 +0.2*g1
  

  out[l,1] <- score_reg(g1,x,y)
  
  ### now compute CIT 
  cit.fit <- cit.cp(g1,x,y)
  out[l,2] <- cit.fit[1]
  print(out[l,])
} 

#mean(out[,1]<0.05)
#mean(out[,2]<0.05)

load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData")
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/imputation_GEME.RData")
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/ge_annohg19.RData")
studyno_GE=colnames(GE)
studyno_ME=colnames(ME)
studyno_GT=colnames(highrisk)[6:ncol(highrisk)]
studyno=intersect(studyno_GE,studyno_ME)
studyno=intersect(studyno,studyno_GT)

highrisk=highrisk[complete.cases(highrisk),]
idx=match(studyno,colnames(highrisk))
highrisk1=highrisk[,idx]
highrisk1_pos=data.frame(snp=highrisk$V1,chr=highrisk$chr,pos=highrisk$V2,stringsAsFactors = F)

idx=match(studyno,colnames(ME))
ME1=ME[,idx]


GE1_pos=read.table("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/GE_POS.txt",header=T,sep="\t",row.names = 1)
GE1=read.table("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/GE1.txt",header=T,sep="\t",row.names = 1)
colnames(GE1)=gsub("X","",colnames(GE1))
sum(rownames(GE1)==rownames(GE1_pos))
idx=match(studyno,colnames(GE1))
GE1=GE1[,idx]
ME1_pos=read.table("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/ME_POS.txt",header=T,sep="\t",row.names = 1)
sum(rownames(ME1)==rownames(ME1_pos))
sum(colnames(ME1)==colnames(GE1))
sum(colnames(ME1)==colnames(highrisk1))

library(GenomicRanges)
gr_ME1=GRanges(seqnames = ME1_pos$chr,ranges=IRanges(start=ME1_pos$s1,end=ME1_pos$s2))
gr_GE1=GRanges(seqnames = GE1_pos$chr,ranges=IRanges(start=GE1_pos$s1,end=GE1_pos$s2))

#first form the probe idx table
testtable=NULL
cutoff=1e6
for (i in 1:nrow(highrisk1))
{
  if (i %%10==0) cat(i,'..')
  gprobe=i
  #find methylation probes
  gr_gt=GRanges(seqnames = highrisk1_pos$chr[i],ranges=IRanges(start=highrisk1_pos$pos[i],width = 1))
  tmp=distance(gr_ME1,gr_gt)
  idx=which(tmp<cutoff) #ME
  tmp1=distance(gr_GE1,gr_gt)
  idx1=which(tmp1<cutoff) #GE
  if (length(idx)>0 & length(idx1)>0)
  {
    numrow=length(idx)*length(idx1)
    tmp2=data.frame(gprobe=rep(gprobe,numrow),xprobe=rep(idx,each=length(idx1)),yprobe=rep(idx1,times=length(idx)))
    testtable=rbind(testtable,tmp2)
  }
}
testtable1=testtable
#testtable1 with 1M cutoff
#testtable with 500K cutoff
#save(testtable,testtable1,highrisk1,highrisk1_pos,ME1,ME1_pos,GE1,GE1_pos,covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_testtable.RData")
compute_pvalues=function(jobn)
{
  gprobe=testtable$gprobe[jobn]
  g=unlist(highrisk1[gprobe,])
  xprobe=testtable$xprobe[jobn]
  x=unlist(ME1[xprobe,])
  yprobe=testtable$yprobe[jobn]
  y=unlist(GE1[yprobe,])
  p1=score_reg(g,x,y)
  cit.fit=cit.cp(g,x,y)
  p2=cit.fit[1]
  return(data.frame(pcit=p2,p=p1))
}

#geneexp->methylation
compute_pvalues1=function(jobn)
{
  gprobe=testtable$gprobe[jobn]
  g=unlist(highrisk1[gprobe,])
  xprobe=testtable$xprobe[jobn]
  y=unlist(ME1[xprobe,])
  yprobe=testtable$yprobe[jobn]
  x=unlist(GE1[yprobe,])
  p1=score_reg(g,x,y)
  cit.fit=cit.cp(g,x,y)
  p2=cit.fit[1]
  return(data.frame(pcit=p2,p=p1))
}

mecovariate=read.table("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/COVA_ME.txt",header=T,sep="\t",row.names = 1)
colnames(mecovariate)=gsub("^X","",colnames(mecovariate))
idx=match(colnames(GE1),colnames(mecovariate))
mecovariate=mecovariate[,idx]
sum(colnames(mecovariate)==colnames(GE1))
covariate=matrix(data=NA,nrow=ncol(GE1),ncol=2)
#age
covariate[,1]=unlist(mecovariate["AGEREF",])
covariate[,2]=unlist(mecovariate["Gleason_score",])
#methylation->geneexp, adjust for covariates
compute_pvalues2=function(jobn)
{
  gprobe=testtable1$gprobe[jobn]
  g=unlist(highrisk1[gprobe,])
  xprobe=testtable1$xprobe[jobn]
  x=unlist(ME1[xprobe,])
  yprobe=testtable1$yprobe[jobn]
  y=unlist(GE1[yprobe,])
  #p1=score_reg(g,x,y)
  cit.fit=cit.cp(g,x,y,covariate)
  return(t(data.frame(cit.fit)))
}

#geneexp->methylation, adjust for covariates
compute_pvalues3=function(jobn)
{
  gprobe=testtable1$gprobe[jobn]
  g=unlist(highrisk1[gprobe,])
  xprobe=testtable1$xprobe[jobn]
  y=unlist(ME1[xprobe,])
  yprobe=testtable1$yprobe[jobn]
  x=unlist(GE1[yprobe,])
  #p1=score_reg(g,x,y)
  cit.fit=cit.cp(g,x,y,covariate)
  return(t(data.frame(cit.fit)))
}
#salloc -t 1-1 -n 100 mpirun -n 1 R --interactive
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_testtable.RData")
mpi.bcast.Robj2slave(testtable)
mpi.bcast.Robj2slave(highrisk1)
mpi.bcast.Robj2slave(ME1)
mpi.bcast.Robj2slave(GE1)
mpi.bcast.Robj2slave(score_reg)
mpi.bcast.Robj2slave(compute_pvalues)
mpi.bcast.Robj2slave(compute_pvalues1)
mpi.bcast.Robj2slave(compute_pvalues2)
mpi.bcast.Robj2slave(compute_pvalues3)
mpi.bcast.Robj2slave(testtable1)
mpi.bcast.Robj2slave(covariate)
mpi.bcast.cmd(library(cit))

result=data.frame(matrix(NA,nrow=nrow(testtable),ncol=5))
colnames(result)=c("genotype_idx","methylation_idx","geneexpression_idx","pcit","p")
result[,1:3]=testtable
num=1000
#load("/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result.RData")
nrun <- ceiling(nrow(testtable)/num)
print(paste0("total num of run: ",nrun))
for (j in 1:nrun){
  cat(j,"..")
  if (j < nrun) cseq <- ((j-1)*num+1):(j*num)  else  cseq <- ((j-1)*num+1):nrow(testtable)
  #res=mpi.parSapply(X=cseq,FUN=compute_pvalues,highrisk1=highrisk1,ME1=ME1,GE1=GE1,job.num=njobs)
  res=mpi.parSapply(X=cseq,FUN=compute_pvalues,job.num=njobs)
  res=unlist(res)
  result$pcit[cseq]=res[seq(1,length(res),2)]
  result$p[cseq]=res[seq(2,length(res),2)]
  #if (j %% 10==0) save(result,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result.RData")
}
#save(result,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result.RData")

compute_pvalues(2)

#            pcit         p
# p_cit 0.9148682 0.6354581
#            pcit         p
# p_cit 0.6616844 0.9176898
# 
#            pcit         p
# p_cit 0.7507759 0.9760668

n=sum(!is.na(result$p))


idx=which(is.na(result$p))
length(idx)
if (length(idx)>0)
{
  nrun <- ceiling(length(idx)/num)
  print(paste0("total num of run: ",nrun))
  for (j in 1:nrun){
    cat(j,"..")
    if (j < nrun) cseq <- ((j-1)*num+1):(j*num)  else  cseq <- ((j-1)*num+1):length(idx)
    #res=mpi.parSapply(X=cseq,FUN=compute_pvalues,highrisk1=highrisk1,ME1=ME1,GE1=GE1,job.num=njobs)
    res=mpi.parSapply(X=idx[cseq],FUN=compute_pvalues,job.num=njobs)
    res=unlist(res)
    result$pcit[idx[cseq]]=res[seq(1,length(res),2)]
    result$p[idx[cseq]]=res[seq(2,length(res),2)]
    #if (j %% 10==0) save(result,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result.RData")
  }
}

png("qqplot.png")
par(mfrow=c(2,1))
plot(-log((1:n)/n,base=10),-log(result$pcit[order(result$pcit)],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)",main="CIT")
abline(0,1)
plot(-log((1:n)/n,base=10),-log(result$p[order(result$p)],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)")
abline(0,1)
dev.off()
min(result$pcit,na.rm=T)
min(result$p,na.rm=T)
which.min(result$pcit)
result[841661,]
which.min(result$p)
pcutoff=0.05/nrow(result)
sum(result$pcit<=pcutoff,na.rm=T)
sum(result$p<=pcutoff,na.rm=T)

load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/ge_annohg19.RData")
#anno
#load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/GE_METHY_all.RData")
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
find_triplet=function(gt_idx=45,me_idx=50756,ge_idx=9899)
{
  gtname=highrisk1_pos$snp[gt_idx]
  
  idx=which(anno$IlmnID==rownames(ME1_pos)[me_idx])
  mename=as.character(anno$UCSC_RefGene_Name[idx])
  meprobe=as.character(anno$IlmnID[idx])
  idx=which(ge_annohg19$Probe_Id==rownames(GE1_pos)[ge_idx])
  # gename=ge_anno$Symbol[idx]
  # geprobe=ge_anno$Probe_Id[idx]
  gename=ge_annohg19$Symbol[idx]
  geprobe=ge_annohg19$Probe_Id[idx]
  #return(list(gtname=gtname,mename=mename,meprobe=meprobe,gename=gename,geprobe=geprobe))
  return(paste0(gtname,"-",meprobe,'(',mename,')',"-",geprobe,'(',gename,')'))
}
idxorder=order(result$p)
head(result[idxorder,])
find_triplet(gt_idx = result$genotype_idx[idxorder][1],me_idx=result$methylation_idx[idxorder][1],ge_idx=result$geneexpression_idx[idxorder][1])
#"rs3129859:32400939:G:C-cg22627029(HLA-DRB6)-ILMN_1715169(HLA-DRB1)"
find_triplet(gt_idx = result$genotype_idx[idxorder][2],me_idx=result$methylation_idx[idxorder][2],ge_idx=result$geneexpression_idx[idxorder][2])
# "rs3096702-cg08265274(HLA-DRB5)-ILMN_1715169(HLA-DRB1)"

idxorder=order(result$pcit)
head(result[idxorder,])
for (i in 1:6)
{
  print(i)
  print(find_triplet(gt_idx = result$genotype_idx[idxorder][i],me_idx=result$methylation_idx[idxorder][i],ge_idx=result$geneexpression_idx[idxorder][i]))
}
# [1] 1
# [1] "rs3129859:32400939:G:C-cg15011943(HLA-DRB5)-ILMN_2159694(HLA-DRB4)"
# [1] 2
# [1] "rs3129859:32400939:G:C-cg24147543(HLA-DRB1)-ILMN_2159694(HLA-DRB4)"
# [1] 3
# [1] "rs3096702-cg22627029(HLA-DRB6)-ILMN_2159694(HLA-DRB4)"
# [1] 4
# [1] "rs3096702-cg15708909(HLA-DRB5)-ILMN_2159694(HLA-DRB4)"
# [1] 5
# [1] "rs3096702-cg24147543(HLA-DRB1)-ILMN_2159694(HLA-DRB4)"
# [1] 6
# [1] "rs3096702-cg05341252(HLA-DQB1)-ILMN_1697499(HLA-DRB5)"

#geneexp->methylation
result1=data.frame(matrix(NA,nrow=nrow(testtable),ncol=5))
colnames(result1)=c("genotype_idx","methylation_idx","geneexpression_idx","pcit","p")
result1[,1:3]=testtable
num=1000
nrun <- ceiling(nrow(testtable)/num)
print(paste0("total num of run: ",nrun))
for (j in 1:nrun){
  cat(j,"..")
  if (j < nrun) cseq <- ((j-1)*num+1):(j*num)  else  cseq <- ((j-1)*num+1):nrow(testtable)
  #res=mpi.parSapply(X=cseq,FUN=compute_pvalues,highrisk1=highrisk1,ME1=ME1,GE1=GE1,job.num=njobs)
  res=mpi.parSapply(X=cseq,FUN=compute_pvalues1,job.num=njobs)
  res=unlist(res)
  result1$pcit[cseq]=res[seq(1,length(res),2)]
  result1$p[cseq]=res[seq(2,length(res),2)]
  if (j %% 50==0) save(result1,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result1.RData")
}
#save(result1,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result1.RData")
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result1.RData")
n=nrow(result1)
par(mfrow=c(2,1))
plot(-log((1:n)/n,base=10),-log(result1$pcit[order(result1$pcit)],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)",main="CIT")
abline(0,1)
plot(-log((1:n)/n,base=10),-log(result1$p[order(result1$p)],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)")
abline(0,1)
pcutoff=0.05/nrow(result1)
sum(result1$pcit<=pcutoff,na.rm=T)
#[1] 16
sum(result1$p<=pcutoff,na.rm=T)
#[1] 8
idxorder=order(result1$p)
head(result1[idxorder,])
find_triplet(gt_idx = result1$genotype_idx[idxorder][1],me_idx=result1$methylation_idx[idxorder][1],ge_idx=result1$geneexpression_idx[idxorder][1])
#"rs3129859:32400939:G:C-cg22627029(HLA-DRB6)-ILMN_1715169(HLA-DRB1)"
find_triplet(gt_idx = result1$genotype_idx[idxorder][2],me_idx=result1$methylation_idx[idxorder][2],ge_idx=result1$geneexpression_idx[idxorder][2])
# "rs3096702-cg08265274(HLA-DRB5)-ILMN_1715169(HLA-DRB1)"

idxorder=order(result1$pcit)
head(result1[idxorder,])
for (i in 1:16)
{
  print(i)
  print(find_triplet(gt_idx = result1$genotype_idx[idxorder][i],me_idx=result1$methylation_idx[idxorder][i],ge_idx=result1$geneexpression_idx[idxorder][i]))
}

result2=data.frame(matrix(NA,nrow=nrow(testtable1),ncol=8))
colnames(result2)=c("genotype_idx","methylation_idx","geneexpression_idx","p1","p2","p3","p4","p5")
result2[,1:3]=testtable1
num=1000
nrun <- ceiling(nrow(testtable1)/num)
print(paste0("total num of run: ",nrun))
for (j in 1:nrun){
  if (j %%100 ==0) cat(j,"..")
  if (j < nrun) cseq <- ((j-1)*num+1):(j*num)  else  cseq <- ((j-1)*num+1):nrow(testtable1)
  #res=mpi.parSapply(X=cseq,FUN=compute_pvalues,highrisk1=highrisk1,ME1=ME1,GE1=GE1,job.num=njobs)
  res=mpi.parSapply(X=cseq,FUN=compute_pvalues2,job.num=njobs)
  result2[cseq,4:8]=t(res)
  #if (j %% 50==0) save(result1,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result1.RData")
}
save(result2,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result2.RData")

result3=data.frame(matrix(NA,nrow=nrow(testtable1),ncol=8))
colnames(result3)=c("genotype_idx","methylation_idx","geneexpression_idx","p1","p2","p3","p4","p5")
result3[,1:3]=testtable1
num=1000
nrun <- ceiling(nrow(testtable1)/num)
print(paste0("total num of run: ",nrun))
for (j in 1:nrun){
  if (j %%100 ==0) cat(j,"..")
  if (j < nrun) cseq <- ((j-1)*num+1):(j*num)  else  cseq <- ((j-1)*num+1):nrow(testtable1)
  #res=mpi.parSapply(X=cseq,FUN=compute_pvalues,highrisk1=highrisk1,ME1=ME1,GE1=GE1,job.num=njobs)
  res=mpi.parSapply(X=cseq,FUN=compute_pvalues3,job.num=njobs)
  result3[cseq,4:8]=t(res)
  if (j %% 500==0) save(result3,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result3.RData")
}
save(result3,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result3.RData")
library(qvalue)
qvalue3=qvalue(result3$p1)
sum(qvalue3$qvalues<0.1,na.rm=T)
#[1] 21
qvalue2=qvalue(result2$p1)
sum(qvalue2$qvalues<0.1,na.rm=T)
#[1] 13
idx1=which(qvalue2$qvalues<0.1)
View(result2[idx1,])
#check result vs result2
idx=which(result$pcit<=0.05/nrow(result))
for (i in 1:length(idx))
{
  print(i)
  print(result[idx[i],])
  idx1=which(result2$genotype_idx==result$genotype_idx[idx[i]] & result2$methylation_idx==result$methylation_idx[idx[i]] & result2$geneexpression_idx==result$geneexpression_idx[idx[i]])
  print(result2[idx1,])
}
sum(result2$p1<=0.05/nrow(result),na.rm = T)

#Nov22---------------------------
#only use those found in e/mQTLs
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/highrisk_mqtlresult.RData")

#put dat1 as qtl of T, dat2 as qtl of G
citqtlpairs=function(dat1=highrisk_cis_hutch_eqtl,dat2=highrisk_cis_hutch_mqtl,ME=hutch_me,GE=hutch_ge,
                     snp_ge=highrisk_hutch_snp_ge,
                     covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
{
  comsamples=intersect(colnames(GE),colnames(ME))
  idx=match(comsamples,colnames(GE))
  GE2=GE[,idx]
  idx=match(comsamples,colnames(ME))
  ME2=ME[,idx]
  idx=match(comsamples,colnames(snp_ge))
  snp2=snp_ge[,idx]
  idx=match(comsamples,colnames(covariate_ge))
  covariate2=t(covariate_ge[,idx])
  idx=match(comsamples,colnames(covariate_me))
  covariate2=cbind.data.frame(covariate2,t(covariate_me[3:nrow(covariate_me),idx]))
  
  tmp=merge(dat1,dat2,by="snp_idx")
  tmp$gene.x=as.character(tmp$gene.x)
  tmp$gene.y=as.character(tmp$gene.y)
  tmp$gene.x=as.character(tmp$gene.x)
  tmp$gene.y=as.character(tmp$gene.y)
  tmp$genename.x=as.character(tmp$genename.x)
  tmp$genename.y=as.character(tmp$genename.y)
  tmp=tmp[tmp$gene.x!=tmp$gene.y,]
  tmp=tmp[sum(is.na(snp2[tmp$snp_idx,]))<10,] #not many NA in the snp (highrisk117)
  if (nrow(tmp)>0)
  {
    citresults=vector('list', nrow(tmp))
    for (i in 1:nrow(tmp))
    {
      if (i %% 100==0) cat(i,"..")
      #L
      g=unlist(snp2[tmp$snp_idx[i],])
      if (grepl("^cg",tmp$gene.y[1])) #trait is methylation
      {
        #G
        idx=match(tmp$gene.y[i],rownames(ME2))
        x=unlist(ME2[idx,])
        #T
        idx=match(tmp$gene.x[i],rownames(GE2))
        y=unlist(GE2[idx,])
      }else #reactive
      {
        #G
        idx=match(tmp$gene.y[i],rownames(GE2))
        x=unlist(GE2[idx,])
        #T
        idx=match(tmp$gene.x[i],rownames(ME2))
        y=unlist(ME2[idx,])
        # fit1=glm(x~y+g)
        # summary(fit1)
        # fit2=glm(y~x+g)
        # summary(fit2)
        
      }
      citresults[[ i ]]=cit.cp(g,x,y,covariate2,n.perm=20)
    }
    fdrresults=fdr.cit(citresults)
    fdrresults=cbind.data.frame(tmp,fdrresults)
  }else
  {
    fdrresults=NA
  }
  
  return(fdrresults)
}

#based on Bonferroni qtl results:
#T_G
hutch_med_ciseqtl_cismqtl=citqtlpairs()
#reactive
hutch_med_cismqtl_ciseqtl=citqtlpairs(dat1=highrisk_cis_hutch_mqtl,dat2=highrisk_cis_hutch_eqtl,ME=hutch_me,GE=hutch_ge,
                                               snp_ge=highrisk_hutch_snp_ge,
                                               covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_ciseqtl_transmqtl=citqtlpairs(dat1=highrisk_cis_hutch_eqtl,dat2=highrisk_trans_hutch_mqtl,ME=hutch_me,GE=hutch_ge,
                                        snp_ge=highrisk_hutch_snp_ge,
                                        covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_transmqtl_ciseqtl=citqtlpairs(dat1=highrisk_trans_hutch_mqtl,dat2=highrisk_cis_hutch_eqtl,ME=hutch_me,GE=hutch_ge,
                                        snp_ge=highrisk_hutch_snp_ge,
                                        covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
#NAs
hutch_med_transeqtl_cismqtl=citqtlpairs(dat1=highrisk_trans_hutch_eqtl,dat2=highrisk_cis_hutch_mqtl,ME=hutch_me,GE=hutch_ge,
                                      snp_ge=highrisk_hutch_snp_ge,
                                      covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_cismqtl_transeqtl=citqtlpairs(dat1=highrisk_cis_hutch_mqtl,dat2=highrisk_trans_hutch_eqtl,ME=hutch_me,GE=hutch_ge,
                                        snp_ge=highrisk_hutch_snp_ge,
                                        covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_transeqtl_transmqtl=citqtlpairs(dat1=highrisk_trans_hutch_eqtl,dat2=highrisk_trans_hutch_mqtl,ME=hutch_me,GE=hutch_ge,
                                        snp_ge=highrisk_hutch_snp_ge,
                                        covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_transmqtl_transeqtl=citqtlpairs(dat1=highrisk_trans_hutch_mqtl,dat2=highrisk_trans_hutch_eqtl,ME=hutch_me,GE=hutch_ge,
                                        snp_ge=highrisk_hutch_snp_ge,
                                        covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
save(hutch_med_ciseqtl_cismqtl,hutch_med_cismqtl_ciseqtl,hutch_med_ciseqtl_transmqtl,hutch_med_transmqtl_ciseqtl,
     hutch_med_transeqtl_cismqtl,hutch_med_cismqtl_transeqtl,hutch_med_transeqtl_transmqtl,hutch_med_transmqtl_transeqtl,
     file="/fh/fast/stanford_j/Xiaoyu/QTL/result/hutch_medresult.RData")

find_meds=function(med1=hutch_med_ciseqtl_cismqtl,med2=hutch_med_cismqtl_ciseqtl,qcutoff=0.05)
{
  meds=NULL
  for (i in 1:nrow(med1))
  {
    #find the other direction
    idx=which(med2$snp_idx==med1$snp_idx[i] & med2$gene.y==med1$gene.x[i] & med2$gene.x==med1$gene.y[i])
    q1=med1$q.cit[i]
    q2=med2$q.cit[idx]
    if (length(q1)*length(q2)>1) print(i)
    if (q1<qcutoff | q2<qcutoff)
    {
      if (q1<qcutoff & q2<qcutoff) #both directions met qcutoff, solve conflict
      {
        #if (med1$p_LindTgvnG[i]<med2$p_LindTgvnG[idx])
        if (med1$q.cit[i]<med2$q.cit[idx])
        {
          meds=rbind.data.frame(meds,med1[i,])
        }else
        {
          meds=rbind.data.frame(meds,med2[idx,])
        }
      }else
      {
        if (q1<qcutoff)
        {
          meds=rbind.data.frame(meds,med1[i,])
        }else
        {
          meds=rbind.data.frame(meds,med2[idx,])
        }
      }
      # if (q1<=qcutoff)
      # {
      #   meds=rbind.data.frame(meds,med1[i,])
      # }
      # if (q2<=qcutoff)
      # {
      #   meds=rbind.data.frame(meds,med2[idx,])
      # }
    }
  }
  return(meds)
}
hutch_meds_ciseqtl_cismqtl=find_meds()
hutch_meds_ciseqtl_transmqtl=find_meds(med1=hutch_med_ciseqtl_transmqtl,med2=hutch_med_transmqtl_ciseqtl,qcutoff=0.05)

validate_meds_innormal=function(meds=hutch_meds_ciseqtl_cismqtl,ME=normal_me,GE=normal_ge,snp_ge=highrisk_normal_snp_ge,
                                covariate_ge=normal_ge_covariate,covariate_me=normal_me_covariate)
{
  comsamples=intersect(colnames(GE),colnames(ME))
  idx=match(comsamples,colnames(GE))
  GE2=GE[,idx]
  idx=match(comsamples,colnames(ME))
  ME2=ME[,idx]
  idx=match(comsamples,colnames(snp_ge))
  snp2=snp_ge[,idx]
  idx=match(comsamples,colnames(covariate_ge))
  covariate2=t(covariate_ge[,idx])
  idx=match(comsamples,colnames(covariate_me))
  covariate2=cbind.data.frame(covariate2,t(covariate_me[3:nrow(covariate_me),idx]))
  
  #check howmany can be verified in normals due to incompleteness of probes
  idx=which(meds$gene.x %in% c(rownames(normal_ge),rownames(normal_me)) & meds$gene.y %in% c(rownames(normal_ge),rownames(normal_me)))
  validres=NULL
  for (i in 1:length(idx))
  {
    g=unlist(snp2[meds$snp_idx[idx[i]],])
    if (grepl("^cg",meds$gene.x[idx[i]])) #T is methylation
    {
      idx1=match(meds$gene.x[idx[i]],rownames(ME2))
      y=unlist(ME2[idx1,])
      idx1=match(meds$gene.y[idx[i]],rownames(GE2))
      x=unlist(GE2[idx1,])
    }else
    {
      idx1=match(meds$gene.x[idx[i]],rownames(GE2))
      y=unlist(GE2[idx1,])
      idx1=match(meds$gene.y[idx[i]],rownames(ME2))
      x=unlist(ME2[idx1,])
    }
    validres=rbind(validres,cit.cp(g,x,y,covariate2))
  }
  meds1=cbind(meds,data.frame(matrix(ncol=ncol(validres),nrow=nrow(meds))))
  colnames(meds1)[(ncol(meds)+1):ncol(meds1)]=paste0("valid_",colnames(validres))
  meds1[idx,(ncol(meds)+1):ncol(meds1)]=validres
  return(meds1)
}

normal_meds_ciseqtl_cismqtl=validate_meds_innormal(meds=hutch_meds_ciseqtl_cismqtl,ME=normal_me,GE=normal_ge,snp_ge=highrisk_normal_snp_ge,
                                                            covariate_ge=normal_ge_covariate,covariate_me=normal_me_covariate)
normal_meds_ciseqtl_transmqtl=validate_meds_innormal(meds=hutch_meds_ciseqtl_transmqtl,ME=normal_me,GE=normal_ge,snp_ge=highrisk_normal_snp_ge,
                                                 covariate_ge=normal_ge_covariate,covariate_me=normal_me_covariate)
#based on fdr qtl results:
hutch_med_ciseqtl_cismqtl_fdr=citqtlpairs(dat1=highrisk_cis_hutch_eqtl_fdr,dat2=highrisk_cis_hutch_mqtl_fdr,ME=hutch_me,GE=hutch_ge,
                                      snp_ge=highrisk_hutch_snp_ge,
                                      covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
#reactive
hutch_med_cismqtl_ciseqtl_fdr=citqtlpairs(dat1=highrisk_cis_hutch_mqtl_fdr,dat2=highrisk_cis_hutch_eqtl_fdr,ME=hutch_me,GE=hutch_ge,
                                      snp_ge=highrisk_hutch_snp_ge,
                                      covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_ciseqtl_transmqtl_fdr=citqtlpairs(dat1=highrisk_cis_hutch_eqtl_fdr,dat2=highrisk_trans_hutch_mqtl_fdr,ME=hutch_me,GE=hutch_ge,
                                        snp_ge=highrisk_hutch_snp_ge,
                                        covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_transmqtl_ciseqtl_fdr=citqtlpairs(dat1=highrisk_trans_hutch_mqtl_fdr,dat2=highrisk_cis_hutch_eqtl_fdr,ME=hutch_me,GE=hutch_ge,
                                        snp_ge=highrisk_hutch_snp_ge,
                                        covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
#NAs
hutch_med_transeqtl_cismqtl_fdr=citqtlpairs(dat1=highrisk_trans_hutch_eqtl_fdr,dat2=highrisk_cis_hutch_mqtl_fdr,ME=hutch_me,GE=hutch_ge,
                                        snp_ge=highrisk_hutch_snp_ge,
                                        covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_cismqtl_transeqtl_fdr=citqtlpairs(dat1=highrisk_cis_hutch_mqtl_fdr,dat2=highrisk_trans_hutch_eqtl_fdr,ME=hutch_me,GE=hutch_ge,
                                        snp_ge=highrisk_hutch_snp_ge,
                                        covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_transeqtl_transmqtl_fdr=citqtlpairs(dat1=highrisk_trans_hutch_eqtl_fdr,dat2=highrisk_trans_hutch_mqtl_fdr,ME=hutch_me,GE=hutch_ge,
                                          snp_ge=highrisk_hutch_snp_ge,
                                          covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_med_transmqtl_transeqtl_fdr=citqtlpairs(dat1=highrisk_trans_hutch_mqtl_fdr,dat2=highrisk_trans_hutch_eqtl_fdr,ME=hutch_me,GE=hutch_ge,
                                          snp_ge=highrisk_hutch_snp_ge,
                                          covariate_ge=hutch_ge_covariate,covariate_me=hutch_me_covariate)
hutch_meds_ciseqtl_cismqtl_fdr=find_meds(med1=hutch_med_ciseqtl_cismqtl_fdr,med2=hutch_med_cismqtl_ciseqtl_fdr,qcutoff=0.05)
hutch_meds_ciseqtl_transmqtl_fdr=find_meds(med1=hutch_med_ciseqtl_transmqtl_fdr,med2=hutch_med_transmqtl_ciseqtl_fdr,qcutoff=0.05)
normal_meds_ciseqtl_cismqtl_fdr=validate_meds_innormal(meds=hutch_meds_ciseqtl_cismqtl_fdr,ME=normal_me,GE=normal_ge,snp_ge=highrisk_normal_snp_ge,
                                                   covariate_ge=normal_ge_covariate,covariate_me=normal_me_covariate)
normal_meds_ciseqtl_transmqtl_fdr=validate_meds_innormal(meds=hutch_meds_ciseqtl_transmqtl_fdr,ME=normal_me,GE=normal_ge,snp_ge=highrisk_normal_snp_ge,
                                                     covariate_ge=normal_ge_covariate,covariate_me=normal_me_covariate)


# save(hutch_med_ciseqtl_cismqtl,hutch_med_cismqtl_ciseqtl,hutch_med_ciseqtl_transmqtl,hutch_med_transmqtl_ciseqtl,
#      hutch_med_transeqtl_cismqtl,hutch_med_cismqtl_transeqtl,hutch_med_transeqtl_transmqtl,hutch_med_transmqtl_transeqtl,
#      hutch_meds_ciseqtl_cismqtl,hutch_meds_ciseqtl_transmqtl,normal_meds_ciseqtl_cismqtl,normal_meds_ciseqtl_transmqtl,
#      hutch_med_ciseqtl_cismqtl_fdr,hutch_med_cismqtl_ciseqtl_fdr,hutch_med_ciseqtl_transmqtl_fdr,hutch_med_transmqtl_ciseqtl_fdr,
#      hutch_med_transeqtl_cismqtl_fdr,hutch_med_cismqtl_transeqtl_fdr,hutch_med_transeqtl_transmqtl_fdr,hutch_med_transmqtl_transeqtl_fdr,
#      hutch_meds_ciseqtl_cismqtl_fdr,hutch_meds_ciseqtl_transmqtl_fdr,normal_meds_ciseqtl_cismqtl_fdr,normal_meds_ciseqtl_transmqtl_fdr,
#      file="/fh/fast/stanford_j/Xiaoyu/QTL/result/hutch_medresult.RData")

#check the overlap of result2 with  hutch_med_ciseqtl_cismqtl 
result22=result2
result22$genotype_idx=highrisk1_pos$snp[result2$genotype_idx]
result22$methylation_idx=rownames(ME1)[result2$methylation_idx]
result22$geneexpression_idx=rownames(GE1)[result2$geneexpression_idx]
idx=which(qvalue2$qvalues<=0.1)
for (i in 1:length(idx))
{
  idx1=which(highrisk_snp_pos$snp[hutch_med_ciseqtl_cismqtl$snp_idx]==result22$genotype_idx[idx[i]] & hutch_med_ciseqtl_cismqtl$gene.x==result22$geneexpression_idx[idx[i]]
             & hutch_med_ciseqtl_cismqtl$gene.y==result22$methylation_idx[idx[i]])
  if (length(idx1)>0)
  {
    print(i)
    print(hutch_med_ciseqtl_cismqtl[idx1,])
  }
}

par(mfrow=c(1,1))
n=nrow(hutch_med_ciseqtl_cismqtl)
plot(-log((1:n)/n,base=10),-log(hutch_med_ciseqtl_cismqtl$p.raw[order(hutch_med_ciseqtl_cismqtl$p.raw)],base=10),xlab="expected p-value (log base 10)",ylab="observed p-value (log base 10)")
abline(0,1)
