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
cutoff=5e5
for (i in 1:nrow(highrisk1))
{
  if (i %%10==0) cat(i,'..')
  gprobe=i
  #find methylation probes
  gr_gt=GRanges(seqnames = highrisk1_pos$chr[i],ranges=IRanges(start=highrisk1_pos$pos[i],width = 1))
  tmp=distance(gr_ME1,gr_gt)
  idx=which(tmp<cutoff)
  tmp1=distance(gr_GE1,gr_gt)
  idx1=which(tmp1<cutoff)
  if (length(idx)>0 & length(idx1)>0)
  {
    for (j in 1:length(idx))
    {
      xprobe=idx[j] #methylation
      testtable=rbind.data.frame(testtable,data.frame(gprobe=rep(gprobe,length(idx1)),xprobe=rep(xprobe,length(idx1)),yprobe=idx1))
        # for (k in 1:length(idx1))
        # {
        #   yprobe=idx1[k] #gene expression
        #   testtable=rbind.data.frame(testtable,data.frame(gprobe=gprobe,xprobe=xprobe,yprobe=yprobe))
        # }
    
    }
  }
}

save(testtable,highrisk1,highrisk1_pos,ME1,ME1_pos,GE1,GE1_pos,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_testtable.RData")
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
mpi.bcast.cmd(library(cit))

result=data.frame(matrix(NA,nrow=nrow(testtable),ncol=5))
colnames(result)=c("genotype_idx","methylation_idx","geneexpression_idx","pcit","p")
result[,1:3]=testtable
num=1000
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/mediation_result.RData")
nrun <- ceiling(nrow(testtable)/num)
print(paste0("total num of run: ",nrun))
for (j in 1450:nrun){
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
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/GE_METHY_all.RData")
find_triplet=function(gt_idx=45,me_idx=50756,ge_idx=9899)
{
  gtname=highrisk1_pos$snp[gt_idx]
  
  idx=which(anno$IlmnID==rownames(ME1_pos)[me_idx])
  mename=as.character(anno$UCSC_RefGene_Name[idx])
  meprobe=as.character(anno$IlmnID[idx])
  idx=which(ge_anno$Probe_Id==rownames(GE1_pos)[ge_idx])
  gename=ge_anno$Symbol[idx]
  geprobe=ge_anno$Probe_Id[idx]
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