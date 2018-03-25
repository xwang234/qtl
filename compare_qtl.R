#!/usr/bin/env Rscript
library(data.table)
library(sas7bdat)
source("functions.R")
plotgenome=function(data,ylab="",main="",qvalues=NULL,qcutoff=0.05)
{
  chr=data$chr
  pos=data$pos
  value=data$value
  chrs=c(1:22,"X")
  dictfile="/fh/fast/dai_j/CancerGenomics/Tools/database/reference/compact/ucsc.hg19.compact.dict"
  chrlen=read.table(file=dictfile,sep="\t",skip=1,stringsAsFactors = F)
  chrlen=chrlen$V3
  chrlen=gsub("LN:","",chrlen,fixed=T)
  chrlen=as.numeric(chrlen)
  names(chrlen)=chrs
  chrstart=rep(0,length(chrlen))
  for (i in 2:length(chrlen))
  {
    tmp=0
    for (j in 1:(i-1))
    {
      tmp=tmp+chrlen[j]
    }
    chrstart[i]=tmp
  }
  names(chrstart)=names(chrlen)
  chrend=rep(0,length(chrlen))
  for (i in 1:length(chrlen))
  {
    tmp=0
    for (j in 1:i)
    {
      tmp=tmp+chrlen[j]
    }
    chrend[i]=tmp
  }
  names(chrend)=names(chrlen)
  
  res=data.frame(chr=chr,pos=pos,posall=rep(NA,length(chr)),value=value)
  res$chr=gsub(23,"X",res$chr)
  res$chr=gsub(24,"Y",res$chr)
  res$chr=factor(res$chr,levels=chrs)
  res=res[order(res$chr,res$pos),]
  res=res[res$chr %in% chrs,]
  for (mychr in chrs)
  {
    idx1=which(res$chr %in% mychr)
    idx2=which(names(chrstart)==mychr)
    res$posall[idx1]=res$pos[idx1]+chrstart[idx2]
  }
  
  ymax=1.1*max(res$value,na.rm=T)
  
  res1=res[complete.cases(res),]
  
  par(mar=c(1.1,4.1,1.1,1.1))
  plot(c(0,res$posall[nrow(res)]),c(-ymax/8,ymax),xaxt="n",yaxt="n",type="n", mgp=c(1,1,0),
       xlab="",ylab=ylab,cex=1.2,cex.lab=1.2,cex.axis=1.2,frame=F,main=main)
  
  # x.poly <- c(res1$posall, res1$posall[length(res1$posall)], res1$posall[1])         # Adjoin two x-coordinates
  # y.poly <- c(res1$value, 0, 0)                     #
  # polygon(x.poly, y.poly, col=gray(0.95), border=NA)          # Show the polygon fill only
  #lines(res1$posall,res1$value)
  segments(res1$posall,0,res1$posall,res1$value)
  axis(side=2, cex.axis=1.2,cex.lab=1.2,pos=0,lwd.ticks=1.1)
  #points(res1$posall,res1$value)
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  for (i in 1:length(chrs))
  {
    mychr=chrs[i]
    # rect(chrstart[i],1.05*max(res$dist2,na.rm=T),chrend[i],1.15*max(res$dist2,na.rm=T),
    #      col=color2[i],border=F)
    # text(0.5*(chrstart[i]+chrend[i]),1.15*max(res$dist2,na.rm=T),labels=chrs[i],col=color2[i],cex=1.3)
    rect(chrstart[i],-ymax/20,chrend[i],0,col=color2[i],border=F)
    #if (i<=15 | (i>15 & i %% 2==1))
    text(0.5*(chrstart[i]+chrend[i]),-ymax/10,labels=chrs[i],col=color2[i],cex=1.6,font=1.2)
    lines(c(chrend[i],chrend[i]),c(0,ymax),lty=2,col="gray",lwd=1)  
  }
  pvaluecutoff=NULL
  if (!is.null(qvalues))
  {
    idx=which(qvalues<=qcutoff)
    if (length(idx)>0) pvaluecutoff=max(10^(-value[idx]))
  }
  if (!is.null(pvaluecutoff))
  {
    segments(0,-log10(pvaluecutoff),par('usr')[2],-log10(pvaluecutoff),col="red")
    #text(sum(chrlen1)/2,-log10(pvaluecutoff)+0.5,paste0("FDR=",qcutoff))
    text(sum(chrlen)*0.75,-log10(pvaluecutoff)+0.1,paste0("FDR=",qcutoff),cex=1.2)
  }else
  {
    if (!is.null(qvalues))
    {
      text(sum(chrlen)*0.75,ymax-0.5,paste0("Minimum FDR=",round(min(qvalues,na.rm=T),2)),cex=1.2)
    }
  }
  return(res)
}

#load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/ge_annohg19.RData")


#CIS-------------------------------------------------------------------------------
#for highrisk HUTCH eqtl
#for number of tests
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_peer_highrisk.RData")
highrisk_HUTCH_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_HUTCH_cis_eqtl=addgenename_GE(highrisk_HUTCH_cis_eqtl)
length(unique(highrisk_HUTCH_cis_eqtl$genename))
highrisk_all_HUTCH_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
highrisk_all_HUTCH_cis_eqtl=addgenename_GE(highrisk_all_HUTCH_cis_eqtl)
length(unique(highrisk_all_HUTCH_cis_eqtl$genename))

#for highrisk NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk.RData")
highrisk_NORMAL_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_NORMAL_cis_eqtl=addgenename_GE(highrisk_NORMAL_cis_eqtl)
length(unique(highrisk_NORMAL_cis_eqtl$genename))
highrisk_all_NORMAL_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_cis",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt")
highrisk_all_NORMAL_cis_eqtl=addgenename_GE(highrisk_all_NORMAL_cis_eqtl)
length(unique(highrisk_all_NORMAL_cis_eqtl$genename))
#use less stringent cutoff pvalue>=0.01
length(unique(highrisk_all_HUTCH_cis_eqtl$snp_idx[highrisk_all_HUTCH_cis_eqtl$value>=2]))
length(unique(highrisk_all_NORMAL_cis_eqtl$snp_idx[highrisk_all_NORMAL_cis_eqtl$value>=2]))
length(intersect(unique(highrisk_all_HUTCH_cis_eqtl$snp_idx[highrisk_all_HUTCH_cis_eqtl$value>=2]),unique(highrisk_all_NORMAL_cis_eqtl$snp_idx[highrisk_all_NORMAL_cis_eqtl$value>=2])))
pvalue_overlapsnp(147,79,43,32)
#for HUTCH eqtl

rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl.RData")
HUTCH_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",bonfadj=qtl$cis$ntests)
HUTCH_cis_eqtl=addgenename_GE(HUTCH_cis_eqtl)
length(unique(HUTCH_cis_eqtl$genename))
#for NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl.RData")
NORMAL_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_cis",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$cis$ntests)
NORMAL_cis_eqtl=addgenename_GE(NORMAL_cis_eqtl)
length(unique(NORMAL_cis_eqtl$genename))
length(unique(intersect(HUTCH_cis_eqtl$genename,NORMAL_cis_eqtl$genename)))
checkoverlapsnps()
pvalue_overlapsnp(6500779,98928,2286,1792)

postscript("/fh/fast/stanford_j/Xiaoyu/QTL//result/cis_eqtl_HUTCH_NORMAL.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
par(mfrow=c(2,1))
tmp1=plotgenome(data=data.frame(chr=HUTCH_cis_eqtl$chr,pos=HUTCH_cis_eqtl$opos2,value=HUTCH_cis_eqtl$value))
tmp2=plotgenome(data=data.frame(chr=NORMAL_cis_eqtl$chr,pos=NORMAL_cis_eqtl$opos2,value=NORMAL_cis_eqtl$value))
dev.off()

#for highrisk HUTCH mqtl--------------------------------------------------
#for number of tests
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk.RData")
highrisk_HUTCH_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_HUTCH_cis_mqtl=addgenename_ME(highrisk_HUTCH_cis_mqtl)
length(unique(highrisk_HUTCH_cis_mqtl$genename))
highrisk_all_HUTCH_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
highrisk_all_HUTCH_cis_mqtl=addgenename_ME(highrisk_all_HUTCH_cis_mqtl)
length(unique(highrisk_all_HUTCH_cis_mqtl$genename))

#for highrisk NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk.RData")
highrisk_NORMAL_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_cis",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_NORMAL_cis_mqtl=addgenename_ME(highrisk_NORMAL_cis_mqtl)
length(unique(highrisk_NORMAL_cis_mqtl$genename))
sum(length(unique(intersect(highrisk_HUTCH_cis_mqtl$snp_idx,highrisk_NORMAL_cis_mqtl$snp_idx))))
highrisk_all_NORMAL_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_cis",
                                          snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                          geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt")
highrisk_all_NORMAL_cis_mqtl=addgenename_ME(highrisk_all_NORMAL_cis_mqtl)
length(unique(highrisk_all_NORMAL_cis_mqtl$genename))
checkoverlapsnps(qtlres1 = highrisk_HUTCH_cis_mqtl,qtlres2 = highrisk_NORMAL_cis_mqtl,
                 snpposfile1 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
pvalue_overlapsnp(147,68,9,7)
#use less stringent cutoff pvalue>=0.001
length(unique(highrisk_all_HUTCH_cis_mqtl$snp_idx[highrisk_all_HUTCH_cis_mqtl$value>=3]))
length(unique(highrisk_all_NORMAL_cis_mqtl$snp_idx[highrisk_all_NORMAL_cis_mqtl$value>=3]))
length(intersect(unique(highrisk_all_HUTCH_cis_mqtl$snp_idx[highrisk_all_HUTCH_cis_mqtl$value>=3]),unique(highrisk_all_NORMAL_cis_mqtl$snp_idx[highrisk_all_NORMAL_cis_mqtl$value>=3])))
#pvalue_overlapsnp(147,120,72,67)
#for HUTCH mqtl
rm(qtl)
#load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl.RData")
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI/mqtl.RData") #cis/trans/num_cis/num_trans
HUTCH_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI/mqtl_cis",
                            snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                            geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=sum(num_cis))
HUTCH_cis_mqtl=addgenename_ME(HUTCH_cis_mqtl)
length(unique(HUTCH_cis_mqtl$genename))

# rm(qtl)
# load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl.RData")
# HUTCH_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_cis",
#                              snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
#                              geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=qtl$cis$ntests)

#for NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl.RData")
NORMAL_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_cis",
                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt",
                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$cis$ntests)
NORMAL_cis_mqtl=addgenename_ME(NORMAL_cis_mqtl)
length(unique(NORMAL_cis_mqtl$genename))
length(unique(intersect(HUTCH_cis_mqtl$genename,NORMAL_cis_mqtl$genename)))
checkoverlapsnps(qtlres1=HUTCH_cis_mqtl,qtlres2=NORMAL_cis_mqtl,
                          snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                          snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt")
pvalue_overlapsnp(6500800,1063155,24447,23208)

postscript("/fh/fast/stanford_j/Xiaoyu/QTL//result/cis_mqtl_HUTCH_NORMAL.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
par(mfrow=c(2,1))
tmp1=plotgenome(data=data.frame(chr=HUTCH_cis_mqtl$chr,pos=HUTCH_cis_mqtl$opos2,value=HUTCH_cis_mqtl$value))
tmp2=plotgenome(data=data.frame(chr=NORMAL_cis_mqtl$chr,pos=NORMAL_cis_mqtl$opos2,value=NORMAL_cis_mqtl$value))
dev.off()

#
#Trans----------------------------------------------------------------------------------------------------------
#for highrisk HUTCH eqtl
#for number of tests
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk.RData")
highrisk_HUTCH_trans_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_HUTCH_trans_eqtl=addgenename_GE(highrisk_HUTCH_trans_eqtl)
length(unique(highrisk_HUTCH_trans_eqtl$genename))
highrisk_all_HUTCH_trans_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
highrisk_all_HUTCH_trans_eqtl=addgenename_GE(highrisk_all_HUTCH_trans_eqtl)
length(unique(highrisk_all_HUTCH_trans_eqtl$genename))

#for highrisk NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk.RData")
highrisk_NORMAL_trans_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_trans",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_NORMAL_trans_eqtl=addgenename_GE(highrisk_NORMAL_trans_eqtl)
length(unique(highrisk_NORMAL_trans_eqtl$genename))
checkoverlapsnps(qtlres1 = highrisk_HUTCH_trans_eqtl,qtlres2 = highrisk_NORMAL_trans_eqtl,
                 snpposfile1 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
highrisk_all_NORMAL_trans_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_trans",
                                          snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                          geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt")
highrisk_all_NORMAL_trans_eqtl=addgenename_GE(highrisk_all_NORMAL_trans_eqtl)
length(unique(highrisk_all_NORMAL_trans_eqtl$genename))
#use less stringent cutoff pvalue>=0.01
length(unique(highrisk_all_HUTCH_trans_eqtl$snp_idx[highrisk_all_HUTCH_trans_eqtl$value>=5]))
length(unique(highrisk_all_NORMAL_trans_eqtl$snp_idx[highrisk_all_NORMAL_trans_eqtl$value>=5]))
length(intersect(unique(highrisk_all_HUTCH_trans_eqtl$snp_idx[highrisk_all_HUTCH_trans_eqtl$value>=5]),unique(highrisk_all_NORMAL_trans_eqtl$snp_idx[highrisk_all_NORMAL_trans_eqtl$value>=5])))
pvalue_overlapsnp(147,35,28,11)
#for HUTCH eqtl

rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl.RData")
HUTCH_trans_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_trans",
                            snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                            geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",bonfadj=qtl$trans$ntests)
HUTCH_trans_eqtl=addgenename_GE(HUTCH_trans_eqtl)
length(unique(HUTCH_trans_eqtl$genename))
#for NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl.RData")
NORMAL_trans_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_trans",
                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt",
                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$trans$ntests)
NORMAL_trans_eqtl=addgenename_GE(NORMAL_trans_eqtl)
length(unique(NORMAL_trans_eqtl$genename))
length(unique(intersect(HUTCH_trans_eqtl$genename,NORMAL_trans_eqtl$genename)))
checkoverlapsnps(qtlres1 = HUTCH_trans_eqtl,qtlres2 = NORMAL_trans_eqtl,
                 snpposfile1 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                 snpposfile2 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt")
pvalue_overlapsnp(6500779,15111,7918,20)

postscript("/fh/fast/stanford_j/Xiaoyu/QTL//result/trans_eqtl_HUTCH_NORMAL.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
par(mfrow=c(2,1))
tmp1=plotgenome(data=data.frame(chr=HUTCH_trans_eqtl$chr,pos=HUTCH_trans_eqtl$opos2,value=HUTCH_trans_eqtl$value))
tmp2=plotgenome(data=data.frame(chr=NORMAL_trans_eqtl$chr,pos=NORMAL_trans_eqtl$opos2,value=NORMAL_trans_eqtl$value))
dev.off()

#for highrisk HUTCH mqtl--------------------------------------------------
#for number of tests
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk.RData")
highrisk_HUTCH_trans_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_HUTCH_trans_mqtl=addgenename_ME(highrisk_HUTCH_trans_mqtl)
length(unique(highrisk_HUTCH_trans_mqtl$genename))
highrisk_all_HUTCH_trans_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
highrisk_all_HUTCH_trans_mqtl=addgenename_ME(highrisk_all_HUTCH_trans_mqtl)
length(unique(highrisk_all_HUTCH_trans_mqtl$genename))

#for highrisk NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk.RData")
highrisk_NORMAL_trans_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_trans",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_NORMAL_trans_mqtl=addgenename_ME(highrisk_NORMAL_trans_mqtl)
length(unique(highrisk_NORMAL_trans_mqtl$genename))
sum(length(unique(intersect(highrisk_HUTCH_trans_mqtl$snp_idx,highrisk_NORMAL_trans_mqtl$snp_idx))))
highrisk_all_NORMAL_trans_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_trans",
                                          snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                          geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt")
highrisk_all_NORMAL_trans_mqtl=addgenename_ME(highrisk_all_NORMAL_trans_mqtl)
length(unique(highrisk_all_NORMAL_trans_mqtl$genename))

checkoverlapsnps(qtlres1=highrisk_HUTCH_trans_mqtl,qtlres2=highrisk_NORMAL_trans_mqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
pvalue_overlapsnp(147,17,9,6)
#use less stringent cutoff pvalue>=0.001
length(unique(highrisk_all_HUTCH_trans_mqtl$snp_idx[highrisk_all_HUTCH_trans_mqtl$value>=7]))
length(unique(highrisk_all_NORMAL_trans_mqtl$snp_idx[highrisk_all_NORMAL_trans_mqtl$value>=7]))
length(intersect(unique(highrisk_all_HUTCH_trans_mqtl$snp_idx[highrisk_all_HUTCH_trans_mqtl$value>=7]),unique(highrisk_all_NORMAL_trans_mqtl$snp_idx[highrisk_all_NORMAL_trans_mqtl$value>=7])))

#for HUTCH mqtl
rm(qtl)
#load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl.RData")
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI/mqtl.RData") #cis/trans/num_cis/num_trans
HUTCH_trans_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI/mqtl_trans",
                            snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                            geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=sum(num_trans))
HUTCH_trans_mqtl=addgenename_ME(HUTCH_trans_mqtl)
length(unique(HUTCH_trans_mqtl$genename))
#for NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl.RData")
NORMAL_trans_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_trans",
                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt",
                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$trans$ntests)
NORMAL_trans_mqtl=addgenename_ME(NORMAL_trans_mqtl)
length(unique(NORMAL_trans_mqtl$genename))
length(unique(intersect(HUTCH_trans_mqtl$genename,NORMAL_trans_mqtl$genename)))
checkoverlapsnps(qtlres1=HUTCH_trans_mqtl,qtlres2=NORMAL_trans_mqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt")
pvalue_overlapsnp(6500800,150311,15504,1764)

postscript("/fh/fast/stanford_j/Xiaoyu/QTL//result/trans_mqtl_HUTCH_NORMAL.ps",
           horizontal=T,width = 18, height = 4,pointsize = 6)
par(mfrow=c(2,1))
tmp1=plotgenome(data=data.frame(chr=HUTCH_trans_mqtl$chr,pos=HUTCH_trans_mqtl$opos2,value=HUTCH_trans_mqtl$value))
tmp2=plotgenome(data=data.frame(chr=NORMAL_trans_mqtl$chr,pos=NORMAL_trans_mqtl$opos2,value=NORMAL_trans_mqtl$value))
dev.off()

#compare with Thibodeau's result
Thibosnps=read.table("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/Thibo.txt",header=T,sep="\t",fill=T,stringsAsFactors = F)
idx=which(Thibosnps$P.value2<=1.96e-7)
length(unique(Thibosnps$Gene[idx]))
#[1] 88 genes
#Thibosnps=Thibosnps[Thibosnps$Risk.Region!="",]
idx=Thibosnps$LD..r2.==1 & !is.na(Thibosnps$position..bp...1) & Thibosnps$position..bp...1==Thibosnps$position..bp.1
Thibosnps1=Thibosnps[idx,c("PC.Risk.SNP","Chr","position..bp...1","Gene","P.value2")]
colnames(Thibosnps1)=c("snp","chr","pos","gene","pvalue")
tmp=Thibosnps[!idx,c("Associated.SNP","Chr","position..bp.1","Gene","P.value2")]
colnames(tmp)=c("snp","chr","pos","gene","pvalue")
Thibosnps1=rbind.data.frame(Thibosnps1,tmp)
Thibosnps1$hutch=0
for (i in 1:nrow(Thibosnps1))
{
  if (sum(highrisk_SNP_POS$chr==Thibosnps1$chr[i] & highrisk_SNP_POS$pos==Thibosnps1$pos[i])>0)
  {
    Thibosnps1$hutch[i]=1
  }
}
Thibosnps1$hutch_pvalue=NA
Thibosnps1$hutch_gene=NA
for (i in 1:nrow(Thibosnps1))
{
  if(Thibosnps1$hutch[i]==1)
  {
    idx1=which(highrisk_all_HUTCH_cis_eqtl$chr==Thibosnps1$chr[i] & highrisk_all_HUTCH_cis_eqtl$opos1==Thibosnps1$pos[i])
    if (length(idx1)>0)
    {
      idx2=which.max(highrisk_all_HUTCH_cis_eqtl$value[idx1])
      Thibosnps1$hutch_pvalue[i]=10^-highrisk_all_HUTCH_cis_eqtl$value[idx1[idx2]]
      Thibosnps1$hutch_gene[i]=as.character(highrisk_all_HUTCH_cis_eqtl$genename[idx1[idx2]])
    }
  }
}

Thibosnps1$normal_pvalue=NA
Thibosnps1$normal_gene=NA
for (i in 1:nrow(Thibosnps1))
{
  if(Thibosnps1$hutch[i]==1)
  {
    idx1=which(highrisk_all_NORMAL_cis_eqtl$chr==Thibosnps1$chr[i] & highrisk_all_NORMAL_cis_eqtl$opos1==Thibosnps1$pos[i])
    if (length(idx1)>0)
    {
      idx2=which.max(highrisk_all_NORMAL_cis_eqtl$value[idx1])
      Thibosnps1$normal_pvalue[i]=10^-highrisk_all_NORMAL_cis_eqtl$value[idx1[idx2]]
      Thibosnps1$normal_gene[i]=as.character(highrisk_all_NORMAL_cis_eqtl$genename[idx1[idx2]])
    }
  }
}
sum(Thibosnps1$hutch==1)
#[1] 12
sum(Thibosnps1$hutch==1 &Thibosnps1$pvalue<1.96e-7)
#[1] 6
idx1=which(Thibosnps1$hutch==1 &Thibosnps1$pvalue<1.96e-7)
View(Thibosnps1[idx1,])
#confirmed by hutch
sum(Thibosnps1$hutch==1 &Thibosnps1$pvalue<1.96e-7 & Thibosnps1$hutch_pvalue<=1.15e-5,na.rm=T)
#[1] 6
quantile(Thibosnps1$hutch_pvalue[Thibosnps1$hutch==1 &Thibosnps1$pvalue<1.96e-7],na.rm=T)
#          0%          25%          50%          75%         100% 
#1.374295e-13 3.056508e-10 3.056508e-10 3.056508e-10 5.728093e-10
#confirmed by normal
sum(Thibosnps1$hutch==1 &Thibosnps1$pvalue<1.96e-7 & Thibosnps1$normal_pvalue<=1.6e-5,na.rm=T)
#[1] 0
quantile(Thibosnps1$normal_pvalue[Thibosnps1$hutch==1 &Thibosnps1$pvalue<1.96e-7],na.rm=T)
#          0%          25%          50%          75%         100% 
#0.001381946 0.001461700 0.001461700 0.001461700 0.020331717

#work on highrisk Nov21---------------
#check the imputation quality of the 147 SNPs in Hutch data and in normal data
#summary about how many SNPs have to be imputed, what are their imputation qualities
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData")
hutch_highrisk=highrisk
hutch_allsnps=allsnps
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_allhighrisksnps_new.RData") 
tcga_highrisk=highrisk
tcga_allsnps=allsnps
tcga_allsnps[130,]
#            snp chr position source info exp_freq_a1
#130 rs138213197  17 46805705  Known   NA          NA
tmp=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation4/SNP6_info_chr17.txt")
tmp=data.frame(tmp)
sum(tmp$position==tcga_allsnps$position[130])
#[1] 0 the snp doesnot appear in imputed data
147-sum(hutch_allsnps$r2_icogs==1,na.rm=T)
#[1] 83
147-sum(hutch_allsnps$r2_onco==1,na.rm=T)
#[1] 64
sum(hutch_allsnps$r2_icogs<1 & hutch_allsnps$r2_onco<1,na.rm=T)
#[1] 54
147-sum(tcga_allsnps$info==1,na.rm=T)
#[1] 119
quantile(hutch_allsnps$r2_icogs,na.rm=T)
#       0%       25%       50%       75%      100% 
#0.4951550 0.9480705 0.9966270 1.0000000 1.0000000 
quantile(hutch_allsnps$r2_onco,na.rm=T)
#      0%       25%       50%       75%      100% 
#0.7724670 0.9943795 1.0000000 1.0000000 1.0000000 
quantile(tcga_allsnps$info,na.rm=T)
#     0%     25%     50%     75%    100% 
#0.48800 0.94850 0.98400 0.99675 1.00000
hutchclinical=read.sas7bdat("/fh/fast/stanford_j/Janet/alldata_2016dec20.sas7bdat")

hutch_geneexp_samples=hutchclinical$studyno[which(hutchclinical$race_==1 & (hutchclinical$in_iCOGS==1 | hutchclinical$in_OncoArray==1) & hutchclinical$GeneExpr==1)]
hutch_methylation_samples=hutchclinical$studyno[which(hutchclinical$race_==1 & (hutchclinical$in_iCOGS==1 | hutchclinical$in_OncoArray==1) & hutchclinical$Methy==1)]
sum(hutchclinical$race_==1 & hutchclinical$in_iCOGS==1 &hutchclinical$GeneExpr==1,na.rm=T)
#[1] 208
sum(hutchclinical$race_==1 & hutchclinical$in_OncoArray==1 &hutchclinical$GeneExpr==1,na.rm=T)
#[1] 147
sum(hutchclinical$race_==1 & hutchclinical$in_iCOGS==1 &hutchclinical$Methy==1,na.rm=T)
#[1] 219
sum(hutchclinical$race_==1 & hutchclinical$in_OncoArray==1 & hutchclinical$Methy==1,na.rm=T)
#[1] 158
sum(hutchclinical$race_==1 & hutchclinical$in_iCOGS==1 & (hutchclinical$GeneExpr==1 | hutchclinical$Methy==1),na.rm=T)
#[1] 234
sum(hutchclinical$race_==1 & hutchclinical$in_OncoArray==1 & (hutchclinical$GeneExpr==1 | hutchclinical$Methy==1),na.rm=T)
#[1] 161
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/allnormaldata.RData")
sum(colnames(hutch_geneexp) %in% hutchclinical$studyno[hutchclinical$in_iCOGS==1])
#[1] 11
sum(colnames(hutch_geneexp) %in% hutchclinical$studyno[hutchclinical$in_OncoArray==1])
#[3]


# Among 147 risk loci, how many of them are eQTL alone, mQTL alone, both eQTL and mQTL, neither; use both FDR<0.05 and FWER <0.05 cut off
extractpairs=function(dat=highrisk_all_HUTCH_cis_eqtl,fdrcutoff=NULL)
{
  if (!is.null(fdrcutoff))
  {
    dat=dat[dat$fdr<=fdrcutoff,]
  }
  res=data.frame(dat[,c("snp_idx","chr","opos2","value","fdr","gene","genename")])
}


hutch_ge_pos=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",header=T,stringsAsFactors = F)
tcga_ge_pos=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",header=T,stringsAsFactors = F)
hutch_me_pos=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",header=T,stringsAsFactors = F)
tcga_me_pos=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",header=T,stringsAsFactors = F)
highrisk_snp_pos=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",header=T,stringsAsFactors = F)
hutch_ge=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",header=T,sep="\t",stringsAsFactors = F,row.names = 1)
colnames(hutch_ge)=gsub("^X","",colnames(hutch_ge))
hutch_me=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",header=T,sep="\t",stringsAsFactors = F,row.names = 1)
colnames(hutch_me)=gsub("^X","",colnames(hutch_me))
normal_ge=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt",header=T,sep="\t",stringsAsFactors = F,row.names = 1)
colnames(normal_ge)=gsub("^X","",colnames(normal_ge))
normal_me=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME.txt",header=T,sep="\t",stringsAsFactors = F,row.names = 1)
colnames(normal_me)=gsub("^X","",colnames(normal_me))
highrisk_hutch_snp_ge=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",header=T,sep="\t",stringsAsFactors = F,row.names = 1)
colnames(highrisk_hutch_snp_ge)=gsub("^X","",colnames(highrisk_hutch_snp_ge))
highrisk_hutch_snp_me=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_ME.txt",header=T,sep="\t",stringsAsFactors = F,row.names = 1)
colnames(highrisk_hutch_snp_me)=gsub("^X","",colnames(highrisk_hutch_snp_me))
highrisk_normal_snp_ge=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_GE.txt",header=T,sep="\t",stringsAsFactors = F,row.names = 1)
colnames(highrisk_normal_snp_ge)=gsub("^X","",colnames(highrisk_normal_snp_ge))
highrisk_normal_snp_me=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_ME.txt",header=T,sep="\t",stringsAsFactors = F,row.names = 1)
colnames(highrisk_normal_snp_me)=gsub("^X","",colnames(highrisk_normal_snp_me))
hutch_ge_covariate=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(hutch_ge_covariate)=gsub("^X","",colnames(hutch_ge_covariate))
hutch_me_covariate=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(hutch_me_covariate)=gsub("^X","",colnames(hutch_me_covariate))
normal_ge_covariate=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_GE.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(normal_ge_covariate)=gsub("^X","",colnames(normal_ge_covariate))
normal_me_covariate=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_ME.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(normal_me_covariate)=gsub("^X","",colnames(normal_me_covariate))

highrisk_cis_hutch_eqtl=extractpairs(dat=highrisk_HUTCH_cis_eqtl)
highrisk_cis_hutch_eqtl_fdr=extractpairs(dat=highrisk_all_HUTCH_cis_eqtl,fdrcutoff = 0.05)
length(unique(highrisk_cis_hutch_eqtl_fdr$snp_idx)) #37
length(unique(highrisk_cis_hutch_eqtl_fdr$gene)) #60
highrisk_cis_hutch_mqtl=extractpairs(dat=highrisk_HUTCH_cis_mqtl) #68
highrisk_cis_hutch_mqtl_fdr=extractpairs(dat=highrisk_all_HUTCH_cis_mqtl,fdrcutoff = 0.05) #1405
length(unique(highrisk_cis_hutch_mqtl_fdr$snp_idx)) #116
length(unique(highrisk_cis_hutch_mqtl_fdr$gene)) #1290

highrisk_cis_normal_eqtl=extractpairs(dat=highrisk_NORMAL_cis_eqtl)
highrisk_cis_normal_eqtl_fdr=extractpairs(dat=highrisk_all_NORMAL_cis_eqtl,fdrcutoff = 0.05) #37
length(highrisk_cis_normal_eqtl_fdr$snp_idx)
highrisk_cis_normal_mqtl=extractpairs(dat=highrisk_NORMAL_cis_mqtl) #68
highrisk_cis_normal_mqtl_fdr=extractpairs(dat=highrisk_all_NORMAL_cis_mqtl,fdrcutoff = 0.05) #116
length(unique(highrisk_cis_normal_mqtl_fdr$snp_idx)) #18
length(unique(highrisk_cis_normal_mqtl_fdr$gene)) #29
length(intersect(unique(highrisk_cis_normal_mqtl_fdr$snp_idx),unique(highrisk_cis_hutch_mqtl_fdr$snp_idx)))

highrisk_trans_hutch_eqtl=extractpairs(dat=highrisk_HUTCH_trans_eqtl)
highrisk_trans_hutch_eqtl_fdr=extractpairs(dat=highrisk_all_HUTCH_trans_eqtl,fdrcutoff = 0.05) #37
length(unique(highrisk_trans_hutch_eqtl_fdr$snp_idx)) #5
length(unique(highrisk_trans_hutch_eqtl_fdr$gene)) #16

highrisk_trans_hutch_mqtl=extractpairs(dat=highrisk_HUTCH_trans_mqtl) #68
highrisk_trans_hutch_mqtl_fdr=extractpairs(dat=highrisk_all_HUTCH_trans_mqtl,fdrcutoff = 0.05) #116
length(unique(highrisk_trans_hutch_mqtl_fdr$snp_idx)) #91
length(unique(highrisk_trans_hutch_mqtl_fdr$gene)) #2302

highrisk_trans_normal_eqtl=extractpairs(dat=highrisk_NORMAL_trans_eqtl)
highrisk_trans_normal_eqtl_fdr=extractpairs(dat=highrisk_all_NORMAL_trans_eqtl,fdrcutoff = 0.05) #37
length(unique(highrisk_trans_normal_eqtl_fdr$snp_idx)) #18
length(unique(highrisk_trans_normal_eqtl_fdr$gene)) #157
length(intersect(unique(highrisk_trans_normal_eqtl_fdr$snp_idx),unique(highrisk_trans_hutch_eqtl_fdr$snp_idx)))

highrisk_trans_normal_mqtl=extractpairs(dat=highrisk_NORMAL_trans_mqtl) #68
highrisk_trans_normal_mqtl_fdr=extractpairs(dat=highrisk_all_NORMAL_trans_mqtl,fdrcutoff = 0.05) #116
length(unique(highrisk_trans_normal_mqtl_fdr$snp_idx)) #97
length(unique(highrisk_trans_normal_mqtl_fdr$gene)) #4253
length(intersect(unique(highrisk_trans_normal_mqtl_fdr$snp_idx),unique(highrisk_trans_hutch_mqtl_fdr$snp_idx)))

pvalue_overlapsnp(147,24,1,1)
pvalue_overlapsnp(147,68,9,7)
pvalue_overlapsnp(147,2,7,1)
pvalue_overlapsnp(147,17,9,6)
pvalue_overlapsnp(147,37,1,1)
pvalue_overlapsnp(147,116,18,17)
pvalue_overlapsnp(147,5,18,3)
pvalue_overlapsnp(147,92,97,66)

save(highrisk_all_HUTCH_cis_eqtl,highrisk_all_NORMAL_cis_eqtl,highrisk_all_HUTCH_cis_mqtl,highrisk_all_NORMAL_cis_mqtl,
     highrisk_all_HUTCH_trans_eqtl,highrisk_all_NORMAL_trans_eqtl,highrisk_all_HUTCH_trans_mqtl,highrisk_all_NORMAL_trans_mqtl,
     highrisk_hutch_snp_ge,highrisk_hutch_snp_me,highrisk_normal_snp_ge,highrisk_normal_snp_me,
     hutch_ge,hutch_me,normal_ge,normal_me,hutch_ge_covariate,hutch_me_covariate,normal_ge_covariate,normal_me_covariate,
     highrisk_snp_pos,hutch_ge_pos,hutch_me_pos,tcga_ge_pos,tcga_me_pos,highrisk_cis_hutch_eqtl,highrisk_cis_hutch_eqtl_fdr,
     highrisk_cis_hutch_mqtl,highrisk_cis_hutch_mqtl_fdr,highrisk_cis_normal_eqtl,highrisk_cis_normal_eqtl_fdr,highrisk_cis_normal_mqtl,
     highrisk_cis_normal_mqtl_fdr,
     highrisk_trans_hutch_eqtl,highrisk_trans_hutch_eqtl_fdr,
     highrisk_trans_hutch_mqtl,highrisk_trans_hutch_mqtl_fdr,highrisk_trans_normal_eqtl,highrisk_trans_normal_eqtl_fdr,highrisk_trans_normal_mqtl,
     highrisk_trans_normal_mqtl_fdr,
     file="/fh/fast/stanford_j/Xiaoyu/QTL/result/highrisk_qtlresult.RData")

olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl=overlappairs(dat1=highrisk_cis_hutch_eqtl,dat2=highrisk_cis_hutch_mqtl)
olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl_fdr=overlappairs(dat1=highrisk_cis_hutch_eqtl_fdr,dat2=highrisk_cis_hutch_mqtl_fdr)
olap_highrisk_cis_normal_eqtl_cis_normal_mqtl=overlappairs(dat1=highrisk_cis_normal_eqtl,dat2=highrisk_cis_normal_mqtl)
olap_highrisk_cis_hutch_eqtl_cis_normal_eqtl=overlappairs(dat1=highrisk_cis_hutch_eqtl,dat2=highrisk_cis_normal_eqtl)
olap_highrisk_cis_hutch_mqtl_cis_normal_mqtl=overlappairs(dat1=highrisk_cis_hutch_mqtl,dat2=highrisk_cis_normal_mqtl)



tcga_me_pos=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_POS.txt",header=T,stringsAsFactors = F) #the order is different from normal
validate_highrisk_cis_hutch_eqtl=validatepairs()
validate_highrisk_cis_hutch_mqtl=validatepairs(dat1=highrisk_cis_hutch_mqtl,dat2=highrisk_all_NORMAL_cis_mqtl,phenotype2=tcga_me_pos)
validate_highrisk_trans_hutch_eqtl=validatepairs(dat1=highrisk_trans_hutch_eqtl,dat2=highrisk_all_NORMAL_trans_eqtl,phenotype2=tcga_ge_pos)
validate_highrisk_trans_hutch_mqtl=validatepairs(dat1=highrisk_trans_hutch_mqtl,dat2=highrisk_all_NORMAL_trans_mqtl,phenotype2=tcga_me_pos)
validtriplet_highrisk_cis_hutch_eqtl_cis_hutch_mqtl=validtriplet()

validate_highrisk_cis_hutch_eqtl_fdr=validatepairs(dat1=highrisk_cis_hutch_eqtl_fdr,dat2=highrisk_all_NORMAL_cis_eqtl,phenotype2=tcga_ge_pos)
validate_highrisk_cis_hutch_mqtl_fdr=validatepairs(dat1=highrisk_cis_hutch_mqtl_fdr,dat2=highrisk_all_NORMAL_cis_mqtl,phenotype2=tcga_me_pos)
validtriplet_highrisk_cis_hutch_eqtl_cis_hutch_mqtl=validtriplet(dat1=olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl_fdr$triplet,validateeqtl=highrisk_all_NORMAL_cis_eqtl,
                                                                 validatemqtl=highrisk_all_NORMAL_cis_mqtl,phenotype1=tcga_ge_pos,phenotype2=tcga_me_pos)
validate_highrisk_trans_hutch_eqtl_fdr=validatepairs(dat1=highrisk_trans_hutch_eqtl_fdr,dat2=highrisk_all_NORMAL_trans_eqtl,phenotype2=tcga_ge_pos)
validate_highrisk_trans_hutch_mqtl_fdr=validatepairs(dat1=highrisk_trans_hutch_mqtl_fdr,dat2=highrisk_all_NORMAL_trans_mqtl,phenotype2=tcga_me_pos)

load("/fh/fast/stanford_j/Xiaoyu/QTL/result/hutch_medresult.RData")
validatemediation=function(dat1=hutch_med_ciseqtl_cismqtl,dat2=normal_meds_ciseqtl_cismqtl,phenotype1=tcga_ge_pos,phenotype2=tcga_me_pos)
{
  print(paste0("number of triplet:",nrow(dat1)))
  print(paste0("number of triplets can be validate:",sum(dat1$gene.x %in% phenotype1[,1] & dat1$gene.y %in% phenotype2[,1],na.rm=T)))
  idxTeqtl=!grepl("^cg",dat2$gene.x)
  print(paste0("number of causal triplet:",sum(idxTeqtl)))
  print(paste0("number of causal triplet can be validated:",sum(!is.na(dat2$valid_p_cit[idxTeqtl]))))
  print(paste0("number of causal triplet validated:",sum(dat2$valid_p_cit[idxTeqtl]<=0.05,na.rm=T)))
  print(paste0("number of reactive triplet:",sum(!idxTeqtl)))
  print(paste0("number of reactive triplet can be validated:",sum(!is.na(dat2$valid_p_cit[!idxTeqtl]))))
  print(paste0("number of reactive triplet validated:",sum(dat2$valid_p_cit[!idxTeqtl]<=0.05,na.rm=T)))
  tmp=dat2[c("snp_idx","gene.x","genename.x","chr.x","opos2.x","gene.y","genename.y","chr.y","opos2.y","p.raw","q.cit","valid_p_cit")]
  tmp=tmp[order(tmp$p.raw),]
  tmp$snp_idx=highrisk_snp_pos$snp[tmp$snp_idx]
  allqtl1=dat1[,c("snp_idx","gene.x")]
  allqtl1$snp_idx=highrisk_snp_pos$snp[allqtl1$snp_idx]
  allqtl1=unique(allqtl1)
  allqtl2=dat1[,c("snp_idx","gene.y")]
  allqtl2$snp_idx=highrisk_snp_pos$snp[allqtl2$snp_idx]
  allqtl2=unique(allqtl2)
  T_eqtl=dat2[idxTeqtl,c("snp_idx","gene.x")]
  T_eqtl$snp_idx=highrisk_snp_pos$snp[T_eqtl$snp_idx]
  T_eqtl=unique(T_eqtl)
  G_mqtl=dat2[idxTeqtl,c("snp_idx","gene.y")]
  G_mqtl$snp_idx=highrisk_snp_pos$snp[G_mqtl$snp_idx]
  G_mqtl=unique(G_mqtl)
  
  T_mqtl=dat2[!idxTeqtl,c("snp_idx","gene.x")]
  T_mqtl$snp_idx=highrisk_snp_pos$snp[T_mqtl$snp_idx]
  T_mqtl=unique(T_mqtl)
  G_eqtl=dat2[!idxTeqtl,c("snp_idx","gene.y")]
  G_eqtl$snp_idx=highrisk_snp_pos$snp[G_eqtl$snp_idx]
  G_eqtl=unique(G_eqtl)
  
  comeqtl=merge(T_eqtl,G_eqtl,by.x=c("snp_idx","gene.x"),by.y=c("snp_idx","gene.y"))
  commqtl=merge(T_mqtl,G_mqtl,by.x=c("snp_idx","gene.x"),by.y=c("snp_idx","gene.y"))
  return(list(alltable=dat1,medtable=tmp,T_eqtl=T_eqtl,G_mqtl=G_mqtl,T_mqtl=T_mqtl,G_eqtl=G_eqtl,
              comeqtl=comeqtl,commqtl=commqtl))
}

validatemediation_hutch_med_ciseqtl_cismqtl=validatemediation()
idxTeqtl=!grepl("^cg",validatemediation_hutch_med_ciseqtl_cismqtl$gene.x)
write.table(validatemediation_hutch_med_ciseqtl_cismqtl$medtable[idxTeqtl,],file="/fh/fast/stanford_j/Xiaoyu/QTL/result/validatemediation_hutch_med_ciseqtl_cismqtl.txt",
            row.names = F,col.names = T,sep="\t",quote=F)
write.table(validatemediation_hutch_med_ciseqtl_cismqtl$medtable[!idxTeqtl,],file="/fh/fast/stanford_j/Xiaoyu/QTL/result/validatemediation_hutch_med_cismqtl_ciseqtl.txt",
            row.names = F,col.names = T,sep="\t",quote=F)
validatemediation_hutch_med_ciseqtl_transmqtl=validatemediation(dat1=hutch_med_ciseqtl_transmqtl,dat2=normal_meds_ciseqtl_transmqtl)
idxTeqtl=!grepl("^cg",validatemediation_hutch_med_ciseqtl_transmqtl$gene.x)
write.table(validatemediation_hutch_med_ciseqtl_transmqtl$medtable[idxTeqtl,],file="/fh/fast/stanford_j/Xiaoyu/QTL/result/validatemediation_hutch_med_ciseqtl_transmqtl.txt",
            row.names = F,col.names = T,sep="\t",quote=F)
write.table(validatemediation_hutch_med_ciseqtl_transmqtl$medtable[!idxTeqtl,],file="/fh/fast/stanford_j/Xiaoyu/QTL/result/validatemediation_hutch_med_transmqtl_ciseqtl.txt",
            row.names = F,col.names = T,sep="\t",quote=F)

classifymediation_pairs=function(validateresult=validatemediation_hutch_med_ciseqtl_cismqtl,
                                 alleqtl=highrisk_HUTCH_cis_eqtl,allmqtl=highrisk_cis_hutch_mqtl)
{
  alleqtl=alleqtl[,c("snp_idx","gene")]
  alleqtl$snp_idx=highrisk_snp_pos$snp[alleqtl$snp_idx]
  alleqtl=unique(alleqtl)
  alleqtl$T=0
  alleqtl$G=0
  idx=which(alleqtl$snp_idx %in% validateresult$T_eqtl$snp_idx & alleqtl$gene %in% validateresult$T_eqtl$gene.x)
  alleqtl$T[idx]=1
  idx=which(alleqtl$snp_idx %in% validateresult$G_eqtl$snp_idx & alleqtl$gene %in% validateresult$G_eqtl$gene.y)
  alleqtl$G[idx]=1
  
  allmqtl=allmqtl[,c("snp_idx","gene")]
  allmqtl$snp_idx=highrisk_snp_pos$snp[allmqtl$snp_idx]
  allmqtl=unique(allmqtl)
  allmqtl$T=0
  allmqtl$G=0
  idx=which(allmqtl$snp_idx %in% validateresult$T_mqtl$snp_idx & allmqtl$gene %in% validateresult$T_mqtl$gene.x)
  allmqtl$T[idx]=1
  idx=which(allmqtl$snp_idx %in% validateresult$G_mqtl$snp_idx & allmqtl$gene %in% validateresult$G_mqtl$gene.y)
  allmqtl$G[idx]=1
  
  print(paste0("number of eqtl pairs, total:",nrow(alleqtl)))
  print(paste0("number of eqtl pairs, worked as T:",sum(alleqtl$T==1)))
  print(paste0("number of eqtl pairs, worked as G:",sum(alleqtl$G==1)))
  print(paste0("number of eqtl pairs, worked as T & G:",sum(alleqtl$G==1 &alleqtl$T==1)))
  print(paste0("number of eqtl pairs, not worked as T | G:",sum(alleqtl$G==0 &alleqtl$T==0)))
  print(paste0("number of mqtl pairs, total:",nrow(allmqtl)))
  print(paste0("number of mqtl pairs, worked as T:",sum(allmqtl$T==1)))
  print(paste0("number of mqtl pairs, worked as G:",sum(allmqtl$G==1)))
  print(paste0("number of mqtl pairs, worked as T & G:",sum(allmqtl$G==1 &allmqtl$T==1)))
  print(paste0("number of mqtl pairs, not worked as T | G:",sum(allmqtl$G==0 &allmqtl$T==0)))
  
  return(list(alleqtl=alleqtl,allmqtl=allmqtl))
}  

classifymediation_hutch_med_ciseqtl_cismqtl=classifymediation_pairs()
classifymediation_hutch_med_ciseqtl_transmqtl=classifymediation_pairs(validateresult=validatemediation_hutch_med_ciseqtl_transmqtl,
                                                                      alleqtl=highrisk_HUTCH_cis_eqtl,allmqtl=highrisk_trans_hutch_mqtl)

classifymediation_all=funcion(classifyresult1=classifymediation_hutch_med_ciseqtl_cismqtl,
                              classifyresult2=classifymediation_hutch_med_ciseqtl_transmqtl)
{
  alleqtl=merge(classifyresult1$alleqtl,classifyresult2$alleqtl,by=c("snp_idx","gene"))
  allmqtl=merge(classifyresult1$allmqtl,classifyresult2$allmqtl,by=c("snp_idx","gene"),all.x=T) #no overlap in gene
  return(alleqtl)
}

library(VennDiagram)
#set allqtl before using it
likes <- function(animals) {
  if ("T.x" %in% colnames(allqtl))
  {
    ppl <- allqtl[,c("T.x","G.x","T.y","G.y")]
    names(ppl) <- c("T.x","G.x","T.y","G.y")
  }else
  {
    ppl <- allqtl[,c("T","G")]
    names(ppl) <- c("T","G")
  }
  
  for (i in 1:length(animals)) {
    ppl <- subset(ppl, ppl[animals[i]] == 1)
  }
  
  nrow(ppl)
}
likes("T.x")
likes(c("T.x","T.y"))
likes("T")
plotAnimals <- function(a, ...) {
  grid.newpage()
  if (length(a) == 1) {
    out <- draw.single.venn(likes(a), ...)
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(likes(a[1]), likes(a[2]), likes(a[1:2]), ...)
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[1:2]), 
                            likes(a[2:3]), likes(a[c(1, 3)]), likes(a), ...)
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[4]), 
                          likes(a[1:2]), likes(a[c(1, 3)]), likes(a[c(1, 4)]), likes(a[2:3]), 
                          likes(a[c(2, 4)]), likes(a[3:4]), likes(a[1:3]), likes(a[c(1, 2, 
                                                                                     4)]), likes(a[c(1, 3, 4)]), likes(a[2:4]), likes(a), ...)
  }
  if (!exists("out")) 
    out <- "Oops"
  return(out)
}
allqtl=alleqtl
library(Cairo)
library(ggplot2)
CairoPDF("eqtl_pairs.pdf")
plotAnimals(c("T.x","G.x","T.y","G.y"), category = c("T_cism-cise", "G_cism-cise", 
                                                     "T_transm-cise", "G_transm-cise"), lty = "blank", fill = c("skyblue", "pink1", 
                                                                                                                "mediumorchid", "orange"))
dev.off()
sum(alleqtl$T.x==0&alleqtl$G.x==0&alleqtl$T.y==0&alleqtl$G.y==0)
#[1] 12
allqtl=classifymediation_hutch_med_ciseqtl_cismqtl$allmqtl
sum(allqtl$T==0 & allqtl$G==0)
#[1] 300
CairoPDF("cis-mqtl_pairs.pdf",9,5)
plotAnimals(c("T","G"), category = c("T_cism-cise", "G_cism-cise"), lty = "blank", fill = c("skyblue", "pink1"))
dev.off()
allqtl=classifymediation_hutch_med_ciseqtl_transmqtl$allmqtl
sum(allqtl$T==0 & allqtl$G==0)
#[1] 436
CairoPDF("trans-mqtl_pairs.pdf",9,5)
plotAnimals(c("T","G"), category = c("T_transm-cise", "G_transm-cise"), lty = "blank", fill = c("skyblue", "pink1"))
dev.off()
idx=which(alleqtl$T.x==1&alleqtl$G.x==1&alleqtl$T.y==1&alleqtl$G.y==1)
eqtlhub=alleqtl[idx,]
eqtlhub$genename=eqtlhub$chr=eqtlhub$pos=NA
for (i in 1:nrow(eqtlhub))
{
  idx1=which(hutch_ge_anno$Probe_Id==eqtlhub$gene[i])
  eqtlhub$genename[i]=hutch_ge_anno$Symbol[idx1]
  eqtlhub$chr[i]=hutch_ge_anno$Chromosome[idx1]
  eqtlhub$pos[i]=hutch_ge_anno$start[idx1]
}
write.table(eqtlhub[,c(1,2,7:9)],file="cis-eqtlhub.txt",col.names = T,row.names = F,sep="\t",quote=F)

rmduplicatedgenename=function(dat1=highrisk_cis_hutch_mqtl)
{
  dat1$genename=as.character(dat1$genename)
  for (i in 1:nrow(dat1))
  {
    if (grepl(";",dat1$genename[i]))
    {
      tmp=unlist(strsplit(dat1$genename[i],";"))
      tmp=unique(tmp)
      tmp=paste0(tmp,collapse = ";")
      dat1$genename[i]=tmp
    }
  }
  return(dat1)
}
highrisk_cis_hutch_mqtl=rmduplicatedgenename()
highrisk_cis_hutch_mqtl_fdr=rmduplicatedgenename(dat1=highrisk_cis_hutch_mqtl_fdr)
highrisk_trans_hutch_mqtl=rmduplicatedgenename(dat1=highrisk_trans_hutch_mqtl)
highrisk_trans_hutch_mqtl_fdr=rmduplicatedgenename(dat1=highrisk_trans_hutch_mqtl_fdr)
highrisk_all_HUTCH_cis_mqtl=rmduplicatedgenename(dat1=highrisk_all_HUTCH_cis_mqtl)
highrisk_all_HUTCH_trans_mqtl=rmduplicatedgenename(dat1=highrisk_all_HUTCH_trans_mqtl)

idx=which(olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl$triplet$genename.x==olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl$triplet$genename.y)
olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl1=olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl$triplet[idx,]
length(unique(olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl1$snp_idx))
as.character(unique(olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl1$genename.x))
as.character(unique(highrisk_cis_hutch_eqtl$genename))
length(unique(highrisk_cis_hutch_mqtl_fdr$snp_idx))
length(unique(highrisk_cis_hutch_eqtl_fdr$snp_idx))
sum(unique(highrisk_cis_hutch_eqtl_fdr$snp_idx) %in% unique(highrisk_cis_hutch_mqtl_fdr$snp_idx))


length(unique(highrisk_cis_hutch_eqtl_fdr$genename))
length(unique(highrisk_cis_hutch_mqtl_fdr$genename))
allgenes=NULL
for (i in 1:length(highrisk_cis_hutch_mqtl_fdr$genename))
{
  tmp=unlist(strsplit(highrisk_cis_hutch_mqtl_fdr$genename[i],";"))
  tmp=unique(tmp)
  allgenes=c(allgenes,tmp)
  allgenes=unique(allgenes)
}
length(allgenes)
#[1] 321

snp_eqtl_mqtl=data.frame(snp=1:147,eqtl=0,mqtl=0)
for (i in 1:nrow(snp_eqtl_mqtl))
{
  if (i %in% highrisk_cis_hutch_eqtl_fdr$snp_idx)
  {
    snp_eqtl_mqtl$eqtl[i]=1
  }
  if (i %in% highrisk_cis_hutch_mqtl_fdr$snp_idx)
  {
    snp_eqtl_mqtl$mqtl[i]=1
  }
}
sum(snp_eqtl_mqtl$eqtl==1) #37
sum(snp_eqtl_mqtl$mqtl==1) #116
sum(snp_eqtl_mqtl$eqtl==1 & snp_eqtl_mqtl$mqtl==1) #35
sum(snp_eqtl_mqtl$eqtl==0 & snp_eqtl_mqtl$mqtl==0) #29

highrisk_cis_hutch_mqtl_fdr$snp_idx=highrisk_snp_pos$snp[highrisk_cis_hutch_mqtl_fdr$snp_idx]
highrisk_cis_hutch_eqtl_fdr$snp_idx=highrisk_snp_pos$snp[highrisk_cis_hutch_eqtl_fdr$snp_idx]
#problematic 450k probes:
library(gdata)
polycpg=read.xls("/fh/fast/dai_j/CancerGenomics/Tools/database/other/48640-polymorphic-CpGs-Illumina450k.xlsx")
snpcpg=read.xls("/fh/fast/dai_j/CancerGenomics/Tools/database/other/48640-polymorphic-CpGs-Illumina450k.xlsx",sheet=2)

#DEC20---
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors_allhighrisksnps_new.RData") 
tcga_highrisk=highrisk
tcga_allsnps=allsnps
tcga_allsnps[130,]
quantile(tcga_allsnps$info,na.rm=T)
# 0%     25%     50%     75%    100% 
# 0.57500 0.94725 0.97900 0.99500 1.00000 
#HUTCH eqtl----------
# #if use Bonferroni
# rm(qtl)
# load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_peer_highrisk.RData")
highrisk_HUTCH_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_cis",
                                           snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                           geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
#in the table, value:-log10(pvalue), gene:probeid, opos1:position of snp, opos2:position of geneexp
# [1] "processed qtl file: /fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_cis"
# [1] "number of SNPs: 147"
# [1] "number of phenotypeprobes: 25019"
# [1] "number of pairs: 4336"
# [1] "number of unique SNPs: 147"
# [1] "number of unique phenotypes: 3736"
# [1] "FDR cutoff: 0.05"
# [1] "pvalue cutoff: 0.000566975243517476"
# [1] "number of pairs (after adjustment)62"
# [1] "number of unique SNPs (after adjustment)40, 0.272108843537415"
# [1] "number of unique phenotypes (after adjustment)59, 0.0023582077621008"
highrisk_HUTCH_cis_eqtl=addgenename_GE(highrisk_HUTCH_cis_eqtl,anno=hutch_ge_anno)
length(unique(highrisk_HUTCH_cis_eqtl$genename))
highrisk_HUTCH_cis_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_cis",
                                   snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                   geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",fdrcutoff = 1)
highrisk_HUTCH_cis_eqtl_all=addgenename_GE(highrisk_HUTCH_cis_eqtl_all,anno=hutch_ge_anno)
#GE covariates use PCA
highrisk_HUTCH_pca_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
highrisk_HUTCH_pca_cis_eqtl=addgenename_GE(highrisk_HUTCH_pca_cis_eqtl,anno=hutch_ge_anno)
length(unique(highrisk_HUTCH_pca_cis_eqtl$genename))
#Old result use gleason score, GE_pca, and no SNP_pca
highrisk_HUTCH_old_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis_old",
                                       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                       geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
# [1] "pvalue cutoff: 0.000698999955987891"
# [1] "number of pairs (after adjustment)64"
# [1] "number of unique SNPs (after adjustment)37, 0.251700680272109"
# [1] "number of unique phenotypes (after adjustment)60, 0.00239817738518726"
highrisk_HUTCH_old_cis_eqtl=addgenename_GE(highrisk_HUTCH_old_cis_eqtl,anno=hutch_ge_anno)
length(unique(highrisk_HUTCH_old_cis_eqtl$genename))
#old result, no gleason score
highrisk_HUTCH_old_nogleason_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis_old_nogleason",
                                       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                       geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
# [1] "pvalue cutoff: 0.000674719235162517"
# [1] "number of pairs (after adjustment)62"
# [1] "number of unique SNPs (after adjustment)37, 0.251700680272109"
# [1] "number of unique phenotypes (after adjustment)58, 0.00231823813901435"
highrisk_HUTCH_old_nogleason_cis_eqtl=addgenename_GE(highrisk_HUTCH_old_nogleason_cis_eqtl,anno=hutch_ge_anno)
length(unique(highrisk_HUTCH_old_nogleason_cis_eqtl$genename))
#old result, no gleason score no age
highrisk_HUTCH_old_nogleasonage_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis_old_nogleasonage",
                                                 snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                                 geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
# [1] "pvalue cutoff: 0.00061858469996702"
# [1] "number of pairs (after adjustment)63"
# [1] "number of unique SNPs (after adjustment)38, 0.258503401360544"
# [1] "number of unique phenotypes (after adjustment)59, 0.0023582077621008"
highrisk_HUTCH_old_nogleasonage_cis_eqtl=addgenename_GE(highrisk_HUTCH_old_nogleasonage_cis_eqtl,anno=hutch_ge_anno)
length(unique(highrisk_HUTCH_old_nogleasonage_cis_eqtl$genename))

#compare peer with pca
checkoverlapsnps(qtlres1=highrisk_HUTCH_cis_eqtl,qtlres2=highrisk_HUTCH_pca_cis_eqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
# [1] "number of unique SNP1:40"
# [1] "number of unique SNP2:39"
# [1] "number of overlapped SNP:38"
# in2 notin2
# in1     38      1
# notin1   2    106
# [1] 1.306219e-31
ovarlappairs_HUTCH_cis_eqtl_peer_pca=checkoverlappairs_sameplatform(qtlres1=highrisk_HUTCH_cis_eqtl,qtlres2=highrisk_HUTCH_pca_cis_eqtl)
# [1] "number of pairs in qtlres1: 62"
# [1] "number of pairs in qtlres2: 61"
# [1] "total number of pairs: 64"
# [1] "number of common pairs: 59"
#compare peer with old result
checkoverlapsnps(qtlres1=highrisk_HUTCH_cis_eqtl,qtlres2=highrisk_HUTCH_old_cis_eqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
# [1] "number of unique SNP1:40"
# [1] "number of unique SNP2:37"
# [1] "number of overlapped SNP:35"
# in2 notin2
# in1     35      2
# notin1   5    105
# [1] 4.735097e-26
ovarlappairs_HUTCH_cis_eqtl_peer_pca=checkoverlappairs_sameplatform(qtlres1=highrisk_HUTCH_cis_eqtl,qtlres2=highrisk_HUTCH_old_cis_eqtl)
# [1] "number of pairs in qtlres1: 62"
# [1] "number of pairs in qtlres2: 64"
# [1] "total number of pairs: 72"
# [1] "number of common pairs: 54"


#TCGA_tumors eqtl
highrisk_TCGA_tumors_cis_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk_peer_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_POS.txt",fdrcutoff = 1)
qqplot(pvalue=10^-highrisk_TCGA_tumors_cis_eqtl_all$value)
#check some pairs manually:
computep()
highrisk_TCGA_tumors_cis_eqtl_all=addgenename_GE(highrisk_TCGA_tumors_cis_eqtl_all,anno=tcga_ge_anno)
length(unique(highrisk_TCGA_tumors_cis_eqtl_all$genename))
validatepairs_HUTCH_ciseqtl=validatepairs_difplatform()
validatepairs_HUTCH_ciseqtl=update_highrisk_snp_idx()
# [1] "number of qtl pairs to validate:62"
# [1] "number of qtl pairs can be validated:57"
# [1] "number of pairs validated:39"
# [1] "number of phenotype (probe) to validate:59"
# [1] "number of phenotype (probe) can be validated:55"
# [1] "number of phenotype (probe) validated:37"
# [1] "number of phenotype (gene) to validate:52"
# [1] "number of phenotype (gene) can be validated:49"
# [1] "number of phenotype (gene) validated:33"
# [1] "number of snp to validate:40"
# [1] "number of snp can be validate within pairs:39"
# [1] "number of snp validated:29"
write_qtl_input(validatepairs_HUTCH_ciseqtl,filename ="/fh/fast/stanford_j/Xiaoyu/QTL/result/validatepairs_HUTCH_ciseqtl.txt")
validatepairs_HUTCH_ciseqtl_fwer=validatepairs_difplatform(opt="fwer")
# [1] "number of qtl pairs to validate:62"
# [1] "number of qtl pairs can be validated:57"
# [1] "number of pairs validated:28"
# [1] "number of phenotype (probe) to validate:59"
# [1] "number of phenotype (probe) can be validated:55"
# [1] "number of phenotype (probe) validated:26"
# [1] "number of phenotype (gene) to validate:52"
# [1] "number of phenotype (gene) can be validated:49"
# [1] "number of phenotype (gene) validated:22"
# [1] "number of snp to validate:40"
# [1] "number of snp can be validate within pairs:39"
# [1] "number of snp validated:22"
validatepairs_HUTCH_ciseqtl_fwer=update_highrisk_snp_idx(dat=validatepairs_HUTCH_ciseqtl_fwer)
write_qtl_input(validatepairs_HUTCH_ciseqtl_fwer,filename ="/fh/fast/stanford_j/Xiaoyu/QTL/result/validatepairs_HUTCH_ciseqtl_fwer.txt")

#mqtl--------------
highrisk_HUTCH_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis",
                                   snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                   geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
# [1] "number of SNPs: 147"
# [1] "number of phenotypeprobes: 478998"
# [1] "number of pairs: 104610"
# [1] "number of unique SNPs: 147"
# [1] "number of unique phenotypes: 80644"
# [1] "FDR cutoff: 0.05"
# [1] "pvalue cutoff: 0.000567510638710477"
# [1] "number of pairs (after adjustment)1189"
# [1] "number of unique SNPs (after adjustment)116, 0.789115646258503"
# [1] "number of unique phenotypes (after adjustment)1136, 0.0023716174180268"
highrisk_HUTCH_cis_mqtl=addgenename_ME(highrisk_HUTCH_cis_mqtl)
length(unique(highrisk_HUTCH_cis_mqtl$genename))

checkoverlapsnps(qtlres1=highrisk_HUTCH_cis_eqtl,qtlres2=highrisk_HUTCH_cis_mqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
# in2 notin2
# in1     38     78
# notin1   2     29
# [1] 0.002741112
#Old result use gleason score and ME_pca
highrisk_HUTCH_old_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis_old",
                                       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                       geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
# [1] "pvalue cutoff: 0.000670786510527841"
# [1] "number of pairs (after adjustment)1405"
# [1] "number of unique SNPs (after adjustment)116, 0.789115646258503"
# [1] "number of unique phenotypes (after adjustment)1290, 0.00269312189194944"
highrisk_HUTCH_old_cis_mqtl=addgenename_ME(highrisk_HUTCH_old_cis_mqtl)
length(unique(highrisk_HUTCH_old_cis_mqtl$genename))
#old result, no gleason score
highrisk_HUTCH_old_nogleason_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis_old_nogleason",
                                                 snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                                 geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
# [1] "pvalue cutoff: 0.000676825062735074"
# [1] "number of pairs (after adjustment)1419"
# [1] "number of unique SNPs (after adjustment)117, 0.795918367346939"
# [1] "number of unique phenotypes (after adjustment)1306, 0.00272652495417517"
highrisk_HUTCH_old_nogleason_cis_mqtl=addgenename_ME(highrisk_HUTCH_old_nogleason_cis_mqtl)
length(unique(highrisk_HUTCH_old_nogleason_cis_mqtl$genename))
#old result, no gleason score no age
highrisk_HUTCH_old_nogleasonage_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis_old_nogleasonage",
                                                 snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                                 geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
# [1] "pvalue cutoff: 0.000672238306805779"
# [1] "number of pairs (after adjustment)1410"
# [1] "number of unique SNPs (after adjustment)117, 0.795918367346939"
# [1] "number of unique phenotypes (after adjustment)1296, 0.00270564804028409"
highrisk_HUTCH_old_nogleasonage_cis_mqtl=addgenename_ME(highrisk_HUTCH_old_nogleasonage_cis_mqtl)
length(unique(highrisk_HUTCH_old_nogleasonage_cis_mqtl$genename))
highrisk_HUTCH_old_nogleasonage_add3snppc_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis_old_nogleasonage_add3snppc",
                                                    snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                                    geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
# [1] "pvalue cutoff: 0.000577958582165979"
# [1] "number of pairs (after adjustment)1211"
# [1] "number of unique SNPs (after adjustment)119, 0.80952380952381"
# [1] "number of unique phenotypes (after adjustment)1151, 0.00240293278886342"
highrisk_HUTCH_old_nogleasonage_add3snppc_cis_mqtl=addgenename_ME(highrisk_HUTCH_old_nogleasonage_add3snppc_cis_mqtl)
length(unique(highrisk_HUTCH_old_nogleasonage_add3snppc_cis_mqtl$genename))
#compare with old result
checkoverlapsnps(qtlres1=highrisk_HUTCH_cis_mqtl,qtlres2=highrisk_HUTCH_old_cis_mqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
# [1] "number of unique SNP1:115"
# [1] "number of unique SNP2:116"
# [1] "number of overlapped SNP:113"
# in2 notin2
# in1    113      3
# notin1   2     29
# [1] 5.261159e-25
checkoverlapsnps(qtlres1=highrisk_HUTCH_old_nogleasonage_cis_mqtl,qtlres2=highrisk_HUTCH_old_nogleasonage_add3snppc_cis_mqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")

ovarlappairs_HUTCH_cis_mqtl=checkoverlappairs_sameplatform(qtlres1=highrisk_HUTCH_cis_mqtl,qtlres2=highrisk_HUTCH_old_cis_mqtl)
# [1] "number of pairs in qtlres1: 1179"
# [1] "number of pairs in qtlres2: 1405"
# [1] "total number of pairs: 1518"
# [1] "number of common pairs: 1066"

ovarlappairs_HUTCH_cis_mqtl=checkoverlappairs_sameplatform(qtlres1=highrisk_HUTCH_old_nogleasonage_cis_mqtl,
                                                           qtlres2=highrisk_HUTCH_old_nogleasonage_add3snppc_cis_mqtl)
# [1] "number of pairs in qtlres1: 1410"
# [1] "number of pairs in qtlres2: 1211"
# [1] "total number of pairs: 1485"
# [1] "number of common pairs: 1136"

#TCGA_tumors mqtl
highrisk_TCGA_tumors_cis_mqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/mqtl_highrisk_cis",
                                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_POS.txt",fdrcutoff = 1)
qqplot(pvalue=10^-highrisk_TCGA_tumors_cis_mqtl_all$value)
#check some pairs manually:
computep(snpidx=c(73,81),phenotypeidx=c(12142,8562),
         snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_ME_updatename.txt",
         phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME.txt",
         covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_ME_used.txt")

idx=which(10^-highrisk_TCGA_tumors_cis_mqtl_all$value<=0.05)
highrisk_TCGA_tumors_cis_mqtl_all=highrisk_TCGA_tumors_cis_mqtl_all[idx,]
highrisk_TCGA_tumors_cis_mqtl_all=addgenename_ME(highrisk_TCGA_tumors_cis_mqtl_all)
length(unique(highrisk_TCGA_tumors_cis_mqtl_all$genename))
validatepairs_HUTCH_cismqtl=validatepairs(qtlres1=highrisk_HUTCH_cis_mqtl,qtlres2=highrisk_TCGA_tumors_cis_mqtl_all,phenotype2=tcga_me_pos)
# [1] "number of qtl pairs to validate:1189"
# [1] "number of qtl pairs can be validated:731"
# [1] "number of pairs validated:579"
# [1] "number of phenotype to validate:1136"
# [1] "number of phenotype can be validated:710"
# [1] "number of phenotype validated:563"
# [1] "number of snp to validate:116"
# [1] "number of snp can be validated within pairs:101"
# [1] "number of snp validated:86"
validatepairs_HUTCH_cismqtl=update_highrisk_snp_idx(dat=validatepairs_HUTCH_cismqtl)
write_qtl_input(validatepairs_HUTCH_cismqtl,filename ="/fh/fast/stanford_j/Xiaoyu/QTL/result/validatepairs_HUTCH_cismqtl.txt")
validatepairs_HUTCH_cismqtl_fwer=validatepairs(qtlres1=highrisk_HUTCH_cis_mqtl,qtlres2=highrisk_TCGA_tumors_cis_mqtl_all,phenotype2=tcga_me_pos,opt="fwer")
# [1] "number of qtl pairs to validate:1189"
# [1] "number of qtl pairs can be validated:731"
# [1] "number of pairs validated:385"
# [1] "number of phenotype to validate:1136"
# [1] "number of phenotype can be validated:710"
# [1] "number of phenotype validated:563"
# [1] "number of snp to validate:116"
# [1] "number of snp can be validated within pairs:101"
# [1] "number of snp validated:70"

overlap_hutch_highrisk_cis_eqtl_hutch_cis_mqtl=overlappairs(dat1=highrisk_HUTCH_cis_eqtl,dat2=highrisk_HUTCH_cis_mqtl)
# [1] "number of overlapped SNP:38"
# [1] "number of phenoytype1:57"
# [1] "number of phenoytype2:834"
# [1] "number of triplets:1949"
validtriplet_HUTCH_ciseqtl_cismqtl=validtriplet_difplatform(dat1=overlap_hutch_highrisk_cis_eqtl_hutch_cis_mqtl$triplet,validateeqtl=highrisk_TCGA_tumors_cis_eqtl_all,
                      validatemqtl=highrisk_TCGA_tumors_cis_mqtl_all,phenotypemap=hutch_ge_anno,phenotype2=tcga_me_pos)
# [1] "number of triplets to validate:1949"
# [1] "number of triplets can be validated considering eqtl genename mapping:1606"
# [1] "number of triplets can be validated:903"
# [1] "number of triplets validated:578"
# [1] "number of snp can be validate:36"
# [1] "number of snp validated:26"
validtriplet_HUTCH_ciseqtl_cismqtl=update_highrisk_snp_idx(validtriplet_HUTCH_ciseqtl_cismqtl)
write_qtl_input(validtriplet_HUTCH_ciseqtl_cismqtl,filename ="/fh/fast/stanford_j/Xiaoyu/QTL/result/validtriplet_HUTCH_ciseqtl_cismqtl.txt")

#Trans
highrisk_HUTCH_trans_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_trans",
                                   snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                   geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
# [1] "number of SNPs: 147"
# [1] "number of phenotypeprobes: 25019"
# Read 3673457 rows and 6 (of 6) columns from 0.357 GB file in 00:00:04
# [1] "number of pairs: 3673457"
# [1] "number of unique SNPs: 147"
# [1] "number of unique phenotypes: 25019"
# [1] "FDR cutoff: 0.05"
# [1] "pvalue cutoff: 1.88004842249937e-07"
# [1] "number of pairs (after adjustment)16"
# [1] "number of unique SNPs (after adjustment)5, 0.0340136054421769"
# [1] "number of unique phenotypes (after adjustment)16, 0.000639513969383269"
highrisk_HUTCH_trans_eqtl=addgenename_GE(highrisk_HUTCH_trans_eqtl,anno=hutch_ge_anno)

highrisk_HUTCH_trans_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_trans",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",fdrcutoff = 1)
#idx=which(10^-highrisk_HUTCH_trans_eqtl_all$value<=0.05)
#highrisk_HUTCH_trans_eqtl_all=highrisk_HUTCH_trans_eqtl_all[idx,]
highrisk_HUTCH_trans_eqtl_all=addgenename_GE(highrisk_HUTCH_trans_eqtl_all,anno=hutch_ge_anno)
highrisk_HUTCH_trans_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
                                   snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                   geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
# [1] "number of SNPs: 147"
# [1] "number of phenotypeprobes: 478998"
# Read 7035745 rows and 6 (of 6) columns from 0.675 GB file in 00:00:09
# [1] "number of pairs: 7035745"
# [1] "number of unique SNPs: 147"
# [1] "number of unique phenotypes: 478998"
# [1] "FDR cutoff: 0.05"
# [1] "pvalue cutoff: 1.80397193943627e-06"
# [1] "number of pairs (after adjustment)2537"
# [1] "number of unique SNPs (after adjustment)92, 0.625850340136054"
# [1] "number of unique phenotypes (after adjustment)2528, 0.00527768383166527"
highrisk_HUTCH_trans_mqtl=addgenename_ME(highrisk_HUTCH_trans_mqtl)

checkoverlapsnps(qtlres1=highrisk_HUTCH_trans_eqtl,qtlres2=highrisk_HUTCH_trans_mqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")

highrisk_TCGA_tumors_trans_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk_peer_trans",
                                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_POS.txt",fdrcutoff = 1)
qqplot(pvalue=10^-highrisk_TCGA_tumors_trans_eqtl_all$value)
highrisk_TCGA_tumors_trans_eqtl_all=addgenename_GE(highrisk_TCGA_tumors_trans_eqtl_all,anno=tcga_ge_anno)
length(unique(highrisk_TCGA_tumors_trans_eqtl_all$genename))
validatepairs_HUTCH_transeqtl=validatepairs_difplatform(qtlres1=highrisk_HUTCH_trans_eqtl,qtlres2=highrisk_TCGA_tumors_trans_eqtl_all)
# [1] "number of qtl pairs to validate:16"
# [1] "number of qtl pairs can be validated:15"
# [1] "number of pairs validated:0"
# [1] "number of phenotype (probe) to validate:16"
# [1] "number of phenotype (probe) can be validated:15"
# [1] "number of phenotype (probe) validated:0"
# [1] "number of phenotype (gene) to validate:16"
# [1] "number of phenotype (gene) can be validated:15"
# [1] "number of phenotype (gene) validated:0"
# [1] "number of snp to validate:5"
# [1] "number of snp can be validate within pairs:5"
# [1] "number of snp validated:0"
highrisk_TCGA_tumors_trans_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk_peer_trans",
                                               snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                               geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_POS.txt")
highrisk_TCGA_tumors_trans_eqtl=addgenename_GE(highrisk_TCGA_tumors_trans_eqtl,anno=tcga_ge_anno)
#use TCGA as discovery set
validatepairs_HUTCH_transeqtl_1=validatepairs_difplatform(qtlres1=highrisk_TCGA_tumors_trans_eqtl,qtlres2=highrisk_HUTCH_trans_eqtl_all,phenotypemap = tcga_ge_anno)
validatepairs_HUTCH_transeqtl_1=update_highrisk_snp_idx(validatepairs_HUTCH_transeqtl_1)

highrisk_TCGA_tumors_trans_mqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/mqtl_highrisk_trans",
                                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_POS.txt",fdrcutoff = 1)
qqplot(pvalue=10^-highrisk_TCGA_tumors_trans_mqtl_all$value)
idx=which(10^-highrisk_TCGA_tumors_trans_mqtl_all$value<=0.05)
highrisk_TCGA_tumors_trans_mqtl_all=highrisk_TCGA_tumors_trans_mqtl_all[idx,]
#highrisk_TCGA_tumors_trans_mqtl_all=addgenename_ME(highrisk_TCGA_tumors_trans_mqtl_all)

validatepairs_HUTCH_transmqtl=validatepairs(qtlres1=highrisk_HUTCH_trans_mqtl,qtlres2=highrisk_TCGA_tumors_trans_mqtl_all,phenotype2=tcga_me_pos)
# [1] "number of qtl pairs to validate:2537"
# [1] "number of qtl pairs can be validated:2188"
# [1] "number of pairs validated:82"
# [1] "number of phenotype to validate:2528"
# [1] "number of phenotype can be validated:2179"
# [1] "number of phenotype validated:90"
# [1] "number of snp to validate:92"
# [1] "number of snp can be validated within pairs:75"
# [1] "number of snp validated:18"
validatepairs_HUTCH_transmqtl_fwer=validatepairs(qtlres1=highrisk_HUTCH_trans_mqtl,qtlres2=highrisk_TCGA_tumors_trans_mqtl_all,phenotype2=tcga_me_pos,opt="fwer")
# [1] "number of qtl pairs to validate:2537"
# [1] "number of qtl pairs can be validated:2188"
# [1] "number of pairs validated:11"
# [1] "number of phenotype to validate:2528"
# [1] "number of phenotype can be validated:2179"
# [1] "number of phenotype validated:90"
# [1] "number of snp to validate:92"
# [1] "number of snp can be validated within pairs:75"
# [1] "number of snp validated:8"
validatepairs_HUTCH_transmqtl_fwer=update_highrisk_snp_idx(validatepairs_HUTCH_transmqtl_fwer)

#Mediation analysis
Sys.time()
HUTCH_med_ciseqtl_cismqtl=citqtlpairs_inputfromfiles()
Sys.time()
HUTCH_med_cismqtl_ciseqtl=citqtlpairs_inputfromfiles(dat1=highrisk_HUTCH_cis_mqtl,dat2=highrisk_HUTCH_cis_eqtl)
Sys.time()
HUTCH_med_ciseqtl_cismqtl=HUTCH_med_ciseqtl_cismqtl[order(HUTCH_med_ciseqtl_cismqtl$snp_idx,HUTCH_med_ciseqtl_cismqtl$gene.x,HUTCH_med_ciseqtl_cismqtl$gene.y),]
HUTCH_med_cismqtl_ciseqtl=HUTCH_med_cismqtl_ciseqtl[order(HUTCH_med_cismqtl_ciseqtl$snp_idx,HUTCH_med_cismqtl_ciseqtl$gene.y,HUTCH_med_cismqtl_ciseqtl$gene.x),]
save(HUTCH_med_ciseqtl_cismqtl,HUTCH_med_cismqtl_ciseqtl,file="med_tmp.RData")
load("../result/tmp_HUTCH_med_ciseqtl_cismqtl.RData")
load("../result/tmp_HUTCH_med_cismqtl_ciseqtl.RData")
HUTCH_med_ciseqtl_cismqtl=HUTCH_med_ciseqtl_cismqtl[order(HUTCH_med_ciseqtl_cismqtl$snp_idx,HUTCH_med_ciseqtl_cismqtl$gene.x,HUTCH_med_ciseqtl_cismqtl$gene.y),]
HUTCH_med_cismqtl_ciseqtl=HUTCH_med_cismqtl_ciseqtl[order(HUTCH_med_cismqtl_ciseqtl$snp_idx,HUTCH_med_cismqtl_ciseqtl$gene.y,HUTCH_med_cismqtl_ciseqtl$gene.x),]
sum(HUTCH_med_ciseqtl_cismqtl$gene.x==HUTCH_med_cismqtl_ciseqtl$gene.y)
sum(HUTCH_med_ciseqtl_cismqtl$gene.y==HUTCH_med_cismqtl_ciseqtl$gene.x)
sum(HUTCH_med_ciseqtl_cismqtl$q.cit<=0.05)
sum(HUTCH_med_cismqtl_ciseqtl$q.cit<=0.05)
sum(HUTCH_med_ciseqtl_cismqtl$q.cit<=0.05 & HUTCH_med_cismqtl_ciseqtl$q.cit<=0.05)
HUTCH_triplet=data.frame(snp_idx=HUTCH_med_ciseqtl_cismqtl$snp_idx,geneexp=HUTCH_med_ciseqtl_cismqtl$gene.x,methylation=HUTCH_med_ciseqtl_cismqtl$gene.y)
HUTCH_triplet$validateinTCGA=NA
for (i in 1:nrow(HUTCH_triplet))
{
  idx=which(validtriplet_HUTCH_ciseqtl_cismqtl$snp_idx==HUTCH_triplet$snp_idx[i] & validtriplet_HUTCH_ciseqtl_cismqtl$gene.x==HUTCH_triplet$geneexp[i]
            & validtriplet_HUTCH_ciseqtl_cismqtl$gene.y==HUTCH_triplet$methylation[i])
  if (length(idx)>0)
  {
    HUTCH_triplet$validateinTCGA[i]=validtriplet_HUTCH_ciseqtl_cismqtl$validate[idx]
  }
}
sum(HUTCH_med_ciseqtl_cismqtl$gene.x==HUTCH_triplet$geneexp)
sum(HUTCH_med_cismqtl_ciseqtl$gene.x==HUTCH_triplet$methylation)
HUTCH_triplet=cbind.data.frame(HUTCH_triplet,HUTCH_med_ciseqtl_cismqtl[,26:35])
colnames(HUTCH_triplet)[5:ncol(HUTCH_triplet)]=paste0("meth2gen_",colnames(HUTCH_triplet)[5:ncol(HUTCH_triplet)])
HUTCH_triplet=cbind.data.frame(HUTCH_triplet,HUTCH_med_cismqtl_ciseqtl[,26:35])
sum(HUTCH_triplet$meth2gen_q.cit<=0.05) #449,new 442
sum(HUTCH_triplet$meth2gen_q.cit<=0.05 & HUTCH_triplet$q.cit>0.05) #118 meth2gen
sum(HUTCH_triplet$q.cit<=0.05) #428, new 429
sum(HUTCH_triplet$q.cit<=0.05 & HUTCH_triplet$meth2gen_q.cit>0.05) #97 gen2meth
HUTCH_triplet$medclass=NA
HUTCH_triplet$medclass[HUTCH_triplet$meth2gen_q.cit<=0.05 & HUTCH_triplet$q.cit>0.05]="meth2gen"
HUTCH_triplet$medclass[HUTCH_triplet$q.cit<=0.05 & HUTCH_triplet$meth2gen_q.cit>0.05]="gen2meth"
idx=which(HUTCH_triplet$meth2gen_q.cit<=0.05 & HUTCH_triplet$q.cit<=0.05) #331, new 330
View(HUTCH_triplet[idx,])
sum(HUTCH_triplet$meth2gen_q.cit[idx]<HUTCH_triplet$q.cit[idx]) #144 meth2gen
sum(HUTCH_triplet$meth2gen_q.cit[idx]>HUTCH_triplet$q.cit[idx]) #187 gen2meth
HUTCH_triplet$medclass[HUTCH_triplet$meth2gen_q.cit<=0.05 & HUTCH_triplet$q.cit<=0.05 & HUTCH_triplet$meth2gen_q.cit<HUTCH_triplet$q.cit]="ambiguous_meth2gen"
HUTCH_triplet$medclass[HUTCH_triplet$meth2gen_q.cit<=0.05 & HUTCH_triplet$q.cit<=0.05 & HUTCH_triplet$meth2gen_q.cit>HUTCH_triplet$q.cit]="ambiguous_gen2meth"
table(HUTCH_triplet$medclass)

# ambiguous_gen2meth ambiguous_meth2gen           gen2meth           meth2gen 
# 164                166                 99                112 
#validation
TCGA_med_ciseqtl_cismqtl=computep_meds_inputfromfiles()
TCGA_med_cismqtl_ciseqtl=computep_meds_inputfromfiles(opt="gene2meth")
save(TCGA_med_ciseqtl_cismqtl,TCGA_med_cismqtl_ciseqtl,file="../result/tmp_TCGA_med_ciseqtl_cismqtl.RData")
HUTCH_triplet$TCGA_meds_meth2gen=HUTCH_triplet$TCGA_meds_gen2meth=NA
for (i in 1:nrow(HUTCH_triplet))
{
  idx=which(TCGA_med_ciseqtl_cismqtl$snp_idx==HUTCH_triplet$snp_idx[i] & TCGA_med_ciseqtl_cismqtl$gene.x==HUTCH_triplet$geneexp[i] & TCGA_med_ciseqtl_cismqtl$gene.y==HUTCH_triplet$methylation[i])
  if (length(idx)>0)
  {
    HUTCH_triplet$TCGA_meds_meth2gen[i]=TCGA_med_ciseqtl_cismqtl$valid_p_cit[idx]
  }
  idx=which(TCGA_med_cismqtl_ciseqtl$snp_idx==HUTCH_triplet$snp_idx[i] & TCGA_med_cismqtl_ciseqtl$gene.x==HUTCH_triplet$geneexp[i] & TCGA_med_cismqtl_ciseqtl$gene.y==HUTCH_triplet$methylation[i])
  if (length(idx)>0)
  {
    HUTCH_triplet$TCGA_meds_gen2meth[i]=TCGA_med_cismqtl_ciseqtl$valid_p_cit[idx]
  }
}
idx=which(HUTCH_triplet$TCGA_meds_gen2meth<=0.05 &HUTCH_triplet$medclass=="gen2meth")
meds_gen2meth=HUTCH_triplet[idx,]
idx=which(meds_gen2meth$TCGA_meds_meth2gen>0.05)
meds_gen2meth=meds_gen2meth[idx,]
meds_gen2meth=meds_gen2meth[order(meds_gen2meth$q.cit),]
meds_gen2meth$gene=meds_gen2meth$geneexp
meds_gen2meth=addgenename_GE(datatable = meds_gen2meth,anno=hutch_ge_anno,colname="GE_genename")
meds_gen2meth$gene=meds_gen2meth$methylation
meds_gen2meth=addgenename_ME(datatable = meds_gen2meth,colname="ME_genename")
meds_gen2meth=update_highrisk_snp_idx(meds_gen2meth)
meds_gen2meth=meds_gen2meth[,c("snp","geneexp","methylation","GE_genename","ME_genename")]
colnames(meds_gen2meth)[4:5]=c("geneexp_genename","methylation_genename")
idx=which(HUTCH_triplet$TCGA_meds_meth2gen<=0.05 &HUTCH_triplet$medclass=="meth2gen")
meds_meth2gen=HUTCH_triplet[idx,]
write.table(meds_gen2meth,file="../result/mediation_gen2meth.txt",col.names = T,row.names = F,sep="\t",quote=F)
idx=which(meds_meth2gen$TCGA_meds_gen2meth>0.05)
if (length(idx)>0) meds_meth2gen=meds_meth2gen[idx,]
meds_meth2gen=meds_meth2gen[order(meds_meth2gen$meth2gen_q.cit),]
meds_meth2gen$gene=meds_meth2gen$geneexp
meds_meth2gen=addgenename_GE(datatable = meds_meth2gen,anno=hutch_ge_anno,colname="GE_genename")
meds_meth2gen$gene=meds_meth2gen$methylation
meds_meth2gen=addgenename_ME(datatable = meds_meth2gen,colname="ME_genename")
meds_meth2gen=update_highrisk_snp_idx(meds_meth2gen)
meds_meth2gen=meds_meth2gen[,c("snp","geneexp","methylation","GE_genename","ME_genename")]
colnames(meds_meth2gen)[4:5]=c("geneexp_genename","methylation_genename")
write.table(meds_meth2gen,file="../result/mediation_meth2gen.txt",col.names = T,row.names = F,sep="\t",quote=F)
#JAN5
load("../data/GTEx/gtex_ge_anno.RData")
highrisk_gtex_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/gtex/eqtl_highrisk_cis",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE_POS.txt")
highrisk_gtex_cis_eqtl=addgenename_GE(highrisk_gtex_cis_eqtl,anno=gtex_ge_anno)
highrisk_gtex_cis_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/gtex/eqtl_highrisk_cis",
                                       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                       geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE_POS.txt",fdrcutoff = 1)
highrisk_gtex_cis_eqtl_all=addgenename_GE(highrisk_gtex_cis_eqtl_all,anno=gtex_ge_anno)
#find the overlapped pairs
# tmp=merge(highrisk_HUTCH_cis_eqtl_all,highrisk_gtex_cis_eqtl_all,by=c("chr","opos1","opos2"))
# nrow(tmp)
# #[1] 7 complete match
# pvalue_diff_HUTCH_normal_eqtl=diff_tumror_normal_eqtl()
# length(pvalue_diff_HUTCH_normal_eqtl)
# #[1] 2221 eqtl-ge pairs
# hist(pvalue_diff_HUTCH_normal_eqtl,main="",xlab="p-value",breaks = 20,cex.lab=1.3,cex.axis=1.3)
# quantile(pvalue_diff_HUTCH_normal_eqtl)
# # 0%        25%        50%        75%       100% 
# # 0.05690559 0.45764924 0.62707903 0.81280973 1.00000000
# qqplot(pvalue_diff_HUTCH_normal_eqtl)
# pvalue_diff_TCGA_tumors_normal_eqtl=diff_tumror_normal_eqtl(dat1 = highrisk_TCGA_tumors_cis_eqtl_all)
# length(pvalue_diff_TCGA_tumors_normal_eqtl)
# #[1] 2858
# quantile(pvalue_diff_TCGA_tumors_normal_eqtl)
# # 0%        25%        50%        75%       100% 
# # 0.05595134 0.46699256 0.63761736 0.82720337 1.00000000 
# qqplot(pvalue_diff_TCGA_tumors_normal_eqtl)

#removed duplicates transcripts in a gene
pvalue_diff_tumor_normal_eqtl=diff_tumror_normal_eqtl()
save(pvalue_diff_tumor_normal_eqtl,file="../result/tmp_pvalue_diff_tumor_normal_eqtl.RData")
#pvalue_diff_tumor_normal_eqtl=update_highrisk_snp_idx(pvalue_diff_tumor_normal_eqtl)
length(pvalue_diff_tumor_normal_eqtl$pvalue)
#[1] 2494
postscript(file="../result/qqplot_diff_tumor_normal_eqtl.ps")
par(mar=c(5.1,5.1,4.1,2.1))
qqplot(pvalue_diff_tumor_normal_eqtl$pvalue_diff)
dev.off()
which.min(pvalue_diff_tumor_normal_eqtl$pvalue_diff)
#[1] 1
pvalue_diff_tumor_normal_eqtl[1,]
sum(pvalue_diff_tumor_normal_eqtl$pvalue_diff<=0.05/length(pvalue_diff_tumor_normal_eqtl$pvalue_diff))
#[1] 3
pvalue_diff_tumor_normal_eqtl[1:3,c("snp_idx","genename","pvalue","pvalue2")]
pvalue_diff_tumor_normal_eqtl$pvalue_diff[1:3]
#check TCGA validated eQTLs in normals
# validate_TCGA_highrisk_cis_eqtl_in_gtex=merge(validatepairs_HUTCH_ciseqtl,highrisk_gtex_cis_eqtl_all,by=c("snp_idx","genename"))
# idx=which(10^-validate_TCGA_highrisk_cis_eqtl_in_gtex$value<=0.05/nrow(validate_TCGA_highrisk_cis_eqtl_in_gtex))
# length(unique(validate_TCGA_highrisk_cis_eqtl_in_gtex$snp_idx[idx]))
validate_TCGA_highrisk_ciseqtl_pear_in_gtex=merge(validatepairs_HUTCH_ciseqtl_fwer,highrisk_gtex_cis_eqtl_all,by=c("snp_idx","genename"))
idx=which(10^-validate_TCGA_highrisk_ciseqtl_pear_in_gtex$value<=0.05/nrow(validate_TCGA_highrisk_ciseqtl_pear_in_gtex)
          &validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta.x*validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta>0)
length(unique(validatepairs_HUTCH_ciseqtl_fwer$snp_idx))
length(unique(validate_TCGA_highrisk_ciseqtl_pear_in_gtex$snp_idx[idx]))
#use p-value<0.05
idx=which(10^-validate_TCGA_highrisk_ciseqtl_pear_in_gtex$value<=0.05
          &validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta.x*validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta>0)
length(unique(validate_TCGA_highrisk_ciseqtl_pear_in_gtex$snp_idx[idx]))
dat_test=merge(pvalue_diff_tumor_normal_eqtl,validatepairs_HUTCH_ciseqtl_fwer)
test=merge(pvalue_diff_tumor_normal_eqtl,validatepairs_HUTCH_ciseqtl_fwer,all.y=T)
#use p-value<0.1
idx=which(10^-validate_TCGA_highrisk_ciseqtl_pear_in_gtex$value<=0.1
          &validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta.x*validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta>0)
length(unique(validate_TCGA_highrisk_ciseqtl_pear_in_gtex$snp_idx[idx]))
#only consider direction
idx=which(validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta.x*validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta>0)
length(unique(validate_TCGA_highrisk_ciseqtl_pear_in_gtex$snp_idx[idx]))
idx=which(validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta.x*validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta<0)
validate_TCGA_highrisk_ciseqtl_pear_in_gtex$value[idx]

#MAR19, check betta difference
source("functions.R")
#for highrisk SNPs
highrisk_tbd_cis_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_highrisk_peer_cis_all",
                        snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                        geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",fdrcutoff = 1)

highrisk_hutch_cis_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_highrisk_peer_cis",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_POS.txt",fdrcutoff = 1)
diff_highrisk_hutch_tbd=diff_2pairs_tumror_normal_eqtl()
#diff_highrisk_hutch_tbd$value.x: -log10(pvalue) of hutch eqtl result
#diff_highrisk_hutch_tbd$value.y: -log10(pvalue) of mayo eqtl result
#pvalue is for the beta difference

qqplot(diff_highrisk_hutch_tbd$pvalue,ylim=c(0,10))
nrow(diff_highrisk_hutch_tbd) #[1] 2154
sum(diff_highrisk_hutch_tbd$pvalue<0.05/nrow(diff_highrisk_hutch_tbd)) #[1] 25
idx=which(diff_highrisk_hutch_tbd$pvalue<0.05/nrow(diff_highrisk_hutch_tbd))
sum(length(unique(diff_highrisk_hutch_tbd$gene[idx]))) #25 genes
sum(length(unique(diff_highrisk_hutch_tbd$pos1[idx]))) #21 snps
#for genome wide snps
hutch_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_cis (copy)",
                        snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                        geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_POS.txt",fdrcutoff = 1)

tbd_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_cis (copy)",
                        snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_POS.txt",
                        geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",fdrcutoff = 1)
diff_hutch_tbd=diff_2pairs_tumror_normal_eqtl(dat1=hutch_cis_eqtl,dat2=tbd_cis_eqtl)

nrow(diff_hutch_tbd) #[1] 172127
sum(diff_hutch_tbd$pvalue<0.05/nrow(diff_hutch_tbd)) #[1] 48877
idx=which(diff_hutch_tbd$pvalue<0.05/nrow(diff_hutch_tbd))
sum(length(unique(diff_hutch_tbd$gene[idx]))) #590 genes
sum(length(unique(diff_hutch_tbd$pos1[idx]))) #45998 snps
diff_hutch_tbd[1:3,c(2,5,6,11,12,13,18,19,20,21,22,23)]
#       chr    gene  value.x   tstat.x     beta.x  value.y   tstat.y    beta.y       var.x       var.y    tscore
# 137903  22 FAM118A 72.09058 -23.49460 -1.4718698 202.9762  55.51274  1.904622 0.003924666 0.001177151 -47.27192
# 137904  22 FAM118A 72.32216  23.55460  1.4783437 201.0987 -54.90545 -1.883243 0.003939120 0.001176473  46.99984
# 82158   11  MRPL21 41.06566 -15.64257 -0.6836446 160.3839  42.85082  1.220973 0.001910049 0.000811884 -36.50641
#        pvalue
# 137903  0.000000e+00
# 137904  0.000000e+00
# 82158  8.775074e-292
diff_highrisk_hutch_tbd[1:3,c(2,5,6,11,12,13,18,19,20,21,22,23)]



hutch_gene=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene.txt",header=T)
tbd_gene=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE.txt",header=T)
hutch_snp=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",header=T)
hutch_snppos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",header=T)
tbd_snp=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_GE.txt",header=T)
tbd_snppos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_POS.txt",header=T)
testid=3
idx1=which(hutch_gene$id==diff_hutch_tbd$gene[testid])
idx2=which(tbd_gene$id==diff_hutch_tbd$gene[testid])
tmp=data.frame(normal=unlist(tbd_gene[idx2,2:ncol(tbd_gene)]))
tmp$tumor=NA
tmp$tumor[1:(ncol(hutch_gene)-1)]=unlist(hutch_gene[idx1,2:ncol(hutch_gene)])
boxplot(tmp)

boxplot(unlist(hutch_snp[diff_hutch_tbd$snp_idx.x[testid],2:ncol(hutch_snp)]))
boxplot(unlist(tbd_snp[diff_hutch_tbd$snp_idx.y[testid],2:ncol(tbd_snp)]))
snp_pvalue=rep(NA,nrow(diff_hutch_tbd))
for (i in 1:nrow(diff_hutch_tbd))
{
  if (i %%10000==0) cat(i,"..")
  snp_pvalue[i]=t.test(unlist(hutch_snp[diff_hutch_tbd$snp_idx.x[i],2:ncol(hutch_snp)]),
                       unlist(tbd_snp[diff_hutch_tbd$snp_idx.y[i],2:ncol(tbd_snp)]))$p.value
}
png("qqplot_topgenotype.png")
qqplot(snp_pvalue,main="genotype (p<1e-4)")
dev.off()
diff_chrs_hutch_tbd=NULL
for (chr in 7:23)
{
  cat(chr,"..")
  hutch_cis_eqtl_chr=readqtlres(qtlresfile=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_cis.chr",chr),
                            snpposfile=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/chr/HUTCH_SNP_POS.txt.chr",chr),
                            geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_POS.txt",fdrcutoff = 1)
  
  tbd_cis_eqtl_chr=readqtlres(qtlresfile=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_cis.chr",chr),
                          snpposfile=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/chr/TBD_SNP_POS.txt.chr",chr),
                          geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",fdrcutoff = 1)
  tmp=diff_2pairs_tumror_normal_eqtl(dat1=hutch_cis_eqtl_chr,dat2=tbd_cis_eqtl_chr)
  diff_chrs_hutch_tbd=rbind(diff_chrs_hutch_tbd,tmp)
}

qqplot(diff_chrs_hutch_tbd$pvalue,main="All pairs")

#compare hutch and tbd data
snppos1=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",header=T)
snppos2=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_POS.txt",header=T)
genepos1=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_POS.txt")
genepos2=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt")
tmp=merge(snppos1,snppos2,by=c("chr","pos"))
nrow(tmp) #[1] 6842398
nrow(snppos1) #[1] 6835832
nrow(snppos2) #[1] 7262931
nrow(genepos1) #[1] 18077
nrow(genepos2) #[1] 16898


