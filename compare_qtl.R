#!/usr/bin/env Rscript
library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library(sas7bdat)
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

chrs=c(1:22,"X","Y")
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
totalpos=function(chr=NULL,pos=NULL)
{
  chr=gsub(23,"X",chr)
  chr=gsub(24,"Y",chr)
  chrs=c(1:22,"X","Y")
  
  idx=which(chr %in% chrs)
  chr=chr[idx]
  pos=pos[idx]
  #chr=factor(chr,chrs)
  
  allchrs=unique(chr)
  posall=rep(NA,length(chr))
  for (mychr in allchrs)
  {
    idx=which(chr %in% mychr)
    idx1=which(names(chrstart)==mychr)
    posall[idx]=pos[idx]+chrstart[idx1]
  }
  return(posall)
}

#extract qtl result
formlevelmat=function(qtlresfile="eqtl_cis",colpvalue=5,snpposfile="eqtl_POS_SNP_GE.txt",geposfile="eqtl_POS_GE.txt",bonfadj=NULL)
{
  print(paste0("processed qtl file: ",qtlresfile))
  snpspos=fread(snpposfile,sep="\t")
  print(paste0("number of SNPs: ",nrow(snpspos)))
  snpspos$chr=gsub(23,"X",snpspos$chr)
  snpspos$chr=gsub(24,"Y",snpspos$chr)
  
  genepos=fread(geposfile,sep="\t")
  print(paste0("number of phenotypeprobes: ",nrow(genepos)))
  genepos$chr=gsub(23,"X", genepos$chr)
  genepos$chr=gsub(24,"Y", genepos$chr)
  
  eqtl_cis=fread(qtlresfile,header=T,sep="\t",fill=T,check.names = F)
  eqtl_cis=as.data.frame(eqtl_cis)
  if (! is.null(bonfadj))
  {
    print(paste0("number of pairs: ",nrow(eqtl_cis)))
    print(paste0("number of unique SNPs: ",length(unique(eqtl_cis$SNP))))
    print(paste0("number of unique phenotypes: ",length(unique(eqtl_cis$gene))))
    print(paste0("Bonferroni cutoff: ",0.05/bonfadj))
    idx=which(eqtl_cis[,colpvalue]<=0.05/bonfadj)
    print(paste0("number of pairs (after adjustment)",length(idx)))
    print(paste0("number of unique SNPs (after adjustment)",length(unique(eqtl_cis$SNP[idx])),", ",length(unique(eqtl_cis$SNP[idx]))/nrow(snpspos)))
    print(paste0("number of unique phenotypes (after adjustment)",length(unique(eqtl_cis$gene[idx])),", ",length(unique(eqtl_cis$gene[idx]))/nrow(genepos)))
  }else
  {
    print(paste0("number of pairs: ",nrow(eqtl_cis)))
    print(paste0("number of unique SNPs: ",length(unique(eqtl_cis$SNP))))
    print(paste0("number of unique phenotypes: ",length(unique(eqtl_cis$gene))))
    fdrcutoff=1
    print(paste0("FDR cutoff: ",fdrcutoff))
    tmp=which(eqtl_cis$FDR<=fdrcutoff)
    print(paste0("pvalue cutoff: ",eqtl_cis[max(tmp),colpvalue]))
    idx=which(eqtl_cis$FDR<=fdrcutoff)
    print(paste0("number of pairs (after adjustment)",length(idx)))
    print(paste0("number of unique SNPs (after adjustment)",length(unique(eqtl_cis$SNP[idx])),", ",length(unique(eqtl_cis$SNP[idx]))/nrow(snpspos)))
    print(paste0("number of unique phenotypes (after adjustment)",length(unique(eqtl_cis$gene[idx])),", ",length(unique(eqtl_cis$gene[idx]))/nrow(genepos)))
  }
  
  eqtl_cis=eqtl_cis[idx,]
  
  idx=match(eqtl_cis$SNP,snpspos$snp)
  eqtl_cis$SNP_chr=snpspos$chr[idx]
  eqtl_cis$SNP_pos=snpspos$pos[idx]
  #which snps were included
  eqtl_cis$SNP_idx=idx
  
  idx=match(eqtl_cis$gene,genepos$geneid)
  eqtl_cis$gene_chr=genepos$chr[idx]
  eqtl_cis$gene_pos=genepos$s1[idx]
  #get the total position
  eqtl_cis$SNP_posall=totalpos(chr=eqtl_cis$SNP_chr,pos=eqtl_cis$SNP_pos)
  eqtl_cis$gene_posall=totalpos(chr=eqtl_cis$gene_chr,pos=eqtl_cis$gene_pos)
  levelmat=data.frame(value=-log10(eqtl_cis[,colpvalue]),pos1=eqtl_cis$SNP_posall,pos2=eqtl_cis$gene_posall,
                      chr=eqtl_cis$SNP_chr,gene_chr=eqtl_cis$gene_chr,snp_idx=eqtl_cis$SNP_idx,
                      opos1=eqtl_cis$SNP_pos,opos2=eqtl_cis$gene_pos,gene=eqtl_cis$gene,fdr=eqtl_cis$FDR)
  levelmat$chr=gsub(23,"X",levelmat$chr)
  levelmat$chr=gsub(24,"Y",levelmat$chr)
  # cc = colorRampPalette( c("green","white","red"))
  # numcol=21
  # colors=cc(numcol)
  # zrng <- range(levelmat$value)       # what's the range of z
  # tol <- 1e-2            # what tolerance is necessary?
  # zquant=quantile(levelmat$value,c(0,0.5,0.9,1))
  # colorBreaks1 <- seq(zquant[1]-0.01,zquant[2],length.out = ceiling(numcol*0.5))
  # colorBreaks2 <- seq(zquant[2]+0.1,zquant[4]+0.01,length.out = numcol-ceiling(numcol*0.5))
  # colorBreaks <- c(colorBreaks1,colorBreaks2)
  # df_breaks=seq(zrng[1],zrng[2],length.out = numcol)
  # # levelplot(cor~ge*me,data=levelmat,xlab="Gene expression",ylab="Methylation")
  # levelmat$color=rep(NA,nrow(levelmat))
  # for (i in 1:(numcol-1))
  # {
  #   idx=which(levelmat$value>=colorBreaks[i] & levelmat$value<colorBreaks[i+1])
  #   levelmat$color[idx]=i
  # }
  # # idx=which(levelmat$color>numcol/2)
  # # levelmat=levelmat[idx,]
  # levelmat$colors=colors[levelmat$color]
  levelmat=levelmat[order(levelmat$value),]
  return(levelmat)
}
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/ge_annohg19.RData")
addgenename_GE=function(datatable=tmp5,anno=ge_annohg19,colname="genename")
{
  idx=match(datatable$gene,anno$Probe_Id)
  datatable=cbind(datatable,anno$Symbol[idx])
  colnames(datatable)[ncol(datatable)]=colname
  return(datatable)
}
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
addgenename_ME=function(datatable=tmp7,anno1=anno,colname="genename")
{
  idx=match(datatable$gene,anno1$IlmnID)
  datatable=cbind(datatable,anno1$UCSC_RefGene_Name[idx])
  colnames(datatable)[ncol(datatable)]=colname
  return(datatable)
}
library(GenomicRanges)
checkoverlapsnps=function(qtlres1=HUTCH_cis_eqtl,qtlres2=NORMAL_cis_eqtl,
                          snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                          snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt")
{
  snppos1=fread(snpposfile1,sep="\t")
  snppos2=fread(snpposfile2,sep="\t")
  gr_res1=GRanges(seqnames = qtlres1$chr,ranges = IRanges(start=snppos1$pos[qtlres1$snp_idx],width=1))
  gr_res1=unique(gr_res1)
  gr_res2=GRanges(seqnames = qtlres2$chr,ranges = IRanges(start=snppos2$pos[qtlres2$snp_idx],width=1))
  gr_res2=unique(gr_res2)
  print(paste0("number of unique SNP1:",length(gr_res1)))
  print(paste0("number of unique SNP2:",length(gr_res2)))
  print(paste0("number of overlapped SNP:",length(intersect(gr_res1,gr_res2))))
}

pvalue_overlapsnp=function(totalnum=147,num1=33,num2=127,num12=32)
{
  m=matrix(c(num12,num1-num12,num2-num12,totalnum-num1-num2+num12),nrow=2,dimnames = list(c("in1","notin1"),c("in2","notin2")))
  print(m)
  print(fisher.test(m)$p.value)
}

#CIS-------------------------------------------------------------------------------
#for highrisk HUTCH eqtl
#for number of tests
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk.RData")
highrisk_HUTCH_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_HUTCH_cis_eqtl=addgenename_GE(highrisk_HUTCH_cis_eqtl)
length(unique(highrisk_HUTCH_cis_eqtl$genename))
highrisk_all_HUTCH_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
highrisk_all_HUTCH_cis_eqtl=addgenename_GE(highrisk_all_HUTCH_cis_eqtl)
length(unique(highrisk_all_HUTCH_cis_eqtl$genename))

#for highrisk NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk.RData")
highrisk_NORMAL_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_NORMAL_cis_eqtl=addgenename_GE(highrisk_NORMAL_cis_eqtl)
length(unique(highrisk_NORMAL_cis_eqtl$genename))
highrisk_all_NORMAL_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_cis",
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
HUTCH_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",bonfadj=qtl$cis$ntests)
HUTCH_cis_eqtl=addgenename_GE(HUTCH_cis_eqtl)
length(unique(HUTCH_cis_eqtl$genename))
#for NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl.RData")
NORMAL_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_cis",
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
highrisk_HUTCH_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_HUTCH_cis_mqtl=addgenename_ME(highrisk_HUTCH_cis_mqtl)
length(unique(highrisk_HUTCH_cis_mqtl$genename))
highrisk_all_HUTCH_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
highrisk_all_HUTCH_cis_mqtl=addgenename_ME(highrisk_all_HUTCH_cis_mqtl)
length(unique(highrisk_all_HUTCH_cis_mqtl$genename))

#for highrisk NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk.RData")
highrisk_NORMAL_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_cis",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_NORMAL_cis_mqtl=addgenename_ME(highrisk_NORMAL_cis_mqtl)
length(unique(highrisk_NORMAL_cis_mqtl$genename))
sum(length(unique(intersect(highrisk_HUTCH_cis_mqtl$snp_idx,highrisk_NORMAL_cis_mqtl$snp_idx))))
highrisk_all_NORMAL_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_cis",
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
HUTCH_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI/mqtl_cis",
                            snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                            geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=sum(num_cis))
HUTCH_cis_mqtl=addgenename_ME(HUTCH_cis_mqtl)
length(unique(HUTCH_cis_mqtl$genename))

# rm(qtl)
# load("/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl.RData")
# HUTCH_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_cis",
#                              snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
#                              geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=qtl$cis$ntests)

#for NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl.RData")
NORMAL_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_cis",
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
highrisk_HUTCH_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_HUTCH_trans_eqtl=addgenename_GE(highrisk_HUTCH_trans_eqtl)
length(unique(highrisk_HUTCH_trans_eqtl$genename))
highrisk_all_HUTCH_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
highrisk_all_HUTCH_trans_eqtl=addgenename_GE(highrisk_all_HUTCH_trans_eqtl)
length(unique(highrisk_all_HUTCH_trans_eqtl$genename))

#for highrisk NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk.RData")
highrisk_NORMAL_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_trans",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_NORMAL_trans_eqtl=addgenename_GE(highrisk_NORMAL_trans_eqtl)
length(unique(highrisk_NORMAL_trans_eqtl$genename))
checkoverlapsnps(qtlres1 = highrisk_HUTCH_trans_eqtl,qtlres2 = highrisk_NORMAL_trans_eqtl,
                 snpposfile1 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
highrisk_all_NORMAL_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_trans",
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
HUTCH_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_trans",
                            snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                            geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",bonfadj=qtl$trans$ntests)
HUTCH_trans_eqtl=addgenename_GE(HUTCH_trans_eqtl)
length(unique(HUTCH_trans_eqtl$genename))
#for NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl.RData")
NORMAL_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_trans",
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
highrisk_HUTCH_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_HUTCH_trans_mqtl=addgenename_ME(highrisk_HUTCH_trans_mqtl)
length(unique(highrisk_HUTCH_trans_mqtl$genename))
highrisk_all_HUTCH_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
highrisk_all_HUTCH_trans_mqtl=addgenename_ME(highrisk_all_HUTCH_trans_mqtl)
length(unique(highrisk_all_HUTCH_trans_mqtl$genename))

#for highrisk NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk.RData")
highrisk_NORMAL_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_trans",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_NORMAL_trans_mqtl=addgenename_ME(highrisk_NORMAL_trans_mqtl)
length(unique(highrisk_NORMAL_trans_mqtl$genename))
sum(length(unique(intersect(highrisk_HUTCH_trans_mqtl$snp_idx,highrisk_NORMAL_trans_mqtl$snp_idx))))
highrisk_all_NORMAL_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_trans",
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
HUTCH_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI/mqtl_trans",
                            snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                            geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",bonfadj=sum(num_trans))
HUTCH_trans_mqtl=addgenename_ME(HUTCH_trans_mqtl)
length(unique(HUTCH_trans_mqtl$genename))
#for NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl.RData")
NORMAL_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_trans",
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
overlappairs=function(dat1=highrisk_cis_hutch_eqtl,dat2=highrisk_cis_hutch_mqtl)
{
  tmp=merge(dat1,dat2,by="snp_idx")
  print(paste0("number of overlapped SNP:",length(unique(tmp$snp_idx))))
  print(paste0("number of phenoytype1:",length(unique(tmp$gene.x))))
  print(paste0("number of phenoytype2:",length(unique(tmp$gene.y))))
  print(paste0("number of triplets:",nrow(tmp)))
  idx1=which(as.character(tmp$gene.x)==as.character(tmp$gene.y))
  print(paste0("number of overapped pairs:",sum(as.character(tmp$gene.x)==as.character(tmp$gene.y),na.rm=0)))
  return(list(snpidx=unique(tmp$snp_idx),pair=tmp[idx1,c("snp_idx","gene.x")],triplet=tmp))
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

validatepairs=function(dat1=highrisk_cis_hutch_eqtl,dat2=highrisk_all_NORMAL_cis_eqtl,phenotype2=tcga_ge_pos)
{
  tmp=merge(dat1,dat2,by=c("snp_idx","gene"))
  idx=which(10^-tmp$value.y<=0.05)
  print(paste0("number of qtl pairs to validate:",nrow(dat1)))
  print(paste0("number of qtl pairs can be validated:",sum(dat1$gene %in% phenotype2[,1],na.rm=T)))
  print(paste0("number of pairs validated:",length(idx)))
  
  print(paste0("number of phenotype to validate:",length(unique(dat1$gene))))
  print(paste0("number of phenotype can be validated:",length(unique(dat1$gene[dat1$gene %in% phenotype2[,1]]))))
  print(paste0("number of phenotype validated:",length(unique(tmp$gene))))
  
  print(paste0("number of snp to validate:",length(unique(dat1$snp_idx))))
  print(paste0("number of snp validated:",length(unique(tmp$snp_idx[idx]))))
  return(tmp)
}
validtriplet=function(dat1=olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl$triplet,validateeqtl=highrisk_all_NORMAL_cis_eqtl,
                      validatemqtl=highrisk_all_NORMAL_cis_mqtl,phenotype1=tcga_ge_pos,phenotype2=tcga_me_pos)
{
  validateeqtl$gene=as.character(validateeqtl$gene)
  validatemqtl$gene=as.character(validatemqtl$gene)
  dat1$gene.x=as.character(dat1$gene.x)
  dat1$gene.y=as.character(dat1$gene.y)
  dat1$validate=0
  dat1$validate_eqtl=NA
  dat1$validate_mqtl=NA
  print(paste0("number of triplets to validate:",nrow(dat1)))
  print(paste0("number of triplets can be validate:",sum(dat1$gene.x %in% phenotype1[,1] & dat1$gene.y %in% phenotype2[,1],na.rm=T)))
  for (i in 1:nrow(dat1))
  {
    idx1=which(validateeqtl$snp_idx==dat1$snp_idx[i] & validateeqtl$gene==dat1$gene.x[i])
    idx2=which(validatemqtl$snp_idx==dat1$snp_idx[i] & validatemqtl$gene==dat1$gene.y[i])
    if (length(idx1)*length(idx2)==1)
    {
      if (10^-validateeqtl$value[idx1]<=0.05 & 10^-validatemqtl$value[idx2]<=0.05)
      {
        dat1$validate[i]=1
        dat1$validate_eqtl[i]=validateeqtl$value[idx1]
        dat1$validate_mqtl[i]=validatemqtl$value[idx2]
      }
    }
  }
  print(paste0("number of triplets validated:", sum(dat1$validate==1,na.rm = T)))
  print(paste0("number of snp to validate:", length(unique(dat1$snp_idx))))
  print(paste0("number of snp validated:", length(unique(dat1$snp_idx[dat1$validate==1]))))
  return(dat1)
}

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
library(gdata)
polycpg=read.xls("/fh/fast/dai_j/CancerGenomics/Tools/database/other/48640-polymorphic-CpGs-Illumina450k.xlsx")
snpcpg=read.xls("/fh/fast/dai_j/CancerGenomics/Tools/database/other/48640-polymorphic-CpGs-Illumina450k.xlsx",sheet=2)
