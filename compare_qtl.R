#!/usr/bin/env Rscript

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
                      opos1=eqtl_cis$SNP_pos,opos2=eqtl_cis$gene_pos,gene=eqtl_cis$gene)
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
addgenename_GE=function(datatable=tmp5,anno=ge_annohg19,colname="GE_gene")
{
  idx=match(datatable$gene,anno$Probe_Id)
  datatable=cbind(datatable,anno$Symbol[idx])
  colnames(datatable)[ncol(datatable)]=colname
  return(datatable)
}
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
addgenename_ME=function(datatable=tmp7,anno1=anno,colname="ME_gene")
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
length(unique(highrisk_HUTCH_cis_eqtl$GE_gene))
highrisk_all_HUTCH_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
highrisk_all_HUTCH_cis_eqtl=addgenename_GE(highrisk_all_HUTCH_cis_eqtl)
length(unique(highrisk_all_HUTCH_cis_eqtl$GE_gene))

#for highrisk NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk.RData")
highrisk_NORMAL_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_cis",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_NORMAL_cis_eqtl=addgenename_GE(highrisk_NORMAL_cis_eqtl)
length(unique(highrisk_NORMAL_cis_eqtl$GE_gene))
highrisk_all_NORMAL_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_cis",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt")
highrisk_all_NORMAL_cis_eqtl=addgenename_GE(highrisk_all_NORMAL_cis_eqtl)
length(unique(highrisk_all_NORMAL_cis_eqtl$GE_gene))
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
length(unique(HUTCH_cis_eqtl$GE_gene))
#for NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl.RData")
NORMAL_cis_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_cis",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$cis$ntests)
NORMAL_cis_eqtl=addgenename_GE(NORMAL_cis_eqtl)
length(unique(NORMAL_cis_eqtl$GE_gene))
length(unique(intersect(HUTCH_cis_eqtl$GE_gene,NORMAL_cis_eqtl$GE_gene)))
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
length(unique(highrisk_HUTCH_cis_mqtl$ME_gene))
highrisk_all_HUTCH_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
highrisk_all_HUTCH_cis_mqtl=addgenename_ME(highrisk_all_HUTCH_cis_mqtl)
length(unique(highrisk_all_HUTCH_cis_mqtl$ME_gene))

#for highrisk NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk.RData")
highrisk_NORMAL_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_cis",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$cis$ntests)
highrisk_NORMAL_cis_mqtl=addgenename_ME(highrisk_NORMAL_cis_mqtl)
length(unique(highrisk_NORMAL_cis_mqtl$ME_gene))
sum(length(unique(intersect(highrisk_HUTCH_cis_mqtl$snp_idx,highrisk_NORMAL_cis_mqtl$snp_idx))))
highrisk_all_NORMAL_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_cis",
                                          snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                          geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt")
highrisk_all_NORMAL_cis_mqtl=addgenename_ME(highrisk_all_NORMAL_cis_mqtl)
length(unique(highrisk_all_NORMAL_cis_mqtl$ME_gene))
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
length(unique(HUTCH_cis_mqtl$ME_gene))
#for NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl.RData")
NORMAL_cis_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_cis",
                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt",
                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$cis$ntests)
NORMAL_cis_mqtl=addgenename_ME(NORMAL_cis_mqtl)
length(unique(NORMAL_cis_mqtl$ME_gene))
length(unique(intersect(HUTCH_cis_mqtl$ME_gene,NORMAL_cis_mqtl$ME_gene)))
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
length(unique(highrisk_HUTCH_trans_eqtl$GE_gene))
highrisk_all_HUTCH_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
highrisk_all_HUTCH_trans_eqtl=addgenename_GE(highrisk_all_HUTCH_trans_eqtl)
length(unique(highrisk_all_HUTCH_trans_eqtl$GE_gene))

#for highrisk NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk.RData")
highrisk_NORMAL_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_trans",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_NORMAL_trans_eqtl=addgenename_GE(highrisk_NORMAL_trans_eqtl)
length(unique(highrisk_NORMAL_trans_eqtl$GE_gene))
checkoverlapsnps(qtlres1 = highrisk_HUTCH_trans_eqtl,qtlres2 = highrisk_NORMAL_trans_eqtl,
                 snpposfile1 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2 = "/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
highrisk_all_NORMAL_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_trans",
                                          snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                          geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt")
highrisk_all_NORMAL_trans_eqtl=addgenename_GE(highrisk_all_NORMAL_trans_eqtl)
length(unique(highrisk_all_NORMAL_trans_eqtl$GE_gene))
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
length(unique(HUTCH_trans_eqtl$GE_gene))
#for NORMAL eqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl.RData")
NORMAL_trans_eqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_trans",
                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt",
                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",bonfadj=qtl$trans$ntests)
NORMAL_trans_eqtl=addgenename_GE(NORMAL_trans_eqtl)
length(unique(NORMAL_trans_eqtl$GE_gene))
length(unique(intersect(HUTCH_trans_eqtl$GE_gene,NORMAL_trans_eqtl$GE_gene)))
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
length(unique(highrisk_HUTCH_trans_mqtl$ME_gene))
highrisk_all_HUTCH_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
highrisk_all_HUTCH_trans_mqtl=addgenename_ME(highrisk_all_HUTCH_trans_mqtl)
length(unique(highrisk_all_HUTCH_trans_mqtl$ME_gene))

#for highrisk NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk.RData")
highrisk_NORMAL_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_trans",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$trans$ntests)
highrisk_NORMAL_trans_mqtl=addgenename_ME(highrisk_NORMAL_trans_mqtl)
length(unique(highrisk_NORMAL_trans_mqtl$ME_gene))
sum(length(unique(intersect(highrisk_HUTCH_trans_mqtl$snp_idx,highrisk_NORMAL_trans_mqtl$snp_idx))))
highrisk_all_NORMAL_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_trans",
                                          snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                          geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt")
highrisk_all_NORMAL_trans_mqtl=addgenename_ME(highrisk_all_NORMAL_trans_mqtl)
length(unique(highrisk_all_NORMAL_trans_mqtl$ME_gene))

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
length(unique(HUTCH_trans_mqtl$ME_gene))
#for NORMAL mqtl
rm(qtl)
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl.RData")
NORMAL_trans_mqtl=formlevelmat(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_trans",
                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt",
                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",bonfadj=qtl$trans$ntests)
NORMAL_trans_mqtl=addgenename_ME(NORMAL_trans_mqtl)
length(unique(NORMAL_trans_mqtl$ME_gene))
length(unique(intersect(HUTCH_trans_mqtl$ME_gene,NORMAL_trans_mqtl$ME_gene)))
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