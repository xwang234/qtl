#!/usr/bin/env Rscript
library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
#generate inputs for impute2
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/TCGAnormals.RData")
snp6_anno <- read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
#the SNP probes:
sum(rownames(genotypedata) %in% snp6_anno$Probe.Set.ID)==nrow(genotypedata)
idx=match(rownames(genotypedata),snp6_anno$Probe.Set.ID)
snp6_anno=snp6_anno[idx,]
snp6_anno$Physical.Position=as.integer(snp6_anno$Physical.Position)
sum(is.na(snp6_anno$Chromosome))
sum(is.na(snp6_anno$Physical.Position))
table(snp6_anno$Chromosome)

chrs=c(1:22,"X")
#generate strand file and genotype file:
for (chr in chrs)
{
  idx=which(snp6_anno$Chromosome==chr & !is.na(snp6_anno$Physical.Position))
  annochrtable=snp6_anno[idx,]
  gtchrtable=genotypedata[idx,]
  if(sum(duplicated(annochrtable$Physical.Position))>0)
  {
    warning(paste0("There duplicated positions...","in chr",chr))
    idx1=!duplicated(annochrtable$Physical.Position)
    annochrtable=annochrtable[idx1,]
    gtchrtable=gtchrtable[idx1,]
  }
  
  if (sum(annochrtable$Probe.Set.ID==rownames(gtchrtable))!=nrow(gtchrtable))
  {
    warning(paste0("The order of SNPs have problem...","in chr",chr))
  }
  #change genotype notation,0->1 0 0, 1->0 1 0, 2-> 0 0 1
  gtchrtable1=data.frame(matrix(NA,nrow=nrow(gtchrtable),ncol=3*ncol(gtchrtable)))
  for (i in 1:ncol(gtchrtable))
  {
    colid=(i-1)*3+1
    idx0=which(gtchrtable[,i]==0)
    if (length(idx0)>0) gtchrtable1[idx0,c(colid,colid+1,colid+2)]=data.frame(rep(1,length(idx0)),rep(0,length(idx0)),rep(0,length(idx0)))
    idx1=which(gtchrtable[,i]==1)
    if (length(idx1)>0) gtchrtable1[idx1,c(colid,colid+1,colid+2)]=data.frame(rep(0,length(idx1)),rep(1,length(idx1)),rep(0,length(idx1)))
    idx2=which(gtchrtable[,i]==2)
    if (length(idx2)>0) gtchrtable1[idx2,c(colid,colid+1,colid+2)]=data.frame(rep(0,length(idx2)),rep(0,length(idx2)),rep(1,length(idx2)))
  }
  
  chrorder=order(annochrtable$Physical.Position)
  strandtable=data.frame(pos=annochrtable$Physical.Position[chrorder],strand=annochrtable$Strand[chrorder], stringsAsFactors = F)
  rownames(strandtable)=annochrtable$Probe.Set.ID[chrorder]
  write.table(strandtable,file=paste0("../result/imputation/","SNP6_strand_chr",chr,".txt"),row.names = F,col.names = F,sep="\t",quote=F)
  gttable=data.frame(probeid=annochrtable$Probe.Set.ID[chrorder],dbsnp=annochrtable$dbSNP.RS.ID[chrorder],pos=annochrtable$Physical.Position[chrorder],
                     ref=annochrtable$Allele.A[chrorder],alt=annochrtable$Allele.B[chrorder], stringsAsFactors = F)
  gttable=cbind.data.frame(gttable,gtchrtable1[chrorder,])
  write.table(gttable,file=paste0("../result/imputation/","SNP6_gtdata_chr",chr,".txt"),row.names = F,col.names = F,sep="\t",quote=F)
}

#to generate sample file used for chrX
samplefile="../result/imputation/SNP6_samplefile.txt"
tmp=data.frame(matrix(NA,nrow=36,ncol=4),stringsAsFactors = F)
colnames(tmp)=c("ID_1","ID_2","missing","sex")
tmp[,1]=tmp[,2]=0:35
tmp[,3]=0
tmp[,4]=1
tmp[1,]=c(0,0,0,"D")
write.table(tmp,file=samplefile,row.names = F,sep=" ",quote=F)

samplefile="../result/imputation3/SNP6_allsamplefile.txt"
tmp=data.frame(matrix(NA,nrow=497,ncol=4),stringsAsFactors = F)
colnames(tmp)=c("ID_1","ID_2","missing","sex")
tmp[,1]=tmp[,2]=0:496
tmp[,3]=0
tmp[,4]=1
tmp[1,]=c(0,0,0,"D")
write.table(tmp,file=samplefile,row.names = F,sep=" ",quote=F)

pca <- prcomp(t(expr[-1]), center=TRUE, scale = TRUE)
pc <- pca$x
covar <- SlicedData$new()
covar$CreateFromMatrix(t(pc[,1:10]))

