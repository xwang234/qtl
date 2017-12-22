#!/usr/bin/env Rscript
#create Rdata for genotype/geneexp/methylation

library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library(sas7bdat)

#Hutch--------------------------
#Methylation data
ME=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",header=T,sep="\t")
ME=as.data.frame(ME)
rownames(ME)=ME$id
ME=ME[,-1]
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
ME_anno=anno
ME_pos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",header=T,sep="\t")
ME_pos=as.data.frame(ME_pos)
rownames(ME_pos)=ME_pos$geneid
ME_pos=ME_pos[,-1]
sum(rownames(ME)==rownames(ME_pos))
idx=match(rownames(ME),ME_anno$IlmnID)
sum(ME_pos$s1==ME_anno$MAPINFO[idx])
#covariate
hutchclinical=read.sas7bdat("/fh/fast/stanford_j/Janet/alldata_2016dec20.sas7bdat")
#PCA
pcME=prcomp(t(ME),scale=T,center=T)
varprop=pcME$sdev^2/sum(pcME$sdev^2)
numpc=sum(varprop>0.01) #9
idx=match(colnames(ME),hutchclinical$studyno)
COVA_ME=data.frame(matrix(NA,nrow=6+numpc,ncol=ncol(ME)+1))
colnames(COVA_ME)=c("variable",colnames(ME))
COVA_ME[,1]=c("Gleason_score","smoke_status","AGEREF","bmi","alcoholr","study",paste0("pc",1:numpc))
rownames(COVA_ME)=c("Gleason_score","smoke_status","AGEREF","bmi","alcoholr","study",paste0("pc",1:numpc))
COVA_ME[1,2:ncol(COVA_ME)]=hutchclinical$Gleason_score[idx]
COVA_ME[2,2:ncol(COVA_ME)]=hutchclinical$smoke_status[idx]
COVA_ME[3,2:ncol(COVA_ME)]=hutchclinical$AGEREF[idx]
COVA_ME[4,2:ncol(COVA_ME)]=hutchclinical$bmi[idx]
COVA_ME[5,2:ncol(COVA_ME)]=hutchclinical$alcoholr[idx]
COVA_ME[6,2:ncol(COVA_ME)]=hutchclinical$study[idx]

idx=match(colnames(ME),rownames(pcME$x))
for (i in 7:(6+numpc))
{
  COVA_ME[i,2:ncol(COVA_ME)]=pcME$x[idx,i-6]
}
COVA_ME=COVA_ME[,-1]


#Gene expression
GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",header=T,sep="\t")
GE=as.data.frame(GE)
rownames(GE)=GE$id
GE=GE[,-1]
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/ge_annohg19.RData")
GE_anno=ge_annohg19
GE_anno=GE_anno[GE_anno$Probe_Id %in% rownames(GE),]
GE_pos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",header=T,sep="\t")
GE_pos=as.data.frame(GE_pos)
rownames(GE_pos)=GE_pos$geneid
GE_pos=GE_pos[,-1]
#covariate
hutchclinical=read.sas7bdat("/fh/fast/stanford_j/Janet/alldata_2016dec20.sas7bdat")
#PCA
pcGE=prcomp(t(GE),scale=T,center=T)
varprop=pcGE$sdev^2/sum(pcGE$sdev^2)
numpc=sum(varprop>0.01) #14
idx=match(colnames(GE),hutchclinical$studyno)
COVA_GE=data.frame(matrix(NA,nrow=6+numpc,ncol=ncol(GE)+1))
colnames(COVA_GE)=c("variable",colnames(GE))
COVA_GE[,1]=c("Gleason_score","smoke_status","AGEREF","bmi","alcoholr","study",paste0("pc",1:numpc))
rownames(COVA_GE)=c("Gleason_score","smoke_status","AGEREF","bmi","alcoholr","study",paste0("pc",1:numpc))
COVA_GE[1,2:ncol(COVA_GE)]=hutchclinical$Gleason_score[idx]
COVA_GE[2,2:ncol(COVA_GE)]=hutchclinical$smoke_status[idx]
COVA_GE[3,2:ncol(COVA_GE)]=hutchclinical$AGEREF[idx]
COVA_GE[4,2:ncol(COVA_GE)]=hutchclinical$bmi[idx]
COVA_GE[5,2:ncol(COVA_GE)]=hutchclinical$alcoholr[idx]
COVA_GE[6,2:ncol(COVA_GE)]=hutchclinical$study[idx]

idx=match(colnames(GE),rownames(pcGE$x))
for (i in 7:(6+numpc))
{
  COVA_GE[i,2:ncol(COVA_GE)]=pcGE$x[idx,i-6]
}
COVA_GE=COVA_GE[,-1]


#Genotype
SNP_ME=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_ME.txt",header=T,sep="\t")
SNP_ME=as.data.frame(SNP_ME)
rownames(SNP_ME)=SNP_ME$id
SNP_ME=SNP_ME[,-1]
sum(colnames(SNP_ME)==colnames(ME))
SNP_GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",header=T,sep="\t")
SNP_GE=as.data.frame(SNP_GE)
rownames(SNP_GE)=SNP_GE$id
SNP_GE=SNP_GE[,-1]
sum(colnames(SNP_GE)==colnames(GE))
#quality score
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData")
snps=read.xls("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/RESUB_Supplementary_Table16_v9.xlsx",skip=1,sep=",")
snps=snps[1:147,c(1:3,5)]
colnames(snps)=c("snp","chr","position","allele")
snps$chr=gsub(23,"X",snps$chr)
checksum=rep(0,147)
for (i in c(1:129,131:147))
{
  idx=which(highrisk$chr==allsnps$chr[i] & highrisk$V2==allsnps$position[i])
  idx1=which(snps$chr==allsnps$chr[i] & snps$position==allsnps$position[i])
  #check allele
  if (snps$allele[idx1] == highrisk$V3[idx] | snps$allele[idx1] == highrisk$V4[idx])
  {
    checksum[i]=1
  }
}
#allels were correct

#map of SNPid
idx=rep(0,147)
for(i in 1:147)
{
  idx[i]=which(highrisk$chr==allsnps$chr[i] & highrisk$V2==allsnps$position[i])
}
allsnps$imputedsnp=highrisk$V1[idx]

SNP_ALL=highrisk
SNP_impscore=allsnps

#form the QTL pairids
library(GenomicRanges)
formmeQTLpairs=function(pos=ME_pos,cutoff=1e6)
{
  res=NULL
  gr_snp=GRanges(seqnames = SNP_ALL$chr,ranges=IRanges(start=SNP_ALL$V2,width = 1))
  gr_pos=GRanges(seqnames = pos$chr,ranges=IRanges(start=pos$s1,end=pos$s2))
  for (i in 1:length(gr_snp))
  {
    
    tmp=distance(gr_pos,gr_snp[i])
    idx=which(tmp<cutoff)
    res=rbind.data.frame(res,data.frame(snpidx=i,yidx=idx))
  }
  for (i in 1:length(gr_snp))
  {
    if (sum(res$snpidx==i)==0) print(i)
  }
  return(res)
}

SNP_ME_pairs=formmeQTLpairs()
colnames(SNP_ME_pairs)[2]="methyidx"
SNP_GE_pairs=formmeQTLpairs(pos=GE_pos)
colnames(SNP_GE_pairs)[2]="genexpidx"

#SNP_GE_PCA,SNP_ME_PCA
SNP_GE_PCA=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE_PCA.txt",header=T,sep="\t")
colnames(SNP_GE_PCA)=gsub("^X","",colnames(SNP_GE_PCA))
rownames(SNP_GE_PCA)=SNP_GE_PCA$id
SNP_GE_PCA=SNP_GE_PCA[,-1]

SNP_ME_PCA=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_ME_PCA.txt",header=T,sep="\t")
colnames(SNP_ME_PCA)=gsub("^X","",colnames(SNP_ME_PCA))
rownames(SNP_ME_PCA)=SNP_ME_PCA$id
SNP_ME_PCA=SNP_ME_PCA[,-1]

#peer
GE_PEER=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_PEER.txt",header=T,sep="\t")
rownames(GE_PEER)=GE_PEER$id
GE_PEER=GE_PEER[,-1]
colnames(GE_PEER)=gsub("^X","",colnames(GE_PEER))
save(GE,ME,COVA_GE,GE_PEER,COVA_ME,GE_anno,ME_anno,SNP_ALL,SNP_GE,SNP_ME,SNP_GE_PCA,SNP_ME_PCA,SNP_impscore,SNP_GE_pairs,SNP_ME_pairs,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/hutch_SNP_GE_ME.RData")


#TCGA tumors-----------------------------------------------
#Methylation data
ME=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME.txt",header=T,sep="\t")
ME=as.data.frame(ME)
rownames(ME)=ME$id
ME=ME[,-1]
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
ME_anno=anno
ME_pos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_POS.txt",header=T,sep="\t")
ME_pos=as.data.frame(ME_pos)
rownames(ME_pos)=ME_pos$geneid
ME_pos=ME_pos[,-1]
sum(rownames(ME)==rownames(ME_pos))
idx=match(rownames(ME),ME_anno$IlmnID)
sum(ME_pos$s1==ME_anno$MAPINFO[idx])
#covariate
tcgaclinical=read.table("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/clinical/dc99e120-0159-4470-ab3c-6581032935e9/nationwidechildrens.org_clinical_patient_prad.txt",
                        skip=1,stringsAsFactors = F,header=T,sep="\t")
tcgaclinical=tcgaclinical[-1,]
for (i in 1:ncol(tcgaclinical))
{
  tcgaclinical[,i]=gsub("\\[Not Applicable\\]",NA,tcgaclinical[,i])
  tcgaclinical[,i]=gsub("\\[Not Available\\]",NA,tcgaclinical[,i])
  tcgaclinical[,i]=gsub("\\[Unknown\\]",NA,tcgaclinical[,i])
}
tcgaclinical$sampleid=tcgaclinical$bcr_patient_barcode
tcgaclinical$gleason=tcgaclinical$gleason_score
tcgaclinical$age=-as.numeric(tcgaclinical$days_to_birth)/365

#PCA
pcME=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_PCA.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(pcME)=gsub(".","-",colnames(pcME),fixed = T)
colnames(pcME)[1]="variable"
COVA_ME=data.frame(matrix(NA,nrow=2,ncol=ncol(pcME)))
colnames(COVA_ME)=colnames(pcME)
COVA_ME[,1]=c("Gleason_score","AGEREF")
idx=match(colnames(pcME)[2:ncol(pcME)],tcgaclinical$sampleid)
COVA_ME[1,2:ncol(COVA_ME)]=tcgaclinical$gleason[idx]
COVA_ME[2,2:ncol(COVA_ME)]=round(tcgaclinical$age[idx],digits = 2)
COVA_ME=rbind.data.frame(COVA_ME,pcME)
rownames(COVA_ME)=COVA_ME$variable
COVA_ME=COVA_ME[,-1]


#Gene expression
GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE.txt",header=T,sep="\t")
GE=as.data.frame(GE)
rownames(GE)=GE$id
GE=GE[,-1]
GE_pos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_POS.txt",header=T,sep="\t")
GE_pos=as.data.frame(GE_pos)
rownames(GE_pos)=GE_pos$geneid
GE_pos=GE_pos[,-1]
#covariate
#PCA
pcGE=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_PCA.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(pcGE)=gsub(".","-",colnames(pcGE),fixed = T)
colnames(pcGE)[1]="variable"
COVA_GE=data.frame(matrix(NA,nrow=2,ncol=ncol(pcGE)))
colnames(COVA_GE)=colnames(pcGE)
COVA_GE[,1]=c("Gleason_score","AGEREF")
idx=match(colnames(pcGE)[2:ncol(pcGE)],tcgaclinical$sampleid)
COVA_GE[1,2:ncol(COVA_GE)]=tcgaclinical$gleason[idx]
COVA_GE[2,2:ncol(COVA_GE)]=round(tcgaclinical$age[idx],digits = 2)
COVA_GE=rbind.data.frame(COVA_GE,pcGE)
rownames(COVA_GE)=COVA_GE$variable
COVA_GE=COVA_GE[,-1]

#Genotype
SNP_ME=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_ME.txt",header=T,sep="\t")
SNP_ME=as.data.frame(SNP_ME)
rownames(SNP_ME)=SNP_ME$id
SNP_ME=SNP_ME[,-1]
sum(colnames(SNP_ME)==colnames(ME))
SNP_GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_GE.txt",header=T,sep="\t")
SNP_GE=as.data.frame(SNP_GE)
rownames(SNP_GE)=SNP_GE$id
SNP_GE=SNP_GE[,-1]
sum(colnames(SNP_GE)==colnames(GE))
#quality score
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors_allhighrisksnps_new.RData")
#snps=read.xls("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/RESUB_Supplementary_Table16_v9.xlsx",skip=1,sep=",")
snps=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/data/RESUB_Supplementary_Table16_v9.txt",sep="\t",header=T,stringsAsFactors = F)

snps=snps[1:147,c(1:3,5)]
colnames(snps)=c("snp","chr","position","allele")
snps$chr=gsub(23,"X",snps$chr)
checksum=rep(0,147)
for (i in c(1:129,131:147))
{
  idx=which(highrisk$chr==allsnps$chr[i] & highrisk$V2==allsnps$position[i])
  idx1=which(snps$chr==allsnps$chr[i] & snps$position==allsnps$position[i])
  #check allele
  if (snps$allele[idx1] == highrisk$V3[idx] | snps$allele[idx1] == highrisk$V4[idx])
  {
    checksum[i]=1
  }
}
#allels were correct

#map of SNPid
idx=rep(0,147)
for(i in 1:147)
{
  idx[i]=which(highrisk$chr==allsnps$chr[i] & highrisk$V2==allsnps$position[i])
}
allsnps$imputedsnp=highrisk$V1[idx]

SNP_ALL=highrisk
SNP_impscore=allsnps

#form the QTL pairids
library(GenomicRanges)
formmeQTLpairs=function(pos=ME_pos,cutoff=1e6)
{
  res=NULL
  gr_snp=GRanges(seqnames = SNP_ALL$chr,ranges=IRanges(start=SNP_ALL$V2,width = 1))
  gr_pos=GRanges(seqnames = pos$chr,ranges=IRanges(start=pos$s1,end=pos$s2))
  for (i in 1:length(gr_snp))
  {
    
    tmp=distance(gr_pos,gr_snp[i])
    idx=which(tmp<cutoff)
    res=rbind.data.frame(res,data.frame(snpidx=i,yidx=idx))
  }
  for (i in 1:length(gr_snp))
  {
    if (sum(res$snpidx==i)==0) print(i)
  }
  return(res)
}

SNP_ME_pairs=formmeQTLpairs()
colnames(SNP_ME_pairs)[2]="methyidx"
SNP_GE_pairs=formmeQTLpairs(pos=GE_pos)
colnames(SNP_GE_pairs)[2]="genexpidx"

#SNP_GE_PCA,SNP_ME_PCA
SNP_GE_PCA=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_SNP_GE_PCA.txt",header=T,sep="\t")
colnames(SNP_GE_PCA)=gsub("^X","",colnames(SNP_GE_PCA))
colnames(SNP_GE_PCA)=gsub(".","-",colnames(SNP_GE_PCA),fixed=T)
rownames(SNP_GE_PCA)=SNP_GE_PCA$id
SNP_GE_PCA=SNP_GE_PCA[,-1]
SNP_ME_PCA=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_SNP_ME_PCA.txt",header=T,sep="\t")
colnames(SNP_ME_PCA)=gsub("^X","",colnames(SNP_ME_PCA))
colnames(SNP_ME_PCA)=gsub(".","-",colnames(SNP_ME_PCA),fixed=T)
rownames(SNP_ME_PCA)=SNP_ME_PCA$id
SNP_ME_PCA=SNP_ME_PCA[,-1]

#peer
GE_PEER=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_PEER.txt",header=T,sep="\t")
rownames(GE_PEER)=GE_PEER$id
GE_PEER=GE_PEER[,-1]
colnames(GE_PEER)=gsub(".","-",colnames(GE_PEER),fixed = T)

save(GE,ME,COVA_GE,GE_PEER,COVA_ME,ME_anno,SNP_ALL,SNP_GE,SNP_GE_PCA,SNP_ME,SNP_ME_PCA,SNP_impscore,SNP_GE_pairs,SNP_ME_pairs,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/tcga_tumors_SNP_GE_ME.RData")

#NORMALS-----------------------------------------------
#Methylation data
ME=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME.txt",header=T,sep="\t")
ME=as.data.frame(ME)
rownames(ME)=ME$id
ME=ME[,-1]
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
ME_anno=anno
ME_pos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",header=T,sep="\t")
ME_pos=as.data.frame(ME_pos)
rownames(ME_pos)=ME_pos$geneid
ME_pos=ME_pos[,-1]
sum(rownames(ME)==rownames(ME_pos))
idx=match(rownames(ME),ME_anno$IlmnID)
sum(ME_pos$s1==ME_anno$MAPINFO[idx])
#covariate
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/clinical.RData")

#PCA
pcME=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_PCA.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(pcME)=gsub(".","-",colnames(pcME),fixed = T)
colnames(pcME)=gsub("^X","",colnames(pcME))
colnames(pcME)[1]="variable"
COVA_ME=data.frame(matrix(NA,nrow=2,ncol=ncol(pcME)))
colnames(COVA_ME)=colnames(pcME)
COVA_ME[,1]=c("Gleason_score","AGEREF")
idx=match(colnames(pcME)[2:ncol(pcME)],clinical$sampleid)
COVA_ME[1,2:ncol(COVA_ME)]=clinical$gleason[idx]
COVA_ME[2,2:ncol(COVA_ME)]=round(clinical$age[idx],digits = 2)
COVA_ME=rbind.data.frame(COVA_ME,pcME)
rownames(COVA_ME)=COVA_ME$variable
COVA_ME=COVA_ME[,-1]


#Gene expression
GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt",header=T,sep="\t")
GE=as.data.frame(GE)
rownames(GE)=GE$id
GE=GE[,-1]
GE_pos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",header=T,sep="\t")
GE_pos=as.data.frame(GE_pos)
rownames(GE_pos)=GE_pos$geneid
GE_pos=GE_pos[,-1]
#covariate
#PCA
pcGE=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_PCA.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(pcGE)=gsub(".","-",colnames(pcGE),fixed = T)
colnames(pcGE)=gsub("^X","",colnames(pcGE))
colnames(pcGE)[1]="variable"
COVA_GE=data.frame(matrix(NA,nrow=2,ncol=ncol(pcGE)))
colnames(COVA_GE)=colnames(pcGE)
COVA_GE[,1]=c("Gleason_score","AGEREF")
idx=match(colnames(pcGE)[2:ncol(pcGE)],clinical$sampleid)
COVA_GE[1,2:ncol(COVA_GE)]=clinical$gleason[idx]
COVA_GE[2,2:ncol(COVA_GE)]=round(clinical$age[idx],digits = 2)
COVA_GE=rbind.data.frame(COVA_GE,pcGE)
rownames(COVA_GE)=COVA_GE$variable
COVA_GE=COVA_GE[,-1]

#Genotype
SNP_ME=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_ME.txt",header=T,sep="\t")
SNP_ME=as.data.frame(SNP_ME)
rownames(SNP_ME)=SNP_ME$id
SNP_ME=SNP_ME[,-1]
sum(colnames(SNP_ME)==colnames(ME))
SNP_GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_GE.txt",header=T,sep="\t")
SNP_GE=as.data.frame(SNP_GE)
rownames(SNP_GE)=SNP_GE$id
SNP_GE=SNP_GE[,-1]
sum(colnames(SNP_GE)==colnames(GE))
#quality score
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData")
hutch_normal_highrisk=highrisk
hutch_normal_allsnps=allsnps
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_allhighrisksnps_new.RData")
tcga_normal_highrisk=highrisk
tcga_normal_allsnps=allsnps
# #read normal info
# load("/fh/fast/stanford_j/Xiaoyu/QTL/data/allnormaldata.RData")
# tcga_geneexp_samples=colnames(tcga_geneexp)
# tcga_methylation_samples=colnames(tcga_methy)
# tcga_genotype_samples=unique(c(colnames(tcga_geneexp),colnames(tcga_methy))) #after removing seminal vescile samples,it is different from colnames(genotypedata)
# hutch_normal_samples=colnames(hutch_geneexp)
# normal_samples=c(hutch_normal_samples,tcga_genotype_samples)
# normal_geneexp_samples=c(hutch_normal_samples,tcga_geneexp_samples)
# normal_methylation_samples=c(hutch_normal_samples,tcga_methylation_samples)

snps=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/data/RESUB_Supplementary_Table16_v9.txt",sep="\t",header=T,stringsAsFactors = F)

# snps=snps[1:147,c(1:3,5)]
# colnames(snps)=c("snp","chr","position","allele")
# snps$chr=gsub(23,"X",snps$chr)
# checksum=rep(0,147)
# for (i in c(1:129,131:147))
# {
#   idx=which(highrisk$chr==allsnps$chr[i] & highrisk$V2==allsnps$position[i])
#   idx1=which(snps$chr==allsnps$chr[i] & snps$position==allsnps$position[i])
#   #check allele
#   if (snps$allele[idx1] == highrisk$V3[idx] | snps$allele[idx1] == highrisk$V4[idx])
#   {
#     checksum[i]=1
#   }
# }
# #allels were correct

#map of SNPid
idx=rep(0,147)
for(i in 1:147)
{
  idx[i]=which(hutch_normal_highrisk$chr==hutch_normal_allsnps$chr[i] & hutch_normal_highrisk$V2==hutch_normal_allsnps$position[i])
}
hutch_normal_allsnps$imputedsnp=hutch_normal_highrisk$V1[idx]
hutch_normal_SNP_impscore=hutch_normal_allsnps

idx=rep(0,147)
for(i in 1:147)
{
  idx[i]=which(tcga_normal_highrisk$chr==tcga_normal_allsnps$chr[i] & tcga_normal_highrisk$V2==tcga_normal_allsnps$position[i])
}
tcga_normal_allsnps$imputedsnp=tcga_normal_highrisk$V1[idx]
tcga_normal_SNP_impscore=tcga_normal_allsnps
SNP_ALL=cbind.data.frame(tcga_normal_highrisk,hutch_normal_highrisk[,6:ncol(hutch_normal_highrisk)])
#form the QTL pairids
library(GenomicRanges)
formmeQTLpairs=function(pos=ME_pos,cutoff=1e6)
{
  res=NULL
  gr_snp=GRanges(seqnames = SNP_ALL$chr,ranges=IRanges(start=SNP_ALL$V2,width = 1))
  gr_pos=GRanges(seqnames = pos$chr,ranges=IRanges(start=pos$s1,end=pos$s2))
  for (i in 1:length(gr_snp))
  {
    
    tmp=distance(gr_pos,gr_snp[i])
    idx=which(tmp<cutoff)
    res=rbind.data.frame(res,data.frame(snpidx=i,yidx=idx))
  }
  for (i in 1:length(gr_snp))
  {
    if (sum(res$snpidx==i)==0) print(i)
  }
  return(res)
}

SNP_ME_pairs=formmeQTLpairs()
colnames(SNP_ME_pairs)[2]="methyidx"
SNP_GE_pairs=formmeQTLpairs(pos=GE_pos)
colnames(SNP_GE_pairs)[2]="genexpidx"

#SNP_GE_PCA,SNP_ME_PCA
SNP_GE_PCA=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_GE_PCA.txt",header=T,sep="\t")
colnames(SNP_GE_PCA)=gsub("^X","",colnames(SNP_GE_PCA))
colnames(SNP_GE_PCA)=gsub(".","-",colnames(SNP_GE_PCA),fixed=T)
rownames(SNP_GE_PCA)=SNP_GE_PCA$id
SNP_GE_PCA=SNP_GE_PCA[,-1]

SNP_ME_PCA=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_ME_PCA.txt",header=T,sep="\t")
colnames(SNP_ME_PCA)=gsub("^X","",colnames(SNP_ME_PCA))
colnames(SNP_ME_PCA)=gsub(".","-",colnames(SNP_ME_PCA),fixed=T)
rownames(SNP_ME_PCA)=SNP_ME_PCA$id
SNP_ME_PCA=SNP_ME_PCA[,-1]

#peer
GE_PEER=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_PEER.txt",header=T,sep="\t")
rownames(GE_PEER)=GE_PEER$id
GE_PEER=GE_PEER[,-1]
colnames(GE_PEER)=gsub("^X","",colnames(GE_PEER))
colnames(GE_PEER)=gsub(".","-",colnames(GE_PEER),fixed = T)

save(GE,ME,COVA_GE,GE_PEER,COVA_ME,ME_anno,SNP_ALL,SNP_GE,SNP_ME,SNP_GE_PCA,SNP_ME_PCA,hutch_normal_SNP_impscore,tcga_normal_SNP_impscore,SNP_GE_pairs,SNP_ME_pairs,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/normal_SNP_GE_ME.RData")

