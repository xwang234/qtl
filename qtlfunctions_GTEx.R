#!/usr/bin/env Rscript
#filter SNP data
library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
source("functions.R")

#GE
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/highrisk_foreqtl.RData")
#match sample names in GE and highrisk
idx=which(gtex_ge_subjects %in% colnames(gtex_highrisk))
gtex_ge_samples1=gtex_ge_samples[idx]
gtex_ge_subjects1=gtex_ge_subjects[idx]
sum(duplicated(gtex_ge_subjects1))
#[1] 2
#remove duplicated GE samples
idx=which(duplicated(gtex_ge_subjects1))
gtex_ge_samples2=gtex_ge_samples1[-idx]
gtex_ge_subjects2=gtex_ge_subjects1[-idx]
idx=match(gtex_ge_samples2,colnames(gtex_GE))
GE=gtex_GE[,idx]
GE=removeconstrows(dat=GE)
#GE=imputemissdat(dat=GE)
#standardize
GE=t(scale(t(GE)))
tmp=cbind.data.frame(id=rownames(GE),GE)
write.table(tmp,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE.txt",col.names = T,row.names = F,sep="\t",quote=F)
#GE position
GE_POS=data.frame(matrix(NA,nrow=nrow(GE),ncol=4))
colnames(GE_POS)=c("geneid","chr","s1","s2")
GE_POS$geneid=rownames(GE)
idx=match(rownames(GE),gtex_ge_anno$Probe_Id)
GE_POS$chr=gtex_ge_anno$Chromosome[idx]
GE_POS$s1=gtex_ge_anno$start[idx]
GE_POS$s2=gtex_ge_anno$end[idx]
write.table(GE_POS,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE_POS.txt",col.names = T,row.names = F,sep="\t",quote=F)


GE_PEER=peer_number(dat=GE)#15,0.354
write.table(GE_PEER,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE_PEER.txt",col.names = T,row.names = F,sep="\t",quote=F)

#highrisk SNPs
idx=match(gtex_ge_subjects2,colnames(gtex_highrisk))
gtex_highrisk_SNP_GE=cbind.data.frame(id=gtex_highrisk$V1,gtex_highrisk[,idx])
gtex_highrisk_SNP_GE$id=as.character(gtex_highrisk_SNP_GE$id)
write.table(gtex_highrisk_SNP_GE,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_highrisk_SNP_GE.txt",row.names = F,col.names = T,sep="\t",quote=F)


