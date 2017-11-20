#!/usr/bin/env Rscript
#combine methylation,geneexp
rm(list=ls())

#methylation
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGAnormals.RData")
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/normtumor_GEME.RData")
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/imputation_GEME.RData") #load GE, ME all hutch samples
hutch_samples=colnames(normal_ge)[colnames(normal_ge) %in% colnames(GE)]
idx=match(hutch_samples,colnames(normal_methy))
hutch_methy=normal_methy[,idx]
rownames(hutch_methy)=rownames(ME)

#read seminal vesicles info for TCGA
normclinical=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/data/nationwidechildrens.org_ssf_normal_controls_prad.txt",skip=1,sep="\t",header=T,stringsAsFactors = F)
normclinical=normclinical[-1,]
seminal_sampleid_bio=normclinical$bcr_patient_barcode[normclinical$normal_tissue_anatomic_site=="Seminal Vesicle"]
tcga_methylation_samples=colnames(methydata)
tcga_methylation_samples=tcga_methylation_samples[!tcga_methylation_samples %in% seminal_sampleid_bio]
length(tcga_methylation_samples)
#[1] 46 methylation sample size
idx=match(tcga_methylation_samples,colnames(methydata))
tcga_methy=methydata[,idx]
sum(rownames(tcga_methy) %in% rownames(hutch_methy))
idx=match(rownames(hutch_methy),rownames(tcga_methy))
tcga_methy=tcga_methy[idx,]
sum(rownames(tcga_methy) == rownames(hutch_methy))
allnormal_methy=cbind.data.frame(tcga_methy,hutch_methy)

#gene expression
rownames(normal_ge)=rownames(GE)
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/ge_annohg19.RData")
#with nonNA hg19 positions
hutch_ge_pos=read.table("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/GE_POS.txt",header=T,sep="\t",stringsAsFactors = F)
idx=match(hutch_ge_pos$geneid,ge_annohg19$Probe_Id)
hutch_ge_anno=ge_annohg19[idx,]
hutch_ge_anno$start=hutch_ge_pos$s1
hutch_ge_anno$end=hutch_ge_pos$s2
idx=match(hutch_ge_anno$Probe_Id,rownames(normal_ge))
normal_ge1=normal_ge[idx,]
idx=match(hutch_samples,colnames(normal_ge1))
normal_ge2=normal_ge1[,idx]

hutch_ge_acc=sapply(hutch_ge_anno$Accession,function(x){
  tmp=unlist(strsplit(x,".",fixed = T))[1]
})
hutch_ge_anno$acc=hutch_ge_acc

rnaseqv2_anno=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/TCGA.hg19.June2011.gaf",header=T,sep="\t",comment.char = "",stringsAsFactors = F)
tcga_ge_acc=sapply(rownames(geneexpdata),function(x){
  idx=which(rnaseqv2_anno$Gene==x)
  tmp=rnaseqv2_anno$FeatureAliases[idx]
  tmp=tmp[tmp!=""]
  tmp=unique(tmp)
  if (length(tmp)>0)
  {
    tmp=paste0(tmp,collapse = "|")
  }
  return(tmp)
})
tcga_ge_anno=data.frame(matrix(NA,nrow=nrow(geneexpdata),ncol=6))
colnames(tcga_ge_anno)=c("gene","chr","start","end","strand","acc")
tcga_ge_anno$acc=tcga_ge_acc
tcga_ge_gene=unlist(strsplit(rownames(geneexpdata),"|",fixed=T))
tcga_ge_gene=tcga_ge_gene[seq(1,length(tcga_ge_gene),2)]
tcga_ge_gene[tcga_ge_gene=="?"]=NA
tcga_ge_anno$gene=tcga_ge_gene
for (i in 1:nrow(tcga_ge_anno))
{
  if (i %% 1000==0) cat(i,'..')
  idx=which(rnaseqv2_anno$Gene==rownames(geneexpdata)[i])
  locus=rnaseqv2_anno$GeneLocus[idx][1]
  tmp=unlist(strsplit(locus,":"))
  tcga_ge_anno$chr[i]=gsub("chr","",tmp[1])
  tcga_ge_anno$strand[i]=tmp[3]
  tmp=unlist(strsplit(tmp[2],"-"))
  tcga_ge_anno$start[i]=as.integer(tmp[1])
  tcga_ge_anno$end[i]=as.integer(tmp[2])
}
tcga_ge_anno$strand=gsub("+\\S+","+",tcga_ge_anno$strand)
tcga_ge_anno$strand=gsub("-\\S+","-",tcga_ge_anno$strand)
tcga_ge_anno$idx=1:nrow(tcga_ge_anno)

tmp=NULL
for (i in 1:nrow(tcga_ge_anno))
{
  tmp=c(tmp,unlist(strsplit(tcga_ge_anno$acc[i],"\\|")))
}
length(unique(tmp))
length(unique(tcga_ge_anno$acc))
#[1] 20531
length(unique(hutch_ge_anno$acc))
#[1] 20725
length(unique(intersect(tcga_ge_anno$acc,hutch_ge_anno$acc)))
#[1] 12906
length(unique(tcga_ge_anno$gene))
#[1] 20502
length(unique(hutch_ge_anno$Symbol))
#[1] 18077
length(unique(intersect(tcga_ge_anno$gene,hutch_ge_anno$Symbol)))
#17136

library(GenomicRanges)
gr_tcga_ge=GRanges(seqnames = tcga_ge_anno$chr,ranges = IRanges(start=tcga_ge_anno$start,end=tcga_ge_anno$end))
gr_hutch_ge=GRanges(seqnames = hutch_ge_anno$Chromosome,ranges = IRanges(start=hutch_ge_anno$start,end=hutch_ge_anno$end))
tcga_ge_anno$hutchidx=NA
tcga_ge_anno$match=NA
distcutoff=100000 #only combine probes within the cutoff

#use genename
for (i in 1:nrow(tcga_ge_anno))
{
  if (i %% 1000==0) cat(i,'..')
  genename=tcga_ge_anno$gene[i]
  if (!is.na(genename))
  {
    idx=which(hutch_ge_anno$Symbol==genename & hutch_ge_anno$Chromosome==tcga_ge_anno$chr[i]) #make sure the genes in the same chr
    if (length(idx)>0)
    {
      tcga_ge_anno$match[i]="gene"
      if (length(idx)==1)
      {
        tcga_ge_anno$hutchidx[i]=idx
      }else #more hits
      {
        dists=distance(gr_tcga_ge[i],gr_hutch_ge[idx])
        idx1=idx[dists<=distcutoff]
        if (length(idx1)>0) #in same chr
        {
          tcga_ge_anno$hutchidx[i]=paste0(idx1,collapse = "|")
        }else #didn't find the match within cutoff
        {
            idx2=idx[which.min(dists)]
            tcga_ge_anno$hutchidx[i]=idx2
         }
      }
    }
  }
}
which(duplicated(tcga_ge_anno$gene))
# [1]     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21    22
#[22]    23    24    25    26    27    28    29 16302
idx=which(tcga_ge_anno$gene==tcga_ge_anno$gene[16302])
idx
#[1] 16301 16302 #one duplicated gene
for (i in 1:length(idx))
{
  tcga_ge_anno$hutchidx[idx[i]]=which(hutch_ge_anno$acc==tcga_ge_anno$acc[idx[i]])
}

tcga_gene_hutchidx=NULL
for (i in 1:nrow(tcga_ge_anno))
{
  if (!is.na(tcga_ge_anno$match[i]) & tcga_ge_anno$match[i]=="gene") tcga_gene_hutchidx=c(tcga_gene_hutchidx,unlist(strsplit(tcga_ge_anno$hutchidx[i],"\\|")))
}
length(unique(tcga_gene_hutchidx))==length(tcga_gene_hutchidx)
#[1] TRUE

#check if can be improved using Accesssion
sum(tcga_ge_anno$acc[is.na(tcga_ge_anno$hutchidx)] %in% hutch_ge_anno$acc)
#[1] 261
idx=which(is.na(tcga_ge_anno$hutchidx))
tcga_all_hutchidx=tcga_gene_hutchidx
for (i in 1:length(idx))
{
  tmpacc=unlist(strsplit(tcga_ge_anno$acc[idx[i]],"|",fixed=T))
  tmpacc=tmpacc[tmpacc %in% hutch_ge_anno$acc]
  if (length(tmpacc)>0)
  {
    tmp1=NULL
    for (j in 1:length(tmpacc))
    {
      idx1=which(hutch_ge_anno$acc==tmpacc[j] & hutch_ge_anno$Chromosome==tcga_ge_anno$chr[idx[i]])
      if (length(idx1)>0)
      {
        for (idx2 in idx1)
        {
          if (!idx2 %in% tcga_all_hutchidx) #make sure no overlap acc
          {
            tmp1=unique(c(tmp1,idx2))
            tcga_all_hutchidx=c(tcga_all_hutchidx,idx2)
          }
        }
      }
    }
    if (length(tmp1)>0)
    {
      tcga_ge_anno$match[idx[i]]="acc"
      dists=distance(gr_hutch_ge[tmp1],gr_tcga_ge[idx[i]])
      tmp2=tmp1[dists<=distcutoff]
      if (length(tmp2)>0)
      {
        tcga_ge_anno$hutchidx[idx[i]]=paste0(tmp1,collapse = "|")
      }else
      {
        tcga_ge_anno$hutchidx[idx[i]]=tmp1[which.min(dists)]
      }
    }
  }
}

#no duplicated hutchidx
tmp=NULL
for (i in 1:nrow(tcga_ge_anno))
{
  if (!is.na(tcga_ge_anno$match[i])) tmp=c(tmp,unlist(strsplit(tcga_ge_anno$hutchidx[i],"\\|")))
}
length(unique(tmp))==length(tmp)
#[1] TRUE

#number of hutchprobes used
length(tcga_all_hutchidx)
#[1] 22655
#number of probes
sum(!is.na(tcga_ge_anno$match))
#[1] 17457

#check missings
idx=which(is.na(tcga_ge_anno$match))
dists=distanceToNearest(gr_tcga_ge[idx],gr_hutch_ge)

for (i in 1:length(idx))
{
  if (dists@elementMetadata$distance[i]<=distcutoff & !dists@to[i] %in% tcga_all_hutchidx)
  {
    tcga_ge_anno$hutchidx[idx[i]]=dists@to[i]
    tcga_all_hutchidx=c(tcga_all_hutchidx,dists@to[i])
    tcga_ge_anno$match[idx[i]]="distance"
  }
}
length(unique(tcga_all_hutchidx))==length(tcga_all_hutchidx)

length(tcga_all_hutchidx)
#[1] 22867
#number of probes
sum(!is.na(tcga_ge_anno$match))
#[1] 17669

table(tcga_ge_anno$match)
#acc distance     gene 
#336      212    17121 

#number of unique matched probe
sum(!grepl("|",tcga_ge_anno$hutchidx,fixed = T))
#[1] 16913

idx=match(tcga_geneexp_samples,colnames(geneexpdata))
tcga_geneexp=geneexpdata[,idx]
idx=!is.na(tcga_ge_anno$hutchidx)
tcga_geneexp=tcga_geneexp[idx,]

idx=which(!is.na(tcga_ge_anno$hutchidx))
hutch_geneexp=data.frame(matrix(NA,nrow=length(idx),ncol=ncol(normal_ge2)))
colnames(hutch_geneexp)=colnames(normal_ge2)
allnormal_geneexp_pos=data.frame(matrix(NA,nrow=length(idx),ncol=4))
colnames(allnormal_geneexp_pos)=c("geneid","chr","s1","s2")
allnormal_geneexp_pos$chr=tcga_ge_anno$chr[idx]
for (i in 1:length(idx))
{
  if (i %% 1000==0) cat(i,"..")
  tmp=as.integer(unlist(strsplit(tcga_ge_anno$hutchidx[idx[i]],"|",fixed=T)))
  if (length(tmp)==1)
  {
    hutch_geneexp[i,]=normal_ge2[tmp,]
    allnormal_geneexp_pos$geneid[i]=hutch_ge_anno$Probe_Id[tmp]
    allnormal_geneexp_pos$s1[i]=hutch_ge_anno$start[tmp]
    allnormal_geneexp_pos$s2[i]=hutch_ge_anno$end[tmp]
  }else
  {
    hutch_geneexp[i,]=colMeans(normal_ge2[tmp,])
    allnormal_geneexp_pos$geneid[i]=hutch_ge_anno$Probe_Id[tmp][1]
    allnormal_geneexp_pos$s1[i]=min(hutch_ge_anno$start[tmp])
    allnormal_geneexp_pos$s2[i]=max(hutch_ge_anno$end[tmp])
  }
}

tcga_geneexp_samples=tcga_geneexp_samples[!tcga_geneexp_samples %in% seminal_sampleid_bio]
length(tcga_geneexp_samples)
#[1] 48 number of gene exp samples
idx=match(tcga_geneexp_samples,colnames(tcga_geneexp))
tcga_geneexp=tcga_geneexp[,idx]
save(hutch_samples,tcga_methy,hutch_methy,allnormal_methy,allnormal_geneexp_pos,hutch_ge_anno,tcga_ge_anno,hutch_geneexp,tcga_geneexp,
     file="/fh/fast/stanford_j/Xiaoyu/QTL/data/allnormaldata.RData")

allnormal_geneexp=cbind.data.frame(tcga_geneexp,hutch_geneexp)

#normalization
library(sva)
# pheno=read.table("test.txt",sep="\t")
# pheno=pheno[,-1]
pheno=data.frame(matrix(NA,nrow=ncol(allnormal_geneexp),ncol=2))
colnames(pheno)=c("sample","batch")
pheno$sample=colnames(allnormal_geneexp)
pheno$batch=c(rep("snp6",48),rep("rnaseq",14))
batch<-pheno$batch
modcombat<-model.matrix(~1, data=pheno)
combat_allnormal_geneexp= ComBat(dat=allnormal_geneexp, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

library(preprocessCore)
mynormalize=function(data=allnormal_geneexp,targetid=1:48,normid=49:ncol(allnormal_geneexp))
{
  res=data
  for (i in 1:nrow(data))
  {
    if (i %% 5000==0) cat(i,'..')
    if(any(data[i,targetid]!=0)) #for GE
    {
      res[i,normid]=normalize.quantiles.use.target(as.matrix(unlist(data[i,normid])),target=unlist(data[i,targetid]))
    }
  }
  return(res)
}

# plot(unlist(allnormal_methy[1,]))
# points(unlist(allnormal_methy_norm[1,]),col="red")
allnormal_methy_norm=mynormalize(data=allnormal_methy,targetid = 1:46,normid=47:60)
allnormal_geneexp_norm=mynormalize(data=allnormal_geneexp)
save(tcga_methy,hutch_methy,allnormal_methy,allnormal_geneexp_pos,hutch_ge_anno,tcga_ge_anno,hutch_geneexp,tcga_geneexp,
     allnormal_methy_norm,allnormal_geneexp_norm,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/allnormaldata.RData")

load("/fh/fast/stanford_j/Xiaoyu/QTL/data/allnormaldata.RData")
#check batch effect
tmp=cbind.data.frame(tcga_geneexp,hutch_geneexp)
tmp1=prcomp(t(tmp),scale=T,center=T)
plot(tmp1)
plot(tmp1$x[,1],tmp1$x[,2],col=c(rep("red",48),rep("blue",14)),xlab="PC1",ylab="PC2",cex.axis=1.3,cex.lab=1.3)
legend("topright",c("TCGA","Hutch"),col=c("red","blue"),pch = 1,cex=1.3)
tmp1=prcomp(t(allnormal_geneexp_norm),scale=T,center=T)
plot(tmp1$x[,1],tmp1$x[,2],col=c(rep("red",48),rep("blue",14)),xlab="PC1",ylab="PC2",cex.axis=1.3,cex.lab=1.3)
legend("topright",c("TCGA","Hutch"),col=c("red","blue"),pch = 1,cex=1.3)
tmp1=prcomp(t(combat_allnormal_geneexp),scale=T,center=T)
plot(tmp1$x[,1],tmp1$x[,2],col=c(rep("red",48),rep("blue",14)),xlab="PC1",ylab="PC2",cex.axis=1.3,cex.lab=1.3)
legend("topright",c("TCGA","Hutch"),col=c("red","blue"),pch = 1,cex=1.3)
