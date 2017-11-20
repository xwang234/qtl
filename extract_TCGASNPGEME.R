#!/usr/bin/env Rscript
rm(list=ls())
library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")


#read seminal vesicles info
#load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/TCGAnormals.RData")
normclinical=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/data/nationwidechildrens.org_ssf_normal_controls_prad.txt",skip=1,sep="\t",header=T,stringsAsFactors = F)
normclinical=normclinical[-1,]
idx=match(colnames(genotypedata),normclinical$bcr_patient_barcode)
table(normclinical$normal_tissue_anatomic_site[idx])
seminal_sampleid_bio=normclinical$bcr_patient_barcode[normclinical$normal_tissue_anatomic_site=="Seminal Vesicle"]
#idx=which(colnames(genotypedata) %in% seminal_sampleid_bio) #seminal vesicle samples to be removed

#process gene expression, RNA-Seq, Illumina HiSeq
#to find the sampleID
#the manifest file used to download data:
manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/gdc_manifest_geneexpression.txt"
manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
#the ID column stores file UUID. With it, we can find its caseID or patient UUID using json file (downloaed from GDC, export table)
metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/metadata_geneexpression.json"
metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
#using the case ID, we can find TCGA sample ID using clinical table
clinicfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/clinical_solidtissue/e94c9009-baf3-4c0b-a736-1c5158a6d4b6/nationwidechildrens.org_biospecimen_portion_prad.txt"
clinicdata=read.table(clinicfile,header = T,sep="\t",stringsAsFactors = F)

getTCGAsampleID=function(manifestfile,metadatafile,clinicfile)
{
  manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
  metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
  clinicdata=read.table(clinicfile,header = T,sep="\t",stringsAsFactors = F)
  colnames(manifest)[1]="fileid"
  manifest$caseid=NA
  manifest$sampleid=NA
  for (i in 1:nrow(manifest))
  {
    fileid=manifest$fileid[i]
    idxrow=which(grepl(fileid,metadata[,1]))
    for (j in idxrow:nrow(metadata))
    {
      if (grepl("case_id",metadata[j,1]))
      {
        caseid=gsub("case_id:","",metadata[j,1])
        caseid=gsub(" ","",caseid)
        manifest$caseid[i]=caseid
        break
      }
    }
    idxrow=which(tolower(clinicdata$bcr_patient_uuid)==tolower(caseid))
    if (length(idxrow)>0)
    {
      tmp=unlist(strsplit(clinicdata$bcr_sample_barcode[idxrow][1],"-",fixed = T))
      manifest$sampleid[i]=paste0(tmp[1:3],collapse = "-")
    }
  }
  return(manifest)
}

geneexp_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile)
length(unique(geneexp_manifest$sampleid))
#[1] 52
length(unique(geneexp_manifest$sampleid)[!unique(geneexp_manifest$sampleid) %in% seminal_sampleid_bio])
#[1] 48
tcga_geneexp_samples=unique(geneexp_manifest$sampleid)

#methylation, Illumina Human Methylation 450
manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/gdc_manifest_methylation.txt"
manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/metadata_methylation.json"
metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
methylation_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile)
length(unique(methylation_manifest$sampleid))
#[1] 50
length(unique(methylation_manifest$sampleid)[!unique(methylation_manifest$sampleid) %in% seminal_sampleid_bio])
#[1] 46
sum(unique(geneexp_manifest$sampleid) %in% unique(methylation_manifest$sampleid))
#[1] 35
comsampleid=unique(geneexp_manifest$sampleid)[unique(geneexp_manifest$sampleid) %in% unique(methylation_manifest$sampleid)]
tcga_methylation_samples=unique(methylation_manifest$sampleid)

#genotype/copynumber,Affymetrix SNP Array 6.0
# manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/gdc_manifest_snp6_solidtissue_estimatecopynumber.txt"
# manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
# metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/metadata_snp6_solidtissue_estimatecopynumber.json"
# metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
# genotype_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile)
# sum(comsampleid %in% genotype_manifest$sampleid)
# #[1] 35
# #35 samples have all 3 types data

#solid normals
manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/gdc_manifest_snv_genotype_solidtissue.txt"
manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/metadata_snv_genotype_solidtissue.json"
metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
genotype_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile)
length(unique(genotype_manifest$sampleid))
#[1] 115
sum(comsampleid %in% genotype_manifest$sampleid)
#[1] 35
#35 samples have all 3 types data

#blood normals:
manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/gdc_manifest_snv_genotype_blood.txt"
manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/metadata_snv_genotype_blood.json"
metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
genotype_blood_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/clinical_blood/e94c9009-baf3-4c0b-a736-1c5158a6d4b6/nationwidechildrens.org_biospecimen_portion_prad.txt")
length(unique(genotype_blood_manifest$sampleid))
#[1] 439

#check overlap of solid normals and blood
length(intersect(genotype_manifest$sampleid,genotype_blood_manifest$sampleid))
#[1] 58

#check overlap of geneexp with solid normals genotype
sum(tcga_geneexp_samples %in% genotype_manifest$sampleid)
#[1] 50
sum(tcga_geneexp_samples %in% genotype_blood_manifest$sampleid)
#[1] 46
sum(tcga_geneexp_samples %in% c(genotype_manifest$sampleid,genotype_blood_manifest$sampleid))
#[1] 52
#46 from blood, 6 from solid normals
tcga_blood_geneexp_samples=tcga_geneexp_samples[tcga_geneexp_samples %in% genotype_blood_manifest$sampleid]
tcga_solid_geneexp_samples=tcga_geneexp_samples[!tcga_geneexp_samples %in% tcga_blood_geneexp_samples]

#check overlap of methylation with solid normals genotype
sum(tcga_methylation_samples %in% genotype_manifest$sampleid)
#[1] 50
sum(tcga_methylation_samples %in% genotype_blood_manifest$sampleid)
#[1] 37
sum(tcga_methylation_samples %in% c(genotype_manifest$sampleid,genotype_blood_manifest$sampleid))
# 37 from blood 13 from solid normals
tcga_blood_methylation_samples=tcga_methylation_samples[tcga_methylation_samples %in% genotype_blood_manifest$sampleid]
tcga_solid_methylation_samples=tcga_methylation_samples[!tcga_methylation_samples %in% tcga_blood_methylation_samples]
length(unique(c(tcga_blood_geneexp_samples,tcga_solid_geneexp_samples,tcga_blood_methylation_samples,tcga_solid_methylation_samples)))
#[1] 67 total number of samples
tcga_blood_genotype_samples=unique(c(tcga_blood_geneexp_samples,tcga_blood_methylation_samples))
length(tcga_blood_genotype_samples)
#[1] 54
length(tcga_blood_genotype_samples[!tcga_blood_genotype_samples %in% seminal_sampleid_bio])
#[1] 51
tcga_solid_genotype_samples=unique(c(tcga_solid_geneexp_samples,tcga_solid_methylation_samples))
length(tcga_solid_genotype_samples)
#[1] 13
length(tcga_solid_genotype_samples[!tcga_solid_genotype_samples %in% seminal_sampleid_bio])
#[1] 11
sum(tcga_geneexp_samples %in% c(tcga_blood_genotype_samples,tcga_solid_genotype_samples))
sum(tcga_methylation_samples %in% c(tcga_blood_genotype_samples,tcga_solid_genotype_samples))



length(intersect(comsampleid,genotype_blood_manifest$sampleid))
#[1] 29
head(which(genotype_blood_manifest$sampleid %in% comsampleid))
genotype_blood_manifest$sampleid[6]

#snp6 data
# rawdata=fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/copynumber_solidtissue/estimate/96d803b0-dbd6-45f6-b717-9a6bf44aaf2e/CLEAR_p_TCGA_b91_SNP_1N_GenomeWideSNP_6_A03_755736.raw.copynumber.data.txt",skip=1,sep="\t",header=T,stringsAsFactors = F)
# rawdata=as.data.frame(rawdata)
# tangentdata=fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/copynumber_solidtissue/estimate/8cc69c7b-2324-4800-a898-762bf5b9250a/CLEAR_p_TCGA_b91_SNP_1N_GenomeWideSNP_6_A03_755736.tangent.copynumber.data.txt",skip=1,sep="\t",header=T,stringsAsFactors = F)
# tangentdata=as.data.frame(tangentdata)
# alleldata=fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/copynumber_solidtissue/estimate/e60fe0e5-1b8f-4222-98c3-3a394859092a/CLEAR_p_TCGA_b91_SNP_1N_GenomeWideSNP_6_A03_755736.byallele.copynumber.data.txt",skip=1,sep="\t",header=T,stringsAsFactors = F)
# alleldata=as.data.frame(alleldata)
# sum(alleldata$`Composite Element REF` %in% snp6_anno$Probe.Set.ID)
genotypedata=fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/snv_genotype/0a8ee016-649c-4b67-95d4-63aa18f597bb/CLADE_p_TCGASNP_184_195_N_GenomeWideSNP_6_C03_1039826.birdseed.data.txt",skip=1,sep="\t",header=T,stringsAsFactors = F)
genotypedata=as.data.frame(genotypedata)
sum(genotypedata$`Composite Element REF` %in% snp6_anno$Probe.Set.ID)

#read data from manifest table
#genotype
readdatafrommanifest=function(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/snv_genotype",manifest=NULL,sampleid=NULL)
{
  res=NULL
  idx=which(manifest$sampleid==sampleid)
  if (length(idx)>0)
  {
    thefile=paste0(loc,"/",manifest$fileid[idx],"/",manifest$filename[idx])
    #res1=read.table(thefile,skip=1,header=T,sep="\t",row.names = 1)
    res=fread(thefile,skip=1,header = T,sep="\t")
    res=as.data.frame(res)
    rownames(res)=res[,1]
    res=res[,2:ncol(res)]
  }
  return(res)
}

#geneexp
readdatafrommanifest1=function(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/snv_genotype",manifest=NULL,sampleid=NULL)
{
  res=NULL
  idx=which(manifest$sampleid==sampleid & grepl("normalized",manifest$filename))
  
  if (length(idx)>0)
  {
    thefile=paste0(loc,"/",manifest$fileid[idx],"/",manifest$filename[idx])
    #res=read.table(thefile,header=T,sep="\t",row.names = 1)
    res=fread(thefile,header = T,sep="\t")
    res=as.data.frame(res)
    rownames(res)=res[,1]
    res=res[,2:ncol(res),drop=F]
  }
  return(res)
}

#snp6 data
#first read genotype solid
genotypedata=NULL
for (i in 1:length(tcga_solid_genotype_samples))
{
  tmp=readdatafrommanifest(manifest = genotype_manifest,sampleid=tcga_solid_genotype_samples[i])
  genotypedata=cbind(genotypedata,tmp[,1])
  rownames(genotypedata)=rownames(tmp)
  colnames(genotypedata)[i]=tcga_solid_genotype_samples[i]
}
#then add genotype blood
for (i in 1:length(tcga_blood_genotype_samples))
{
  tmp=readdatafrommanifest(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/genotype_blood",manifest = genotype_blood_manifest,sampleid=tcga_blood_genotype_samples[i])
  genotypedata=cbind(genotypedata,tmp[,1])
  rownames(genotypedata)=rownames(tmp)
  colnames(genotypedata)[i+length(tcga_solid_genotype_samples)]=tcga_blood_genotype_samples[i]
}
genotypedata=as.data.frame(genotypedata)


#add other the blood normals

other_bloodnormals=genotype_blood_manifest$sampleid[! genotype_blood_manifest$sampleid %in% c(tcga_blood_genotype_samples,tcga_solid_genotype_samples)]
other_bloodgenotypedata=NULL
for (i in 1:length(other_bloodnormals))
{
  if (i %% 100==0) cat(i,'..')
  tmp=readdatafrommanifest(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/genotype_blood",manifest = genotype_blood_manifest,sampleid=other_bloodnormals[i])
  other_bloodgenotypedata=cbind(other_bloodgenotypedata,tmp[,1])
  rownames(other_bloodgenotypedata)=rownames(tmp)
  colnames(other_bloodgenotypedata)[i]=other_bloodnormals[i]
}
other_bloodgenotypedata=as.data.frame(other_bloodgenotypedata)
sum(colnames(genotypedata) %in% colnames(other_bloodgenotypedata))
allnormalgenotypedata=cbind(genotypedata,other_bloodgenotypedata)

save(genotypedata,allnormalgenotypedata,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGAnormals.RData")

#geneexpression data
geneexpdata=NULL
for (i in 1:length(c(tcga_solid_geneexp_samples,tcga_blood_geneexp_samples)))
{
  tmp=readdatafrommanifest1(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/geneexpression",manifest = geneexp_manifest,sampleid=c(tcga_solid_geneexp_samples,tcga_blood_geneexp_samples)[i])
  geneexpdata=cbind(geneexpdata,tmp[,1])
  rownames(geneexpdata)=rownames(tmp)
  colnames(geneexpdata)[i]=c(tcga_solid_geneexp_samples,tcga_blood_geneexp_samples)[i]
}
geneexpdata=as.data.frame(geneexpdata)

#methylation data
methydata=NULL
for (i in 1:length(c(tcga_solid_methylation_samples,tcga_blood_methylation_samples)))
{
  tmp=readdatafrommanifest(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/methylation",manifest =methylation_manifest,sampleid=c(tcga_solid_methylation_samples,tcga_blood_methylation_samples)[i])
  methydata=cbind(methydata,tmp[,1])
  rownames(methydata)=rownames(tmp)
  colnames(methydata)[i]=c(tcga_solid_methylation_samples,tcga_blood_methylation_samples)[i]
}
methydata=as.data.frame(methydata)
save(genotypedata,geneexpdata,methydata,allnormalgenotypedata,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGAnormals.RData")
#add probe position info

#genotype
idx=match(rownames(genotypedata),snp6_anno$Probe.Set.ID)
genotypedata=cbind.data.frame(chr=snp6_anno$Chromosome[idx],pos=snp6_anno$Physical.Position[idx],genotypedata)
#geneexpression
rnaseqv2_anno=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/TCGA.hg19.June2011.gaf",header=T,sep="\t",comment.char = "",stringsAsFactors = F)
idx=match(rownames(geneexpdata),rnaseqv2_anno$FeatureID)
#compare with array
load("ge_annohg19.RData")
load("imputation/imputation_GEME.RData")
idx=match(rownames(GE),ge_annohg19$Probe_Id)
ge_array_anno=ge_annohg19[idx,]
ge_rnaseq_genes=unlist(strsplit(rownames(geneexpdata),"|",fixed=T))
ge_rnaseq_genes=ge_rnaseq_genes[seq(1,length(ge_rnaseq_genes),2)]
sum(unique(ge_rnaseq_genes) %in% ge_array_anno$Symbol)
sum(unique(ge_rnaseq_genes) %in% ge_array_anno$Symbol)/length(unique(ge_array_anno$Symbol))
#[1] 0.922678 overlap of RNA-seq and array geneexp on gene level
length(unique(ge_rnaseq_genes))-1 #-"?"
#[1] 20501 number of genes on rnaseq
length(unique(ge_array_anno$Symbol))
#[1] 18895 number of genes on array geneexp


