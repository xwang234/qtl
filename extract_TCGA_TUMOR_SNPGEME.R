#!/usr/bin/env Rscript
library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
#using the case ID, we can find TCGA sample ID using clinical table
clinicfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/clinical/1f86efca-b015-4aec-83f9-fbcedc5bcd0a/nationwidechildrens.org_biospecimen_aliquot_prad.txt"

getTCGAsampleID=function(manifestfile,metadatafile,clinicfile,colfileid=2)
{
  manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
  metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
  clinicdata=read.table(clinicfile,header = T,sep="\t",stringsAsFactors = F)
  #colnames(manifest)[1]="fileid"
  colnames(manifest)[colfileid]="fileid"
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

manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/gdc_manifest_geneexpression_tumor.txt"
#manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
#the ID column stores file UUID. With it, we can find its caseID or patient UUID using json file (downloaed from GDC, export table)
metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/metadata_geneexpression_tumor.json"
#metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)

geneexp_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile)
length(unique(geneexp_manifest$sampleid))
#[1] 497

#methylation, Illumina Human Methylation 450
manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/gdc_manifest_methylation_tumor.txt"
#manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/metadata_methylation_tumor.json"
#metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
methylation_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile)
length(unique(methylation_manifest$sampleid))
#[1] 498
#there are 4 duplicated samples -01A- and -01B-, pick -01A-
idx=rep(T,length(unique(methylation_manifest$sampleid)))
idx1=which(duplicated(methylation_manifest$sampleid))
for (i in idx1)
{
  idx2=which(methylation_manifest$sampleid==methylation_manifest$sampleid[i])
  portionA=grepl("-01A-",methylation_manifest$fileid[idx2])
  idx[idx2[!portionA]]=F
}
methylation_manifest=methylation_manifest[idx,]

sum(unique(geneexp_manifest$sampleid) %in% unique(methylation_manifest$sampleid))
#[1] 497
comsampleid=unique(geneexp_manifest$sampleid)[unique(geneexp_manifest$sampleid) %in% unique(methylation_manifest$sampleid)]

#genotype
manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/gdc_manifest_snv_genotype_tumor.txt"
#manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/metadata_snv_genotype_tumor.json"
#metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
genotype_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile)
idx=!grepl("FFPE",genotype_manifest$fileid)
genotype_manifest=genotype_manifest[idx,]
length(unique(genotype_manifest$sampleid))

manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/gdc_manifest_snv_genotype_blood.txt"
#manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/metadata_snv_genotype_blood.json"
#metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
genotype_blood_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile)
length(unique(genotype_blood_manifest$sampleid))

#read data from manifest table
#genotype
readdatafrommanifest=function(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/snv_genotype",manifest=NULL,sampleid=NULL)
{
  res=NULL
  idx=which(manifest$sampleid==sampleid)
  if (length(idx)>0)
  {
    thefile=paste0(loc,"/",manifest[idx,1],"/",manifest[idx,2])
    #res1=read.table(thefile,skip=1,header=T,sep="\t",row.names = 1)
    res=fread(thefile,skip=1,header = T,sep="\t")
    res=as.data.frame(res)
    rownames(res)=res[,1]
    res=res[,2:ncol(res)]
  }
  return(res)
}

#geneexp
readdatafrommanifest1=function(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/snv_genotype",manifest=NULL,sampleid=NULL)
{
  res=NULL
  idx=which(manifest$sampleid==sampleid & grepl("normalized",manifest[,2]))
  
  if (length(idx)>0)
  {
    thefile=paste0(loc,"/",manifest[idx,1],"/",manifest[idx,2])
    #res=read.table(thefile,header=T,sep="\t",row.names = 1)
    res=fread(thefile,header = T,sep="\t")
    res=as.data.frame(res)
    rownames(res)=res[,1]
    res=res[,2:ncol(res),drop=F]
  }
  return(res)
}

#snp6 data
tcga_tumor_genotype_samples=genotype_manifest$sampleid
genotypedata=NULL
for (i in 1:length(tcga_tumor_genotype_samples))
{
  if (i %%50==0) cat(i,"..")
  tmp=readdatafrommanifest(manifest = genotype_manifest,sampleid=tcga_tumor_genotype_samples[i])
  genotypedata=cbind(genotypedata,tmp[,1])
  rownames(genotypedata)=rownames(tmp)
  colnames(genotypedata)[i]=tcga_tumor_genotype_samples[i]
}
tcga_genotype_blood_samples=genotype_blood_manifest$sampleid
blood_genotypedata=NULL
for (i in 1:length(tcga_genotype_blood_samples))
{
  if (i %%50==0) cat(i,"..")
  tmp=readdatafrommanifest(loc = "/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/genotype_blood",
                           manifest = genotype_blood_manifest,sampleid=tcga_genotype_blood_samples[i])
  blood_genotypedata=cbind(blood_genotypedata,tmp[,1])
  rownames(blood_genotypedata)=rownames(tmp)
  colnames(blood_genotypedata)[i]=tcga_genotype_blood_samples[i]
}
blood_genotypedata=as.data.frame(blood_genotypedata)



#geneexpression data
geneexpdata=NULL
tcga_tumor_geneexp_samples=unique(geneexp_manifest$sampleid)

for (i in 1:length(tcga_tumor_geneexp_samples))
{
  if (i %%50==0) cat(i,"..")
  tmp=readdatafrommanifest1(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/geneexpression",manifest = geneexp_manifest,sampleid=tcga_tumor_geneexp_samples[i])
  geneexpdata=cbind(geneexpdata,tmp[,1])
  rownames(geneexpdata)=rownames(tmp)
  colnames(geneexpdata)[i]=tcga_tumor_geneexp_samples[i]
}
geneexpdata=as.data.frame(geneexpdata)

#methylation data
methydata=NULL
tcga_tumor_methylation_samples=methylation_manifest$sampleid
for (i in 1:length(tcga_tumor_methylation_samples))
{
  if (i %%50==0) cat(i,"..")
  tmp=readdatafrommanifest(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/methylation",manifest =methylation_manifest,sampleid=tcga_tumor_methylation_samples[i])
  methydata=cbind(methydata,tmp[,1])
  rownames(methydata)=rownames(tmp)
  colnames(methydata)[i]=tcga_tumor_methylation_samples[i]
}
methydata=as.data.frame(methydata)
save(genotypedata,blood_genotypedata,geneexpdata,methydata,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGAtumors.RData")

idx=match(colnames(blood_genotypedata),colnames(genotypedata))
tumor_genotypedata=genotypedata[,idx]
for (i in 1:10)
{
  print(sum(tumor_genotypedata[,i]==blood_genotypedata[,i])/nrow(tumor_genotypedata))
  #print(cor(tumor_genotypedata[,i],blood_genotypedata[,i]))
}

