#the file was copied from /fh/fast/daij/.../Tools/wang/prostate/process_TCGAnormals.R to generate genotypes from all 496 samples 
setwd("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate")
library(GenomicRanges)

#process SNP6
load("../../../Tools/wang/prostate/allhighrisksnps_new.RData")
snp6_anno <- read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
snp6_anno$Physical.Position=as.integer(snp6_anno$Physical.Position)
sum(allsnps$snp %in% snp6_anno$dbSNP.RS.ID)
overlapsnps=allsnps$snp[allsnps$snp %in% snp6_anno$dbSNP.RS.ID]
idx1=match(overlapsnps,allsnps$snp)
View(allsnps[idx1,])
idx1=match(overlapsnps,snp6_anno$dbSNP.RS.ID)
View(snp6_anno[idx1,])
gr_allsnps=GRanges(seqnames = allsnps$chr,ranges=IRanges(start=allsnps$position,width = 1))
idxNoNA=which(!is.na(snp6_anno$Physical.Position))
gr_snp6_anno=GRanges(seqnames = snp6_anno$Chromosome[idxNoNA],ranges = IRanges(start=snp6_anno$Physical.Position[idxNoNA],width=1))
length(intersect(gr_allsnps,gr_snp6_anno))
tmp=distanceToNearest(gr_allsnps,gr_snp6_anno)
quantile(tmp@elementMetadata$distance[tmp@elementMetadata$distance!=0])
hist(tmp@elementMetadata$distance)


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

#methylation, Illumina Human Methylation 450
manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/gdc_manifest_methylation.txt"
manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/metadata_methylation.json"
metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
methylation_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile)
length(unique(methylation_manifest$sampleid))
#[1] 50
sum(unique(geneexp_manifest$sampleid) %in% unique(methylation_manifest$sampleid))
#[1] 35
comsampleid=unique(geneexp_manifest$sampleid)[unique(geneexp_manifest$sampleid) %in% unique(methylation_manifest$sampleid)]

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

length(intersect(genotype_manifest$sampleid,genotype_blood_manifest$sampleid))
#[1] 58
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
    res=read.table(thefile,skip=1,header=T,sep="\t",row.names = 1)
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
    res=read.table(thefile,header=T,sep="\t",row.names = 1)
  }
  return(res)
}

#snp6 data
genotypedata=NULL
for (i in 1:length(comsampleid))
{
  tmp=readdatafrommanifest(manifest = genotype_manifest,sampleid=comsampleid[i])
  genotypedata=cbind(genotypedata,tmp[,1])
  rownames(genotypedata)=rownames(tmp)
  colnames(genotypedata)[i]=comsampleid[i]
}
genotypedata=as.data.frame(genotypedata)

#for all the solid normals
solidnormals=c(genotype_manifest$sampleid[genotype_manifest$sampleid %in% comsampleid],genotype_manifest$sampleid[!genotype_manifest$sampleid %in% comsampleid])
solidgenotypedata=NULL
for (i in 1:length(solidnormals))
{
  tmp=readdatafrommanifest(manifest = genotype_manifest,sampleid=solidnormals[i])
  solidgenotypedata=cbind(solidgenotypedata,tmp[,1])
  rownames(solidgenotypedata)=rownames(tmp)
  colnames(solidgenotypedata)[i]=solidnormals[i]
}
solidgenotypedata=as.data.frame(solidgenotypedata)
idx=match(comsampleid,colnames(solidgenotypedata))
solidgenotypedata=cbind(solidgenotypedata[,idx],solidgenotypedata[,36:ncol(solidgenotypedata)])

#for all the blood normals

bloodnormals=c(genotype_blood_manifest$sampleid[genotype_blood_manifest$sampleid %in% comsampleid],genotype_blood_manifest$sampleid[!genotype_blood_manifest$sampleid %in% comsampleid])
bloodgenotypedata=NULL
for (i in 1:length(bloodnormals))
{
  if (i %% 10==0) cat(i,'..')
  tmp=readdatafrommanifest(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/genotype_blood",manifest = genotype_blood_manifest,sampleid=bloodnormals[i])
  bloodgenotypedata=cbind(bloodgenotypedata,tmp[,1])
  rownames(bloodgenotypedata)=rownames(tmp)
  colnames(bloodgenotypedata)[i]=bloodnormals[i]
}
bloodgenotypedata=as.data.frame(bloodgenotypedata)
length(colnames(solidgenotypedata))
#[1] 115
#firt 115 columns are from solid tissue normals, first 35 are from the comsamples
allnormalgenotypedata=cbind(solidgenotypedata,bloodgenotypedata[,!colnames(bloodgenotypedata) %in% colnames(solidgenotypedata)])


#geneexpression data
geneexpdata=NULL
for (i in 1:length(comsampleid))
{
  tmp=readdatafrommanifest1(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/geneexpression",manifest = geneexp_manifest,sampleid=comsampleid[i])
  geneexpdata=cbind(geneexpdata,tmp[,1])
  rownames(geneexpdata)=rownames(tmp)
  colnames(geneexpdata)[i]=comsampleid[i]
}
geneexpdata=as.data.frame(geneexpdata)

#methylation data
methydata=NULL
for (i in 1:length(comsampleid))
{
  tmp=readdatafrommanifest(loc="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/methylation",manifest =methylation_manifest,sampleid=comsampleid[i])
  methydata=cbind(methydata,tmp[,1])
  rownames(methydata)=rownames(tmp)
  colnames(methydata)[i]=comsampleid[i]
}
methydata=as.data.frame(methydata)
sum(colnames(allnormalgenotypedata)[1:35]==colnames(geneexpdata))
save(genotypedata,geneexpdata,methydata,solidgenotypedata,bloodgenotypedata,allnormalgenotypedata,file="TCGAnormals.RData")
save(all496normalgenotypedata,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGA496normalgenotype.RData")
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

#read seminal vesicles info
#load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/TCGAnormals.RData")
normclinical=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/data/nationwidechildrens.org_ssf_normal_controls_prad.txt",skip=1,sep="\t",header=T,stringsAsFactors = F)
normclinical=normclinical[-1,]
idx=match(colnames(genotypedata),normclinical$bcr_patient_barcode)
table(normclinical$normal_tissue_anatomic_site[idx])
seminal_sampleid_bio=normclinical$bcr_patient_barcode[normclinical$normal_tissue_anatomic_site=="Seminal Vesicle"]
idx=which(colnames(genotypedata) %in% seminal_sampleid_bio) #seminal vesicle samples to be removed
#[1]  6 15 25

# manifestfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/gdc_manifest_clinical_xml.txt"
# manifest=read.table(file=manifestfile,header=T,sep="\t",stringsAsFactors=F)
# metadatafile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/metadata_clinical_xml.json"
# metadata=read.table(metadatafile,header=F,sep="\t",stringsAsFactors=F)
# clinical_xml_manifest=getTCGAsampleID(manifestfile,metadatafile,clinicfile="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/clinical_blood/e94c9009-baf3-4c0b-a736-1c5158a6d4b6/nationwidechildrens.org_biospecimen_portion_prad.txt")
# length(unique(genotype_blood_manifest$sampleid))
# seminal_sampleid_bio=normclinical$bcr_patient_barcode[normclinical$normal_tissue_anatomic_site=="Seminal Vesicle"]
# length(seminal_sampleid_bio)
# length(unique(seminal_sampleid_bio))

#using xml files
folder="/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/normals/clinical_SSF"
idxseminal=NULL
xmlsampleid=NULL
tmp=system(paste0("grep -i -c Seminal ",folder,"/","*/*.xml"),intern = TRUE)
for (i in 1:length(tmp))
{
    tmp1=unlist(strsplit(tmp[i],":"))
    if (tmp1[2]!="0") idxseminal=c(idxseminal,i)
    tmp1=unlist(strsplit(tmp[i],".",fixed = T))
    xmlsampleid=c(xmlsampleid,tmp1[3])
}
length(idxseminal)
length(unique(xmlsampleid[idxseminal]))
seminal_sampleid_xml=unique(xmlsampleid[idxseminal])
sum(seminal_sampleid_xml %in% seminal_sampleid_bio)

