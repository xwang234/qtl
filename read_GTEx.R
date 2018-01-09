rm(list=ls())
library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
test=fread("../data/GTEx/59597/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000219.v4.GTEx_Pilot_Imputation.genotype-imputed-data.c1/allSNPs/chr1.dosage")
#dosage data
test1=fread("../data/GTEx/59597/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/test.txt",sep="\t")
gtex_ge_anno=fread("../data/GTEx/59597/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/ExpressionFiles/phe000020.v1.GTEx_RNAseq.marker-info.MULTI/gencode.v19.genes.v7.patched_contigs.gtf",sep="\t")
gtex_ge_anno=as.data.frame(gtex_ge_anno)
tmp=strsplit(gtex_ge_anno$V9,";",fixed = T)
gtex_ge_anno$Probe_Id=NA
gtex_ge_anno$Symbol=NA
for (i in 1:nrow(gtex_ge_anno))
{
  if (i %% 10000==0) cat(i,"..")
  tmp1=tmp[[i]]
  tmp2=unlist(strsplit(tmp1[1],'\"',fixed=T))
  tmp3=unlist(strsplit(tmp1[5],'\"',fixed=T))
  gtex_ge_anno$Probe_Id[i]=tmp2[2]
  gtex_ge_anno$Symbol[i]=tmp3[2]
}
colnames(gtex_ge_anno)[1]="Chromosome"
colnames(gtex_ge_anno)[4]="start"
colnames(gtex_ge_anno)[5]="end"
save(gtex_ge_anno,file="../data/GTEx/gtex_ge_anno.RData")

#Prostate samples---------------
gtex_samples=fread("../data/GTEx_v7_Annotations_SampleAttributesDS.txt",sep="\t")
gtex_samples=as.data.frame(gtex_samples)
gtex_samples=gtex_samples$SAMPID[gtex_samples$SMTS=="Prostate"]
gtex_subjects=rep(NA,length(gtex_samples))
for (i in 1:length(gtex_samples))
{
  tmp=unlist(strsplit(gtex_samples[i],"-"))
  gtex_subjects[i]=paste0(tmp[1:2],collapse = "-")
}

#expression data----------------
test=fread("../data/GTEx/59597/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/test.txt",skip=2,sep="\t")
sum(colnames(test) %in% gtex_samples)
#[1] 158
gtex_ge_samples=gtex_samples[gtex_samples %in% colnames(test)]
gtex_ge_subjects=rep(NA,length(gtex_ge_samples))
for (i in 1:length(gtex_ge_samples))
{
  tmp=unlist(strsplit(gtex_ge_samples[i],"-"))
  gtex_ge_subjects[i]=paste0(tmp[1:2],collapse = "-")
}

gtex_GE=fread("../data/GTEx/59597/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct",
               skip=2,sep="\t",select=c(colnames(test)[1:2],gtex_ge_samples))
gtex_GE=as.data.frame(gtex_GE)
rownames(gtex_GE)=gtex_GE$Name
gtex_GE=gtex_GE[,-1]
colnames(gtex_GE)[1]="Symbol"
rmconstrows=function(dat)
{
  idxconst=rep(F,nrow(dat))
  for (i in 1:nrow(dat))
  {
    if (var(unlist(dat[i,]))==0) idxconst[i]=T
  }
  dat=dat[!idxconst,]
  return(res=list(dat=dat,idxconst=idxconst))
}
tmp=rmconstrows(dat=gtex_GE[,2:ncol(gtex_GE)])
gtex_GE=gtex_GE[!tmp$idxconst,]

# #standardize GE
# tmp=t(scale(t(gtex_GE[,2:ncol(gtex_GE)])))
# sum(colnames(tmp)==colnames(gtex_GE)[2:ncol(gtex_GE)])
# #[1] 158
# gtex_GE[,2:ncol(gtex_GE)]=tmp

#genotype data--------------------
#mid_point
test=fread("../data/GTEx/59597/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/test.txt",sep="\t")
midpoint_samples=colnames(test)[10:ncol(test)]
midpoint_subjects=rep(NA,length(midpoint_samples))
for (i in 1:length(midpoint_samples))
{
  tmp=unlist(strsplit(midpoint_samples[i],"-"))
  midpoint_subjects[i]=paste0(tmp[1:2],collapse = "-")
}
sum(gtex_subjects %in% midpoint_subjects)
sum(gtex_ge_subjects %in% midpoint_subjects)
#[1] 98
#pilot
test1=read.table("../data/GTEx/59597/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000219.v4.GTEx_Pilot_Imputation.genotype-imputed-data.c1/data_used_for_imputation/data_hg19_1KGZ.fam",stringsAsFactors = F)
pilot_samples=test1$V1
pilot_subjects=rep(NA,length(pilot_samples))
for (i in 1:length(pilot_samples))
{
  tmp=unlist(strsplit(pilot_samples[i],"-"))
  pilot_subjects[i]=paste0(tmp[1:2],collapse = "-")
}
#WGS
test2=fread("../data/GTEx/59597/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/test.txt",skip=224)
wgs_subjects=colnames(test2)[10:ncol(test2)]
sum(wgs_subjects %in% gtex_ge_subjects)
#[1] 133
sum(midpoint_subjects %in% wgs_subjects)
sum(unique(c(wgs_subjects,midpoint_subjects)) %in% gtex_ge_subjects)
#[1] 141
#read genotype data from WGS
snpfile="../data/GTEx/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1.GRU_gtex_highrisk.vcf"
snpdata=read.table(snpfile,sep="\t",stringsAsFactors = F)
snpdata$V1=gsub("X",23,snpdata$V1)
tmp=read.table("../data/RESUB_Supplementary_Table16_v9.txt",header=T,sep="\t",stringsAsFactors = F)
load("../result/TCGA_tumors_allhighrisksnps_new.RData")
gtex_highrisk=data.frame(matrix(NA,nrow=nrow(tmp),ncol=5+ncol(test2)-9))
colnames(gtex_highrisk)[6:ncol(gtex_highrisk)]=wgs_subjects
rownames(gtex_highrisk)=tmp$SNP
gtex_highrisk[,1]=tmp$Chr
gtex_highrisk[,3]=tmp$Position
for (i in 1:nrow(tmp))
{
  idx=which(snpdata$V1==tmp$Chr[i] & snpdata$V2==tmp$Position[i])
  if (length(idx)>0)
  {
    if (length(idx)==1)
    {
      gtex_highrisk[i,2]=snpdata$V3[idx]
      gtex_highrisk[i,4]=snpdata$V4[idx]
      gtex_highrisk[i,5]=snpdata$V5[idx]
    }else
    {
      idx1=which(highrisk$chr==tmp$Chr[i] & highrisk$V2==tmp$Position[i])
      idx2=which(snpdata$V4[idx]==highrisk$V3[idx1]& snpdata$V5[idx]==highrisk$V4[idx1])
      idx=idx[idx2]
      if (length(idx)>0)
      {
        gtex_highrisk[i,2]=snpdata$V3[idx]
        gtex_highrisk[i,4]=snpdata$V4[idx]
        gtex_highrisk[i,5]=snpdata$V5[idx]
      }
    }
    if (length(idx)==1)
    {
      for (j in 6:ncol(gtex_highrisk))
      {
        tmp1=snpdata[idx,j+4]
        tmp2=unlist(strsplit(tmp1,":"))
        tmp2=as.numeric(unlist(strsplit(tmp2[5],",")))
        #https://software.broadinstitute.org/gatk/documentation/article?id=5913
        tmp2=10^-(tmp2/10)
        gtprob=tmp2/sum(tmp2)
        dosage=round(gtprob[2]+gtprob[3]*2,digits = 4)
        gtex_highrisk[i,j]=dosage
      }
    }
  }
}
gtex_highrisk$X1=gsub(23,"X",gtex_highrisk$X1)
#make the order of 147 snp consistent with the previous data
idx=rep(NA,nrow(highrisk))
for (i in 1:length(idx))
{
  idx[i]=which(gtex_highrisk$X1==highrisk$chr[i] & gtex_highrisk$X3==highrisk$V2[i])
}
gtex_highrisk=gtex_highrisk[idx,]
colnames(gtex_highrisk)[1:5]=c("chr","V1","V2","V3","V4")
save(gtex_ge_anno,gtex_samples,gtex_subjects,gtex_ge_samples,gtex_ge_subjects,gtex_GE,gtex_highrisk,file="../data/GTEx/highrisk_foreqtl.RData")
sum(gtex_subjects %in% colnames(gtex_highrisk))
#[1] 136
length(unique(gtex_subjects[gtex_subjects %in% colnames(gtex_highrisk)]))
