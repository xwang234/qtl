#!/usr/bin/env Rscript
#generate sample file used for shapeit on chrX
samplefile="../result/imputation_tbd/tbd_allsamplefile.txt"
tmp=data.frame(matrix(NA,nrow=471+1,ncol=4),stringsAsFactors = F)
colnames(tmp)=c("ID_1","ID_2","missing","sex")
tmp[,1]=tmp[,2]=0:(ncol(genotypedata))
tmp[,3]=0
tmp[,4]=2
tmp[1,]=c(0,0,0,"D")
for (i in 2:nrow(tmp))
{
  idx=which(clinicdata$bcr_patient_barcode==colnames(genotypedata)[i-1])
  if (clinicdata$gender[idx]=="MALE")
  {
    tmp[i,4]=1
  }
}
write.table(tmp,file=samplefile,row.names = F,sep=" ",quote=F)


library(data.table)
library(snpStats)
bedfile="/fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/PhenoGenotypeFiles/RootStudyConsentSet_phs000985.ProstateCancer_RiskSNPs.v1.p1.c1.DS-PC-PUB-MDS/GenotypeFiles/phg000725.v1.NCI_ProstateCancer_RiskSNPs.genotype-calls-matrixfmt.HumanOmni2_5-8v1_A.c1.DS-PC-PUB-MDS/genotypes.bed"
bimfile="/fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/PhenoGenotypeFiles/RootStudyConsentSet_phs000985.ProstateCancer_RiskSNPs.v1.p1.c1.DS-PC-PUB-MDS/GenotypeFiles/phg000725.v1.NCI_ProstateCancer_RiskSNPs.genotype-calls-matrixfmt.HumanOmni2_5-8v1_A.c1.DS-PC-PUB-MDS/genotypes.bim"
famfile="/fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/PhenoGenotypeFiles/RootStudyConsentSet_phs000985.ProstateCancer_RiskSNPs.v1.p1.c1.DS-PC-PUB-MDS/GenotypeFiles/phg000725.v1.NCI_ProstateCancer_RiskSNPs.genotype-calls-matrixfmt.HumanOmni2_5-8v1_A.c1.DS-PC-PUB-MDS/genotypes.fam"
tmp=read.plink(bed=bedfile,bim=bimfile,fam=famfile)
genotypedata=t(as.data.frame(tmp$genotypes@.Data))
geneexpdata=fread("/fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/PhenoGenotypeFiles/RootStudyConsentSet_phs000985.ProstateCancer_RiskSNPs.v1.p1.c1.DS-PC-PUB-MDS/ExpressionFiles/phe000016.v1.NCI_ProstateCancer_RiskSNPs.expression-data-matrixfmt.RNAseq_probe_set_grc37.c1.DS-PC-PUB-MDS/normalizedcounts.txt",header=T)
geneexpdata=as.data.frame(geneexpdata)
rownames(geneexpdata)=geneexpdata[,1]
geneexpdata=geneexpdata[,-1]
geneexp_sample=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/PhenoGenotypeFiles/RootStudyConsentSet_phs000985.ProstateCancer_RiskSNPs.v1.p1.c1.DS-PC-PUB-MDS/ExpressionFiles/phe000016.v1.NCI_ProstateCancer_RiskSNPs.sample-info.MULTI/_fam_sub_mot_fat_sam_sex_con_con_con_use_sra_twn-471_expression",stringsAsFactors = F)
clinicdata=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/PhenoGenotypeFiles/RootStudyConsentSet_phs000985.ProstateCancer_RiskSNPs.v1.p1.c1.DS-PC-PUB-MDS/PhenotypeFiles/phs000985.v1.pht005103.v1.p1.c1.PrCa_Risk_SNPs_Subject_Phenotypes.DS-PC-PUB-MDS.txt",sep="\t",skip=10,header=T,stringsAsFactors = F)
genotype_sample=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/PhenoGenotypeFiles/RootStudyConsentSet_phs000985.ProstateCancer_RiskSNPs.v1.p1.c1.DS-PC-PUB-MDS/GenotypeFiles/phg000725.v1.NCI_ProstateCancer_RiskSNPs.sample-info.MULTI/_fam_sub_mot_fat_sam_sex_con_con_con_use_sra_twn-471_genotypes",stringsAsFactors = F)
save(geneexpdata,geneexp_sample,clinicdata,genotypedata,genotype_sample,file="../data/TBD.RData")

