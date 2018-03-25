#!/usr/bin/env Rscript

source("functions.R")

#for HUTCH highrisk eQTL
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_PCA.txt")
sum(colnames(snp_pca)!=colnames(pheno_pca))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_pca)
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_used.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk.RData")

#use peer factor
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_peer=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_PEER.txt")
sum(colnames(snp_pca)!=colnames(pheno_peer))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_peer)
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_PEER_used.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_trans",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_peer_highrisk.RData")


#for HUTCH highrisk mQTL
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_ME_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_PCA.txt")
sum(colnames(snp_pca)!=colnames(pheno_pca))
#includ age
covariate=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME.txt")
sum(colnames(snp_pca)!=colnames(covariate))
#use top 3 snp_pca
covariate=rbind(covariate[1,],rbind(snp_pca[1:3,],pheno_pca))
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME_used.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_ME.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
       cutoff_cis=1,cutoff_trans=0.1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk.RData")


#for TCGA_tumors highrisk eQTL
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_PCA.txt")
sum(colnames(snp_pca)!=colnames(pheno_pca))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_pca)
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_GE_used.txt")
#udpate the 147snp name
snpspos = read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt", header = TRUE, stringsAsFactors = FALSE)
snp_pheno=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_GE.txt")
snp_pheno$id=snpspos$snp
write_qtl_input(snp_pheno,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_GE_updatename.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_GE_updatename.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_GE_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk_trans",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk.RData")

#use peer factor
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_peer=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_PEER.txt")
sum(colnames(snp_pca)!=colnames(pheno_peer))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_peer)
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_GE_PEER_used.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_GE_updatename.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_GE_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk_peer_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk_peer_trans",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_peer_highrisk.RData")

#for TCGA_tumors highrisk mQTL
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_SNP_ME_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_PCA.txt")
sum(colnames(snp_pca)!=colnames(pheno_pca))
#includ age
load("../data/clinical.RData")
tcgaclinical$age=-as.numeric(tcgaclinical$days_to_birth)/365
idx=match(colnames(snp_pca)[2:ncol(snp_pca)],tcgaclinical$bcr_patient_barcode)
#use top 3 snp_pca
covariate=rbind(c("age",tcgaclinical$age[idx]),rbind(snp_pca[1:3,],pheno_pca))
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_ME_used.txt")
snp_pheno=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_ME.txt")
snp_pheno$id=snpspos$snp
write_qtl_input(snp_pheno,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_ME_updatename.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_ME_updatename.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_ME_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/mqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/mqtl_highrisk_trans",
       cutoff_cis=1,cutoff_trans=0.1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/mqtl_highrisk.RData")

#for NORMAL highrisk eQTL
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_PCA.txt")
sum(colnames(snp_pca)!=colnames(pheno_pca))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_pca[1:15,])
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_GE_used.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_GE_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_trans",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk.RData")
#use peer factor
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_peer=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_PEER.txt")
sum(colnames(snp_pca)!=colnames(pheno_peer))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_peer)
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_GE_PEER_used.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_GE_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_peer_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_peer_trans",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_peer_highrisk.RData")


#for NORMAL highrisk mQTL
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_ME_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_PCA.txt")
sum(colnames(snp_pca)!=colnames(pheno_pca))
#includ age
covariate=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_ME.txt")
sum(colnames(snp_pca)!=colnames(covariate))
#use top 3 snp_pca
covariate=rbind(covariate[1,],rbind(snp_pca[1:3,],pheno_pca[1:15,]))
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_ME_used.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_ME.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_ME_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_trans",
       cutoff_cis=1,cutoff_trans=0.1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk.RData")

#GTEx
#udpate the 147snp name
snpspos = read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt", header = TRUE, stringsAsFactors = FALSE)
snp_pheno=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_highrisk_SNP_GE.txt")
snp_pheno$id=snpspos$snp
write_qtl_input(snp_pheno,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_highrisk_SNP_GE_updatename.txt")

do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_highrisk_SNP_GE_updatename.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE_PEER.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/gtex/eqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/gtex/eqtl_highrisk_trans",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/gtex/eqtl_highrisk.RData")


#TBD
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_peer=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_PEER.txt")
sum(colnames(snp_pca)!=colnames(pheno_peer))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_peer)
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_COVA_GE_PEER_used.txt")
snpspos = read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt", header = TRUE, stringsAsFactors = FALSE)
snp_pheno=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_highrisk_SNP_GE.txt")
snp_pheno$id=snpspos$snp
write_qtl_input(snp_pheno,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_highrisk_SNP_GE_updatename.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_highrisk_SNP_GE_updatename.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_COVA_GE_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_highrisk_peer_cis_all",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_highrisk_peer_trans_all",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_highrisk_all.RData")

do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_COVA_GE_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_trans",
       cutoff_cis=1e-4,cutoff_trans=1e-8,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer.RData")
# do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_GE.txt",
#        snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_POS.txt",
#        phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE.txt",
#        phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",
#        covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_COVA_GE_PEER_used.txt",
#        output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_cis",
#        output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_trans",
#        cutoff_cis=1,cutoff_trans=1e-7,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer.RData")
#salloc -t 1-1 -n 24 mpirun -n 1 /app/easybuild/software/R/3.3.3-foss-2016b-fh1/bin/R --interactive
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
res=mpi_do_qtl_chr(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_GE.txt",
                        snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_POS.txt",
                        phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE.txt",
                        phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",
                        covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_COVA_GE_PEER_used.txt",
                        output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_cis",
                        output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_trans",
                        cutoff_cis=1,cutoff_trans=1e-6,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer.RData")
do_qtl_chrs(chrs=10:23,snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_GE.txt",
               snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_POS.txt",
               phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE.txt",
               phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",
               covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_COVA_GE_PEER_used.txt",
               output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_cis",
               output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer_trans",
               cutoff_cis=1,cutoff_trans=1e-6,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD/eqtl_peer.RData")



#HUTCH GE gene
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_peer=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_PEER.txt")
sum(colnames(snp_pca)!=colnames(pheno_peer))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_peer)
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_gene_PEER_used.txt")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_gene_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_trans",
       cutoff_cis=1e-4,cutoff_trans=1e-8,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer.RData")
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_gene_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_trans",
       cutoff_cis=1,cutoff_trans=1e-7,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer.RData")

do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_gene_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_highrisk_peer_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_highrisk_peer_trans",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_highrisk.RData")
res=mpi_do_qtl_chr(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",
                   snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                   phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene.txt",
                   phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_POS.txt",
                   covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_gene_PEER_used.txt",
                   output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_cis",
                   output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_trans",
                   cutoff_cis=1,cutoff_trans=1e-6,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer.RData")
do_qtl_chrs(chrs=1:23,snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",
                   snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                   phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene.txt",
                   phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_gene_POS.txt",
                   covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_gene_PEER_used.txt",
                   output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_cis",
                   output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer_trans",
                   cutoff_cis=1,cutoff_trans=1e-6,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_gene_peer.RData")

#TCGA normals
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_normalsonly_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_peer=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_normalsonly_GE_PEER.txt")
sum(colnames(snp_pca)!=colnames(pheno_peer))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_peer)
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_normalsonly_COVA_GE_PEER_used.txt")
do_qtl_chrs(chrs=1:23,snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_normalsonly_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_normalsonly_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_normalsonly_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_normalsonly_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_normalsonly_COVA_GE_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_normals/eqtl_normalsonly_peer_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_normals/eqtl_normalsonly_peer_trans",
       cutoff_cis=1,cutoff_trans=1e-6,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_normals/eqtl_normalsonly_peer.RData")

#TCGA tumors
snp_pca=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumorsonly_SNP_GE_PCA.txt")
snp_pca$id=paste0("snp_",snp_pca$id)
pheno_peer=read_qtl_input("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumorsonly_GE_PEER.txt")
sum(colnames(snp_pca)!=colnames(pheno_peer))
#use top 3 snp_pca
covariate=rbind(snp_pca[1:3,],pheno_peer)
write_qtl_input(covariate,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumorsonly_COVA_GE_PEER_used.txt")
do_qtl_chrs(chrs=1:23,snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumorsonly_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumorsonly_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumorsonly_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumorsonly_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumorsonly_COVA_GE_PEER_used.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_tumorsonly_peer_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_tumorsonly_peer_trans",
       cutoff_cis=1,cutoff_trans=1e-6,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_tumorsonly_peer.RData")
