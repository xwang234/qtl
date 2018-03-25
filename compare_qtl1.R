library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library(sas7bdat)
source("functions.R")
highrisk_HUTCH_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_cis",
                                   snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                   geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
highrisk_HUTCH_cis_eqtl=addgenename_GE(highrisk_HUTCH_cis_eqtl,anno=hutch_ge_anno)
#eqtl-gene result
highrisk_HUTCH_cis_eqtl=unique_eqtlpairs(dat1=highrisk_HUTCH_cis_eqtl)
# [1] "number of pairs: 55"
# [1] "number of unique SNPs: 40"
# [1] "number of unique genes: 52"
# [1] "snp: rs12500426,PDLIM5|BMPR1B"
# [1] "snp: rs6062509,ZGPAT|LIME1"
# [1] "snp: rs3850699,NT5C2|TMEM180|C10orf32"
# [1] "snp: rs3129859,TAPBP|HLA-DRB3|HLA-DRB5|HLA-DRB4"
# [1] "snp: rs9364554,SLC22A3|SLC22A2"
# [1] "snp: rs6465657,TECPR1|LMTK2|BHLHA15"
# [1] "snp: rs3096702,HLA-DQB2|HLA-DRB4|HLA-DQB1|HLA-DRB5"
# [1] "snp: rs7767188,HLA-A|HCG4"
# [1] "snp: rs684232,FAM57A|VPS53"

length(unique(hutch_ge_anno$Symbol))
#[1] 18077

highrisk_HUTCH_cis_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_cis",
                                       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                       geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",fdrcutoff = 1)
highrisk_HUTCH_cis_eqtl_all=addgenename_GE(highrisk_HUTCH_cis_eqtl_all,anno=hutch_ge_anno)
highrisk_HUTCH_cis_eqtl_all=unique_eqtlpairs(dat1=highrisk_HUTCH_cis_eqtl_all,printflag = F)
# [1] "number of pairs: 3201"
# [1] "number of unique SNPs: 147"
# [1] "number of unique genes: 2755"

#TCGA_tumors eqtl
highrisk_TCGA_tumors_cis_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk_peer_cis",
                                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_POS.txt",fdrcutoff = 1)
highrisk_TCGA_tumors_cis_eqtl_all=addgenename_GE(highrisk_TCGA_tumors_cis_eqtl_all,anno=tcga_ge_anno)
highrisk_TCGA_tumors_cis_eqtl_all=unique_eqtlpairs(dat1=highrisk_TCGA_tumors_cis_eqtl_all,printflag = F)
# [1] "number of pairs: 3665"
# [1] "number of unique SNPs: 147"
# [1] "number of unique genes: 3133"
validatepairs_HUTCH_ciseqtl_fwer=validatepairs_difplatform(opt="fwer")
# [1] "number of qtl pairs to validate:55"
# [1] "number of qtl pairs can be validated:51"
# [1] "number of pairs validated:24"
# [1] "number of phenotype (probe) to validate:52"
# [1] "number of phenotype (probe) can be validated:49"
# [1] "number of phenotype (probe) validated:22"
# [1] "number of phenotype (gene) to validate:52"
# [1] "number of phenotype (gene) can be validated:49"
# [1] "number of phenotype (gene) validated:22"
# [1] "number of snp to validate:40"
# [1] "number of snp can be validate within pairs:39"
# [1] "number of snp validated:22"
validatepairs_HUTCH_ciseqtl_fwer=update_highrisk_snp_idx(dat=validatepairs_HUTCH_ciseqtl_fwer)

validatepairs_HUTCH_ciseqtl=validatepairs_difplatform()
# [1] "number of qtl pairs to validate:55"
# [1] "number of qtl pairs can be validated:51"
# [1] "number of pairs validated:35"
# [1] "number of phenotype (probe) to validate:52"
# [1] "number of phenotype (probe) can be validated:49"
# [1] "number of phenotype (probe) validated:33"
# [1] "number of phenotype (gene) to validate:52"
# [1] "number of phenotype (gene) can be validated:49"
# [1] "number of phenotype (gene) validated:33"
# [1] "number of snp to validate:40"
# [1] "number of snp can be validate within pairs:39"
# [1] "number of snp validated:29"
validatepairs_HUTCH_ciseqtl=update_highrisk_snp_idx(dat=validatepairs_HUTCH_ciseqtl)

#mqtl---
highrisk_HUTCH_cis_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis",
                                   snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                   geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
# [1] "number of SNPs: 147"
# [1] "number of phenotypeprobes: 478998"
# [1] "number of pairs: 104610"
# [1] "number of unique SNPs: 147"
# [1] "number of unique phenotypes: 80644"
# [1] "FDR cutoff: 0.05"
# [1] "pvalue cutoff: 0.000567510638710477"
# [1] "number of pairs (after adjustment)1189"
# [1] "number of unique SNPs (after adjustment)116, 0.789115646258503"
# [1] "number of unique phenotypes (after adjustment)1136, 0.0023716174180268"
highrisk_HUTCH_cis_mqtl=addgenename_ME(highrisk_HUTCH_cis_mqtl)
checkoverlapsnps(qtlres1=highrisk_HUTCH_cis_eqtl,qtlres2=highrisk_HUTCH_cis_mqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")
# [1] "number of unique SNP1:40"
# [1] "number of unique SNP2:116"
# [1] "number of overlapped SNP:38"
# in2 notin2
# in1     38     78
# notin1   2     29
# [1] 0.002741112
highrisk_TCGA_tumors_cis_mqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/mqtl_highrisk_cis",
                                             snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                             geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_POS.txt",fdrcutoff = 1)
idx=which(10^-highrisk_TCGA_tumors_cis_mqtl_all$value<=0.05)
highrisk_TCGA_tumors_cis_mqtl_all=highrisk_TCGA_tumors_cis_mqtl_all[idx,]
highrisk_TCGA_tumors_cis_mqtl_all=addgenename_ME(highrisk_TCGA_tumors_cis_mqtl_all)
tcga_me_pos=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_POS.txt",header=T,stringsAsFactors = F) #the order is different from normal
validatepairs_HUTCH_cismqtl_fwer=validatepairs(qtlres1=highrisk_HUTCH_cis_mqtl,qtlres2=highrisk_TCGA_tumors_cis_mqtl_all,phenotype2=tcga_me_pos,opt="fwer")
# [1] "number of qtl pairs to validate:1189"
# [1] "number of qtl pairs can be validated:731"
# [1] "number of pairs validated:390"
# [1] "number of phenotype to validate:1136"
# [1] "number of phenotype can be validated:710"
# [1] "number of phenotype validated:380"
# [1] "number of snp to validate:116"
# [1] "number of snp can be validated within pairs:101"
# [1] "number of snp validated:70"
validatepairs_HUTCH_cismqtl_fwer=update_highrisk_snp_idx(validatepairs_HUTCH_cismqtl_fwer)
#draw the SNP location barplot
tmp=refsnps$region
tmp=table(tmp)
tmp=tmp[order(tmp)]
dat=data.frame(matrix(0,nrow=length(tmp),ncol=2))
colnames(dat)=c("Validated cis-eQTL","Validated cis-mQTL")
rownames(dat)=names(tmp)
tmp=unique(validatepairs_HUTCH_ciseqtl_fwer$snp)
idx=which(refsnps$SNP %in% tmp)
tmp=refsnps$region[idx]
tmp=table(tmp)
tmp=tmp[order(tmp)]
for (i in 1:length(tmp))
{
  idx=which(rownames(dat)==names(tmp)[i])
  dat$`Validated cis-eQTL`[idx]=tmp[i]
}
tmp=unique(validatepairs_HUTCH_cismqtl_fwer$snp)
idx=which(refsnps$SNP %in% tmp)
tmp=refsnps$region[idx]
tmp=table(tmp)
tmp=tmp[order(tmp)]
for (i in 1:length(tmp))
{
  idx=which(rownames(dat)==names(tmp)[i])
  dat$`Validated cis-mQTL`[idx]=tmp[i]
}
postscript(file="../result/cis_qtl_validated_barplot.ps")
barplot(as.matrix(dat), beside = TRUE,
        col = c("lightblue", "mistyrose", "lightcyan",
                "lavender", "cornsilk","cyan4","green","cyan"),
        ylim = c(0, 40),cex.axis = 1.5,cex.names =1.5)
legend("topleft",rownames(dat),col = c("lightblue", "mistyrose", "lightcyan",
                                        "lavender", "cornsilk","cyan4","green","cyan")
       ,fill=c("lightblue", "mistyrose", "lightcyan",
               "lavender", "cornsilk","cyan4","green","cyan"),cex=1.5,bty = "n")
dev.off()




validatepairs_HUTCH_cismqtl=validatepairs(qtlres1=highrisk_HUTCH_cis_mqtl,qtlres2=highrisk_TCGA_tumors_cis_mqtl_all,phenotype2=tcga_me_pos)
# [1] "number of qtl pairs to validate:1189"
# [1] "number of qtl pairs can be validated:731"
# [1] "number of pairs validated:579"
# [1] "number of phenotype to validate:1136"
# [1] "number of phenotype can be validated:710"
# [1] "number of phenotype validated:561"
# [1] "number of snp to validate:116"
# [1] "number of snp can be validated within pairs:101"
# [1] "number of snp validated:86"
validatepairs_HUTCH_cismqtl=update_highrisk_snp_idx(validatepairs_HUTCH_cismqtl)

overlap_hutch_highrisk_cis_eqtl_hutch_cis_mqtl=overlappairs(dat1=highrisk_HUTCH_cis_eqtl,dat2=highrisk_HUTCH_cis_mqtl)

#add validation info to the SNP table---
refsnps=read.table("../data/RESUB_Supplementary_Table16_v9.txt",header=T,sep="\t",stringsAsFactors = F)
refsnps$Chr=gsub(23,"X",refsnps$Chr)
refsnps_anno=read.csv("../data/SNPannotation1.csv",stringsAsFactors = F)
refsnps_anno$original_147snpfile_annoted_rs.=gsub(" ","",refsnps_anno$original_147snpfile_annoted_rs.)
refsnps_anno$original_147snpfile_annoted_rs.=gsub("43694598-positionwrong","",refsnps_anno$original_147snpfile_annoted_rs.)
sum(refsnps$SNP == refsnps_anno$original_147snpfile_annoted_rs.)
refsnps$region=refsnps_anno$Func.refGene
refsnps$gene=refsnps_anno$Gene.refGene
refsnps=add_validatedgenes_snp(dat=validatepairs_HUTCH_ciseqtl_fwer,dat1=refsnps,colname="Validated_in_ciseQTL")
validatepairs_HUTCH_cismqtl_fwer$genename=validatepairs_HUTCH_cismqtl_fwer$genename.x
refsnps=add_validatedgenes_snp(dat=validatepairs_HUTCH_cismqtl_fwer,dat1=refsnps,colname="Validated_in_cismQTL")
write.table(refsnps,file="../result/snp_TCGAvalidated_cisQTL.txt",row.names = F,col.names = T,sep="\t")
validatepairs_HUTCH_cismqtl_fwer_CpG_table=table_mqtl()
write.table(validatepairs_HUTCH_cismqtl_fwer_CpG_table,file="../result/validatepairs_cismqtl_fwer_CpG_table.txt",row.names = F,col.names = T,sep="\t",quote=F)

#trans--------------
highrisk_HUTCH_trans_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_trans",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
highrisk_HUTCH_trans_eqtl=addgenename_GE(highrisk_HUTCH_trans_eqtl,anno=hutch_ge_anno)
highrisk_HUTCH_trans_eqtl=unique_eqtlpairs(dat1=highrisk_HUTCH_trans_eqtl)
# [1] "number of pairs: 16"
# [1] "number of unique SNPs: 5"
# [1] "number of unique genes: 16"
# [1] "snp: rs9625483,SEPT5|SCYL1"
# [1] "snp: rs138466039,KLK8|FMR1NB"
# [1] "snp: rs138213197,POTED|POTEC"
# [1] "snp: rs76551843,CHGB|LGI1|STXBP2|OR2T5|PTPRN|VIP|A1CF|SST|HMP19"

highrisk_HUTCH_trans_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_peer_trans",
                                         snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                         geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",fdrcutoff = 1)
highrisk_HUTCH_trans_eqtl_all=addgenename_GE(highrisk_HUTCH_trans_eqtl_all,anno=hutch_ge_anno)
idx=which(10^-highrisk_HUTCH_trans_eqtl_all$value<=0.05)
highrisk_HUTCH_trans_eqtl_all=highrisk_HUTCH_trans_eqtl_all[idx,]
highrisk_HUTCH_trans_eqtl_all=unique_eqtlpairs(dat1=highrisk_HUTCH_trans_eqtl_all,printflag = F)
# [1] "number of pairs: 181375"
# [1] "number of unique SNPs: 147"
# [1] "number of unique genes: 18069"
save(highrisk_HUTCH_trans_eqtl_all,file="../result/tmp_highrisk_HUTCH_trans_eqtl_all.RData")

highrisk_HUTCH_trans_mqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
# [1] "number of SNPs: 147"
# [1] "number of phenotypeprobes: 478998"
# Read 7035745 rows and 6 (of 6) columns from 0.675 GB file in 00:00:09
# [1] "number of pairs: 7035745"
# [1] "number of unique SNPs: 147"
# [1] "number of unique phenotypes: 478998"
# [1] "FDR cutoff: 0.05"
# [1] "pvalue cutoff: 1.80397193943627e-06"
# [1] "number of pairs (after adjustment)2537"
# [1] "number of unique SNPs (after adjustment)92, 0.625850340136054"
# [1] "number of unique phenotypes (after adjustment)2528, 0.00527768383166527"
highrisk_HUTCH_trans_mqtl=addgenename_ME(highrisk_HUTCH_trans_mqtl)

highrisk_HUTCH_trans_mqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
                                     snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                     geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",fdrcutoff = 1)

checkoverlapsnps(qtlres1=highrisk_HUTCH_trans_eqtl,qtlres2=highrisk_HUTCH_trans_mqtl,
                 snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                 snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt")

highrisk_TCGA_tumors_trans_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/eqtl_highrisk_peer_trans",
                                               snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                               geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_POS.txt",fdrcutoff = 1)
idx=which(10^-highrisk_TCGA_tumors_trans_eqtl_all$value<=0.05)
highrisk_TCGA_tumors_trans_eqtl_all=highrisk_TCGA_tumors_trans_eqtl_all[idx,]
highrisk_TCGA_tumors_trans_eqtl_all=addgenename_GE(highrisk_TCGA_tumors_trans_eqtl_all,anno=tcga_ge_anno)
highrisk_TCGA_tumors_trans_eqtl_all=unique_eqtlpairs(dat1=highrisk_TCGA_tumors_trans_eqtl_all,printflag = F)
# [1] "number of pairs: 148134"
# [1] "number of unique SNPs: 146"
# [1] "number of unique genes: 20194"
save(highrisk_TCGA_tumors_trans_eqtl_all,file="../result/tmp_highrisk_TCGA_tumors_trans_eqtl_all.RData")

highrisk_TCGA_tumors_trans_mqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors/mqtl_highrisk_trans",
                                               snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                               geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_POS.txt",fdrcutoff = 1)
idx=which(10^-highrisk_TCGA_tumors_trans_mqtl_all$value<=0.05)
highrisk_TCGA_tumors_trans_mqtl_all=highrisk_TCGA_tumors_trans_mqtl_all[idx,]

validatepairs_HUTCH_transeqtl=validatepairs_difplatform(qtlres1=highrisk_HUTCH_trans_eqtl,qtlres2=highrisk_TCGA_tumors_trans_eqtl_all)
# [1] "number of qtl pairs to validate:16"
# [1] "number of qtl pairs can be validated:15"
# [1] "number of pairs validated:0"
# [1] "number of phenotype (probe) to validate:16"
# [1] "number of phenotype (probe) can be validated:15"
# [1] "number of phenotype (probe) validated:0"
# [1] "number of phenotype (gene) to validate:16"
# [1] "number of phenotype (gene) can be validated:15"
# [1] "number of phenotype (gene) validated:0"
# [1] "number of snp to validate:5"
# [1] "number of snp can be validate within pairs:5"
# [1] "number of snp validated:0"
validatepairs_HUTCH_transeqtl_fwer=validatepairs_difplatform(qtlres1=highrisk_HUTCH_trans_eqtl,qtlres2=highrisk_TCGA_tumors_trans_eqtl_all,opt="fwer")

validatepairs_HUTCH_transmqtl=validatepairs(qtlres1=highrisk_HUTCH_trans_mqtl,qtlres2=highrisk_TCGA_tumors_trans_mqtl_all,phenotype2=tcga_me_pos)
validatepairs_HUTCH_transmqtl_fwer=validatepairs(qtlres1=highrisk_HUTCH_trans_mqtl,qtlres2=highrisk_TCGA_tumors_trans_mqtl_all,phenotype2=tcga_me_pos,opt="fwer")
# [1] "number of qtl pairs to validate:2537"
# [1] "number of qtl pairs can be validated:2188"
# [1] "number of pairs validated:18"
# [1] "number of phenotype to validate:2528"
# [1] "number of phenotype can be validated:2179"
# [1] "number of phenotype validated:90"
# [1] "number of snp to validate:92"
# [1] "number of snp can be validated within pairs:75"
# [1] "number of snp validated:9"
validatepairs_HUTCH_transmqtl_fwer=update_highrisk_snp_idx(validatepairs_HUTCH_transmqtl_fwer)
refsnps=add_validatedgenes_snp(dat=validatepairs_HUTCH_transmqtl_fwer,dat1=refsnps,colname="Validated_in_transmQTL")
write.table(refsnps,file="../result/snp_TCGAvalidated_cisQTL.txt",row.names = F,col.names = T,sep="\t")


#gtex
load("../data/GTEx/gtex_ge_anno.RData")
highrisk_gtex_cis_eqtl=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/gtex/eqtl_highrisk_cis",
                                  snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                  geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE_POS.txt")
highrisk_gtex_cis_eqtl=addgenename_GE(highrisk_gtex_cis_eqtl,anno=gtex_ge_anno)
highrisk_gtex_cis_eqtl=unique_eqtlpairs(dat1=highrisk_gtex_cis_eqtl)

highrisk_gtex_cis_eqtl_all=readqtlres(qtlresfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/gtex/eqtl_highrisk_cis",
                                      snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                                      geposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/gtex_GE_POS.txt",fdrcutoff = 1)
highrisk_gtex_cis_eqtl_all=addgenename_GE(highrisk_gtex_cis_eqtl_all,anno=gtex_ge_anno)
idx=which(10^-highrisk_gtex_cis_eqtl_all$value<=0.05)
highrisk_gtex_cis_eqtl_all=highrisk_gtex_cis_eqtl_all[idx,]
validate_TCGA_highrisk_ciseqtl_pear_in_gtex=merge(validatepairs_HUTCH_ciseqtl_fwer,highrisk_gtex_cis_eqtl_all,by=c("snp_idx","genename"))
idx=which(10^-validate_TCGA_highrisk_ciseqtl_pear_in_gtex$value<=0.05/nrow(validate_TCGA_highrisk_ciseqtl_pear_in_gtex)
          &validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta.x*validate_TCGA_highrisk_ciseqtl_pear_in_gtex$beta>0)
validate_TCGA_highrisk_ciseqtl_pear_in_gtex=validate_TCGA_highrisk_ciseqtl_pear_in_gtex[idx,]
validate_TCGA_highrisk_ciseqtl_pear_in_gtex=update_highrisk_snp_idx(validate_TCGA_highrisk_ciseqtl_pear_in_gtex)
refsnps=add_validatedgenes_snp(dat=validate_TCGA_highrisk_ciseqtl_pear_in_gtex,dat1=refsnps,colname="Validated_in_eQTL_GTEx")
write.table(refsnps,file="../result/snp_TCGAvalidated_cisQTL.txt",row.names = F,col.names = T,sep="\t")

TCGA_GE=read_qtl_input_matrix("../result/qtl_input/TCGA_tumors_GE.txt")
HUTCH_GE=read_qtl_input_matrix("../result/qtl_input/HUTCH_GE.txt")
gtex_GE=read_qtl_input_matrix("../result/qtl_input/gtex_GE.txt")
# idx1=which(rownames(TCGA_GE)=="FAM156A|29057") #SNP145, C/T
# idx2=which(rownames(HUTCH_GE)=="ILMN_1658160")
# idx3=which(rownames(gtex_GE)=="ENSG00000182646.12")

TCGA_GE_SNP=read_qtl_input_matrix("../result/qtl_input/TCGA_tumors_highrisk_SNP_GE_updatename.txt")
HUTCH_GE_SNP=read_qtl_input_matrix("../result/qtl_input/HUTCH_highrisk_SNP_GE.txt")
gtex_GE_SNP=read_qtl_input_matrix("../result/qtl_input/gtex_highrisk_SNP_GE_updatename.txt")
sum(colnames(TCGA_GE)==colnames(TCGA_GE_SNP))
sum(colnames(HUTCH_GE)==colnames(HUTCH_GE_SNP))
sum(grepl(colnames(gtex_GE_SNP),colnames(gtex_GE)))

# postscript(file="../result/diff_tumor_normal_eqtl.ps")
# par(mfrow=c(2,1))
# plot(c(TCGA_GE_SNP[145,],HUTCH_GE_SNP[145,]),c(TCGA_GE[idx1,],HUTCH_GE[idx2,]),
#      xlab="Genotype (dosage)",ylab="Gene expression",main="Tumor samples",cex.lab=1.5,cex.axis=1.5)
# plot(c(gtex_GE_SNP[145,]),c(gtex_GE[idx3,]), xlab="Genotype (dosage)",ylab="Gene expression",
#      main="Normal samples", cex.lab=1.5,cex.axis=1.5)
# dev.off()
#plot the most significant validated ciseqtl in TCGA, SNP144,C/T
#rs5945619 (chrX:51241672, intergenic) | NUDT11
idx1=which(rownames(TCGA_GE)=="NUDT11|55190")
par(mfrow=c(1,1))

postscript(file="../result/top_validated_ciseqtl_in_TCGA.ps")
par(mar=c(5.1,5.1,4.1,2.1))
plot(c(TCGA_GE_SNP[144,]),c(TCGA_GE[idx1,]), xlab="Genotype (dosage)",ylab="Gene expression",
           cex.lab=1.5,cex.axis=1.5)
dev.off()
#most significant validated cismqtl in TCGA, SNP73, G/GCGCA
#rs141536087 (chr10:854691,LARP4B) | cg26597838, chr10:835615, ENHANCER,DHS,North_Shelf
TCGA_ME=read_qtl_input_matrix("../result/qtl_input/TCGA_tumors_ME.txt")
TCGA_ME_SNP=read_qtl_input_matrix("../result/qtl_input/TCGA_tumors_highrisk_SNP_ME_updatename.txt")
idx1=which(rownames(TCGA_ME)=="cg26597838")
postscript(file="../result/top_validated_cismqtl_in_TCGA.ps")
par(mar=c(5.1,5.1,4.1,2.1))
plot(c(TCGA_ME_SNP[73,]),c(TCGA_ME[idx1,]), xlab="Genotype (dosage)",ylab="Methylation",
     cex.lab=1.5,cex.axis=1.5)
dev.off()

#check distribution of CpGs in validated mqtl
dat1=merge(validatepairs_HUTCH_cismqtl_fwer,anno,by.x="gene",by.y="IlmnID")
idx=duplicated(dat1$gene)
dat1=dat1[!idx,]
dat=data.frame(chr=dat1$chr.x,pos=dat1$opos2.x,value=dat1$value.y)
tmp=plotgenome(dat)
tmp1=diff(tmp$posall)
quantile(tmp1)
#          0%         25%         50%         75%        100% 
# 2.0       190.5      2340.0     64981.0 175977933.0 

#whether the  top mediation triplets , the CpG and the gene transcript are correlated in TCGA data
idx1=which(rownames(HUTCH_GE)=="ILMN_1697499")
idx1=which(hutch_ge_anno$Probe_Id=="ILMN_1697499")
genename=hutch_ge_anno$Symbol[idx1]
idx1=which(grepl(genename,rownames(TCGA_GE)))
idx2=which(rownames(TCGA_ME)=="cg15708909")
tcga_comsamples=intersect(colnames(TCGA_GE),colnames(TCGA_ME))
tmp1=match(tcga_comsamples,colnames(TCGA_GE))
tmp2=match(tcga_comsamples,colnames(TCGA_ME))
cor(unlist(TCGA_GE[idx1,tmp1]),unlist(TCGA_ME[idx2,tmp2]))
#[1] -0.5408656
idx1=which(rownames(HUTCH_GE)=="ILMN_1697499")
HUTCH_ME=read_qtl_input_matrix("../result/qtl_input/HUTCH_ME.txt")
idx2=which(rownames(HUTCH_ME)=="cg15708909")
hutch_comsamples=intersect(colnames(HUTCH_GE),colnames(HUTCH_ME))
tmp1=match(hutch_comsamples,colnames(HUTCH_GE))
tmp2=match(hutch_comsamples,colnames(HUTCH_ME))
cor(unlist(HUTCH_GE[idx1,tmp1]),unlist(HUTCH_ME[idx2,tmp2]))
# [1] -0.7100514

