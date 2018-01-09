#!/usr/bin/env Rscript
#filter SNP data
library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library(sas7bdat)
source("functions.R")
filterimput=function(chr,colscore=10,colmaf=9,scorecutoff=0.3,mafcutoff=0.05,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation",
                      prefix="icogs")
{
  library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
  #   #check the case of chr=3
  #   tmp=anno_icogs[anno_icogs$Chrom==3,]
  #   dim(tmp)
  #   #[1] 14454    17, total 14454 SNPs
  #   idx1=which(tmp$IlluminaName %in% infotable$rs_id[idxtype2])
  #   idx2=which(tmp$IlluminaName %in% infotable$rs_id)
  #   head(idx2[!idx2 %in% idx1])
  #   head(tmp$IlluminaName[! tmp$IlluminaName %in% infotable$rs_id])
  #   which(colnames(FHCRC_extra_demo.csv)=="rs2729054")
  #   sum(is.na(FHCRC_extra_demo.csv[1:100,109426]))
  #   missing=tmp$IlluminaName[! tmp$IlluminaName %in% infotable$rs_id]
  #   missingidx=which(colnames(FHCRC_extra_demo.csv) %in% missing)
  #   missingtable=FHCRC_extra_demo.csv[,c(1:2,missingidx)]
  #   View(anno_icogs[anno_icogs$IlluminaName %in% missing,])
  #   #it was found that 2% of genotyped SNPs were not shown in imputed table, check missingtable, may because MAF is small
  #   idx=which(anno_icogs$IlluminaName %in% missing)
  #   missingpos=data.frame(chr=anno_icogs$Chrom[idx],pos=anno_icogs$build37position[idx],rsid=anno_icogs$IlluminaName[idx])
  #   head(infotable$position[infotable$position %in% missingpos$pos])
  # #[1]  620311 1820234 1820234 7848731 8389201 8389201
  #   #the part only have imputed data, not the genotyped data. They were included in the missing table, but we got the imputed data
  #   
  #   dim(missingpos)
  #   #[1] 214   3, 214 seems missing when using probeid to map
  #   idx=!missingpos$pos %in% infotable$position
  #   sum(idx)
  #   #[1] 168 not in infotable, they are completely missing
  #   idx1=which(colnames(FHCRC_extra_demo.csv) %in% missingpos$rsid[idx])
  #   missingtable=FHCRC_extra_demo.csv[,c(1:2,idx1)]
  #   
  #   idx1=which(colnames(FHCRC_extra_demo.csv) %in% missingpos$rsid[!idx])
  #   length(idx1)
  #   #[1] 46, they only have imputed data.
  #   imputedonly=FHCRC_extra_demo.csv[,c(1:2,idx1)]
  #   idx2=which(infotable$position %in% missingpos$pos[!idx])
  #   #save the example file
  #   save(missingtable,imputedonly,FHCRC_extra_demo.csv,file="FHCRC_extra_demo_chr3_missing_in_imputed.RData")
  #   
    print(chr) 
    infotable=fread(infofiles[chr],fill = T)
    infotable=as.data.frame(infotable)
    imputetable=fread(imputfiles[chr])
    imputetable=cbind(chr=chr,imputetable)
    idxindel=grepl("INDEL",imputetable$V1)
    if (sum(idxindel)>0) 
    {
      imputetable=imputetable[!idxindel,]
      infotable=infotable[!idxindel,]
    }
    idxr2=which(infotable[,colscore]>scorecutoff)
    idxmaf=which(infotable[,colmaf]>=mafcutoff & infotable[,colmaf]<=(1-mafcutoff))
    idx=intersect(idxr2,idxmaf)
    print(paste0("number of SNPs: ",length(idx)))
    imputetable=imputetable[idx,]
    infotable=infotable[idx,]
    # #remove duplicates
    # idx=which(duplicated(infotable$position))
    # if (length(idx)>0)
    # {
    #   print(paste0("number of duplicates: ",length(idx)))
    #   idx2rm=NULL
    #   for (i in 1:length(idx))
    #   {
    #     if (i %% 1000==0) cat(i,"..")
    #     idx1=which(infotable$position==infotable$position[idx[i]])
    #     if ("imputation_r2" %in% colnames(infotable))
    #     {
    #       idx2=which.max(infotable$imputation_r2[idx1])
    #     }else
    #     {
    #       idx2=which.max(infotable$imputation_r2_euro[idx1])
    #     }
    #     idx3=c(1:length(idx1))[!1:length(idx1) %in% idx2]
    #     idx2rm=c(idx2rm,idx1[idx3])
    #   }
    # }
    # idx2rm=unique(idx2rm)
    # imputetable=imputetable[-idx2rm,]
    # infotable=infotable[-idx2rm,]
    # print(paste0("number of SNPs, after rm duplicates: ",nrow(infotable)))
    save(imputetable,infotable,file=paste0(outfolder,"/",prefix,"_chr",chr,".RData")) 
    #idxtype2=which(infotable$type==2) 0:only in reference; 2: in both reference and data; 3: only in data
    
    # quantile(infotable$exp_freq_a1[idxtype2],c(0,0.1,0.5,0.9,1))
    return(0)
}

#test=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation3/SNP6_info_chr22.txt",header=T,sep=" ")
folder="/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data"
pre="FHCRC_icogs_imputed_dosages_"
folder1="/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data/INFO"
pre1="icogs_practical_info_"
imputfiles=paste0(folder,"/",pre,"chr",1:23,".txt")
infofiles=paste0(folder1,"/",pre1,"chr",1:23,"_varid.txt")
for (chr in 1:23)
{
  cat(chr,"..")
  filterimput(chr,colscore=10,colmaf=9,scorecutoff=0.3,mafcutoff=0.05,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation",
                       prefix="icogs")
}

#the second batch data
folder="/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data/FHCRC_icogs_missing"
pre="FHCRC_missing_icogs_imputed_dosages_"
imputfiles=paste0(folder,"/",pre,"chr",1:23,".txt")
infofiles=paste0(folder1,"/",pre1,"chr",1:23,"_varid.txt")
for (chr in 1:23)
{
  cat(chr,"..")
  filterimput(chr,colscore=10,colmaf=9,scorecutoff=0.3,mafcutoff=0.05,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation",
              prefix="icogs1")
}

#oncoarray
folder="/fh/fast/stanford_j/Oncoarray_Imputed_2016Jan25/Data"
pre="FHCRC_oncoarray_imputed_dosages_"
folder1="/fh/fast/stanford_j/Oncoarray_Imputed_2016Jan25/Data/INFO"
pre1="onco_practical_info_"
imputfiles=paste0(folder,"/",pre,"chr",1:23,".txt")
infofiles=paste0(folder1,"/",pre1,"chr",1:23,"_varid.txt")
for (chr in 1:23)
{
  cat(chr,"..")
  filterimput(chr,colscore=10,colmaf=9,scorecutoff=0.3,mafcutoff=0.05,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation",
              prefix="onco")
}


#For TCGA normals
folder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation4"
pre="SNP6_info_"
pre1="SNP6_imputed_dosages_"
imputfiles=paste0(folder,"/",pre1,"chr",1:23,".txt")
file.exists(imputfiles)
infofiles=paste0(folder,"/",pre,"chr",1:23,".txt")
file.exists(infofiles)
for (chr in 1:23)
{
  cat(chr,"..")
  filterimput(chr,colscore=7,colmaf=6,scorecutoff=0.3,mafcutoff=0.05,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation",
              prefix="TCGA")
}

#combine icogs and oncoarray
#snp1:infotable_icogs,imputetable_icogs2
#colimpstart first column for dosage in imputation file
combine2imputedsnps=function(infotable_icogs,imputetable_icogs2,infotable_onco,imputetable_onco,colimpstart1=6,colimpstart2=6,outfolder,prefix)
{
  infotable_icogs$overlap=infotable_onco$overlap=0
  imp_ref_icogs=imputetable_icogs2[which(infotable_icogs$snp_id == "---"),]
  imp_imp_icogs=imputetable_icogs2[which(infotable_icogs$snp_id != "---"),]
  imp_ref_onco=imputetable_onco[which(infotable_onco$snp_id == "---"),]
  imp_imp_onco=imputetable_onco[which(infotable_onco$snp_id != "---"),]
  info_ref_icogs=infotable_icogs[which(infotable_icogs$snp_id == "---"),]
  info_imp_icogs=infotable_icogs[which(infotable_icogs$snp_id != "---"),]
  info_ref_onco=infotable_onco[which(infotable_onco$snp_id == "---"),]
  info_imp_onco=infotable_onco[which(infotable_onco$snp_id != "---"),]
  #sum(info_ref_icogs$rs_id %in% info_ref_onco$rs_id)
  
  #SNPS were divided into two blocks in each data (refSNPs and nonrefSNPs) in order to deal with SNPs with the same positions
  #for SNPs imputed from reference, ref1 and ref2
  ref_rs_id=intersect(info_ref_icogs$rs_id,info_ref_onco$rs_id)
  idx=match(ref_rs_id,info_ref_icogs$rs_id)
  imp_block1_icogs=imp_ref_icogs[idx,]
  info_ref_icogs$overlap[idx]=1 # find overlap
  info_block1_icogs=info_ref_icogs[idx,]
  idx=match(ref_rs_id,info_ref_onco$rs_id)
  imp_block1_onco=imp_ref_onco[idx,]
  info_ref_onco$overlap[idx]=1 # find overlap
  info_block1_onco=info_ref_onco[idx,]
  
  #for SNPs imputed not from reference, imp1 and imp2
  #sum(info_imp_icogs$rs_id %in% info_imp_onco$rs_id)
  sum(info_imp_icogs$position %in% info_imp_onco$position)
  length(unique(info_imp_icogs$position))==nrow(info_imp_icogs)
  length(unique(info_imp_onco$position))==nrow(info_imp_onco)
  imp_com_pos=intersect(info_imp_icogs$position,info_imp_onco$position)
  idx=match(imp_com_pos,info_imp_icogs$position)
  imp_block2_icogs=imp_imp_icogs[idx,]
  info_block2_icogs=info_imp_icogs[idx,]
  info_imp_icogs$overlap[idx]=2
  idx=match(imp_com_pos,info_imp_onco$position)
  imp_block2_onco=imp_imp_onco[idx,]
  info_block2_onco=info_imp_onco[idx,]
  info_imp_onco$overlap[idx]=2
  
  #ref1 and imp2
  idx1=which(info_ref_icogs$overlap==0)
  idx2=which(info_imp_onco$overlap==0)
  sum(info_ref_icogs$position[idx1] %in% info_imp_onco$position[idx2])
  ref_imp_com_pos=intersect(info_ref_icogs$position[idx1],info_imp_onco$position[idx2])
  idx=match(ref_imp_com_pos,info_ref_icogs$position[idx1])
  imp_block3_icogs=imp_ref_icogs[idx1[idx],]
  info_block3_icogs=info_ref_icogs[idx1[idx],]
  info_ref_icogs$overlap[idx1[idx]]=3
  idx=match(ref_imp_com_pos,info_imp_onco$position[idx2])
  imp_block3_onco=imp_imp_onco[idx2[idx],]
  info_block3_onco=info_imp_onco[idx2[idx],]
  info_imp_onco$overlap[idx2[idx]]=3
  
  #imp1 and ref2
  idx1=which(info_imp_icogs$overlap==0)
  idx2=which(info_ref_onco$overlap==0)
  sum(info_imp_icogs$position[idx1] %in% info_ref_onco$position[idx2])
  imp_ref_com_pos=intersect(info_imp_icogs$position[idx1],info_ref_onco$position[idx2])
  idx=match(imp_ref_com_pos,info_imp_icogs$position[idx1])
  imp_block4_icogs=imp_imp_icogs[idx1[idx],]
  info_block4_icogs=info_imp_icogs[idx1[idx],]
  info_imp_icogs$overlap[idx1[idx]]=4
  idx=match(imp_ref_com_pos,info_ref_onco$position[idx2])
  imp_block4_onco=imp_ref_onco[idx2[idx],]
  info_block4_onco=info_ref_onco[idx2[idx],]
  info_ref_onco$overlap[idx2[idx]]=4
  
  imputetable=cbind(rbind(imp_block1_icogs,imp_block2_icogs,imp_block3_icogs,imp_block4_icogs),
                    rbind(imp_block1_onco[,colimpstart2:ncol(imp_block1_onco)],
                          imp_block2_onco[,colimpstart2:ncol(imp_block2_onco)],
                          imp_block3_onco[,colimpstart2:ncol(imp_block3_onco)],
                          imp_block4_onco[,colimpstart2:ncol(imp_block4_onco)]))
  #colnames(imputetable)[colimpstart1:(colimpstart+length(samplenames)-1)]=samplenames #samplenames in order
  infotable=rbind(info_block1_icogs,info_block2_icogs,info_block3_icogs,info_block4_icogs)
  save(imputetable,infotable,info_ref_icogs,info_ref_onco,info_imp_icogs,info_imp_onco,file=paste0(outfolder,"/",prefix,"_chr",chr,".RData")) 
  return(nrow(imputetable))
}

combineicogsonco=function(infolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation",colimpstart1=6,colimpstart2=6,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/merged_imputation",prefix="HUTCH")
{
  #process sample studyno
  hutchclinical=read.sas7bdat("/fh/fast/stanford_j/Janet/alldata_2016dec20.sas7bdat")
  FHCRC_onco_map1=read.csv("/fh/fast/stanford_j/Oncoarray_2015Dec16/LinktoIDs.csv")
  ICOGID_studyno <- read.table("/fh/fast/stanford_j/COGS/iCOGID_Studyno.txt", header=TRUE,stringsAsFactors = F)
  imputeicogs_samples=read.table("/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data/FHCRC_icogs_sample_order.txt") #use icogsid
  idx=match(imputeicogs_samples[,1],ICOGID_studyno$iCOGS_ID)
  which(is.na(idx)) #not have studyno
  idx1_icogs=which(is.na(idx)) #remove the column in data
  imputeicogs_samples=imputeicogs_samples[-idx1_icogs,]
  idx=match(imputeicogs_samples,ICOGID_studyno$iCOGS_ID)
  studyno_imputeicogs=ICOGID_studyno$studyno[idx]
  
  imputeicogs1_samples=read.table("/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data/FHCRC_icogs_missing/FHCRC_missing_icogs_sample_order.txt")
  idx=match(imputeicogs1_samples[,1],ICOGID_studyno$iCOGS_ID)
  which(is.na(idx)) #not have studyno
  studyno_imputeicogs1=ICOGID_studyno$studyno[idx]
  idx2_icogs=which(!studyno_imputeicogs1 %in% hutchclinical$studyno[hutchclinical$race_==1]) #one African
  studyno_imputeicogs1=studyno_imputeicogs1[-idx2_icogs]
  
  imputeonco_samples=read.table("/fh/fast/stanford_j/Oncoarray_Imputed_2016Jan25/Data/FHCRC_sample_order.txt")
  idx=match(imputeonco_samples[,1],FHCRC_onco_map1$Onc_ID)
  which(is.na(idx)) #not have studyno
  studyno_onco=FHCRC_onco_map1$STUDYNO[idx]
  
  
  chrs=1:23
  numprobes=0
  for (chr in chrs)
  {
    print(chr)
    load(paste0(infolder,"/","icogs","_chr",chr,".RData"))
    imputetable=as.data.frame(imputetable)
    infotable=as.data.frame(infotable)
    imputetable_icogs=imputetable
    imputetable_icogs=imputetable_icogs[,-(idx1_icogs+5)]
    #colnames were missing in the RData
    colnames(imputetable_icogs)[colimpstart1:ncol(imputetable_icogs)]=studyno_imputeicogs
    infotable_icogs=infotable
    load(paste0(infolder,"/","icogs1","_chr",chr,".RData"))
    imputetable=as.data.frame(imputetable)
    infotable=as.data.frame(infotable)
    imputetable_icogs1=imputetable
    imputetable_icogs1=imputetable_icogs1[,-(5+idx2_icogs)]
    colnames(imputetable_icogs1)[colimpstart1:ncol(imputetable_icogs1)]=studyno_imputeicogs1
    infotable_icogs1=infotable
    if (sum(infotable_icogs$position==infotable_icogs1$position)<nrow(infotable_icogs)) warning("icogs files not consistent!")
    load(paste0(infolder,"/","onco","_chr",chr,".RData"))
    imputetable=as.data.frame(imputetable)
    infotable=as.data.frame(infotable)
    imputetable_onco=imputetable
    colnames(imputetable_onco)[colimpstart2:ncol(imputetable_onco)]=studyno_onco
    infotable_onco=infotable
    imputetable_icogs2=cbind(imputetable_icogs,imputetable_icogs1[,colimpstart1:ncol(imputetable_icogs1)])
    #samplenames=c(studyno_imputeicogs,studyno_imputeicogs1,studyno_onco)
    tmp=combine2imputedsnps(infotable_icogs,imputetable_icogs2,infotable_onco,imputetable_onco,colimpstart1=6,colimpstart2=6,outfolder,prefix)
    numprobes=numprobes+tmp  
  }
  print(paste0("number of probes: ",numprobes))
}

#SNP for normal samples
#for TCGA sample names, hutch normal sample names
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/allnormaldata.RData")
tcga_geneexp_samples=colnames(tcga_geneexp)
tcga_methylation_samples=colnames(tcga_methy)
tcga_genotype_samples=unique(c(colnames(tcga_geneexp),colnames(tcga_methy))) #after removing seminal vescile samples,it is different from colnames(genotypedata)
hutch_normal_samples=colnames(hutch_geneexp)
normal_samples=c(hutch_normal_samples,tcga_genotype_samples)
normal_geneexp_samples=c(hutch_normal_samples,tcga_geneexp_samples)
normal_methylation_samples=c(hutch_normal_samples,tcga_methylation_samples)
#load("/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGAnormals.RData") #the ge/me data were not normalized #it includes seminal vescile samples! the first 67 columns of genotypedata in dosage files
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/allnormaldata.RData")
combinenormalsnps=function(infolder1="/fh/fast/stanford_j/Xiaoyu/QTL/result/merged_imputation",
                           infolder2="/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation",
                           colimpstart=2,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/merged_imputation",prefix="NORMAL")
{
  chrs=1:23
  numprobes=0
  for (chr in chrs)
  {
    print(chr)
    load(paste0(infolder1,"/","HUTCH","_chr",chr,".RData"))
    imputetable=as.data.frame(imputetable)
    infotable=as.data.frame(infotable)
    #hutch normals
    idx=match(hutch_normal_samples,colnames(imputetable))
    imputetable_hutch=imputetable[,c(1:5,idx)]
    infotable_hutch=infotable
    load(paste0(infolder2,"/","TCGA","_chr",chr,".RData"))
    imputetable=as.data.frame(imputetable)
    colnames(imputetable)[2:(1+ncol(genotypedata))]=colnames(genotypedata)
    imputetable_tcga=imputetable
    infotable_tcga=infotable
    #samplenames=c(hutch_normal_samples,colnames(genotypedata))
    #colimpstart = 2 the data in second imp starts from Column2
    tmp=combine2imputedsnps(infotable_icogs=infotable_hutch,imputetable_icogs2=imputetable_hutch,infotable_onco=infotable_tcga,imputetable_onco=imputetable_tcga,colimpstart1 = 6,colimpstart2 = 2,outfolder,prefix)
    numprobes=numprobes+tmp  
  }
  print(paste0("number of probes: ",numprobes))
}
  

  
hutchclinical=read.sas7bdat("/fh/fast/stanford_j/Janet/alldata_2016dec20.sas7bdat")
hutch_geneexp_samples=hutchclinical$studyno[which(hutchclinical$race_==1 & (hutchclinical$in_iCOGS==1 | hutchclinical$in_OncoArray==1) & hutchclinical$GeneExpr==1)]
hutch_methylation_samples=hutchclinical$studyno[which(hutchclinical$race_==1 & (hutchclinical$in_iCOGS==1 | hutchclinical$in_OncoArray==1) & hutchclinical$Methy==1)]
sum(hutch_geneexp_samples %in% hutch_methylation_samples)
#[1] 337
#SNP_GE,SNP_ME and SNP_POS

genereate_SNPfiles=function(genexp_samples=hutch_geneexp_samples,methylation_samples=hutch_methylation_samples,rdfiles,
                            colrsid=5,colpos=6,outfoler="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input",prefix="HUTCH")
{
  # SNP_GE=NULL
  # SNP_ME=NULL
  # SNP_POS=NULL
  chrs=1:23
  # for (chr in chrs)
  # {
  #   print(chr)
  #   load(rdfiles[chr])
  #   idx=match(genexp_samples,colnames(imputetable))
  #   SNP_GE=rbind(SNP_GE,cbind.data.frame(id=infotable[,colrsid],imputetable[,idx]))
  #   idx=match(methylation_samples,colnames(imputetable))
  #   SNP_ME=rbind(SNP_ME,cbind.data.frame(id=infotable[,colrsid],imputetable[,idx]))
  #   SNP_POS=rbind(SNP_POS,data.frame(snp=infotable[,colrsid],chr=chr,pos=infotable[,colpos]))
  # }
  # write.table(SNP_GE,file=paste0(outfoler,"/",prefix,"_SNP_GE.txt"),col.names = T,row.names = F,sep="\t",quote=F)
  # write.table(SNP_ME,file=paste0(outfoler,"/",prefix,"_SNP_ME.txt"),col.names = T,row.names = F,sep="\t",quote=F)
  # write.table(SNP_POS,file=paste0(outfoler,"/",prefix,"_SNP_POS.txt"),col.names = T,row.names = F,sep="\t",quote=F)
  for (chr in chrs)
  {
    print(chr)
    load(rdfiles[chr])
    idx1=match(genexp_samples,colnames(imputetable))
    idx2=match(methylation_samples,colnames(imputetable))
    SNP_GE=cbind.data.frame(id=infotable[,colrsid],imputetable[,idx1])
    SNP_ME=cbind.data.frame(id=infotable[,colrsid],imputetable[,idx2])
    SNP_POS=data.frame(snp=infotable[,colrsid],chr=chr,pos=infotable[,colpos])
    if (chr==chrs[1])
    {
      fwrite(SNP_GE,file=paste0(outfoler,"/",prefix,"_SNP_GE.txt"),col.names = T,row.names = F,sep="\t")
      fwrite(NP_ME,file=paste0(outfoler,"/",prefix,"_SNP_ME.txt"),col.names = T,row.names = F,sep="\t")
      fwrite(SNP_POS,file=paste0(outfoler,"/",prefix,"_SNP_POS.txt"),col.names = T,row.names = F,sep="\t")
    }else
    {
      fwrite(SNP_GE,file=paste0(outfoler,"/",prefix,"_SNP_GE.txt"),col.names = F,row.names = F,sep="\t",append = T)
      fwrite(NP_ME,file=paste0(outfoler,"/",prefix,"_SNP_ME.txt"),col.names = F,row.names = F,sep="\t",append = T)
      fwrite(SNP_POS,file=paste0(outfoler,"/",prefix,"_SNP_POS.txt"),col.names = F,row.names = F,sep="\t",append = T)
    }
  }
  
}
#for HUTCH
rdfiles=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/result/merged_imputation/HUTCH_chr",1:23,".RData")
#SNP_PCA
# bigpca=function(dat=HUTCH_SNP_GE[,2:ncol(HUTCH_SNP_GE)],prefix=1)
# {
#   if (class(dat[1,1])=="character")
#   {
#     for (i in 1:ncol(dat))
#     {
#       dat[,i]=as.numeric(dat[,i])
#     }
#   }
#   
#   library(bigpca)
#   orig.dir <- getwd(); #setwd(tempdir()); # move to temporary dir
#   filebck=paste0("bigpca",prefix,".bck")
#   filedsc=paste0("bigpca",prefix,".dsc")
#   if(file.exists(filebck)) { unlink(c(filebck,filedsc)) }
#   bM <- filebacked.big.matrix(nrow(dat), ncol(dat),
#                               dimnames = list(rownames(dat), colnames(dat)),
#                               backingfile = filebck,  backingpath = getwd(), descriptorfile = filedsc)
#   for (i in 1:ncol(dat))
#   {
#     bM[1:nrow(dat),i] <-dat[,i]
#   }
#   #bM <- get.big.matrix(filedsc)
#   #bM[1:nrow(dat),1] <-tmp[,1]
#   res=big.PCA(bM,pcs.to.keep = ncol(dat),thin = T,return.loadings = T)
#   #res=big.PCA(bM,pcs.to.keep=15,thin = T,return.loadings = T)
#   varprop=res$Evalues/sum(res$Evalues)
#   # numpc=sum(varprop>0.01)
#   # print(paste0("number of PCs: ",numpc))
#   # print(paste0("proportion of variance explained: ",round(sum(varprop[1:numpc]),digits = 3)))
#   numpc=15
#   res1=data.frame(matrix(NA,nrow=numpc,ncol=1+ncol(dat)))
#   colnames(res1)=c("id",colnames(dat))
#   res1$id=paste0("pc",1:numpc)
#   idx=match(colnames(dat),rownames(res$PCs))
#   res1[,2:ncol(res1)]=t(res$PCs[idx,])[1:numpc,]
#   unlink(c(filebck,filedsc))
#   return(res1)
# }
HUTCH_SNP_GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",header=T,sep="\t",drop="id")
#HUTCH_SNP_GE=as.data.frame(HUTCH_SNP_GE)
#HUTCH_SNP_GE_PCA=bigpca(dat=t(HUTCH_SNP_GE))
HUTCH_SNP_GE_PCA=pca(dat=HUTCH_SNP_GE,check=T,bigpca = T,numpc=15)
write.table(HUTCH_SNP_GE_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)
rm(HUTCH_SNP_GE)
HUTCH_SNP_ME=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_ME.txt",header=T,sep="\t",drop="id")
#HUTCH_SNP_ME=as.data.frame(HUTCH_SNP_ME)
#HUTCH_SNP_ME_PCA=bigpca(dat=t(HUTCH_SNP_ME),prefix=2)
HUTCH_SNP_ME_PCA=pca(dat=HUTCH_SNP_ME,check=T,bigpca = T,numpc=15)
write.table(HUTCH_SNP_ME_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_ME_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)
rm(HUTCH_SNP_ME)


#for NORMAL
rdfiles=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/result/merged_imputation/NORMAL_chr",1:23,".RData")
genereate_SNPfiles(genexp_samples=normal_geneexp_samples,methylation_samples=normal_methylation_samples,rdfiles,
                            colrsid=5,colpos=6,outfoler="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input",prefix="NORMAL")
#SNP_PCA
NORMAL_SNP_GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_GE.txt",header=T,sep="\t",drop="id")
# NORMAL_SNP_GE_PCA=bigpca(dat=t(NORMAL_SNP_GE)])
NORMAL_SNP_GE_PCA=pca(dat=NORMAL_SNP_GE,check=T,bigpca = T,numpc=15)
write.table(NORMAL_SNP_GE_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_GE_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)
rm(NORMAL_SNP_GE)
NORMAL_SNP_ME=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_ME.txt",header=T,sep="\t",drop="id")
#NORMAL_SNP_ME=as.data.frame(NORMAL_SNP_ME)
#NORMAL_SNP_ME_PCA=bigpca(dat=t(NORMAL_SNP_ME))
NORMAL_SNP_ME_PCA=pca(dat=NORMAL_SNP_ME,check=T,bigpca = T,numpc=15)
write.table(NORMAL_SNP_ME_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_ME_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)
rm(NORMAL_SNP_ME)

tcgaclinical=read.table("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/clinical/dc99e120-0159-4470-ab3c-6581032935e9/nationwidechildrens.org_clinical_patient_prad.txt",
                        skip=1,stringsAsFactors = F,header=T,sep="\t")
tcgaclinical=tcgaclinical[-1,]
clinical=data.frame(matrix(NA,nrow=nrow(hutchclinical)+nrow(tcgaclinical),ncol=3))
colnames(clinical)=c("sampleid","gleason","age")
clinical$sampleid=c(hutchclinical$studyno,tcgaclinical$bcr_patient_barcode)
clinical$gleason=c(hutchclinical$Gleason_score,tcgaclinical$gleason_score)
clinical$age=c(hutchclinical$AGEREF,-as.numeric(tcgaclinical$days_to_birth)/365)
#save(tcgaclinical,hutchclinical,clinical,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/clinical.RData")
# load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/TCGAnormals.RData")
# load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/normtumor_GEME.RData")
# load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData")

#GE/ME Hutch
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/imputation_GEME.RData") #ME and GE were processed before
GE1=read.table("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/GE1.txt",header=T,sep="\t",stringsAsFactors = F)
colnames(GE1)=gsub("^X","",colnames(GE1))
idx=match(hutch_geneexp_samples,colnames(GE1))
HUTCH_GE=GE1[,c(1,idx)]
write.table(HUTCH_GE,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",row.names = F,sep="\t",quote=F)
#Jan5 update---standardize GE
tmp=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",header=T,sep="\t",stringsAsFactors = F)
rownames(tmp)=tmp$id
colnames(tmp)=gsub("^X","",colnames(tmp))
tmp=tmp[,-1]
tmp1=t(scale(t(tmp)))
tmp1=cbind.data.frame(id=rownames(tmp1),tmp1)
write.table(tmp1,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",row.names = F,sep="\t",quote=F)
#------

HUTCH_GE_norm=normalizeGE(dat=HUTCH_GE[,2:ncol(HUTCH_GE)])
HUTCH_GE_norm=cbind.data.frame(id=HUTCH_GE$id,HUTCH_GE_norm)
write.table(HUTCH_GE_norm,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_norm.txt",row.names = F,sep="\t",quote=F)
tmp=HUTCH_GE
rownames(tmp)=tmp$id
HUTCH_GE_PCA=pca(dat=tmp[,2:ncol(tmp)])#14,0.508
write.table(HUTCH_GE_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)
tmp=HUTCH_GE_norm
rownames(tmp)=tmp$id
HUTCH_GE_norm_PCA=pca(dat=tmp[,2:ncol(tmp)])#14,0.51
write.table(HUTCH_GE_norm_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_norm_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)
#peer factors
HUTCH_GE_PEER=peer(dat=HUTCH_GE[,2:ncol(HUTCH_GE)])
write.table(HUTCH_GE_PEER,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_PEER.txt",col.names = T,row.names = F,sep="\t",quote=F)
HUTCH_GE_norm_PEER=peer(dat=HUTCH_GE_norm[,2:ncol(HUTCH_GE_norm)])#15,0.517
write.table(HUTCH_GE_norm_PEER,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_norm_PEER.txt",col.names = T,row.names = F,sep="\t",quote=F)


idx=match(hutch_methylation_samples,colnames(ME))
HUTCH_ME=ME[,idx]
HUTCH_ME=cbind.data.frame(id=rownames(ME),HUTCH_ME)
write.table(HUTCH_ME,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",row.names = F,sep="\t",quote=F)
HUTCH_ME_PCA=pca(dat=HUTCH_ME[,2:ncol(HUTCH_ME)])#9,0.449
write.table(HUTCH_ME_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)

#GE/ME POS
system("cp /fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/GE_POS.txt /fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
system("cp /fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/ME_POS.txt /fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")


#GE/ME NORMAL
idx=match(normal_geneexp_samples,colnames(allnormal_geneexp_norm))
tmp=cbind.data.frame(id=allnormal_geneexp_pos$geneid,allnormal_geneexp_norm[,idx])
write.table(tmp,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt",col.names = T,row.names = F,sep="\t",quote=F)
write.table(allnormal_geneexp_pos,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",col.names = T,row.names = F,sep="\t",quote=F)
NORMAL_GE=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt",header=T,sep="\t")
NORMAL_GE_PCA=pca(dat=NORMAL_GE[,2:ncol(NORMAL_GE)])#26,0.831
write.table(NORMAL_GE_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_PCA.txt",row.names = F,sep="\t",quote=F)
NORMAL_GE_norm=normalizeGE(dat=NORMAL_GE[,2:ncol(NORMAL_GE)])
NORMAL_GE_norm=cbind.data.frame(id=NORMAL_GE$id,NORMAL_GE_norm)
write.table(NORMAL_GE_norm,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_norm.txt",row.names = F,sep="\t",quote=F)
NORMAL_GE_norm_PCA=pca(dat=NORMAL_GE_norm[,2:ncol(NORMAL_GE_norm)])#27,0.825
write.table(NORMAL_GE_norm_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_norm_PCA.txt",row.names = F,sep="\t",quote=F)
#peer factors
NORMAL_GE_PEER=peer(dat=NORMAL_GE[,2:ncol(NORMAL_GE)])#14,0.573
write.table(NORMAL_GE_PEER,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_PEER.txt",col.names = T,row.names = F,sep="\t",quote=F)
NORMAL_GE_norm_PEER=peer(dat=NORMAL_GE_norm[,2:ncol(NORMAL_GE_norm)])#14,0.577
write.table(NORMAL_GE_norm_PEER,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_norm_PEER.txt",col.names = T,row.names = F,sep="\t",quote=F)


#remove all NA rows of me
numNA=rep(0,nrow(allnormal_methy_norm))
for (i in 1:nrow(allnormal_methy_norm))
{
  if (i %% 10000==0) cat(i,'..')
  numNA[i]=sum(is.na(allnormal_methy_norm[i,1:46]))
}
idx=which(numNA<46)
allnormal_methy_norm1=allnormal_methy_norm[idx,]
sum(is.na(allnormal_methy_norm1))
#fill NA
library(impute)
tmp=impute.knn(as.matrix(allnormal_methy_norm1))
allnormal_methy_norm1=as.data.frame(tmp$data)
sum(is.na(allnormal_methy_norm1))
idx1=match(normal_methylation_samples,colnames(allnormal_methy_norm1))
tmp=cbind.data.frame(id=rownames(allnormal_methy_norm1),allnormal_methy_norm1[,idx1])
write.table(tmp,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME.txt",col.names = T,row.names = F,sep="\t",quote=F)
NORMAL_ME=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME.txt",header=T,sep="\t")
NORMAL_ME_PCA=pca(dat=NORMAL_ME[,2:ncol(NORMAL_ME)])#26,0.775
write.table(NORMAL_ME_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)

tmp1=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
tmp1=as.data.frame(tmp1)
sum(rownames(allnormal_methy_norm)==tmp1$geneid)
#[1] 478998 they are the same HUTCH_ME_POS is NORMAL_ME_POS
write.table(tmp1[idx,],file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",col.names = T,row.names = F,sep="\t",quote=F)



formcovariates=function(ge=normal_ge,me=normal_methy,clinical=clinical,sampleid_ge=NULL,sampleid_me=NULL,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input",prefix="FHnormals")
{
  if (!dir.exists(outfolder)) dir.create(outfolder)
  #remove the id column
  idx=match(sampleid_ge,colnames(ge))
  ge=ge[,idx]
  idx=match(sampleid_me,colnames(me))
  me=me[,idx]
  #Covariates
  #Covariate ge
  #the maximum number of covariates: number of samples -3
  idx=match(sampleid_ge,clinical$sampleid)
  npca=11
  nvar=2 #age and gleason
  COVA_ge=data.frame(matrix(NA,nrow=nvar+npca,ncol=ncol(ge)+1))
  COVA_ge[,1]=c("age","gleason",paste0("pc",1:(npca)))
  rownames(COVA_ge)=c("age","gleason",paste0("pc",1:(npca)))
  COVA_ge[1,2:ncol(COVA_ge)]=clinical$age[idx]
  COVA_ge[2,2:ncol(COVA_ge)]=clinical$gleason[idx]
  
  #PCA
  pcge=prcomp(t(ge),scale=T,center=T)
  colnames(COVA_ge)=c("id",colnames(ge))
  idx=match(colnames(ge),rownames(pcge$x))
  for (i in 1:npca)
  {
    COVA_ge[i+nvar,2:ncol(COVA_ge)]=pcge$x[idx,i]
  }
  write.table(COVA_ge,file=paste0(outfolder,"/",prefix,"_COVA_GE.txt"),row.names = F,col.names = T,sep="\t",quote=F)
  
  #Covariate methy
  idx=match(sampleid_me,clinical$sampleid)
  npca=13
  nvar=2 #age and gleason
  COVA_me=data.frame(matrix(NA,nrow=nvar+npca,ncol=ncol(me)+1))
  COVA_me[,1]=c("age","gleason",paste0("pc",1:(npca)))
  rownames(COVA_me)=c("age","gleason",paste0("pc",1:(npca)))
  COVA_me[1,2:ncol(COVA_me)]=clinical$age[idx]
  COVA_me[2,2:ncol(COVA_me)]=clinical$gleason[idx]
 
  #PCA
  pcme=prcomp(t(me),scale=T,center=T)
  colnames(COVA_me)=c("id",colnames(me))
  idx=match(colnames(me),rownames(pcme$x))
  for (i in 1:npca)
  {
    COVA_me[i+nvar,2:ncol(COVA_me)]=pcme$x[idx,i]
  }
  write.table(COVA_me,file=paste0(outfolder,"/",prefix,"_COVA_ME.txt"),row.names = F,col.names = T,sep="\t",quote=F)
  
  # ge=cbind.data.frame(id=rownames(ge),ge) #rownames of ge should be ID
  # write.table(ge,file=paste0(outfolder,"/",prefix,"_ge.txt"),row.names = F,col.names = T,sep="\t",quote=F)
  # methy=cbind.data.frame(id=rownames(methy),methy) #rownames of methy should be ID
  # write.table(methy,file=paste0(outfolder,"/",prefix,"_methy.txt"),row.names = F,col.names = T,sep="\t",quote=F)
  # colnames(snp)[1]="id"
  # snp_ge=snp_methy=snp
  # write.table(snp_ge,file=paste0(outfolder,"/",prefix,"_snp_ge.txt"),row.names = F,col.names = T,sep="\t",quote=F)
  # write.table(snp_methy,file=paste0(outfolder,"/",prefix,"_snp_methy.txt"),row.names = F,col.names = T,sep="\t",quote=F)
}
#for HUTCH
formcovariates(ge=HUTCH_GE,me=HUTCH_ME,clinical=clinical,sampleid_ge=hutch_geneexp_samples,sampleid_me=hutch_methylation_samples,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input",prefix="HUTCH")
#the file was updated in createdata.R

#for NORMAL
formcovariates(ge=allnormal_geneexp_norm,me=allnormal_methy_norm1,clinical=clinical,sampleid_ge=normal_geneexp_samples,sampleid_me=normal_methylation_samples,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input",prefix="NORMAL")


#HUTCH got SNPfiles of highrisks
extract_highrisk_hutch=function(folder1="/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data/INFO",pre1="icogs_practical_info_",
                           folder11="/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data",pre11="FHCRC_icogs_imputed_dosages_",
                           folder111="/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data/FHCRC_icogs_missing/",pre111="FHCRC_missing_icogs_imputed_dosages_",
                           folder2="/fh/fast/stanford_j/Oncoarray_Imputed_2016Jan25/Data/INFO",pre2="onco_practical_info_",
                           folder22="/fh/fast/stanford_j/Oncoarray_Imputed_2016Jan25/Data",pre22="FHCRC_oncoarray_imputed_dosages_")
{
  alldatamap=read.sas7bdat("/fh/fast/stanford_j/Janet/alldata_2016dec20.sas7bdat")
  FHCRC_onco_map1=read.csv("/fh/fast/stanford_j/Oncoarray_2015Dec16/LinktoIDs.csv")
  ICOGID_studyno <- read.table("/fh/fast/stanford_j/COGS/iCOGID_Studyno.txt", header=TRUE,stringsAsFactors = F)
  imputeicogs_samples=read.table("/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data/FHCRC_icogs_sample_order.txt")
  idx=match(imputeicogs_samples[,1],ICOGID_studyno$iCOGS_ID)
  which(is.na(idx)) #not have studyno
  idx1_icogs=which(is.na(idx)) #remove the column in data
  imputeicogs_samples=imputeicogs_samples[-idx1_icogs,]
  idx=match(imputeicogs_samples,ICOGID_studyno$iCOGS_ID)
  studyno_imputeicogs=ICOGID_studyno$studyno[idx]
  
  imputeicogs1_samples=read.table("/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data/FHCRC_icogs_missing/FHCRC_missing_icogs_sample_order.txt")
  idx=match(imputeicogs1_samples[,1],ICOGID_studyno$iCOGS_ID)
  which(is.na(idx)) #not have studyno
  studyno_imputeicogs1=ICOGID_studyno$studyno[idx]
  idx2_icogs=which(!studyno_imputeicogs1 %in% alldatamap$studyno[alldatamap$race_==1]) #one African
  studyno_imputeicogs1=studyno_imputeicogs1[-idx2_icogs]
  
  imputeonco_samples=read.table("/fh/fast/stanford_j/Oncoarray_Imputed_2016Jan25/Data/FHCRC_sample_order.txt")
  idx=match(imputeonco_samples[,1],FHCRC_onco_map1$Onc_ID)
  which(is.na(idx)) #not have studyno
  studyno_onco=FHCRC_onco_map1$STUDYNO[idx]
  
  
  anno_onco=fread("oncoarray_with_consortium_annotation.csv",stringsAsFactors = F)
  anno_onco=as.data.frame(anno_onco)
  newsnps=read.xls("RESUB_Supplementary_Table16_v9.xlsx",skip=1,sep=",")
  newsnps=newsnps[1:147,1:4]
  colnames(newsnps)=c("snp","chr","position","source")
  
  allsnps=cbind(newsnps,r2_icogs=NA,freq_a1_icogs=NA,r2_onco=NA,freq_a1_onco=NA)
  
  highrisk=NULL
  chrs=1:23
  
  for (chr in chrs)
  {
    idxs=which(allsnps$chr==chr)
    if (length(idxs)>0)
    {
      icogs=NULL
      icogs1=NULL
      onco=NULL
      print(chr)
      infotable_icogs=fread(paste0(folder1,"/",pre1,"chr",chr,"_varid.txt"),fill = T)
      infotable_onco=fread(paste0(folder2,"/",pre2,"chr",chr,"_varid.txt"),fill = T)
      imputetable=fread(paste0(folder11,"/",pre11,"chr",chr,".txt"))
      imputetable=cbind(chr=chr,imputetable)
      #imputetable=imputetable[,-413] #no studyno
      imputetable=imputetable[,-(idx1_icogs+5)] #no studyno
      for (i in idxs)
      {
        idx=which(infotable_icogs$position==allsnps$pos[i])
        idx1=which.max(infotable_icogs$imputation_r2[idx])
        if (length(idx1)>0)
        {
          icogs=rbind.data.frame(as.data.frame(icogs),as.data.frame(imputetable[idx[idx1],]))
        }else
        {
          tmp=data.frame(matrix(NA,nrow=1,ncol=ncol(imputetable)))
          colnames(tmp)=colnames(icogs)
          tmp$chr=chr
          tmp$V1=allsnps$snp[i]
          tmp$V2=allsnps$position[i]
          icogs=rbind.data.frame(as.data.frame(icogs),tmp)
        }
        
      }
      colnames(icogs)[6:ncol(icogs)]=studyno_imputeicogs
      imputetable=fread(paste0(folder111,"/",pre111,"chr",chr,".txt"))
      imputetable=imputetable[,-(idx2_icogs+4)] #no chr in column1
      #imputetable=imputetable[,-138] #no chr in column1
      for (i in idxs)
      {
        idx=which(infotable_icogs$position==allsnps$pos[i])
        idx1=which.max(infotable_icogs$imputation_r2[idx])
        if (length(idx1)>0)
        {
          icogs1=rbind.data.frame(as.data.frame(icogs1),as.data.frame(imputetable[idx[idx1],]))
        }else
        {
          print(paste0("icogs allsnps idx:",i))
          tmp=data.frame(matrix(NA,nrow=1,ncol=ncol(imputetable)))
          colnames(tmp)=colnames(icogs1)
          icogs1=rbind.data.frame(as.data.frame(icogs1),tmp)
        }
        
      }
      colnames(icogs1)[5:ncol(icogs1)]=studyno_imputeicogs1
      imputetable=fread(paste0(folder22,"/",pre22,"chr",chr,".txt"))
      #process onco
      for (i in idxs)
      {
        idx=which(infotable_onco$position==allsnps$pos[i])
        idx1=which.max(infotable_onco$imputation_r2_euro[idx])
        if (length(idx1)>0)
        {
          onco=rbind.data.frame(as.data.frame(onco),as.data.frame(imputetable[idx[idx1],]))
        }else
        {
          print(paste0("onco allsnps idx:",i))
          tmp=data.frame(matrix(NA,nrow=1,ncol=ncol(imputetable)))
          colnames(tmp)=colnames(onco)
          onco=rbind.data.frame(as.data.frame(onco),tmp)
        }
        
      }
      colnames(onco)[5:ncol(onco)]=studyno_onco
      highrisk=rbind.data.frame(highrisk,cbind.data.frame(icogs,icogs1[,5:ncol(icogs1)],onco[,5:ncol(onco)]))
      
      for (i in idxs)
      {
        idx=which(infotable_icogs$position==allsnps$pos[i])
        idx1=which.max(infotable_icogs$imputation_r2[idx])
        if (length(idx1)>0)
        {
          allsnps$r2_icogs[i]=infotable_icogs$imputation_r2[idx[idx1]]
          allsnps$freq_a1_icogs[i]=infotable_icogs$exp_freq_a1[idx[idx1]]
        }
        
        #process onco
        idx=which(infotable_onco$position==allsnps$pos[i])
        idx1=which.max(infotable_onco$imputation_r2_euro[idx])
        if (length(idx1)>0)
        {
          allsnps$r2_onco[i]=infotable_onco$imputation_r2_euro[idx[idx1]]
          allsnps$freq_a1_onco[i]=infotable_onco$exp_freq_a1[idx[idx1]]
        }
      }
    }
  }
  highrisk$chr=as.integer(highrisk$chr)
  highrisk$V2=as.integer(highrisk$V2)
  highrisk=highrisk[order(highrisk$chr,highrisk$V2),]
  allsnps$chr=gsub(23,"X",allsnps$chr)
  allsnps$chr=gsub(24,"Y",allsnps$chr)
  highrisk$chr=gsub(23,"X",highrisk$chr)
  highrisk$chr=gsub(24,"Y",highrisk$chr)
  save(allsnps,highrisk,file="allhighrisksnps_new.RData")
  return(allsnps)
}

load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData")
hutch_highrisk=highrisk
idx=match(hutch_geneexp_samples,colnames(hutch_highrisk))
hutch_highrisk_SNP_GE=cbind.data.frame(id=hutch_highrisk$V1,hutch_highrisk[,idx])
hutch_highrisk_SNP_GE$id=as.character(hutch_highrisk_SNP_GE$id)
idx=match(hutch_methylation_samples,colnames(hutch_highrisk))
hutch_highrisk_SNP_ME=cbind.data.frame(id=hutch_highrisk$V1,hutch_highrisk[,idx])
#hutch_highrisk_SNP_POS=data.frame(snp=hutch_highrisk$V1,chr=hutch_highrisk$chr,pos=hutch_highrisk$V2) #there is an NA
hutch_highrisk_SNP_POS=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/SNP_POSselsnps_new.txt",sep="\t",header=T,stringsAsFactors = F)
hutch_highrisk_SNP_POS$chr=factor(hutch_highrisk_SNP_POS$chr,levels = c(1:22,"X"))
hutch_highrisk_SNP_POS=hutch_highrisk_SNP_POS[order(hutch_highrisk_SNP_POS$chr,hutch_highrisk_SNP_POS$pos),]
sum(hutch_highrisk_SNP_POS$snp==hutch_highrisk_SNP_GE$id,na.rm=T)
write.table(hutch_highrisk_SNP_GE,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",row.names = F,col.names = T,sep="\t",quote=F)
write.table(hutch_highrisk_SNP_ME,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_ME.txt",row.names = F,col.names = T,sep="\t",quote=F)
write.table(hutch_highrisk_SNP_POS,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",row.names = F,col.names = T,sep="\t",quote=F)


#Normal SNPfiles of highrisks
extract_highrisk_normal=function(pre1="SNP6_info_",pre11="SNP6_imputed_dosages_")
{
  library(gdata)
  load("/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGAnormals.RData")
  load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData") #the order of highrisk and allsnps are not the same!
  impfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation4"
  newsnps=read.xls("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/RESUB_Supplementary_Table16_v9.xlsx",skip=1,sep=",")
  newsnps=newsnps[1:147,1:4]
  colnames(newsnps)=c("snp","chr","position","source")
  
  allsnps=cbind(newsnps,info=NA,exp_freq_a1=NA)
  
  highrisk=data.frame(matrix(NA,nrow=nrow(newsnps),ncol=5))
  colnames(highrisk)=c("chr","V1","V2","V3","V4")
  highrisk$chr=newsnps$chr
  highrisk$V1=newsnps$snp
  highrisk$V1=as.character(highrisk$V1)
  highrisk$V2=newsnps$position
  highrisk=cbind.data.frame(highrisk,data.frame(matrix(nrow=nrow(highrisk),ncol=ncol(genotypedata))))
  colnames(highrisk)[6:ncol(highrisk)]=colnames(genotypedata)
  chrs=1:23
  
  for (chr in chrs)
  {
    idxs=which(allsnps$chr==chr)
    if (length(idxs)>0)
    {
      #icogs=NULL
      #icogs1=NULL
      #onco=NULL
      print(chr)
      infotable=fread(paste0(impfolder,"/",pre1,"chr",chr,".txt"),fill = T)
      imputetable=fread(paste0(impfolder,"/",pre11,"chr",chr,".txt"))
      imputetable=imputetable[,1:ncol(genotypedata)]
      colnames(imputetable)=colnames(genotypedata) #from load
      
      for (i in idxs)
      {
        idx=which(infotable$position==allsnps$pos[i])
        if (length(idx)>0)
        {
          if (length(idx)>1) print(paste0(i," has multype snps"))
          idx1=which.max(infotable$info[idx])
          if (length(idx1)>0)
          {
            allsnps$info[i]=infotable$info[idx[idx1]]
            allsnps$exp_freq_a1[i]=infotable$exp_freq_a1[idx[idx1]]
            highrisk[i,6:ncol(highrisk)]=as.numeric(imputetable[idx[idx1],])
            highrisk$V1[i]=infotable$rs_id[idx[idx1]]
            highrisk$V3[i]=infotable$a0[idx[idx1]]
            highrisk$V4[i]=infotable$a1[idx[idx1]]
          }
        }
      }  
    }
  }
  highrisk$chr=as.integer(highrisk$chr)
  highrisk$V2=as.integer(highrisk$V2)
  highrisk=highrisk[order(highrisk$chr,highrisk$V2),]
  allsnps$chr=gsub(23,"X",allsnps$chr)
  allsnps$chr=gsub(24,"Y",allsnps$chr)
  highrisk$chr=gsub(23,"X",highrisk$chr)
  highrisk$chr=gsub(24,"Y",highrisk$chr)
  print(sum(is.na(allsnps$info)))
  print(mean(allsnps$info,na.rm=T)) #0.954
  print(median(allsnps$info,na.rm=T)) #0.984
  save(allsnps,highrisk,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_allhighrisksnps_new.RData")
  return(allsnps)
}
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_allhighrisksnps_new.RData") 
tcga_highrisk=highrisk
sum(hutch_highrisk$V2==tcga_highrisk$V2,na.rm = T)
#[1] 147 they have the same order
#the combined data
tmp=cbind.data.frame(hutch_highrisk,tcga_highrisk[,6:ncol(tcga_highrisk)])
idx=match(normal_geneexp_samples,colnames(tmp))
normal_highrisk_SNP_GE=cbind.data.frame(id=tmp$V1,tmp[,idx])
idx=match(normal_methylation_samples,colnames(tmp))
normal_highrisk_SNP_ME=cbind.data.frame(id=tmp$V1,tmp[,idx])
write.table(normal_highrisk_SNP_GE,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_GE.txt",row.names = F,col.names = T,sep="\t",quote=F)
write.table(normal_highrisk_SNP_ME,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_ME.txt",row.names = F,col.names = T,sep="\t",quote=F)


#For HUTCH mQTL
snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_ME.txt"
snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt"
phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt"
phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt"
covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME.txt"
output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_cis"
output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_trans"
cutoff_cis=1e-9
cutoff_trans=1e-11
recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl.RData"

#for HUTCH highrisk eQTL
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",
                snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
                phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",
                covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE.txt",
                output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis_old",
                output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans_old",
                cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_old_highrisk.RData")

do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_nogleason.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis_old_nogleason",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans_old_nogleason",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_old_nogleason_highrisk.RData")

do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_nogleasonage.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis_old_nogleasonage",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans_old_nogleasonage",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_old_nogleasonage_highrisk.RData")

#for HUTCH highrisk mQTL
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_ME.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis_old",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans_old",
       cutoff_cis=1,cutoff_trans=0.1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_old_highrisk.RData")

do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_ME.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME_nogleason.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis_old_nogleason",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans_old_nogleason",
       cutoff_cis=1,cutoff_trans=0.1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_old_nogleason_highrisk.RData")

do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_ME.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME_nogleasonage.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis_old_nogleasonage",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans_old_nogleasonage",
       cutoff_cis=1,cutoff_trans=0.1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_old_nogleasonage_highrisk.RData")

do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_ME.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME_nogleasonage_add3snppc.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis_old_nogleasonage_add3snppc",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans_old_nogleasonage_add3snppc",
       cutoff_cis=1,cutoff_trans=0.1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_old_nogleasonage_add3snppc_highrisk.RData")


#for NORMAL highrisk eQTL
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_GE.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_cis_old",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_trans_old",
       cutoff_cis=1,cutoff_trans=1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_old_highrisk.RData")

#for NORMAL highrisk mQTL
do_qtl(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_ME.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_ME.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_cis_old",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_trans_old",
       cutoff_cis=1,cutoff_trans=0.1,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_old_highrisk.RData")

#For NORMAL eQTL
basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL"
snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_GE.txt"
snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt"
phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt"
phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt"
covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_GE.txt"
output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_cis"
output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_trans"
cutoff_cis=1e-6
cutoff_trans=1e-8
recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl.RData"

#For NORMAL mQTL
basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL"
snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_ME.txt"
snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt"
phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME.txt"
phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt"
covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_ME.txt"
output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_cis"
output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_trans"
cutoff_cis=1e-9
cutoff_trans=1e-11
recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl.RData"

#test for a subset of SNPS
# system(paste0("head -1000 /fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_GE.txt > /fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/test_SNP_GE.txt"))
# system(paste0("head -1000 /fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt > /fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/test_SNP_POS.txt"))
do_qtl(basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result",snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/test_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/test_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_GE.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/test_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/test_trans",
       cutoff_cis=5e-2,cutoff_trans=1e-3,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/eqtl_test.RData")


#Parallel do_qtl for mQTLs
#salloc -t 2-1 -n 46 mpirun -n 1 R --interactive
#not read big table methylationpos when using Rmpi
do_qtl_loadgenotypepos=function(basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH",snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",
                snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
                genepos,
                covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE.txt",
                output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_cis",
                output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_trans",
                cutoff_cis=1e-6,cutoff_trans=1e-8,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl.RData")
{
  library("MatrixEQTL")
  base.dir = basedir
  
  ## Settings
  
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Genotype file name
  SNP_file_name = snpfile
  snps_location_file_name = snpposfile
  
  # Gene expression file name
  expression_file_name = phenotypefile
  #gene_location_file_name = phenotypeposfile
  
  # Covariates file name
  # Set to character() for no covariates
  covariates_file_name = covariatefile
  
  # Output file name
  output_file_name_cis = output_cis
  output_file_name_tra = output_trans
  
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = cutoff_cis
  pvOutputThreshold_tra = cutoff_trans
  
  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance = numeric();
  # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
  
  # Distance for local gene-SNP pairs
  cisDist = 1e6;
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 10000;      # read file in slices of 10,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  ## Run the analysis, not read big table when using Rmpi
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  #genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
  qtl = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  
  #plot(qtl)
  save(qtl,file=recordfile)
  
}  

library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(do_qtl_loadgenotypepos)
#for HUTCH mQTL
mpi_do_qtl(n=500, basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI",
           allsnpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_ME.txt",
           allsnpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
           phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",
           phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",
           covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME.txt",
           alloutput_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI/mqtl_cis",
           alloutput_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI/mqtl_trans",
           cutoff_cis=1e-9,cutoff_trans=1e-11,
           allrecordfiles="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH_MPI/mqtl.RData")
{
  genepos = read.table(phenotypeposfile, header = TRUE, stringsAsFactors = FALSE)
  mpi.bcast.Robj2slave(genepos)
  mpi.bcast.Robj2slave(basedir)
  # mpi.bcast.Robj2slave(allsnpfile)
  # mpi.bcast.Robj2slave(allsnpposfile)
  mpi.bcast.Robj2slave(phenotypefile)
  mpi.bcast.Robj2slave(phenotypeposfile)
  mpi.bcast.Robj2slave(covariatefile)
  mpi.bcast.Robj2slave(alloutput_cis)
  mpi.bcast.Robj2slave(alloutput_trans)
  mpi.bcast.Robj2slave(cutoff_cis)
  mpi.bcast.Robj2slave(cutoff_trans)
  #mpi.bcast.Robj2slave(allrecordfiles)
  #split SNP files
  #if (!file.exists(paste0(basedir,"/",basename(allsnpfile),".",n)))
  #{
    library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
    allsnp=fread(allsnpfile,header = T)
    allsnp=as.data.frame(allsnp)
    allsnppos=fread(allsnpposfile,header = T)
    allsnppos=as.data.frame(allsnppos)
    idx=seq(1,nrow(allsnp),by=ceiling(nrow(allsnp)/n))
    idx=c(idx,nrow(allsnp))
    for (i in 1:n)
    {
      if (i %%50==0) cat(i,"..")
      tmp=allsnp[idx[i]:(idx[i+1]-1),]
      write.table(tmp,file=paste0(basedir,"/",basename(allsnpfile),".",i),row.names = F,col.names = T,sep="\t",quote=F)
      tmp=allsnppos[idx[i]:(idx[i+1]-1),]
      write.table(tmp,file=paste0(basedir,"/",basename(allsnpposfile),".",i),row.names = F,col.names = T,sep="\t",quote=F)
    }
  #}
  snpfiles=paste0(basedir,"/",basename(allsnpfile),".",1:n)
  snpposfiles=paste0(basedir,"/",basename(allsnpposfile),".",1:n)
  output_cises=paste0(basedir,"/",basename(alloutput_cis),".",1:n)
  output_transes=paste0(basedir,"/",basename(alloutput_trans),".",1:n)
  recordfiles=paste0(basedir,"/",basename(allrecordfiles),".",1:n)
  mpi.bcast.Robj2slave(snpfiles)
  mpi.bcast.Robj2slave(snpposfiles)
  mpi.bcast.Robj2slave(output_cises)
  mpi.bcast.Robj2slave(output_transes)
  mpi.bcast.Robj2slave(recordfiles)
  do_qtl_job=function(i)
  {
    do_qtl_loadgenotypepos(basedir,snpfile=snpfiles[i],
           snpposfile=snpposfiles[i],
           phenotypefile,
           genepos,
           covariatefile,
           output_cis=output_cises[i],
           output_trans=output_transes[i],
           cutoff_cis,cutoff_trans,recordfile=recordfiles[i])
    return(0)
  }
  mpi.bcast.Robj2slave(do_qtl_job)
  res=mpi.parSapply(X=1:n,FUN=do_qtl_job,job.num=njobs)
  cis=NULL
  trans=NULL
  num_cis=rep(0,n)
  num_trans=rep(0,n)
  for (i in 1:n)
  {
    if(i %% 50==0) cat(i,'..')
    tmp=read.table(output_cises[i],header=T,sep="\t",stringsAsFactors = F,fill=T)
    tmp=tmp[!grepl("No",tmp$SNP,ignore.case=T),]
    if (nrow(tmp)>0) cis=rbind.data.frame(cis,tmp)
    tmp=read.table(output_transes[i],header=T,sep="\t",stringsAsFactors = F)
    if (nrow(tmp)>0) trans=rbind.data.frame(trans,tmp)
    load(recordfiles[i])
    num_cis[i]=qtl$cis$ntests
    num_trans[i]=qtl$trans$ntests
  }
  cis=cis[order(cis$p.value),]
  trans=trans[order(trans$p.value),]
  save(cis,trans,num_cis,num_trans,file=allrecordfiles)
  write.table(cis,file=alloutput_cis,sep="\t",row.names = F,col.names = T,quote=F)
  write.table(trans,file=alloutput_trans,sep="\t",row.names = F,col.names = T,quote=F)
}

mpi_do_qtl(n=500, basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL_MPI",
           allsnpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_ME.txt",
           allsnpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt",
           phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME.txt",
           phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",
           covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_ME.txt",
           alloutput_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL_MPI/mqtl_cis",
           alloutput_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL_MPI/mqtl_trans",
           cutoff_cis=1e-9,cutoff_trans=1e-11,
           allrecordfiles="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL_MPI/mqtl.RData")
