#!/usr/bin/env Rscript
#filter SNP data
library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library(sas7bdat)
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



#For TCGA tumorss
folder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation_tumor"
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
              prefix="TCGA_tumors")
}


#SNP_GE,SNP_ME and SNP_POS

genereate_SNPfiles1=function(genexp_samples=hutch_geneexp_samples,methylation_samples=hutch_methylation_samples,rdfiles,
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
  #   colnames(imputetable)[2:ncol(imputetable)]=colnames(genotypedata) #add this
  #   imputetable=as.data.frame(imputetable)
  #   idx=match(genexp_samples,colnames(imputetable))
  #   tmp=cbind.data.frame(id=infotable[,colrsid],imputetable[,idx])
  #   SNP_GE$id=as.character(SNP_GE$id)
  #   tmp$id=as.character(tmp$id)
  #   SNP_GE=rbind(SNP_GE,tmp)
  #   idx=match(methylation_samples,colnames(imputetable))
  #   tmp=cbind.data.frame(id=infotable[,colrsid],imputetable[,idx])
  #   SNP_ME$id=as.character(SNP_ME$id)
  #   tmp$id=as.character(tmp$id)
  #   SNP_ME=rbind(SNP_ME,tmp)
  #   SNP_POS=rbind(SNP_POS,data.frame(snp=infotable[,colrsid],chr=chr,pos=infotable[,colpos]))
  # }
  # write.table(SNP_GE,file=paste0(outfoler,"/",prefix,"_SNP_GE.txt"),col.names = T,row.names = F,sep="\t",quote=F)
  # write.table(SNP_ME,file=paste0(outfoler,"/",prefix,"_SNP_ME.txt"),col.names = T,row.names = F,sep="\t",quote=F)
  # write.table(SNP_POS,file=paste0(outfoler,"/",prefix,"_SNP_POS.txt"),col.names = T,row.names = F,sep="\t",quote=F)
  for (chr in chrs)
  {
    print(chr)
    load(rdfiles[chr])
    imputetable=as.data.frame(imputetable)
    colnames(imputetable)[2:ncol(imputetable)]=colnames(genotypedata) #add this
    idx1=match(genexp_samples,colnames(imputetable))
    idx2=match(methylation_samples,colnames(imputetable))
    SNP_GE=cbind.data.frame(id=infotable[,colrsid],imputetable[,idx1])
    SNP_ME=cbind.data.frame(id=infotable[,colrsid],imputetable[,idx2])
    SNP_POS=data.frame(snp=infotable[,colrsid],chr=chr,pos=infotable[,colpos])
    if (chr==chrs[1])
    {
      fwrite(SNP_GE,file=paste0(outfoler,"/",prefix,"_SNP_GE.txt"),col.names = T,row.names = F,sep="\t")
      fwrite(SNP_ME,file=paste0(outfoler,"/",prefix,"_SNP_ME.txt"),col.names = T,row.names = F,sep="\t")
      fwrite(SNP_POS,file=paste0(outfoler,"/",prefix,"_SNP_POS.txt"),col.names = T,row.names = F,sep="\t")
    }else
    {
      fwrite(SNP_GE,file=paste0(outfoler,"/",prefix,"_SNP_GE.txt"),col.names = F,row.names = F,sep="\t",append = T)
      fwrite(SNP_ME,file=paste0(outfoler,"/",prefix,"_SNP_ME.txt"),col.names = F,row.names = F,append = T,sep="\t")
      fwrite(SNP_POS,file=paste0(outfoler,"/",prefix,"_SNP_POS.txt"),col.names = F,row.names = F,sep="\t",append = T)
    }
  }
}

#GE
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGAtumors.RData")
tcga_tumors_geneexp_samples=intersect(colnames(genotypedata),colnames(geneexpdata))
tcga_tumors_methy_samples=intersect(colnames(genotypedata),colnames(methydata))
idx=match(tcga_tumors_geneexp_samples,colnames(geneexpdata))
GE=geneexpdata[,idx]
removeconstrows=function(dat=GE)
{
  idxconst=rep(F,nrow(dat))
  for (i in 1:nrow(dat))
  {
    if (i %%2000==0) cat(i,"...")
    if (var(unlist(dat[i,]))==0) idxconst[i]=T
  }
  dat=dat[!idxconst,]
  return(dat)
}
GE=removeconstrows(dat=GE)
tmp=cbind.data.frame(id=rownames(GE),GE)
write.table(tmp,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE.txt",col.names = T,row.names = F,sep="\t",quote=F)

#the distribution of RPKMs in each sample was quantile-transformed using the average empirical distribution observed across all samples. Expression measurements for each gene in each tissue were subsequently transformed to the quantiles of the standard normal distribution
normalizeGE=function(dat=GE)
{
  library(preprocessCore)
  avgGE=rowMeans(dat)
  GE_norm=data.frame(matrix(NA,nrow=nrow(dat),ncol=ncol(dat)))
  colnames(GE_norm)=colnames(dat)
  rownames(GE_norm)=rownames(dat)
  for (i in 1:ncol(dat))
  {
    GE_norm[,i]=normalize.quantiles.use.target(as.matrix(dat[,i]),target=avgGE)
  }
  GE_norm=t(scale(t(GE_norm)))
  return(GE_norm)
}
GE_norm=normalizeGE(dat=GE)
tmp=cbind.data.frame(id=rownames(GE_norm),GE_norm)
write.table(tmp,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_norm.txt",col.names = T,row.names = F,sep="\t",quote=F)
#GE position
rnaseqv2_anno=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/TCGA.hg19.June2011.gaf",header=T,sep="\t",comment.char = "",stringsAsFactors = F)
sum(rownames(GE) %in% rnaseqv2_anno$Gene)
GE_POS=data.frame(matrix(NA,nrow=nrow(GE),ncol=4))
colnames(GE_POS)=c("geneid","chr","s1","s2")
GE_POS$geneid=rownames(GE)
numdup=rep(0,nrow(GE))
for (i in 1:nrow(GE))
{
  if (i %% 2000 ==0) cat(i,"..")
  idx=which(rnaseqv2_anno$Gene==rownames(GE)[i])
  locus=rnaseqv2_anno$GeneLocus[idx][1]
  tmp=unlist(strsplit(locus,"chr"))
  tmp=tmp[tmp!=""]
  numdup[i]=length(tmp)
  tmp1=unlist(strsplit(tmp,":"))
  GE_POS$chr[i]=tmp1[1]
  startp=endp=rep(0,length(tmp))
  
  #pick the longer one
  for (j in 1:length(tmp))
  {
    tmp1=unlist(strsplit(tmp[j],":"))
    tmp2=as.integer(unlist(strsplit(tmp1[2],"-")))
    startp[j]=tmp2[1]
    endp[j]=tmp2[2]
  }
  idx1=which.max(endp-startp)
  GE_POS$s1[i]=startp[idx1]
  GE_POS$s2[i]=endp[idx1]
}
write.table(GE_POS,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_POS.txt",col.names = T,row.names = F,sep="\t",quote=F)

#ME
idx=match(tcga_tumors_methy_samples,colnames(methydata))
ME=methydata[,idx]
length(intersect(tcga_tumors_geneexp_samples,tcga_tumors_methy_samples))
#[1] 494

#remove all NA rows of me
numNA=rep(0,nrow(ME))
for (i in 1:nrow(ME))
{
  if (i %% 10000==0) cat(i,'..')
  numNA[i]=sum(is.na(ME[i,1:ncol(ME)]))
}
idx=which(numNA<ncol(ME))
ME1=ME[idx,]
sum(is.na(ME1))
#fill NA
library(impute)
tmp=impute.knn(as.matrix(ME1))
ME1=as.data.frame(tmp$data)
sum(is.na(ME1))
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
#some probes are not in anno
idx=rownames(ME1) %in% anno$IlmnID
ME1=ME1[idx,]
idx1=match(tcga_tumors_methy_samples,colnames(ME1))
tmp=cbind.data.frame(id=rownames(ME1),ME1[,idx1])
write.table(tmp,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME.txt",col.names = T,row.names = F,sep="\t",quote=F)

tmp1=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")
tmp1=as.data.frame(tmp1)
idx1=match(rownames(ME1),tmp1$geneid)
write.table(tmp1[idx1,],file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_POS.txt",col.names = T,row.names = F,sep="\t",quote=F)

#SNP_GE,SNP_ME,SNP_POS
rdfiles=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation/TCGA_tumors_chr",1:23,".RData")
genereate_SNPfiles1(genexp_samples=tcga_tumors_geneexp_samples,methylation_samples=tcga_tumors_methy_samples,rdfiles,
                            colrsid=2,colpos=3,outfoler="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input",prefix="TCGA_tumors")

GE_PCA=pca(dat=GE)#14,0.44
write.table(GE_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)
GE_norm_PCA=pca(dat=GE_norm)#13,0.438
write.table(GE_norm_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_norm_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)

#SNP_PCA
TCGA_tumors_SNP_GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_SNP_GE.txt",header=T,drop="id",sep="\t",stringsAsFactors = F)
#TCGA_tumors_SNP_GE=as.data.frame(TCGA_tumors_SNP_GE)
#TCGA_tumors_SNP_GE_PCA=bigpca(dat=TCGA_tumors_SNP_GE[,2:ncol(TCGA_tumors_SNP_GE)])
TCGA_tumors_SNP_GE_PCA=pca(dat=TCGA_tumors_SNP_GE,check=T,bigpca = T,numpc=15)
write.table(TCGA_tumors_SNP_GE_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_SNP_GE_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)
rm(TCGA_tumors_SNP_GE)
TCGA_tumors_SNP_ME=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_SNP_ME.txt",header=T,drop="id",sep="\t",stringsAsFactors = F)
#TCGA_tumors_SNP_ME=as.data.frame(TCGA_tumors_SNP_ME)
#TCGA_tumors_SNP_ME_PCA=bigpca(dat=TCGA_tumors_SNP_ME,prefix=6)
TCGA_tumors_SNP_ME_PCA=pca(dat=TCGA_tumors_SNP_ME,check=T,bigpca = T,numpc=15)
write.table(TCGA_tumors_SNP_ME_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_SNP_ME_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)
rm(TCGA_tumors_SNP_ME)

ME1=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME.txt",header = T,sep="\t")
ME1=as.data.frame(ME1)
rownames(ME1)=ME1$id
ME1=ME1[,-1]
ME_PCA=pca(dat=ME1) #12,0.452
write.table(ME_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)


peer=function(dat=GE)
{
  
  #scale the data first
  dat=t(scale(t(dat)))
  totvar=0
  for (i in 1:nrow(dat)) totvar=totvar+var(unlist(dat[i,]))
  
  library(peer)
  maxfactors=40
  varexplain=rep(0,maxfactors)
  for (i in 1:maxfactors)
  {
    print(i)
    model=PEER()
    PEER_setNk(model,i)
    PEER_setPhenoMean(model, as.matrix(t(dat)))
    PEER_setNmax_iterations(model, 1000)
    PEER_update(model)
    resi=PEER_getResiduals(model)
    resivar=0
    for (j in 1:ncol(resi)) resivar=resivar+var(resi[,j])
    print((totvar-resivar)/totvar)
    varexplain[i]=(totvar-resivar)/totvar
    if (i>1)
    {
      if (abs(varexplain[i]-varexplain[i-1])<0.01) #difference cutoff
      {
        break
      }
    }
  }
  numfactors=i
  print(paste0("number of factors: ",numfactors))
  print(paste0("variance explained: ",round(varexplain[numfactors],digits = 3)))
  #the factors
  x=PEER_getX(model)
  x=as.data.frame(t(x))
  colnames(x)=colnames(dat)
  res=cbind.data.frame(id=paste0("factor",1:numfactors),x)
  return(res)
}
GE_PEER=peer(dat=GE)#15,0.45
write.table(GE_PEER,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_PEER.txt",col.names = T,row.names = F,sep="\t",quote=F)
GE_norm_PEER=peer(dat=GE_norm)#15,0.455
write.table(GE_norm_PEER,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE_norm_PEER.txt",col.names = T,row.names = F,sep="\t",quote=F)


tcgaclinical=read.table("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/clinical/dc99e120-0159-4470-ab3c-6581032935e9/nationwidechildrens.org_clinical_patient_prad.txt",
                        skip=1,stringsAsFactors = F,header=T,sep="\t")
tcgaclinical=tcgaclinical[-1,]
for (i in 1:ncol(tcgaclinical))
{
  tcgaclinical[,i]=gsub("\\[Not Applicable\\]",NA,tcgaclinical[,i])
  tcgaclinical[,i]=gsub("\\[Not Available\\]",NA,tcgaclinical[,i])
  tcgaclinical[,i]=gsub("\\[Unknown\\]",NA,tcgaclinical[,i])
  
}




clinical=data.frame(matrix(NA,nrow=nrow(hutchclinical)+nrow(tcgaclinical),ncol=3))
colnames(clinical)=c("sampleid","gleason","age")
clinical$sampleid=c(hutchclinical$studyno,tcgaclinical$bcr_patient_barcode)
clinical$gleason=c(hutchclinical$Gleason_score,tcgaclinical$gleason_score)
clinical$age=c(hutchclinical$AGEREF,-as.numeric(tcgaclinical$days_to_birth)/365)
tcgaclinical1=read.table("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/tumors/clinical/aac3407e-e62c-4056-a414-9e8e88cbac74/nationwidechildrens.org_biospecimen_slide_prad.txt",
                         header = T,stringsAsFactors = F,sep="\t")
#save(tcgaclinical,hutchclinical,clinical,file="/fh/fast/stanford_j/Xiaoyu/QTL/data/clinical.RData")
# load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/TCGAnormals.RData")
# load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/normtumor_GEME.RData")
# load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData")

# #GE/ME Hutch
# load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/imputation_GEME.RData") #ME and GE were processed before
# GE1=read.table("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/GE1.txt",header=T,sep="\t",stringsAsFactors = F)
# colnames(GE1)=gsub("^X","",colnames(GE1))
# idx=match(hutch_geneexp_samples,colnames(GE1))
# HUTCH_GE=GE1[,c(1,idx)]
# write.table(HUTCH_GE,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",row.names = F,sep="\t",quote=F)
# idx=match(hutch_methylation_samples,colnames(ME))
# HUTCH_ME=ME[,idx]
# HUTCH_ME=cbind.data.frame(id=rownames(ME),HUTCH_ME)
# write.table(HUTCH_ME,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",row.names = F,sep="\t",quote=F)
# #GE/ME POS
# system("cp /fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/GE_POS.txt /fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt")
# system("cp /fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/imputation/ME_POS.txt /fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt")



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


#Tumor SNPfiles of highrisks
extract_highrisk_tumor=function(pre1="SNP6_info_",pre11="SNP6_imputed_dosages_")
{
  library(gdata)
  load("/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGAtumors.RData")
  load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData") #the order of highrisk and allsnps are not the same!
  impfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation_tumor"
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
  print(mean(allsnps$info,na.rm=T)) #0.952
  print(median(allsnps$info,na.rm=T)) #0.979
  save(allsnps,highrisk,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors_allhighrisksnps_new.RData")
  return(allsnps)
}
quantile(allsnps$info,na.rm=T)
#     0%     25%     50%     75%    100% 
#0.57500 0.94725 0.97900 0.99500 1.00000 
load("/fh/fast/stanford_j/Xiaoyu/QTL/result/TCGA_tumors_allhighrisksnps_new.RData") 
TCGA_tumors_allsnps=allsnps
TCGA_tumors_highrisk=highrisk
sum(is.na(TCGA_tumors_allsnps$info))
#[1] 1
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData")
sum(is.na(allsnps$r2_icogs))
par(mfrow=c(1,1))
png("Highrisk_score1.png")
boxplot(TCGA_tumors_allsnps$info[TCGA_tumors_allsnps$info>0],allsnps$r2_icogs,names=c("SNP6","iCOGS"),cex.names=1.4,cex.axis=1.4,ylab="Imputation score",cex.lab=1.4)
dev.off()
TCGA_tumors_allsnps[is.na(TCGA_tumors_allsnps$info),]
#            snp chr position source info exp_freq_a1
#130 rs138213197  17 46805705  Known   NA          NA

idx=match(tcga_tumors_geneexp_samples,colnames(highrisk))
TCGA_tumors_highrisk_SNP_GE=cbind.data.frame(id=highrisk$V1,highrisk[,idx])
idx=match(tcga_tumors_methy_samples,colnames(highrisk))
TCGA_tumors_highrisk_SNP_ME=cbind.data.frame(id=highrisk$V1,highrisk[,idx])
write.table(TCGA_tumors_highrisk_SNP_GE,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_GE.txt",row.names = F,col.names = T,sep="\t",quote=F)
write.table(TCGA_tumors_highrisk_SNP_ME,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_ME.txt",row.names = F,col.names = T,sep="\t",quote=F)


do_qtl=function(basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH",snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",
                snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
                phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",
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
  gene_location_file_name = phenotypeposfile
  
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
  
  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
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
do_qtl(basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH",snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",
                snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
                phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
                phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",
                covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE.txt",
                output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_cis",
                output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk_trans",
                cutoff_cis=5e-2,cutoff_trans=1e-3,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_highrisk.RData")

#for HUTCH highrisk mQTL
do_qtl(basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH",snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_ME.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk_trans",
       cutoff_cis=5e-2,cutoff_trans=1e-3,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/mqtl_highrisk.RData")

#for NORMAL highrisk eQTL
do_qtl(basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL",snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_GE.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_GE_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_GE.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk_trans",
       cutoff_cis=5e-1,cutoff_trans=5e-2,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/eqtl_highrisk.RData")

#for NORMAL highrisk mQTL
do_qtl(basedir="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL",snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_highrisk_SNP_ME.txt",
       snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",
       phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME.txt",
       phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_ME_POS.txt",
       covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_COVA_ME.txt",
       output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_cis",
       output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk_trans",
       cutoff_cis=5e-1,cutoff_trans=5e-2,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/NORMAL/mqtl_highrisk.RData")

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
