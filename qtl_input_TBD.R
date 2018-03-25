#!/usr/bin/env Rscript

#GE
load("../data/TBD.RData")
samples_geneexpdata=geneexp_sample$V2[match(colnames(geneexpdata),geneexp_sample$V5)] #make it as the order of samples
samples_genotype=genotype_sample$V2[match(colnames(genotypedata),genotype_sample$V5)]
source("functions.R")
#GE=removeconstrows(geneexpdata)
GE=t(scale(t(geneexpdata)))
colnames(GE)=samples_geneexpdata
sum(duplicated(rownames(GE))) #0 at gene level
tmp=cbind.data.frame(id=rownames(GE),GE)
write.table(tmp,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE.txt",col.names = T,row.names = F,sep="\t",quote=F)


#GE position
tbd_geneinfo=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/data/TBD/60547/PhenoGenotypeFiles/RootStudyConsentSet_phs000985.ProstateCancer_RiskSNPs.v1.p1.c1.DS-PC-PUB-MDS/GenotypeFiles/phg000725.v1.NCI_ProstateCancer_RiskSNPs.genotype-auxiliary-files.MULTI/geneinfo.txt",header=T)
sum(rownames(geneexpdata)==tbd_geneinfo$GeneID) #[1] 17233

# rnaseqv2_anno=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/TCGA.hg19.June2011.gaf",header=T,sep="\t",comment.char = "",stringsAsFactors = F)
# GE_POS_rnaseq2=data.frame(matrix(NA,nrow=nrow(GE),ncol=4))
# colnames(GE_POS_rnaseq2)=c("geneid","chr","s1","s2")
# GE_POS_rnaseq2$geneid=rownames(GE)
# numdup=rep(0,nrow(GE))
# for (i in 1:nrow(GE))
# {
#   if (i %% 2000 ==0) cat(i,"..")
#   #idx=which(rnaseqv2_anno$Gene==rownames(GE)[i])
#   idx=which(grepl(paste0(rownames(GE)[i],"|"),rnaseqv2_anno$Gene,fixed=T)==T)
#   locus=rnaseqv2_anno$GeneLocus[idx][1]
#   tmp=unlist(strsplit(locus,"chr"))
#   tmp=tmp[tmp!=""]
#   numdup[i]=length(tmp)
#   tmp1=unlist(strsplit(tmp,":"))
#   GE_POS_rnaseq2$chr[i]=tmp1[1]
#   startp=endp=rep(0,length(tmp))
#   
#   #pick the longer one
#   for (j in 1:length(tmp))
#   {
#     tmp1=unlist(strsplit(tmp[j],":"))
#     tmp2=as.integer(unlist(strsplit(tmp1[2],"-")))
#     startp[j]=tmp2[1]
#     endp[j]=tmp2[2]
#   }
#   idx1=which.max(endp-startp)
#   GE_POS_rnaseq2$s1[i]=startp[idx1]
#   GE_POS_rnaseq2$s2[i]=endp[idx1]
# }
# 
# refseq=fread("/fh/fast/stanford_j/Xiaoyu/Tools/reference/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/refFlat.txt",header=F,sep="\t",stringsAsFactors = F)
# refseq=as.data.frame(refseq)
# #sum(rownames(GE) %in% refseq$V1)
# GE_POS=data.frame(matrix(NA,nrow=nrow(GE),ncol=4))
# colnames(GE_POS)=c("geneid","chr","s1","s2")
# GE_POS$geneid=rownames(GE)
# numdup=rep(0,nrow(GE))
# for (i in 1:nrow(GE))
# {
#   if (i %% 2000 ==0) cat(i,"..")
#   #idx=which(grepl(paste0("gene_name ",'"',rownames(GE)[i],'"'),refseq$V9,fixed=T)==T)
#   idx=which(refseq$V1==rownames(GE)[i])
#   if (length(idx)>0)
#   {
#     # GE_POS$chr[i]=refseq$V1[idx[1]]
#     # GE_POS$s1[i]=min(refseq$V4[idx],na.rm=T)
#     # GE_POS$s2[i]=max(refseq$V5[idx],na.rm=T)
#     GE_POS$chr[i]=refseq$V3[idx[1]]
#     idx1=which.max(refseq$V6[idx]-refseq$V5[idx])
#     GE_POS$s1[i]=refseq$V5[idx[idx1]]
#     GE_POS$s2[i]=refseq$V6[idx[idx1]]
#   }
# }
# idx=which(is.na(GE_POS$chr) & !is.na(GE_POS_rnaseq2$chr))
# GE_POS[idx,]=GE_POS_rnaseq2[idx,]
# refseq1=fread("/fh/fast/stanford_j/Xiaoyu/Tools/reference/Homo_sapiens/GCF_000001405.25_GRCh37.p13_genomic.gff",header=F,skip=9,fill=T)
# idx=which(is.na(GE_POS$chr))
# chridx=which(grepl("chromosome=",refseq1$V9))
# for (i in idx)
# {
#   idx1=which(grepl(rownames(GE)[i],refseq1$V9))
#   if (length(idx1)>0)
#   {
#     idx2=chridx[chridx<=idx1[1]]
#     tmp=refseq1$V9[idx2[length(idx2)]]
#     tmp=unlist(strsplit(tmp,"chromosome=",fixed = T))[2]
#     GE_POS$chr[i]=unlist(strsplit(tmp,";"))[1]
#     idx2=which.max(refseq1$V5[idx1]-refseq1$V4[idx1])
#     GE_POS$s1[i]=refseq1$V4[idx1[idx2]]
#     GE_POS$s2[i]=refseq1$V5[idx1[idx2]]
#   }
# }
# 
# idx=which(!is.na(GE_POS$chr))
# write.table(GE_POS[idx,],file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",col.names = T,row.names = F,sep="\t",quote=F)
# tmp=cbind.data.frame(id=rownames(GE[idx,]),GE[idx,])
# write.table(tmp,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE.txt",col.names = T,row.names = F,sep="\t",quote=F)

GE_POS=data.frame(matrix(NA,nrow=nrow(geneexpdata),ncol=4))
colnames(GE_POS)=c("geneid","chr","s1","s2")
GE_POS$geneid=tbd_geneinfo$GeneID
GE_POS$chr=gsub("chr","",tbd_geneinfo$Chr)
GE_POS$s1=tbd_geneinfo$Start
GE_POS$s2=tbd_geneinfo$Stop
write.table(GE_POS,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_POS.txt",col.names = T,row.names = F,sep="\t",quote=F)


GE_PEER=peer_number(dat=GE)#15,0.594
colnames(GE_PEER)[2:ncol(GE_PEER)]=samples_geneexpdata
write.table(GE_PEER,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_GE_PEER.txt",col.names = T,row.names = F,sep="\t",quote=F)



extract_highrisk_tbd=function(pre1="genotypes_info_",pre11="genotypes_imputed_dosages_")
{
  library(gdata)
  load("/fh/fast/stanford_j/Xiaoyu/QTL/data/TBD.RData") #genotypedata
  load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData") #the order of highrisk and allsnps are not the same!
  impfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation_tbd"
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
          if (length(idx)>1) warning(paste0(i," has multype snps"))
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
  #samples_highrisk=genotype_sample$V2[match(colnames(highrisk)[6:ncol(highrisk)],genotype_sample$V5)]
  idx=match(samples_geneexpdata,samples_genotype)
  highrisk1=cbind(highrisk[,1:5],highrisk[,6:ncol(highrisk)][,idx])
  colnames(highrisk1)[6:ncol(highrisk1)]=samples_geneexpdata
  highrisk=highrisk1
  save(allsnps,highrisk,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/TBD_allhighrisksnps_new.RData")
  return(allsnps)
}

TBD_highrisk_SNP_GE=cbind.data.frame(id=highrisk$V1,highrisk[,6:ncol(highrisk)])
write.table(TBD_highrisk_SNP_GE,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_highrisk_SNP_GE.txt",row.names = F,col.names = T,sep="\t",quote=F)

#SNP data
folder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation_tbd"
pre="genotypes_info_"
pre1="genotypes_imputed_dosages_"
imputfiles=paste0(folder,"/",pre1,"chr",1:23,".txt")
file.exists(imputfiles)
infofiles=paste0(folder,"/",pre,"chr",1:23,".txt")
file.exists(infofiles)
for (chr in 1:23)
{
  cat(chr,"..")
  filterimput(chr,imputfiles,infofiles,colscore=7,colmaf=6,scorecutoff=0.3,mafcutoff=0.05,outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation",
              prefix="TBD")
}


#SNP_GE and SNP_POS
genereate_SNPfiles2=function(chrs=1:23,rdfiles,
                             colrsid=2,colpos=3,outfoler="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input",prefix="TBD")
{
  # SNP_GE=NULL
  # SNP_ME=NULL
  # SNP_POS=NULL
  #chrs=1:23
  for (chr in chrs)
  {
    print(chr)
    load(rdfiles[chr])
    imputetable=as.data.frame(imputetable)
    colnames(imputetable)[2:ncol(imputetable)]=samples_genotype #add this
    
    idx1=match(samples_geneexpdata,colnames(imputetable))
    #idx2=match(methylation_samples,colnames(imputetable))
    SNP_GE=cbind.data.frame(id=infotable[,colrsid],imputetable[,idx1])
    #SNP_ME=cbind.data.frame(id=infotable[,colrsid],imputetable[,idx2])
    SNP_POS=data.frame(snp=infotable[,colrsid],chr=chr,pos=infotable[,colpos])
    if (chr==chrs[1])
    {
      fwrite(SNP_GE,file=paste0(outfoler,"/",prefix,"_SNP_GE.txt"),col.names = T,row.names = F,sep="\t")
      #fwrite(SNP_ME,file=paste0(outfoler,"/",prefix,"_SNP_ME.txt"),col.names = T,row.names = F,sep="\t")
      fwrite(SNP_POS,file=paste0(outfoler,"/",prefix,"_SNP_POS.txt"),col.names = T,row.names = F,sep="\t")
    }else
    {
      fwrite(SNP_GE,file=paste0(outfoler,"/",prefix,"_SNP_GE.txt"),col.names = F,row.names = F,sep="\t",append = T)
      #fwrite(SNP_ME,file=paste0(outfoler,"/",prefix,"_SNP_ME.txt"),col.names = F,row.names = F,append = T,sep="\t")
      fwrite(SNP_POS,file=paste0(outfoler,"/",prefix,"_SNP_POS.txt"),col.names = F,row.names = F,sep="\t",append = T)
    }
  }
}

#SNP_GE,SNP_POS
rdfiles=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/result/filtered_imputation/TBD_chr",1:23,".RData")
genereate_SNPfiles2(chrs=1:23,rdfiles,colrsid=2,colpos=3,outfoler="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input",prefix="TBD")
tbd_snppos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_POS.txt",header=T)

TBD_SNP_GE=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_GE.txt",header=T,drop="id",sep="\t",stringsAsFactors = F)
TBD_SNP_GE_PCA=pca(dat=TBD_SNP_GE,check=T,bigpca = T,numpc=15)
write.table(TBD_SNP_GE_PCA,file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_GE_PCA.txt",col.names = T,row.names = F,sep="\t",quote=F)

