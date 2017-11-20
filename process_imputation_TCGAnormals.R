library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library(GenomicRanges)
library(sas7bdat)
library(gdata)

impfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation4"
#first combine ChrX imputations
impchrX_PAR1=read.table(paste0(impfolder,"/SNP6_imp_chrX_PAR1.txt"),stringsAsFactors = F)
doschrX_PAR1=read.table(paste0(impfolder,"/SNP6_imputed_dosages_chrX_PAR1.txt"),stringsAsFactors = F)
infchrX_PAR1=read.table(paste0(impfolder,"/SNP6_info_chrX_PAR1.txt"),header=T,stringsAsFactors = F)
sum(impchrX_PAR1$V2==infchrX_PAR1$rs_id)
impchrX_PAR2=read.table(paste0(impfolder,"/SNP6_imp_chrX_PAR2.txt"),stringsAsFactors = F)
doschrX_PAR2=read.table(paste0(impfolder,"/SNP6_imputed_dosages_chrX_PAR2.txt"),stringsAsFactors = F)
infchrX_PAR2=read.table(paste0(impfolder,"/SNP6_info_chrX_PAR2.txt"),header=T,stringsAsFactors = F)
sum(impchrX_PAR2$V2==infchrX_PAR2$rs_id)
impchrX_nonPAR=fread(paste0(impfolder,"/SNP6_imp_chrX_nonPAR.txt"),stringsAsFactors = F)
impchrX_nonPAR=as.data.frame(impchrX_nonPAR)
doschrX_nonPAR=fread(paste0(impfolder,"/SNP6_imputed_dosages_chrX_nonPAR.txt"),stringsAsFactors = F)
doschrX_nonPAR=as.data.frame(doschrX_nonPAR)
infchrX_nonPAR=read.table(paste0(impfolder,"/SNP6_info_chrX_nonPAR.txt"),header=T,stringsAsFactors = F)
nrow(infchrX_nonPAR)
#[1] 1482545
nrow(doschrX_nonPAR)
#[1] 1482545
nrow(impchrX_nonPAR)
#[1] 1482545
idx=match(impchrX_nonPAR$V2,infchrX_nonPAR$rs_id)
infchrX_nonPAR=infchrX_nonPAR[idx,]
sum(infchrX_nonPAR$rs_id==impchrX_nonPAR$V2)
doschrX=rbind(doschrX_PAR1,doschrX_nonPAR,doschrX_PAR2)
infchrX=rbind(infchrX_PAR1,infchrX_nonPAR,infchrX_PAR2)
write.table(doschrX,file=paste0(impfolder,"/SNP6_imputed_dosages_chr23.txt"),row.names = F,col.names = F,sep="\t",quote=F)
write.table(infchrX,file=paste0(impfolder,"/SNP6_info_chr23.txt"),row.names = F,col.names = T,sep=" ",quote=F)

#the headers of info file were removed
# process_imput=function(dosfile,impfile,inffile,inffile1)
# {
#   dos=fread(dosfile,stringsAsFactors = F)
#   dos=as.data.frame(dos)
#   imp=fread(impfile,stringsAsFactors = F)
#   imp=as.data.frame(imp)
#   inf=fread(inffile,header=T,stringsAsFactors = F) #contains multiple headers
#   inf=as.data.frame(inf)
#   idx=match(imp$V2,inf$rs_id)
#   inf1=inf[idx,]
#   print(sum(inf1$rs_id==imp$V2)==nrow(imp))
#   print(nrow(dos))
#   print(nrow(inf))
#   write.table(inf1,inffile1,row.names = F,col.names = T,sep=" ",quote=F)
# }
# for (i in 22:22)
# {
#   print(i)
#   dosfile=paste0(impfolder,"/SNP6_imputed_dosages_chr",i,".txt")
#   impfile=paste0(impfolder,"/SNP6_imp_chr",i,".txt")
#   inffile=paste0(impfolder,"/SNP6_info_chr",i,".txt")
#   inffile1=paste0(impfolder,"/SNP6_info1_chr",i,".txt")
#   process_imput(dosfile,impfile,inffile,inffile1)
# }

# impfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation2"
# #shapeit prephasing
# impfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation1"
# for (i in 22:22)
# {
#   print(i)
#   dosfile=paste0(impfolder,"/SNP6_imputed_dosages_chr",i,".txt")
#   impfile=paste0(impfolder,"/SNP6_imp_chr",i,".txt")
#   inffile=paste0(impfolder,"/SNP6_info_chr",i,".txt")
#   inffile1=paste0(impfolder,"/SNP6_info1_chr",i,".txt")
#   process_imput(dosfile,impfile,inffile,inffile1)
# }

#compare imputation scores

comp_imput=function(inffile1,inffile2,chr,inf1col=7,inf2col=10)
{
  
  inf1=fread(inffile1,header=T,stringsAsFactors = F) #contains multiple headers
  inf1=as.data.frame(inf1)
  inf1[,inf1col]=as.numeric(inf1[,inf1col])
  inf2=fread(inffile2,header=T,stringsAsFactors = F) #contains multiple headers
  inf2=as.data.frame(inf2)
  inf2[,inf2col]=as.numeric(inf2[,inf2col])
  # inf3=fread(inffile3,header=T,stringsAsFactors = F,fill=T) #contains multiple headers
  # inf3=as.data.frame(inf3)
  
  idx1=which(inf1$snp_id!="---")
  #number of genotyped snps
  ngtsnp1=length(idx1)
  inf1$genotyped=0
  inf1$genotyped[idx1]=1
  idx2=which(inf2$snp_id!="---")
  ngtsnp2=length(idx2)
  inf2$genotyped=0
  inf2$genotyped[idx2]=1

  #number of common genotyped snps
  ngtsnp12=sum(inf1$position[idx1] %in% inf2$position[idx2])
  #number of total snps
  nsnp1=nrow(inf1)
  nsnp2=nrow(inf2)
  compos=intersect(inf1$position,inf2$position)
  
  # #check referenced SNPs
  # idx1=which(inf1$snp_id=="---")
  # refinf1=inf1[idx1,]
  # idx2=which(inf2$snp_id=="---")
  # refinf2=inf2[idx2,]
  # #legend file
  # legendtable=fread(paste0("/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr",chr,".legend"))
  # legendtable=as.data.frame(legendtable)
  # sum(refinf1$rs_id %in% legendtable$id)
  # sum(refinf2$rs_id %in% legendtable$id)
  # #which(!refinf2$rs_id %in% legendtable$id)
  # refinf1=merge(refinf1,legendtable,by.x="rs_id",by.y="id")
  # refinf2=merge(refinf2,legendtable,by.x="rs_id",by.y="id")
  # sum(refinf2$rs_id %in% refinf1$rs_id)
  # #only appear in refinf1
  # idx=which(!refinf1$rs_id %in% refinf2$rs_id)
  # test=refinf1[idx,]
  # sum(test$SAS<=0.001)
  # #it turned out they only use EUR and EAS, and didn't use SAS
  
  idx1=match(compos,inf1$position)
  inf11=inf1[idx1,]
  idx2=match(compos,inf2$position)
  inf22=inf2[idx2,]
  #number of common snps
  nsnp12=length(compos)
  #proprotion of common snps
  nsnp12prop=nrow(inf22)/min(nrow(inf1),nrow(inf2))
  #mean score
  snp12_mscore1=round(mean(inf11[,inf1col],na.rm=T),digits = 2)
  snp12_mscore2=round(mean(inf22[,inf2col],na.rm=T),digits = 2)
  #median score
  snp12_mescore1=round(median(inf11[,inf1col],na.rm=T),digits = 2)
  snp12_mescore2=round(median(inf22[,inf2col],na.rm=T),digits = 2)
  #correlation of score
  snp12_mcor=cor(inf11[,inf1col],inf22[,inf2col])
  
  png(paste0("imputationscore_chr",chr,".png"))
  par(mfrow=c(2,1))
  hist(inf11[,inf1col],xlab="info_score",main="")
  hist(inf22[inf22[,inf2col]>=0,inf2col],xlab="r2_score",main="")
  dev.off()
  
  res=data.frame(ngtsnp1=ngtsnp1,ngtsnp2=ngtsnp2,ngtsnp12=ngtsnp12,nsnp1=nsnp1,nsnp2=nsnp2,nsnp12=nsnp12,
                 nsnp12prop=nsnp12prop,snp12_mscore1=snp12_mscore1,snp12_mscore2=snp12_mscore2,
                 snp12_mescore1=snp12_mescore1,snp12_mescore2=snp12_mescore2,snp12_mcor=snp12_mcor)
}

#impfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation1"
impfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation3"
icogfolder="/fh/fast/stanford_j/iCOGS_imp_1KGv3_FHCRC/Data/INFO"
oncofolder="/fh/fast/stanford_j/Oncoarray_Imputed_2016Jan25/Data/INFO"
compres=NULL
chrlen=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,
         107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
for (i in 1:22)
{
  cat(i,'..')
  #inffile1=paste0(impfolder,"/SNP6_info1_chr",i,".txt")
  inffile1=paste0(impfolder,"/SNP6_info_chr",i,".txt")
  inffile2=paste0(icogfolder,"/icogs_practical_info_chr",i,"_varid.txt")
  inffile3=paste0(oncofolder,"/onco_practical_info_chr",i,"_varid.txt")
  tmp=comp_imput(inffile1,inffile2,chr=i,inf1col=7,inf2col=10)
  compres=rbind(compres,tmp)
}
compres=cbind.data.frame(compres,len=chrlen[1:nrow(compres)])

sum(compres$ngtsnp1)
sum(compres$ngtsnp2)
png("numberofgenotypedSNPs.png")
par(mar=c(5.1,5.1,4.1,2.1))
boxplot(compres$ngtsnp1/compres$len*1000000,compres$ngtsnp2/compres$len*1000000,names=c("SNP6","iCOGS"),
        cex.names=1.4,cex.axis=1.4,ylab="#genotped SNPs/Mb",cex.lab=1.4)
dev.off()

sum(compres$nsnp1)
sum(compres$nsnp2)
png("numberofimputedSNPs.png")
par(mar=c(5.1,5.1,4.1,2.1))
boxplot(compres$nsnp1/compres$len*1000000,compres$nsnp2/compres$len*1000000,names=c("SNP6","iCOGS"),
        cex.names=1.4,cex.axis=1.4,ylab="#imputed SNPs/Mb",cex.lab=1.4)
dev.off()

png("propofoverlapedgenotypedSNPs.png")
par(mar=c(5.1,5.1,4.1,2.1))
boxplot(compres$ngtsnp12/compres$ngtsnp2, cex.names=1.4,cex.axis=1.4,ylab="Proportion of overlap",cex.lab=1.4)
dev.off()

png("propofoverlapedimputedSNPs.png")
par(mar=c(5.1,5.1,4.1,2.1))
boxplot(compres$nsnp12/compres$nsnp2, cex.names=1.4,cex.axis=1.4,ylab="Proportion of overlap",cex.lab=1.4)
dev.off()

extract_highrisk2=function(pre1="SNP6_info_",pre11="SNP6_imputed_dosages_")
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
TCGA_allsnps=allsnps
TCGA_highrisk=highrisk
sum(is.na(TCGA_allsnps$info))
#15
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/allhighrisksnps_new.RData")
par(mfrow=c(1,1))
png("Highrisk_score.png")
boxplot(TCGA_allsnps$info[TCGA_allsnps$info>0],allsnps$r2_icogs,names=c("SNP6","iCOGS"),cex.names=1.4,cex.axis=1.4,ylab="Imputation score",cex.lab=1.4)
dev.off()
t.test(TCGA_allsnps$info,allsnps$r2_icogs)

legendchr1_2012=fread("/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr1_impute.legend",header=T)
legendchr1_2012=as.data.frame(legendchr1_2012)

#2014 OCT
legendchr1=fread("/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr1.legend",header=T)
legendchr1=as.data.frame(legendchr1)

#filter impute
filter_imput=function(dosfile,inffile,infcol=7)
{
  dos=fread(dosfile,stringsAsFactors = F)
  dos=as.data.frame(dos)
  inf=fread(inffile,header=T,stringsAsFactors = F) #contains multiple headers
  inf=as.data.frame(inf)
  if (nrow(dos)!=nrow(inf)) warning("dosasge/inf0 files have problems!")
  quantile()
}


