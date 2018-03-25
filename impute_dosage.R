#!/usr/bin/env Rscript
#it is faster than pandos

library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")

infolder=outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation_tumor"
infolder=outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation_tbd"

todosage=function(impfile,outfile)
{
  dat=fread(impfile,sep=" ")
  dat=as.data.frame(dat)
  nsample=(ncol(dat)-5)/3
  outdat=data.frame(matrix(NA,nrow(dat),ncol=nsample))
  for (j in 1:nsample)
  {
    nc=(j-1)*3+1+5
    tmp=dat[,nc+1]+dat[,nc+2]*2
    tmp=round(tmp,digits = 4)
    idxna=which(dat[,nc]==0 & dat[,nc+1]==0 & dat[,nc+2]==0)
    tmp[idxna]=NA
    outdat[,j]=tmp
  }
  fwrite(outdat,file=outfile,sep="\t",col.names = F)
}

for (i in 1:22)
{
  print(i)
  print(Sys.time())
  impfile=paste0(infolder,"/","SNP6_imp_chr",i,".txt")
  outfile=paste0(outfolder,"/","SNP6_imputed_dosages_chr",i,".txt")
  todosage(impfile,outfile)
}

impfile=paste0(infolder,"/","SNP6_imp_chrX_PAR1.txt")
outfile=paste0(outfolder,"/","SNP6_imputed_dosages_chrX_PAR1.txt")
todosage(impfile,outfile)
impfile=paste0(infolder,"/","SNP6_imp_chrX_PAR2.txt")
outfile=paste0(outfolder,"/","SNP6_imputed_dosages_chrX_PAR2.txt")
todosage(impfile,outfile)
impfile=paste0(infolder,"/","SNP6_imp_chrX_nonPAR.txt")
outfile=paste0(outfolder,"/","SNP6_imputed_dosages_chrX_nonPAR.txt")
todosage(impfile,outfile)

#TBD
for (i in 1:22)
{
  print(i)
  print(Sys.time())
  impfile=paste0(infolder,"/","genotypes_imp_chr",i,".txt")
  outfile=paste0(outfolder,"/","genotypes_imputed_dosages_chr",i,".txt")
  todosage(impfile,outfile)
}

impfile=paste0(infolder,"/","genotypes_imp_chrX_PAR1.txt")
outfile=paste0(outfolder,"/","genotypes_imputed_dosages_chrX_PAR1.txt")
todosage(impfile,outfile)
impfile=paste0(infolder,"/","genotypes_imp_chrX_PAR2.txt")
outfile=paste0(outfolder,"/","genotypes_imputed_dosages_chrX_PAR2.txt")
todosage(impfile,outfile)
impfile=paste0(infolder,"/","genotypes_imp_chrX_nonPAR.txt")
outfile=paste0(outfolder,"/","genotypes_imputed_dosages_chrX_nonPAR.txt")
todosage(impfile,outfile)
#combine ChrX imputations
impchrX_PAR1=read.table(paste0(impfolder,"/genotypes_imp_chrX_PAR1.txt"),stringsAsFactors = F)
doschrX_PAR1=read.table(paste0(impfolder,"/genotypes_imputed_dosages_chrX_PAR1.txt"),stringsAsFactors = F)
infchrX_PAR1=read.table(paste0(impfolder,"/genotypes_info_chrX_PAR1.txt"),header=T,stringsAsFactors = F)
sum(impchrX_PAR1$V2==infchrX_PAR1$rs_id)
impchrX_PAR2=read.table(paste0(impfolder,"/genotypes_imp_chrX_PAR2.txt"),stringsAsFactors = F)
doschrX_PAR2=read.table(paste0(impfolder,"/genotypes_imputed_dosages_chrX_PAR2.txt"),stringsAsFactors = F)
infchrX_PAR2=read.table(paste0(impfolder,"/genotypes_info_chrX_PAR2.txt"),header=T,stringsAsFactors = F)
sum(impchrX_PAR2$V2==infchrX_PAR2$rs_id)
impchrX_nonPAR=fread(paste0(impfolder,"/genotypes_imp_chrX_nonPAR.txt"),stringsAsFactors = F)
impchrX_nonPAR=as.data.frame(impchrX_nonPAR)
doschrX_nonPAR=fread(paste0(impfolder,"/genotypes_imputed_dosages_chrX_nonPAR.txt"),stringsAsFactors = F)
doschrX_nonPAR=as.data.frame(doschrX_nonPAR)
infchrX_nonPAR=read.table(paste0(impfolder,"/genotypes_info_chrX_nonPAR.txt"),header=T,stringsAsFactors = F)
nrow(infchrX_nonPAR)
#[1] 1478348
nrow(doschrX_nonPAR)
#[1] 1478348
nrow(impchrX_nonPAR)
#[1] 1478348
idx=match(impchrX_nonPAR$V2,infchrX_nonPAR$rs_id)
infchrX_nonPAR=infchrX_nonPAR[idx,]
sum(infchrX_nonPAR$rs_id==impchrX_nonPAR$V2)
doschrX=rbind(doschrX_PAR1,doschrX_nonPAR,doschrX_PAR2)
infchrX=rbind(infchrX_PAR1,infchrX_nonPAR,infchrX_PAR2)
write.table(doschrX,file=paste0(impfolder,"/genotypes_imputed_dosages_chr23.txt"),row.names = F,col.names = F,sep="\t",quote=F)
write.table(infchrX,file=paste0(impfolder,"/genotypes_info_chr23.txt"),row.names = F,col.names = T,sep=" ",quote=F)

