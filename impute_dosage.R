#!/usr/bin/env Rscript
#it is faster than pandos

library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")

infolder=outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation_tumor"

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
