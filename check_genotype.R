source("functions.R")
library(data.table)

hutch_snp=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",header=T)
hutch_snppos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",header=T)
tbd_snp=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_GE.txt",header=T)
tbd_snppos=fread("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TBD_SNP_POS.txt",header=T)
hutch_snppos$allpos=paste0(hutch_snppos$chr,":",hutch_snppos$pos)
tbd_snppos$allpos=paste0(tbd_snppos$chr,":",tbd_snppos$pos)
sum(duplicated(hutch_snp$id))
tmp=which(duplicated(hutch_snppos$allpos))
idx_duplicatedpos_hutch=which(hutch_snppos$allpos %in% hutch_snppos$allpos[tmp])
View(hutch_snppos[idx_duplicatedpos_hutch,])
duplicatedsnp=tbd_snp$id[which(duplicated(tbd_snp$id))]
tmp=which(tbd_snp$id %in% duplicatedsnp)

tmp=which(!grepl(":",hutch_snppos$snp))

tmp=which(!grepl(":",tbd_snppos$snp))

sum(hutch_snp$id %in% tbd_snp$id)
snpid=intersect(hutch_snp$id,tbd_snp$id)
idx_inhutch=which(!hutch_snp$id %in% snpid)
idx_intbd=which(!tbd_snp$id %in% snpid)
View(hutch_snppos[idx_inhutch,])
View(tbd_snppos[idx_intbd,])
idx=match(snpid,hutch_snp$id)
hutch_snp1=hutch_snp[idx,]
idx=match(snpid,tbd_snp$id)
tbd_snp1=tbd_snp[idx,]
sum(tbd_snp1$id==hutch_snp1$id)
#hutch_snp1=hutch_snp1[,-1]
#tbd_snp1=tbd_snp1[,-1]
rowmean_hutch=rowMeans(hutch_snp1[,2:ncol(hutch_snp1)])
rowmean_tbd=rowMeans(tbd_snp1[,2:ncol(tbd_snp1)])
which.max(abs(rowmean_hutch-rowmean_tbd)) #[1] 5345895
testid=snpid[5345895]

idx1=which(tbd_snp$id %in% duplicatedsnp)
View(tbd_snppos[idx1,])
plot(rowmean_hutch[1:10000],rowmean_tbd[1:10000])
snp_pvalue=rep(NA,nrow(hutch_snp1))
for (i in 1:nrow(hutch_snp1))
{
  if (i %%10000==0) cat(i,"..")
  snp_pvalue[i]=t.test(unlist(hutch_snp1[i,]),
                       unlist(tbd_snp1[i,]))$p.value
}
qqplot(snp_pvalue[1:i])
#check those without impute
noimp_hutch=NULL
for (i in 1:nrow(hutch_snp1))
{
  if (i %% 2000==0) cat(i,'..')
  if (sum(hutch_snp1[i,]>0 & hutch_snp1[i,]<1)==0 & sum(hutch_snp1[i,]>1 & hutch_snp1[i,]<2)==0 )
  {
    noimp_hutch=c(noimp_hutch,i)
    if (length(noimp_hutch) %% 10==0) cat(length(noimp_hutch),'..')
  }
}

load("../data/TBD.RData")
