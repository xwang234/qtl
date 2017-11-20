#!/usr/bin/env Rscript
#generate pad file and map file for plink
#load("/fh/fast/dai_j/CancerGenomics/Tools/wang/prostate/TCGAnormals.RData")
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/TCGAnormals.RData")
snp6_anno <- read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
#the SNP probes:
sum(rownames(genotypedata) %in% snp6_anno$Probe.Set.ID)==nrow(genotypedata)
idx=match(rownames(genotypedata),snp6_anno$Probe.Set.ID)
snp6_anno=snp6_anno[idx,]
snp6_anno$Physical.Position=as.integer(snp6_anno$Physical.Position)

#the input files of 35 samples
#outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation1/plink"

#the input files including all normals 496
#outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation3/plink"
#67 +385 other blood
outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation4/plink"
#create ped files for each chr
#http://zzz.bwh.harvard.edu/plink/data.shtml#ped
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype
#map file:
# chromosome (1-22, X, Y or 0 if unplaced)
# rs# or snp identifier
# Genetic distance (morgans)
# Base-pair position (bp units)
generateped_map=function(chr,genotypedata)
{
  if (chr==23) chr="X"
  idx=which(snp6_anno$Chromosome==chr & !is.na(snp6_anno$Physical.Position))
  annochrtable=snp6_anno[idx,]
  gtchrtable=genotypedata[idx,]
  if(sum(duplicated(annochrtable$Physical.Position))>0)
  {
    warning(paste0("There duplicated positions...","in chr",chr))
    idx1=!duplicated(annochrtable$Physical.Position)
    annochrtable=annochrtable[idx1,]
    gtchrtable=gtchrtable[idx1,]
  }
  
  #sort
  ordidx=order(annochrtable$Physical.Position)
  annochrtable=annochrtable[ordidx,]
  gtchrtable=gtchrtable[ordidx,]
  if (sum(annochrtable$Probe.Set.ID==rownames(gtchrtable))!=nrow(gtchrtable))
  {
    warning(paste0("The order of SNPs have problem...","in chr",chr))
  }
  pedfile=paste0(outfolder,"/TCGAnormals_chr",chr,".ped")
  ped=data.frame(matrix(NA,nrow=ncol(gtchrtable),ncol=6+nrow(gtchrtable)*2),stringsAsFactors = F)
  ped[,1]=colnames(gtchrtable)
  ped[,2]=1:nrow(ped) #individula id should be unique to use shapeit
  ped[,3]=ped[,4]=0
  ped[,5]=1
  ped[,6]=0
  #0->1 1;1->1 2;2->2 2
  for (i in 1:ncol(gtchrtable))
  {
    
    idx0=which(gtchrtable[,i]==0)
    if (length(idx0)>0) 
    {
      colid1=6+(idx0-1)*2+1
      colid2=colid1+1
      #ped[i,colid]=c(rep(1,length(idx0)),rep(1,length(idx0)))
      ped[i,colid1]=c(annochrtable$Allele.A[idx0])
      ped[i,colid2]=c(annochrtable$Allele.A[idx0])
    }
      
    idx1=which(gtchrtable[,i]==1)
    if (length(idx1)>0) 
    {
      colid1=6+(idx1-1)*2+1
      colid2=colid1+1
      ped[i,colid1]=c(annochrtable$Allele.A[idx1])
      ped[i,colid2]=c(annochrtable$Allele.B[idx1])
    }
  
    idx2=which(gtchrtable[,i]==2)
    if (length(idx2)>0)
    {
      colid1=6+(idx2-1)*2+1
      colid2=colid1+1
      ped[i,colid1]=c(annochrtable$Allele.B[idx2])
      ped[i,colid2]=c(annochrtable$Allele.B[idx2])
    }
  }
  write.table(ped,file=pedfile,col.names = F,row.names = F,sep="\t",quote=F)
  
  mapfile=paste0(outfolder,"/TCGAnormals_chr",chr,".map")
  map=data.frame(matrix(NA,nrow=nrow(gtchrtable),ncol=4),stringsAsFactors = F)
  map[,1]=chr
  map[,2]=annochrtable$Probe.Set.ID
  idxrs=which(annochrtable$dbSNP.RS.ID!="---")
  map[idxrs,2]=annochrtable$dbSNP.RS.ID[idxrs]
  map[,3]=0
  map[,4]=annochrtable$Physical.Position
  write.table(map,file=mapfile,col.names = F,row.names = F,sep="\t",quote=F)
  
  fliplistfile=paste0(outfolder,"/TCGAnormals_chr",chr,".fliplist")
  idxflip=which(annochrtable$Strand=="-")
  flip=data.frame(snpid=map[idxflip,2],stringsAsFactors = F)
  write.table(flip,file=fliplistfile,col.names = F,row.names = F,sep="\t",quote=F)
  return(0)
}
# for (chr in 1:23)
# {
#   cat(chr,'..')
#   generateped_map(snp6_anno,genotypedata,chr)
# }
#salloc -t 1-1 -n 50 mpirun -n 1 R --interactive

# outfolder
# #[1] "/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation4/plink"
# for (chr in 1:23)
# {
#   cat(chr,'..')
#   generateped_map(snp6_anno,genotypedata=allnormalgenotypedata,chr)
# }

library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(snp6_anno)
mpi.bcast.Robj2slave(allnormalgenotypedata)
mpi.bcast.Robj2slave(outfolder)
mpi.bcast.Robj2slave(generateped_map)
res=mpi.parSapply(X=1:23,FUN=generateped_map,genotypedata=allnormalgenotypedata,job.num=njobs)



