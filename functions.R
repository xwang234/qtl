#/usr/bin/env Rscript
library(data.table,lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library(GenomicRanges)

removeconstrows=function(dat)
{
  idxconst=rep(F,nrow(dat))
  for (i in 1:nrow(dat))
  {
    if (i %% 10000==0) cat(i,"..")
    if (var(unlist(dat[i,]))==0 | is.na(var(unlist(dat[i,])))) idxconst[i]=T
  }
  dat=dat[!idxconst,]
}

library(impute)
imputemissdat=function(dat)
{
  tmp=impute.knn(as.matrix(dat))
  res=as.data.frame(tmp$data)
}

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

#PCA
pca=function(dat=GE,check=F,bigpca=F,numpc=NULL)
{
  if (bigpca==T)
  {
    #use random rows
    idx=sample(1:nrow(dat))[1:floor(nrow(dat)*0.1)]
    dat=as.data.frame(dat[idx,])
    # #use means of 10 rows
    # idx=seq(1,nrow(dat),by=10)
    # if (idx[length(idx)]!=nrow(dat))
    # {
    #   idx=c(idx,nrow(dat))
    # }
    # dat1=data.frame(matrix(NA,nrow=length(idx)-1,ncol=ncol(dat)))
    # colnames(dat1)=colnames(dat)
    # for (i in 1:(length(idx)-1))
    # {
    #   dat1[i,]=colMeans(dat[idx[i]:idx[i+1],],na.rm = T)
    # }
    # dat=dat1
  }
  if (check==T)
  {
    #remove constant rows
    idxconst=rep(F,nrow(dat))
    for (i in 1:nrow(dat))
    {
      if (var(unlist(dat[i,]))==0) idxconst[i]=T
    }
    dat=dat[!idxconst,]
  }
  pc=prcomp(t(dat),scale=T,center=T)
  if (is.null(numpc))
  {
    varprop=pc$sdev^2/sum(pc$sdev^2)
    numpc=sum(varprop>0.01) #9
    print(paste0("number of PCs: ",numpc))
    print(paste0("proportion of variance explained: ",round(sum(varprop[1:numpc]),digits = 3)))
  }
  res=data.frame(matrix(NA,nrow=numpc,ncol=1+ncol(dat)))
  colnames(res)=c("id",colnames(dat))
  res$id=paste0("pc",1:numpc)
  idx=match(colnames(dat),rownames(pc$x))
  res[,2:ncol(res)]=t(pc$x[idx,])[1:numpc,]
  return(res)
}

bigpca=function(dat=t(HUTCH_SNP_GE),prefix=1)
{
  if (class(dat[1,1])=="character")
  {
    if (ncol(dat)<nrow(dat))
    {
      for (i in 1:ncol(dat))
      {
        dat[,i]=as.numeric(dat[,i])
      }
    }else
    {
      for (i in 1:nrow(dat))
      {
        dat[i,]=as.numeric(dat[i,])
      }
    }
  }
  
  library(bigpca)
  orig.dir <- getwd(); #setwd(tempdir()); # move to temporary dir
  filebck=paste0("bigpca",prefix,".bck")
  filedsc=paste0("bigpca",prefix,".dsc")
  if(file.exists(filebck)) { unlink(c(filebck,filedsc)) }
  bM <- as.big.matrix(dat,backingfile = filebck,  backingpath = getwd(), descriptorfile = filedsc)
  # for (i in 1:ncol(dat))
  # {
  #   bM[1:nrow(dat),i] <-dat[,i]
  # }
  #bM <- get.big.matrix(filedsc)
  lmat <- thin(bM,.10,rows=F,pref="unif")
  res=big.PCA(lmat,pcs.to.keep = ncol(dat),thin = F,return.loadings = T)
  #res=big.PCA(bM,pcs.to.keep=15,thin = T,return.loadings = T)
  varprop=res$Evalues/sum(res$Evalues)
  numpc=sum(varprop>0.01)
  print(paste0("number of PCs: ",numpc))
  print(paste0("proportion of variance explained: ",round(sum(varprop[1:numpc]),digits = 3)))
  numpc=15
  res1=data.frame(matrix(NA,nrow=numpc,ncol=1+nrow(dat)))
  colnames(res1)=c("id",rownames(dat))
  res1$id=paste0("pc",1:numpc)
  idx=match(rownames(dat),rownames(res$loadings))
  res1[,2:ncol(res1)]=t(res$PCs[idx,])[1:numpc,]
  unlink(c(filebck,filedsc))
  return(res1)
}

#peer factor, need to impute missing data first (scale constant array would introduce NA)
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
    #PEER_setNmax_iterations(model, 1000)
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

#get peer factors by setting the number of factors first
peer_number=function(dat=GE,numfactors=15)
{
  
  #scale the data first
  dat=t(scale(t(dat)))
  totvar=0
  for (i in 1:nrow(dat)) totvar=totvar+var(unlist(dat[i,]))

  library(peer)
  
  model=PEER()
  PEER_setNk(model,numfactors)
  PEER_setPhenoMean(model, as.matrix(t(dat)))
  #PEER_setNmax_iterations(model, 1000)
  PEER_update(model)
  resi=PEER_getResiduals(model)
  resivar=0
  for (j in 1:ncol(resi)) resivar=resivar+var(resi[,j])
  print((totvar-resivar)/totvar)
  varexplain=(totvar-resivar)/totvar
    
  print(paste0("number of factors: ",numfactors))
  print(paste0("variance explained: ",round(varexplain,digits = 3)))
  #the factors
  x=PEER_getX(model)
  x=as.data.frame(t(x))
  colnames(x)=colnames(dat)
  res=cbind.data.frame(id=paste0("factor",1:numfactors),x)
  return(res)
}

#covariates need to be determined, use a temporary file as input
do_qtl=function(snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_GE.txt",
                snpposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
                phenotypeposfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE_POS.txt",
                covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE.txt",
                output_cis="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_cis",
                output_trans="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl_trans",
                cutoff_cis=1e-6,cutoff_trans=1e-8,recordfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/HUTCH/eqtl.RData")
{
  library("MatrixEQTL")
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
  gene$fileSliceSize = 10000;      # read file in slices of 10,000 rows
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

read_qtl_input=function(filename)
{
  dat=read.table(filename,header=T,sep="\t",stringsAsFactors = F)
  colnames(dat)=gsub("^X","",colnames(dat))
  colnames(dat)=gsub(".","-",colnames(dat),fixed = T)
  return(dat)
}
write_qtl_input=function(dat,filename)
{
  write.table(dat,file=filename,col.names = T,row.names = F,sep="\t",quote=F)
}

totalpos=function(chr=NULL,pos=NULL)
{
  chrs=c(1:22,"X","Y")
  dictfile="/fh/fast/dai_j/CancerGenomics/Tools/database/reference/compact/ucsc.hg19.compact.dict"
  chrlen=read.table(file=dictfile,sep="\t",skip=1,stringsAsFactors = F)
  chrlen=chrlen$V3
  chrlen=gsub("LN:","",chrlen,fixed=T)
  chrlen=as.numeric(chrlen)
  names(chrlen)=chrs
  chrstart=rep(0,length(chrlen))
  for (i in 2:length(chrlen))
  {
    tmp=0
    for (j in 1:(i-1))
    {
      tmp=tmp+chrlen[j]
    }
    chrstart[i]=tmp
  }
  names(chrstart)=names(chrlen)
  chr=gsub(23,"X",chr)
  chr=gsub(24,"Y",chr)
  
  idx=which(chr %in% chrs)
  chr=chr[idx]
  pos=pos[idx]
  #chr=factor(chr,chrs)
  
  allchrs=unique(chr)
  posall=rep(NA,length(chr))
  for (mychr in allchrs)
  {
    idx=which(chr %in% mychr)
    idx1=which(names(chrstart)==mychr)
    posall[idx]=pos[idx]+chrstart[idx1]
  }
  return(posall)
}

#extract qtl result
readqtlres=function(qtlresfile="eqtl_cis",colpvalue=5,snpposfile="eqtl_POS_SNP_GE.txt",geposfile="eqtl_POS_GE.txt",bonfadj=NULL,fdrcutoff=0.05)
{
  print(paste0("processed qtl file: ",qtlresfile))
  snpspos=fread(snpposfile,sep="\t")
  print(paste0("number of SNPs: ",nrow(snpspos)))
  snpspos$chr=gsub(23,"X",snpspos$chr)
  snpspos$chr=gsub(24,"Y",snpspos$chr)
  
  genepos=fread(geposfile,sep="\t")
  print(paste0("number of phenotypeprobes: ",nrow(genepos)))
  genepos$chr=gsub(23,"X", genepos$chr)
  genepos$chr=gsub(24,"Y", genepos$chr)
  
  eqtl_cis=fread(qtlresfile,header=T,sep="\t",fill=T,check.names = F)
  eqtl_cis=as.data.frame(eqtl_cis)
  if (! is.null(bonfadj))
  {
    print(paste0("number of pairs: ",nrow(eqtl_cis)))
    print(paste0("number of unique SNPs: ",length(unique(eqtl_cis$SNP))))
    print(paste0("number of unique phenotypes: ",length(unique(eqtl_cis$gene))))
    print(paste0("Bonferroni cutoff: ",0.05/bonfadj))
    idx=which(eqtl_cis[,colpvalue]<=0.05/bonfadj)
    print(paste0("number of pairs (after adjustment)",length(idx)))
    print(paste0("number of unique SNPs (after adjustment)",length(unique(eqtl_cis$SNP[idx])),", ",length(unique(eqtl_cis$SNP[idx]))/nrow(snpspos)))
    print(paste0("number of unique phenotypes (after adjustment)",length(unique(eqtl_cis$gene[idx])),", ",length(unique(eqtl_cis$gene[idx]))/nrow(genepos)))
  }else
  {
    print(paste0("number of pairs: ",nrow(eqtl_cis)))
    print(paste0("number of unique SNPs: ",length(unique(eqtl_cis$SNP))))
    print(paste0("number of unique phenotypes: ",length(unique(eqtl_cis$gene))))
    #fdrcutoff=0.05
    print(paste0("FDR cutoff: ",fdrcutoff))
    tmp=which(eqtl_cis$FDR<=fdrcutoff)
    print(paste0("pvalue cutoff: ",eqtl_cis[max(tmp),colpvalue]))
    idx=which(eqtl_cis$FDR<=fdrcutoff)
    print(paste0("number of pairs (after adjustment)",length(idx)))
    print(paste0("number of unique SNPs (after adjustment)",length(unique(eqtl_cis$SNP[idx])),", ",length(unique(eqtl_cis$SNP[idx]))/nrow(snpspos)))
    print(paste0("number of unique phenotypes (after adjustment)",length(unique(eqtl_cis$gene[idx])),", ",length(unique(eqtl_cis$gene[idx]))/nrow(genepos)))
  }
  
  eqtl_cis=eqtl_cis[idx,]
  
  idx=match(eqtl_cis$SNP,snpspos$snp)
  eqtl_cis$SNP_chr=snpspos$chr[idx]
  eqtl_cis$SNP_pos=snpspos$pos[idx]
  #which snps were included
  eqtl_cis$SNP_idx=idx
  
  idx=match(eqtl_cis$gene,genepos$geneid)
  eqtl_cis$gene_chr=genepos$chr[idx]
  eqtl_cis$gene_pos=genepos$s1[idx]
  #get the total position
  eqtl_cis$SNP_posall=totalpos(chr=eqtl_cis$SNP_chr,pos=eqtl_cis$SNP_pos)
  eqtl_cis$gene_posall=totalpos(chr=eqtl_cis$gene_chr,pos=eqtl_cis$gene_pos)
  levelmat=data.frame(value=-log10(eqtl_cis[,colpvalue]),pos1=eqtl_cis$SNP_posall,pos2=eqtl_cis$gene_posall,
                      chr=eqtl_cis$SNP_chr,gene_chr=eqtl_cis$gene_chr,snp_idx=eqtl_cis$SNP_idx,
                      opos1=eqtl_cis$SNP_pos,opos2=eqtl_cis$gene_pos,gene=eqtl_cis$gene,fdr=eqtl_cis$FDR,tstat=eqtl_cis$`t-stat`,beta=eqtl_cis$beta)
  levelmat$chr=gsub(23,"X",levelmat$chr)
  levelmat$chr=gsub(24,"Y",levelmat$chr)
  # cc = colorRampPalette( c("green","white","red"))
  # numcol=21
  # colors=cc(numcol)
  # zrng <- range(levelmat$value)       # what's the range of z
  # tol <- 1e-2            # what tolerance is necessary?
  # zquant=quantile(levelmat$value,c(0,0.5,0.9,1))
  # colorBreaks1 <- seq(zquant[1]-0.01,zquant[2],length.out = ceiling(numcol*0.5))
  # colorBreaks2 <- seq(zquant[2]+0.1,zquant[4]+0.01,length.out = numcol-ceiling(numcol*0.5))
  # colorBreaks <- c(colorBreaks1,colorBreaks2)
  # df_breaks=seq(zrng[1],zrng[2],length.out = numcol)
  # # levelplot(cor~ge*me,data=levelmat,xlab="Gene expression",ylab="Methylation")
  # levelmat$color=rep(NA,nrow(levelmat))
  # for (i in 1:(numcol-1))
  # {
  #   idx=which(levelmat$value>=colorBreaks[i] & levelmat$value<colorBreaks[i+1])
  #   levelmat$color[idx]=i
  # }
  # # idx=which(levelmat$color>numcol/2)
  # # levelmat=levelmat[idx,]
  # levelmat$colors=colors[levelmat$color]
  levelmat=levelmat[order(levelmat$value),]
  return(levelmat)
}

#Add gene names to qtl result
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/Geneexp_map_TCGA_HUTCH.RData")
addgenename_GE=function(datatable=tmp5,anno=ge_annohg19,colname="genename")
{
  if (!colname %in% colnames(datatable))
  {
    idx=match(datatable$gene,anno$Probe_Id)
    datatable=cbind(datatable,anno$Symbol[idx])
    colnames(datatable)[ncol(datatable)]=colname
  }
  return(datatable)
}

load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
addgenename_ME=function(datatable=tmp7,anno1=anno,colname="genename",splitstr=";")
{
  if (!colname %in% colnames(datatable))
  {
    idx=match(datatable$gene,anno1$IlmnID)
    datatable=cbind(datatable,anno1$UCSC_RefGene_Name[idx])
    colnames(datatable)[ncol(datatable)]=colname
    datatable[,colname]=as.character(datatable[,colname])
    idx=which(datatable[,colname]=="")
    datatable[idx,colname]=NA
    for (i in 1:nrow(datatable))
    {
      if (i%%10000==0) cat(i,"..")
      if (!is.na(datatable[i,colname]))
      {
        tmp=unlist(strsplit(datatable[i,colname],split=splitstr,fixed = T))
        tmp=unique(tmp)
        datatable[i,colname]=paste0(tmp,collapse = "|")
      }
    }
  }
  return(datatable)
}

checkoverlapsnps=function(qtlres1=HUTCH_cis_eqtl,qtlres2=NORMAL_cis_eqtl,
                          snpposfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_SNP_POS.txt",
                          snpposfile2="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/NORMAL_SNP_POS.txt")
{
  snppos1=fread(snpposfile1,sep="\t")
  snppos2=fread(snpposfile2,sep="\t")
  gr_res1=GRanges(seqnames = qtlres1$chr,ranges = IRanges(start=snppos1$pos[qtlres1$snp_idx],width=1))
  gr_res1=unique(gr_res1)
  gr_res2=GRanges(seqnames = qtlres2$chr,ranges = IRanges(start=snppos2$pos[qtlres2$snp_idx],width=1))
  gr_res2=unique(gr_res2)
  print(paste0("number of unique SNP1:",length(gr_res1)))
  print(paste0("number of unique SNP2:",length(gr_res2)))
  print(paste0("number of overlapped SNP:",length(intersect(gr_res1,gr_res2))))
  gr1=GRanges(seqnames = snppos1$chr,ranges=IRanges(start=snppos1$pos,width=1))
  gr2=GRanges(seqnames = snppos2$chr,ranges=IRanges(start=snppos2$pos,width=1))
  gr=unique(c(gr1,gr2))
  pvalue_overlapsnp(length(gr),length(gr_res1),length(gr_res2),length(intersect(gr_res1,gr_res2)))
}

pvalue_overlapsnp=function(totalnum=147,num1=33,num2=127,num12=32)
{
  m=matrix(c(num12,num1-num12,num2-num12,totalnum-num1-num2+num12),nrow=2,dimnames = list(c("in1","notin1"),c("in2","notin2")))
  print(m)
  print(fisher.test(m)$p.value)
}

checkoverlappairs_sameplatform=function(qtlres1=highrisk_HUTCH_cis_eqtl,qtlres2=highrisk_HUTCH_pca_cis_eqtl)
{
  print(paste0("number of pairs in qtlres1: ",nrow(qtlres1)))
  print(paste0("number of pairs in qtlres2: ",nrow(qtlres2)))
  pairs1=qtlres1[,c("snp_idx","gene")]
  pairs2=qtlres2[,c("snp_idx","gene")]
  allpairs=rbind(pairs1,pairs2)
  allpairs=unique(allpairs)
  print(paste0("total number of pairs: ",nrow(allpairs)))
  tmp=merge(qtlres1,qtlres2,by=c("snp_idx","gene"))
  print(paste0("number of common pairs: ",nrow(tmp)))
  return(tmp)
}

#use gene name for validating pairs
validatepairs_difplatform=function(qtlres1=highrisk_HUTCH_cis_eqtl,qtlres2=highrisk_TCGA_tumors_cis_eqtl_all,phenotypemap=hutch_ge_anno,opt="pvalue")
{
  tmp=merge(qtlres1,qtlres2,by=c("snp_idx","genename"))
  if (opt=="pvalue")
  {
    idx=10^-tmp$value.y<=0.05 & tmp$beta.x*tmp$beta.y>0
  }else
  {
    idx=10^-tmp$value.y<=0.05/nrow(tmp) & tmp$beta.x*tmp$beta.y>0 #FWER
  }
  print(paste0("number of qtl pairs to validate:",nrow(qtlres1)))
  tmp1=merge(qtlres1,phenotypemap,by.x="gene",by.y="Probe_Id")
  print(paste0("number of qtl pairs can be validated:",sum(!is.na(tmp1$TCGA_Probe_Id))))
  print(paste0("number of pairs validated:",sum(idx)))
  
  print(paste0("number of phenotype (probe) to validate:",length(unique(qtlres1$gene))))
  print(paste0("number of phenotype (probe) can be validated:",length(unique(tmp1$gene[!is.na(tmp1$TCGA_Probe_Id)]))))
  print(paste0("number of phenotype (probe) validated:",length(unique(tmp$gene.x[idx]))))
  print(paste0("number of phenotype (gene) to validate:",length(unique(qtlres1$genename))))
  print(paste0("number of phenotype (gene) can be validated:",length(unique(tmp1$genename[!is.na(tmp1$TCGA_Probe_Id)]))))
  print(paste0("number of phenotype (gene) validated:",length(unique(tmp$genename[idx]))))
  
  print(paste0("number of snp to validate:",length(unique(qtlres1$snp_idx))))
  print(paste0("number of snp can be validate within pairs:",length(unique(tmp1$snp_idx[!is.na(tmp1$TCGA_Probe_Id)]))))
  print(paste0("number of snp validated:",length(unique(tmp$snp_idx[idx]))))
  if (sum(idx)>0)
  {
    # tmp1=t.test(tmp$value.x[idx],tmp$value.x[!idx])
    # print(tmp1$p.value)
    
    return(tmp[idx,])
  }else
  {
    return(NA)
  }
}

#in the same platforms
validatepairs=function(qtlres1=highrisk_cis_hutch_eqtl,qtlres2=highrisk_all_NORMAL_cis_eqtl,phenotype2=tcga_ge_pos,opt="pvalue")
{
  tmp=merge(qtlres1,qtlres2,by=c("snp_idx","gene"))
  if (opt=="pvalue")
  {
    idx=10^-tmp$value.y<=0.05 & tmp$beta.x*tmp$beta.y>0
  }else
  {
    idx=10^-tmp$value.y<=0.05/nrow(tmp) & tmp$beta.x*tmp$beta.y>0
  }
  print(paste0("number of qtl pairs to validate:",nrow(qtlres1)))
  print(paste0("number of qtl pairs can be validated:",sum(qtlres1$gene %in% phenotype2[,1],na.rm=T)))
  print(paste0("number of pairs validated:",sum(idx)))
  
  print(paste0("number of phenotype to validate:",length(unique(qtlres1$gene))))
  print(paste0("number of phenotype can be validated:",length(unique(qtlres1$gene[qtlres1$gene %in% phenotype2[,1]]))))
  print(paste0("number of phenotype validated:",length(unique(tmp$gene[idx]))))
  
  print(paste0("number of snp to validate:",length(unique(qtlres1$snp_idx))))
  print(paste0("number of snp can be validated within pairs:",length(unique(qtlres1$snp_idx[qtlres1$gene %in% phenotype2[,1]]))))
  print(paste0("number of snp validated:",length(unique(tmp$snp_idx[idx]))))
  if (length(idx)>0)
  {
    return(tmp[idx,])
  }else
  {
    return(NA)
  }
}

#transform snp_idx to snp names
update_highrisk_snp_idx=function(dat=validatepairs_HUTCH_ciseqtl)
{
  if ("snp_idx" %in% colnames(dat))
  {
    refsnps=read.table("../data/RESUB_Supplementary_Table16_v9.txt",header=T,sep="\t",stringsAsFactors = F)
    refsnps$Chr=gsub(23,"X",refsnps$Chr)
    highrisk_snp=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",header=T,stringsAsFactors = F)
    highrisk_snp$name=NA
    for (i in 1:nrow(highrisk_snp))
    {
      idx=which(refsnps$Chr==highrisk_snp$chr[i] & refsnps$Position==highrisk_snp$pos[i])
      highrisk_snp$name[i]=refsnps$SNP[idx]
    }
    dat$snp_idx=highrisk_snp$name[dat$snp_idx]
    colnames(dat)=gsub("snp_idx","snp",colnames(dat))
  }
  return(dat)
}

qqplot=function(pvalue=NULL,main="")
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="expected p-value (log base 10)",
       ylab="observed p-value (log base 10)",main=main,cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
}

#mannually compute the pvalue,t,beta; snpidx phenotypeidx store the pairs to compute
computep=function(snpidx=c(144,33),phenotypeidx=c(12142,8562),
                  snpfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_GE_updatename.txt",
                  phenotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE.txt",
                  covariatefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_GE_PEER_used.txt")
{
  snp=read_qtl_input(snpfile)
  rownames(snp)=snp[,1]
  snp=snp[,-1]
  phenotype=fread(phenotypefile,header=T,sep="\t")
  phenotype=as.data.frame(phenotype)
  rownames(phenotype)=phenotype[,1]
  phenotype=phenotype[,-1]
  covariate=read_qtl_input(covariatefile)
  rownames(covariate)=covariate[,1]
  covariate=covariate[,-1]
  covariate=t(covariate)
  covariate=as.data.frame(covariate)
  res=data.frame(pvalue=rep(NA,length(snpidx)),tstat=rep(NA,length(snpidx)),beta=rep(NA,length(snpidx)))
  for (i in 1:length(snpidx))
  {
    snpidx1=snpidx[i]
    phenotypeidx1=phenotypeidx[i]
    fm=paste0("unlist(phenotype[phenotypeidx1,])~unlist(snp[snpidx1,])+",
              paste0(colnames(covariate),collapse = "+"))
    fit=glm(as.formula(fm),data=covariate)
    #plot(fit)
    res$pvalue[i]=summary(fit)$coefficients[2,4]
    res$tstat[i]=summary(fit)$coefficients[2,3]
    res$beta[i]=summary(fit)$coefficients[2,1]
  }
  return(res)
}

#overlap of eqtl and mqtl
overlappairs=function(dat1=highrisk_cis_hutch_eqtl,dat2=highrisk_cis_hutch_mqtl)
{
  tmp=merge(dat1,dat2,by="snp_idx")
  print(paste0("number of overlapped SNP:",length(unique(tmp$snp_idx))))
  print(paste0("number of phenoytype1:",length(unique(tmp$gene.x))))
  print(paste0("number of phenoytype2:",length(unique(tmp$gene.y))))
  print(paste0("number of triplets:",nrow(tmp)))
  idx1=which(as.character(tmp$gene.x)==as.character(tmp$gene.y))
  print(paste0("number of overapped pairs:",sum(as.character(tmp$gene.x)==as.character(tmp$gene.y),na.rm=0)))
  return(list(snpidx=unique(tmp$snp_idx),pair=tmp[idx1,c("snp_idx","gene.x")],triplet=tmp))
}

validtriplet=function(dat1=olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl$triplet,validateeqtl=highrisk_all_NORMAL_cis_eqtl,
                      validatemqtl=highrisk_all_NORMAL_cis_mqtl,phenotype1=tcga_ge_pos,phenotype2=tcga_me_pos)
{
  validateeqtl$gene=as.character(validateeqtl$gene)
  validatemqtl$gene=as.character(validatemqtl$gene)
  dat1$gene.x=as.character(dat1$gene.x)
  dat1$gene.y=as.character(dat1$gene.y)
  dat1$validate=0
  dat1$validate_eqtl_value=NA
  dat1$validate_mqtl_value=NA
  print(paste0("number of triplets to validate:",nrow(dat1)))
  print(paste0("number of triplets can be validate:",sum(dat1$gene.x %in% phenotype1[,1] & dat1$gene.y %in% phenotype2[,1],na.rm=T)))
  for (i in 1:nrow(dat1))
  {
    idx1=which(validateeqtl$snp_idx==dat1$snp_idx[i] & validateeqtl$gene==dat1$gene.x[i])
    idx2=which(validatemqtl$snp_idx==dat1$snp_idx[i] & validatemqtl$gene==dat1$gene.y[i])
    if (length(idx1)*length(idx2)==1)
    {
      if (10^-validateeqtl$value[idx1]<=0.05 & 10^-validatemqtl$value[idx2]<=0.05)
      {
        dat1$validate[i]=1
        dat1$validate_eqtl_value[i]=validateeqtl$value[idx1]
        dat1$validate_mqtl_value[i]=validatemqtl$value[idx2]
      }
    }
  }
  print(paste0("number of triplets validated:", sum(dat1$validate==1,na.rm = T)))
  print(paste0("number of snp to validate:", length(unique(dat1$snp_idx))))
  print(paste0("number of snp validated:", length(unique(dat1$snp_idx[dat1$validate==1]))))
  return(dat1)
}
#validate at gene level
validtriplet_difplatform=function(dat1=olap_highrisk_cis_hutch_eqtl_cis_hutch_mqtl$triplet,validateeqtl=highrisk_all_NORMAL_cis_eqtl,
                      validatemqtl=highrisk_all_NORMAL_cis_mqtl,phenotypemap=hutch_ge_anno,phenotype2=tcga_me_pos)
{
  validateeqtl$gene=as.character(validateeqtl$gene)
  validatemqtl$gene=as.character(validatemqtl$gene)
  validateeqtl$genename=as.character(validateeqtl$genename)
  validatemqtl$genename=as.character(validatemqtl$genename)
  dat1$gene.x=as.character(dat1$gene.x)
  dat1$gene.y=as.character(dat1$gene.y)
  dat1$genename.x=as.character(dat1$genename.x)
  dat1$genename.y=as.character(dat1$genename.y)
  idx=which(!is.na(dat1$genename.x))
  dat1=dat1[idx,]
  dat1$validate=0
  dat1$validate_eqtl_value=NA
  dat1$validate_mqtl_value=NA
  print(paste0("number of triplets to validate:",nrow(dat1)))
  phenotypemap$gene.x=phenotypemap$Probe_Id
  dat1=dat1[order(dat1$gene.x),] #need to be ordered first
  tmp1=merge(dat1,phenotypemap,by="gene.x",sort=F)
  dat2=dat1[!is.na(tmp1$TCGA_Probe_Id),]
  print(paste0("number of triplets can be validated considering eqtl genename mapping:",nrow(dat2)))
  idx=which(dat2$gene.y %in% phenotype2[,1])
  dat2=dat2[idx,]
  print(paste0("number of triplets can be validated:",nrow(dat2)))
  for (i in 1:nrow(dat2))
  {
    idx1=which(validateeqtl$snp_idx==dat2$snp_idx[i] & validateeqtl$genename==dat2$genename.x[i])
    idx2=which(validatemqtl$snp_idx==dat2$snp_idx[i] & validatemqtl$gene==dat2$gene.y[i])
    if (length(idx1)*length(idx2)==1)
    {
      if (10^-validateeqtl$value[idx1]<=0.05 & 10^-validatemqtl$value[idx2]<=0.05 &dat2$beta.x[i]*validateeqtl$beta[idx1]>0 &dat2$beta.y[i]*validatemqtl$beta[idx2]>0)
      {
        dat2$validate[i]=1
        dat2$validate_eqtl_value[i]=validateeqtl$value[idx1]
        dat2$validate_mqtl_value[i]=validatemqtl$value[idx2]
      }
    }
  }
  print(paste0("number of triplets validated:", sum(dat2$validate==1,na.rm = T)))
  print(paste0("number of snp can be validate:", length(unique(dat2$snp_idx))))
  print(paste0("number of snp validated:", length(unique(dat2$snp_idx[dat2$validate==1]))))
  return(dat2)
}

#put dat1 as qtl of T, dat2 as qtl of G
citqtlpairs_inputfromfiles=function(dat1=highrisk_HUTCH_cis_eqtl,dat2=highrisk_HUTCH_cis_mqtl,
                                    GE_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_GE.txt",
                                    ME_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_ME.txt",
                                    SNP_GE_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_GE.txt",
                                    covariate_ge_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_GE_PEER_used.txt",
                                    covariate_me_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_COVA_ME_used.txt")
{
  library(cit)
  GE=fread(GE_file,header=T,sep="\t")
  GE=as.data.frame(GE)
  rownames(GE)=GE[,1]
  GE=GE[,-1]
  ME=fread(ME_file,header=T,sep="\t")
  ME=as.data.frame(ME)
  rownames(ME)=ME[,1]
  ME=ME[,-1]
  SNP_GE=fread(SNP_GE_file,header=T,sep="\t")
  SNP_GE=as.data.frame(SNP_GE)
  rownames(SNP_GE)=SNP_GE[,1]
  SNP_GE=SNP_GE[,-1]
  covariate_ge=read_qtl_input(covariate_ge_file)
  rownames(covariate_ge)=covariate_ge[,1]
  covariate_ge=covariate_ge[,-1]
  covariate_me=read_qtl_input(covariate_me_file)
  rownames(covariate_me)=covariate_me[,1]
  covariate_me=covariate_me[,-1]
  rownames(covariate_me)=paste0("me_",rownames(covariate_me))
  comsamples=intersect(colnames(GE),colnames(ME))
  idx=match(comsamples,colnames(GE))
  GE2=GE[,idx]
  idx=match(comsamples,colnames(ME))
  ME2=ME[,idx]
  idx=match(comsamples,colnames(SNP_GE))
  SNP2=SNP_GE[,idx]
  idx=match(comsamples,colnames(covariate_ge))
  covariate2=t(covariate_ge[,idx])
  idx=match(comsamples,colnames(covariate_me))
  covariate2=cbind.data.frame(covariate2,t(covariate_me[c(1,5:nrow(covariate_me)),idx])) #keep age
  
  tmp=merge(dat1,dat2,by="snp_idx")
  tmp$gene.x=as.character(tmp$gene.x)
  tmp$gene.y=as.character(tmp$gene.y)
  tmp$gene.x=as.character(tmp$gene.x)
  tmp$gene.y=as.character(tmp$gene.y)
  tmp$genename.x=as.character(tmp$genename.x)
  tmp$genename.y=as.character(tmp$genename.y)
  tmp=tmp[tmp$gene.x!=tmp$gene.y,]
  tmp=tmp[sum(is.na(SNP2[tmp$snp_idx,]))<10,] #not many NA in the snp (highrisk117)
  if (nrow(tmp)>0)
  {
    citresults=vector('list', nrow(tmp))
    for (i in 1:nrow(tmp))
    {
      if (i %% 100==0) cat(i,"..")
      #L
      g=unlist(SNP2[tmp$snp_idx[i],])
      if (grepl("^cg",tmp$gene.y[1])) #trait is methylation
      {
        #G
        idx=match(tmp$gene.y[i],rownames(ME2))
        x=unlist(ME2[idx,])
        #T
        idx=match(tmp$gene.x[i],rownames(GE2))
        y=unlist(GE2[idx,])
      }else #reactive
      {
        #G
        idx=match(tmp$gene.y[i],rownames(GE2))
        x=unlist(GE2[idx,])
        #T
        idx=match(tmp$gene.x[i],rownames(ME2))
        y=unlist(ME2[idx,])
        # fit1=glm(x~y+g)
        # summary(fit1)
        # fit2=glm(y~x+g)
        # summary(fit2)
        
      }
      citresults[[ i ]]=cit.cp(g,x,y,covariate2,n.perm=100,n.resampl=100,rseed =10000)
    }
    fdrresults=fdr.cit(citresults)
    fdrresults=cbind.data.frame(tmp,fdrresults)
  }else
  {
    fdrresults=NA
  }
  
  return(fdrresults)
}

computep_meds_inputfromfiles=function(meds=validtriplet_HUTCH_ciseqtl_cismqtl, #should use snp_idx here
                                      opt="meth2gene",
                       phenotypemap=hutch_ge_anno,
                       GE_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_GE.txt",
                       ME_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_ME.txt",
                       SNP_GE_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_highrisk_SNP_GE_updatename.txt",
                       covariate_ge_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_GE_PEER_used.txt",
                       covariate_me_file="/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/TCGA_tumors_COVA_ME_used.txt")
{
  library(cit)
  GE=fread(GE_file,header=T,sep="\t")
  GE=as.data.frame(GE)
  rownames(GE)=GE[,1]
  GE=GE[,-1]
  ME=fread(ME_file,header=T,sep="\t")
  ME=as.data.frame(ME)
  rownames(ME)=ME[,1]
  ME=ME[,-1]
  SNP_GE=fread(SNP_GE_file,header=T,sep="\t")
  SNP_GE=as.data.frame(SNP_GE)
  rownames(SNP_GE)=SNP_GE[,1]
  SNP_GE=SNP_GE[,-1]
  covariate_ge=read_qtl_input(covariate_ge_file)
  rownames(covariate_ge)=covariate_ge[,1]
  covariate_ge=covariate_ge[,-1]
  covariate_me=read_qtl_input(covariate_me_file)
  rownames(covariate_me)=covariate_me[,1]
  covariate_me=covariate_me[,-1]
  rownames(covariate_me)=paste0("me_",rownames(covariate_me))
  comsamples=intersect(colnames(GE),colnames(ME))
  idx=match(comsamples,colnames(GE))
  GE2=GE[,idx]
  idx=match(comsamples,colnames(ME))
  ME2=ME[,idx]
  idx=match(comsamples,colnames(SNP_GE))
  SNP2=SNP_GE[,idx]
  idx=match(comsamples,colnames(covariate_ge))
  covariate2=t(covariate_ge[,idx])
  idx=match(comsamples,colnames(covariate_me))
  covariate2=cbind.data.frame(covariate2,t(covariate_me[c(1,5:nrow(covariate_me)),idx])) #keep age
  
  validres=NULL
  for (i in 1:nrow(meds))
  {
    if (i %%100==0) cat(i,"..")
    g=unlist(SNP2[meds$snp_idx[i],])
    idx1=match(meds$gene.x[i],phenotypemap$Probe_Id)
    idx2=match(phenotypemap$TCGA_Probe_Id[idx1],rownames(GE2))
    idx3=match(meds$gene.y[i],rownames(ME2))
    if (opt=="meth2gene")
    {
      y=unlist(GE2[idx2,])
      x=unlist(ME2[idx3,])
    }else
    {
      x=unlist(GE2[idx2,])
      y=unlist(ME2[idx3,])
    }
    validres=rbind(validres,cit.cp(g,x,y,covariate2,n.resampl = 100))
  }
  meds1=cbind(meds,data.frame(matrix(ncol=ncol(validres),nrow=nrow(meds))))
  colnames(meds1)[(ncol(meds)+1):ncol(meds1)]=paste0("valid_",colnames(validres))
  meds1[,(ncol(meds)+1):ncol(meds1)]=validres
  return(meds1)
}

# diff_tumror_normal_eqtl=function(dat1=highrisk_HUTCH_cis_eqtl_all,dat2=highrisk_gtex_cis_eqtl_all)
# {
#   dat1=dat1[!duplicated(dat1$genename),]
#   dat2=dat2[!duplicated(dat2$genename),]
#   tmp=merge(dat1,dat2,by=c("snp_idx","genename"))
#   pvalue=rep(NA,nrow(tmp))
#   for (i in 1:nrow(tmp))
#   {
#     if (tmp$tstat.x[i]*tmp$tstat.y[i]!=0)
#     {
#       tscore=(tmp$beta.x[i]-tmp$beta.y[i])/sqrt((tmp$beta.x[i]/tmp$tstat.x[i])^2+(tmp$beta.y[i]/tmp$tstat.y[i])^2)
#     }else
#     {
#       tscore=0
#     }
#     pvalue[i]=2*pnorm(abs(tscore),lower.tail = F)
#   }
#   return(pvalue)
# }

diff_tumror_normal_eqtl=function(dat1=highrisk_HUTCH_cis_eqtl_all,dat2=highrisk_TCGA_tumors_cis_eqtl_all,dat3=highrisk_gtex_cis_eqtl_all)
{
  #combine analysis of dat1 and dat2 at snp-gene level
  dat1=dat1[order(dat1$snp_idx,dat1$genename),]
  #remove duplicated pairs based on p-value
  dat11=NULL
  for (i in 1:nrow(dat1))
  {
    idx1=which(dat1$snp_idx==dat1$snp_idx[i] & dat1$genename==dat1$genename[i])
    if (length(idx1)==1)
    {
      dat11=rbind.data.frame(dat11,dat1[idx1,])
    }else
    {
      idx2=which.max(dat1$value[idx1])
      idx3=which(dat11$snp_idx==dat1$snp_idx[idx1[idx2]]&dat11$genename==dat1$genename[idx1[idx2]])
      if (length(idx3)==0)
      {
        dat11=rbind(dat11,dat1[idx1[idx2],])
      }
    }
  }
 
  dat2=dat2[order(dat2$snp_idx,dat2$genename),]
  #remove duplicated pairs based on p-value
  dat22=NULL
  for (i in 1:nrow(dat2))
  {
    idx1=which(dat2$snp_idx==dat2$snp_idx[i] & dat2$genename==dat2$genename[i])
    if (length(idx1)==1)
    {
      dat22=rbind.data.frame(dat22,dat2[idx1,])
    }else
    {
      idx2=which.max(dat2$value[idx1])
      idx3=which(dat22$snp_idx==dat2$snp_idx[idx1[idx2]]&dat22$genename==dat2$genename[idx1[idx2]])
      if (length(idx3)==0)
      {
        dat22=rbind(dat22,dat2[idx1[idx2],])
      }
    }
  }
  
  dat3=dat3[order(dat3$snp_idx,dat3$genename),]
  #remove duplicated pairs based on p-value
  dat33=NULL
  for (i in 1:nrow(dat3))
  {
    idx1=which(dat3$snp_idx==dat3$snp_idx[i] & dat3$genename==dat3$genename[i])
    if (length(idx1)==1)
    {
      dat33=rbind.data.frame(dat33,dat3[idx1,])
    }else
    {
      idx2=which.max(dat3$value[idx1])
      idx3=which(dat33$snp_idx==dat3$snp_idx[idx1[idx2]]&dat33$genename==dat3$genename[idx1[idx2]])
      if (length(idx3)==0)
      {
        dat33=rbind(dat33,dat3[idx1[idx2],])
      }
    }
  }
  
  dat=merge(dat11,dat22,by=c("snp_idx","genename"))
  dat=merge(dat,dat33,by=c("snp_idx","genename"))
  dat$var=(dat$beta/dat$tstat)^2
  dat$pvalue=10^-dat$value
  dat$beta2=NA
  dat$var2=NA
  dat$pvalue2=NA
  dat$pvalue_diff=NA
  pvalue=rep(1,nrow(dat))
  for (i in 1:nrow(dat))
  {
    if (dat$tstat.x[i]*dat$tstat.y[i]*dat$tstat[i]!=0)
    {
      v1=(dat$beta.x[i]/dat$tstat.x[i])^2
      v2=(dat$beta.y[i]/dat$tstat.y[i])^2
      w1=1/v1
      w2=1/v2
      beta2=(w1*dat$beta.x[i]+w2*dat$beta.y[i])/(w1+w2)
      dat$beta2[i]=beta2
      var2=1/(w1+w2)
      dat$var2[i]=var2
      tscore2=beta2/sqrt(var2)
      dat$pvalue2[i]=2*pnorm(abs(tscore2),lower.tail = F)
      tscore=(beta2-dat$beta[i])/sqrt(var2+(dat$beta[i]/dat$tstat[i])^2)
      pvalue[i]=2*pnorm(abs(tscore),lower.tail = F)
    }
  }
  idx=order(pvalue)
  pvalue=pvalue[idx]
  dat=dat[idx,]
  dat$pvalue_diff=pvalue
  return(dat)
}

unique_eqtlpairs=function(dat1=highrisk_HUTCH_cis_eqtl,printflag=T)
{
  if ("snp_idx" %in% colnames(dat1) & "genename" %in% colnames(dat1))
  {
    dat2=NULL
    eqtlpairs=unique(dat1[,c("snp_idx","genename")])
    if (nrow(eqtlpairs)>10000) print(paste0("number of pairs:",nrow(eqtlpairs)))
    for (i in 1:nrow(eqtlpairs))
    {
      if (i %% 10000==0) cat(i,"..")
      idx=which(dat1$snp_idx==eqtlpairs$snp_idx[i] & dat1$genename==eqtlpairs$genename[i])
      idx1=which.max(dat1$value[idx])
      dat2=rbind.data.frame(dat2,dat1[idx[idx1],])
    }
  }else
  {
    dat2=dat1
  }
  print(paste0("number of pairs: ",nrow(dat2)))
  print(paste0("number of unique SNPs: ",length(unique(dat2$snp_idx))))
  print(paste0("number of unique genes: ",length(unique(dat2$genename))))
  if (length(unique(dat2$snp_idx))<length(unique(dat2$genename)) & printflag==T)
  {
    refsnps=read.table("../data/RESUB_Supplementary_Table16_v9.txt",header=T,sep="\t",stringsAsFactors = F)
    refsnps$Chr=gsub(23,"X",refsnps$Chr)
    highrisk_snp=read.table("/fh/fast/stanford_j/Xiaoyu/QTL/result/qtl_input/HUTCH_highrisk_SNP_POS.txt",header=T,stringsAsFactors = F)
    uniqgenes=as.character(unique(dat2$genename))
    uniqsnps=as.character(unique(dat2$snp_idx))
    for (i in 1:length(uniqsnps))
    {
      idx=which(dat2$snp_idx==uniqsnps[i])
      if (length(idx)>1)
      {
        idx3=which(refsnps$Chr==highrisk_snp$chr[as.numeric(uniqsnps[i])] &refsnps$Position==highrisk_snp$pos[as.numeric(uniqsnps[i])])
        print(paste0("snp: ",refsnps$SNP[idx3],",",paste0(dat2$genename[idx],collapse = "|")))
      }
    }
  }
  return(dat2)
}

add_validatedgenes_snp=function(dat=validatepairs_HUTCH_ciseqtl_fwer,dat1=refsnps,colname="Validated_in_eQTL")
{
  dat$genename=as.character(dat$genename)
  dat1$tmp=NA
  for (i in 1:nrow(dat1))
  {
    idx=which(dat$snp==dat1$SNP[i])
    if (length(idx)>0)
    {
      numNA=sum(is.na(dat$genename[idx]))
      if (numNA>0)
      {
        genes=paste0(numNA,"intergenic")
      }else
      {
        genes=NULL
      }
      for (j in 1:length(idx))
      {
        if (!is.na(dat$genename[idx[j]]))
        {
          tmp1=unlist(strsplit(dat$genename[idx[j]],"|",fixed = T)) #for mqtls using | in genename
          genes=c(genes,tmp1)
          genes=unique(genes)
        }
      }
      dat1$tmp[i]=paste0(genes,collapse = ";")
    }
  }
  colnames(dat1)=gsub("tmp",colname,colnames(dat1))
  return(dat1)
}

table_mqtl=function(dat=validatepairs_HUTCH_cismqtl_fwer)
{
  res=data.frame(matrix(0,nrow=13,ncol=3))
  colnames(res)=c("annotation","number","percent")
  res$annotation=c("CpG island","North shore","South shore","North shelf","South shelf","TSS1500","TSS200","5'UTR","1st Exon","Body","3'UTR","Enhancer","DHS")
  dat1=merge(dat,anno,by.x="gene",by.y="IlmnID")
  idx=duplicated(dat1$gene)
  dat1=dat1[!idx,]
  tmp=table(dat1$Relation_to_UCSC_CpG_Island)
  res$number[1]=tmp[which(names(tmp)=="Island")]
  res$number[2]=tmp[which(names(tmp)=="N_Shore")]
  res$number[3]=tmp[which(names(tmp)=="S_Shore")]
  res$number[4]=tmp[which(names(tmp)=="N_Shelf")]
  res$number[5]=tmp[which(names(tmp)=="S_Shelf")]
  res$percent[1:5]=floor(100*res$number[1:5]/nrow(dat1))
  dat1$UCSC_RefGene_Group=as.character(dat1$UCSC_RefGene_Group)
  for (i in 1:nrow(dat1))
  {
    if (grepl("TSS1500",dat1$UCSC_RefGene_Group[i]))
      res$number[6]=res$number[6]+1
    if (grepl("TSS200",dat1$UCSC_RefGene_Group[i]))
      res$number[7]=res$number[7]+1
    if (grepl("5'UTR",dat1$UCSC_RefGene_Group[i]))
      res$number[8]=res$number[8]+1
    if (grepl("1stExon",dat1$UCSC_RefGene_Group[i]))
      res$number[9]=res$number[9]+1
    if (grepl("Body",dat1$UCSC_RefGene_Group[i]))
      res$number[10]=res$number[10]+1
    if (grepl("3'UTR",dat1$UCSC_RefGene_Group[i]))
      res$number[11]=res$number[11]+1
  }
  res$percent[6:11]=floor(100*res$number[6:11]/nrow(dat1))
  res$number[12]=sum(dat1$Enhancer==T,na.rm=T)
  res$percent[12]=floor(100*res$number[12]/nrow(dat1))
  res$number[13]=sum(dat1$DHS==T,na.rm=T)
  res$percent[13]=floor(100*res$number[13]/nrow(dat1))
  return(res)
}

read_qtl_input_matrix=function(filename)
{
  dat=fread(filename,header=T,sep="\t",stringsAsFactors = F)
  dat=as.data.frame(dat)
  colnames(dat)=gsub("^X","",colnames(dat))
  colnames(dat)=gsub(".","-",colnames(dat),fixed = T)
  rownames(dat)=dat[,1]
  dat=dat[,-1]
  return(dat)
}