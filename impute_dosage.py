#!/usr/bin/python
#SBATCH -t 0-10
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

#compute dosage for the imputed data
import numpy as np
import pandas as pd

#infolder=outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation"
#infolder=outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation1"
infolder=outfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation3"

def todosage(impfile,outfile):
    data=pd.read_csv(impfile,delimiter=" ",header=None)
    nsample=(len(data.columns)-5)/3
    outdata=pd.DataFrame(index=data.index)
    for j in range(1,nsample+1):
        nc=(j-1)*3+1+4
        tmp=data.iloc[:,nc+1]+data.iloc[:,nc+2]*2
        idxna=(data.iloc[:,nc]==0)&(data.iloc[:,nc+1]==0)&(data.iloc[:,nc+2]==0)
        tmp[idxna]=None
        outdata=pd.concat([outdata,tmp],axis=1)
    np.savetxt(outfile,outdata.values,fmt="%.5g",delimiter="\t")

for i in range(1,23):
    print str(i)+".."
    impfile=infolder+"/"+"SNP6_imp_chr"+str(i)+".txt"
    outfile=outfolder+"/"+"SNP6_imputed_dosages_chr"+str(i)+".txt"
    todosage(impfile,outfile)
    
#impfile=infolder+"/"+"SNP6_imp_chrX_par1.txt"
#outfile=outfolder+"/"+"SNP6_imputed_dosages_chrX_par1.txt"
#todosage(impfile,outfile)
#impfile=infolder+"/"+"SNP6_imp_chrX_par2.txt"
#outfile=outfolder+"/"+"SNP6_imputed_dosages_chrX_par2.txt"
#todosage(impfile,outfile)
#impfile=infolder+"/"+"SNP6_imp_chrX_nonpar.txt"
#outfile=outfolder+"/"+"SNP6_imputed_dosages_chrX_nonpar.txt"
#todosage(impfile,outfile)

impfile=infolder+"/"+"SNP6_imp_chrX_PAR1.txt"
outfile=outfolder+"/"+"SNP6_imputed_dosages_chrX_PAR1.txt"
todosage(impfile,outfile)
impfile=infolder+"/"+"SNP6_imp_chrX_PAR2.txt"
outfile=outfolder+"/"+"SNP6_imputed_dosages_chrX_PAR2.txt"
todosage(impfile,outfile)
impfile=infolder+"/"+"SNP6_imp_chrX_nonPAR.txt"
outfile=outfolder+"/"+"SNP6_imputed_dosages_chrX_nonPAR.txt"
todosage(impfile,outfile)
        