#!/usr/bin/python
#SBATCH -t 0-10
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

#https://genome.sph.umich.edu/wiki/IMPUTE2:_1000_Genomes_Imputation_Cookbook#Pre-Phasing_using_IMPUTE2
import subprocess
from subprocess import call
import math
import sys

impute="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/impute2"
impfolder="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation2"
#i: chr 1:22
i=int(sys.argv[1])
print i
#check the length of chr
#legendfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr"+str(i)+".legend.gz"
#mapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/genetic_map_chr"+str(i)+"_combined_b37.txt"
#hapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr"+str(i)+".hap.gz"
legendfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr"+str(i)+"_impute.legend.gz"
mapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr"+str(i)+"_combined_b37.txt"
hapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr"+str(i)+"_impute.hap.gz"
genotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_gtdata_chr"+str(i)+".txt"
strandfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_strand_chr"+str(i)+".txt"

zcatcmd=['zcat',legendfile]
lastlinecmd=["tail","-n 1"]
zcat=subprocess.Popen(zcatcmd, stdout=subprocess.PIPE)
lastline=subprocess.check_output(lastlinecmd,stdin=zcat.stdout)
zcat.stdout.close()
chrlen=int(lastline.split(" ")[1])
#split chrlen into chunnks of 5MB
nchunks=int(math.ceil(chrlen/5000000.0))

outfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_imp_chr"+str(i)+".txt"
outfiles=[]
outfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_info_chr"+str(i)+".txt"
outfiles1=[]
for chunk in range(1,nchunks+1):
# Step 1: Pre-phasing 
	preout="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_preout_chr"+str(i)+"_chunk"+str(chunk)+".txt"
	start=1+(chunk-1)*5000000
	end=chunk*5000000
	precmd=[impute,"-prephase_g","-m",mapfile,"-g",genotypefile,"-int",str(start),str(end),"-Ne","20000", "-o",preout]
	#call(precmd)
#Step 2: Imputation into pre-phased haplotypes
	impout="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_imp_chr"+str(i)+"_chunk"+str(chunk)+".txt"
	impcmd=[impute,"-use_prephased_g","-m",mapfile,"-h",hapfile,"-l",legendfile,"-known_haps_g",preout+"_haps","-strand_g",strandfile,"-int",str(start),str(end),"-Ne","20000", "-o",impout,"-phase"]
	#call(impcmd)
	outfiles.append(impout)
	infofile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/chunkresult/"+"SNP6_imp_chr"+str(i)+"_chunk"+str(chunk)+".txt_info"
	outfiles1.append(infofile)

catcmd=['cat']+outfiles
with open(outfile,"w") as out:
	call(catcmd,stdout=out)

catcmd=['cat']+outfiles1
with open(outfile1,"w") as out1:
	call(catcmd,stdout=out1)	



	

