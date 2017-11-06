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
#work on PAR1
legendfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_PAR1_impute.legend.gz"
mapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/genetic_map_chrX_PAR1_combined_b37.txt"
hapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_PAR1_impute.hap.gz"
genotypefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_gtdata_chrX.txt"
strandfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_strand_chrX.txt"

zcatcmd=['zcat',legendfile]
lastlinecmd=["tail","-n 1"]
zcat=subprocess.Popen(zcatcmd, stdout=subprocess.PIPE)
lastline=subprocess.check_output(lastlinecmd,stdin=zcat.stdout)
zcat.stdout.close()
lastpos=int(lastline.split(" ")[1])
firstlinecmd=["head","-n 2"]
zcat=subprocess.Popen(zcatcmd, stdout=subprocess.PIPE)
firstline=subprocess.check_output(firstlinecmd,stdin=zcat.stdout)
zcat.stdout.close()
firstpos=int(firstline.split(" ")[12])
preout="/fh/fast/stanford_j/Xiaoyu/QTL/result/"+"SNP6_preout_chrX_par1.txt"
precmd=[impute,"-prephase_g","-m",mapfile,"-g",genotypefile,"-int",str(firstpos),str(lastpos),"-Ne","20000", "-o",preout]
call(precmd)
impout="/fh/fast/stanford_j/Xiaoyu/QTL/result/"+"SNP6_imp_chrX_par1.txt"
impcmd=[impute,"-use_prephased_g","-m",mapfile,"-h",hapfile,"-l",legendfile,"-known_haps_g",preout+"_haps","-strand_g",strandfile,"-int",str(firstpos),str(lastpos),"-Ne","20000", "-o",impout,"-phase"]
call(impcmd)

#work on PAR2
legendfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_PAR2_impute.legend.gz"
mapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/genetic_map_chrX_PAR2_combined_b37.txt"
hapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_PAR2_impute.hap.gz"

zcatcmd=['zcat',legendfile]
lastlinecmd=["tail","-n 1"]
zcat=subprocess.Popen(zcatcmd, stdout=subprocess.PIPE)
lastline=subprocess.check_output(lastlinecmd,stdin=zcat.stdout)
zcat.stdout.close()
lastpos=int(lastline.split(" ")[1])
firstlinecmd=["head","-n 2"]
zcat=subprocess.Popen(zcatcmd, stdout=subprocess.PIPE)
firstline=subprocess.check_output(firstlinecmd,stdin=zcat.stdout)
zcat.stdout.close()
firstpos=int(firstline.split(" ")[12])
preout="/fh/fast/stanford_j/Xiaoyu/QTL/result/"+"SNP6_preout_chrX_par2.txt"
precmd=[impute,"-prephase_g","-m",mapfile,"-g",genotypefile,"-int",str(firstpos),str(lastpos),"-Ne","20000", "-o",preout]
call(precmd)
impout="/fh/fast/stanford_j/Xiaoyu/QTL/result/"+"SNP6_imp_chrX_par2.txt"
impcmd=[impute,"-use_prephased_g","-m",mapfile,"-h",hapfile,"-l",legendfile,"-known_haps_g",preout+"_haps","-strand_g",strandfile,"-int",str(firstpos),str(lastpos),"-Ne","20000", "-o",impout,"-phase"]
call(impcmd)

#work on NONPAR
legendfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend.gz"
mapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/genetic_map_chrX_nonPAR_combined_b37.txt"
hapfile="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap.gz"
zcatcmd=['zcat',legendfile]
lastlinecmd=["tail","-n 1"]
zcat=subprocess.Popen(zcatcmd, stdout=subprocess.PIPE)
lastline=subprocess.check_output(lastlinecmd,stdin=zcat.stdout)
zcat.stdout.close()
lastpos=int(lastline.split(" ")[1])
firstlinecmd=["head","-n 2"]
zcat=subprocess.Popen(zcatcmd, stdout=subprocess.PIPE)
firstline=subprocess.check_output(firstlinecmd,stdin=zcat.stdout)
zcat.stdout.close()
firstpos=int(firstline.split(" ")[12])

#split chrlen into chunnks of 5MB
nchunks=int(math.ceil((lastpos-firstpos)/5000000.0))

outfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_imp_chrX_nonpar.txt"
outfiles=[]
outfile1="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_info_chrX_nonpar.txt"
outfiles1=[]
samplefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/SNP6_samplefile.txt"
for chunk in range(1,nchunks+1):
# Step 1: Pre-phasing 
	impout="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_imp_chrX_nonpar"+"_chunk"+str(chunk)+".txt"
	infofile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/"+"SNP6_imp_chrX_nonpar"+"_chunk"+str(chunk)+".txt_info"
	start=firstpos+(chunk-1)*5000000
	end=firstpos+chunk*5000000
	impcmd=[impute,"-chrX","-prephase_g","-m",mapfile,"-h",hapfile,"-g",genotypefile,"-l",legendfile,"-sample_g",samplefile,"-int",str(start),str(end),"-Ne","20000", "-o",impout]
	#call(impcmd)
	outfiles.append(impout)
	outfiles1.append(infofile)

catcmd=['cat']+outfiles
with open(outfile,"w") as out:
	call(catcmd,stdout=out)

catcmd=['cat']+outfiles1
with open(outfile1,"w") as out1:
	call(catcmd,stdout=out1)	



	

