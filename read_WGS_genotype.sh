#!/usr/bin/env bash
wgsfile="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/59597/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_AllVar_QC_metrics.vcf.gz"
highriskfile="/fh/fast/stanford_j/Xiaoyu/QTL/data/RESUB_Supplementary_Table16_v9.txt"
outfile="../data/GTEx/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1.GRU_highrisk.vcf"
highrisksnps=($(awk '{if (NR>1) print $1}' $highriskfile))
chr=($(awk '{if (NR>1) print $2}' $highriskfile))
pos=($(awk '{if (NR>1) print $3}' $highriskfile))
if [[ -f $outfile ]]; then rm $outfile;fi
for ((i=0;i<${#chr[@]};i++))
do
  if [[ ${chr[$i]} -eq 23 ]]; then chr[$i]="X";fi
  tabix $wgsfile ${chr[$i]}:${pos[$i]}-${pos[$i]} >>$outfile
done
#tabix $wgsfile 1:10109-10140 
