#!/usr/bin/env bash

#generate plink files for each chromosomes for TBD data
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink

basefilename="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation_tbd/plink/genotypes"
chr=${1?"chr"}
do_split ()
{
  chr=$1
  echo $chr
  output="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation_tbd/plink/genotypes_chr$chr"
  $plink --noweb --bfile $basefilename --chr $chr --make-bed --out $output
  $plink --noweb --bfile $basefilename --chr $chr --recode --tab --out $output #generate map file
}

for ((i=1;i<=22;i++))
do
  do_split $i
done
do_split "X"
