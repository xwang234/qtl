#!/usr/bin/env bash
#SBATCH -t 1-10
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

#set -e
set -u
set -o pipefail

#reference version
imp=2014OCT

plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink
impute=/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/impute2
shapeit=/fh/fast/stanford_j/Xiaoyu/Tools/shapeit/bin/shapeit
#imputation1: 2012 MAR, imputation2:2014 OCT
impfolder=/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation4
outfolder=${impfolder}/plink
#ped,map,flip files:
infolder=/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation4/plink

chr=X
#prefix
inf=$infolder/TCGAnormals_chr$chr
f=$outfolder/TCGAnormals_chr$chr
#generate plink file
$plink --noweb --file $inf --flip $inf.fliplist --recode --out ${f}_flip
$plink --noweb --file ${f}_flip --make-bed --out ${f}_flip


do_shapeitX () {
 # parameters
 local THREAT=1
 local EFFECTIVESIZE=11418
 local B=$1 
 local GENMAP_FILE=$2
 local OUTPUT_HAPS=$3
 local OUTPUT_SAMPLE=$4
 local OUTPUT_LOG=$5
 $shapeit -B $B --input-map $GENMAP_FILE --thread $THREAT --effective-size $EFFECTIVESIZE --output-max $OUTPUT_HAPS $OUTPUT_SAMPLE --output-log $OUTPUT_LOG --states 200 --chrX
}

do_impute2X () {

  NE=20000
  #output of prephasing  
  #GWAS_HAPS_FILE=$OUTPUT_HAPS
  GWAS_HAPS_FILE=$1
  
  #ALL_OUTPUTimp="${impfolder}/SNP6_imp_chr${chr}.txt"
  #ALL_OUTPUTinfo="${impfolder}/SNP6_info_chr${chr}.txt"
  GENMAP_FILE=$2
  HAPS_FILE=$3
  LEGEND_FILE=$4
  chunk_OUTPUTimp=$5 #"${impfolder}/chunkresult/SNP6_imp_chr${chr}_chunk"
  ALL_OUTPUTimp=$6
  ALL_OUTPUTinfo=$7
  samplefile=$8

  lastpos1=$(cat $inf.map |tail -n 1 |awk '{print $4}')
  lastpos2=$(cat $GENMAP_FILE |tail -n 1 |awk '{print $1}')
  lastpos=0
  if [[ $lastpos1 -le $lastpos2 ]];then lastpos=$lastpos1;else lastpos=$lastpos2;fi
  ((lastpos=lastpos+500000))
  firstpos1=$(cat $inf.map |head -n 2 |tail -n 1|awk '{print $4}')
  firstpos2=$(cat $GENMAP_FILE |head -n 2 |tail -n 1|awk '{print $1}')
  firstpos=0
  if [[ $firstpos1 -le $firstpos2 ]];then firstpos=$firstpos2;else firstpos=$firstpos1;fi
  ((firstpos=firstpos-500000))
  if [[ $firstpos -lt 1 ]];then firstpos=1;fi
  interval=$((lastpos-$firstpos))  
  nchunks=$(( ($interval+5000000-1)/5000000 ))

  if [[ -f $ALL_OUTPUTimp ]];then rm $ALL_OUTPUTimp;fi
  if [[ -f $ALL_OUTPUTinfo ]];then rm $ALL_OUTPUTinfo;fi
  echo "number of chunks: $nchunks"
  for ((chunk=1;chunk<=nchunks;chunk++)) 
  do
    (( CHUNK_START = firstpos+(chunk-1)*5000000))
    (( CHUNK_END = firstpos+chunk*5000000-1))
    OUTPUT_FILE=$chunk_OUTPUTimp"${chunk}.txt"
    $impute -chrX\
       -use_prephased_g\
       -m $GENMAP_FILE \
       -known_haps_g $GWAS_HAPS_FILE \
       -h $HAPS_FILE \
       -l $LEGEND_FILE \
       -pgs_miss \
       -buffer 500 \
       -filt_rules_l 'filter==0' \
       -Ne $NE \
       -int $CHUNK_START $CHUNK_END \
       -o $OUTPUT_FILE \
       -allow_large_regions \
       -sample_g $samplefile\
       -phase

    cat $OUTPUT_FILE >> $ALL_OUTPUTimp
    if [[ $chunk -eq 1 ]]
    then
      cat ${OUTPUT_FILE}_info >> $ALL_OUTPUTinfo
    else
      awk '{if (NR>1) print}' ${OUTPUT_FILE}_info >> $ALL_OUTPUTinfo #not include the header again
    fi
  done
}


# reference data files

#2014 OCT
GENMAP_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/genetic_map_chr${chr}_PAR1_combined_b37.txt"
HAPS_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr${chr}_PAR1.hap.gz"
LEGEND_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr${chr}_PAR1.legend_filter1.gz" #use EUR and EAS

OUTPUT_HAPS=${outfolder}/TCGAnormals_chr${chr}_PAR1.haps
OUTPUT_SAMPLE=${outfolder}/TCGAnormals_chr${chr}_PAR1.sample
OUTPUT_LOG=${outfolder}/TCGAnormals_chr${chr}_PAR1.log
#do_shapeitX ${f}_flip $GENMAP_FILE $OUTPUT_HAPS $OUTPUT_SAMPLE $OUTPUT_LOG
samplefile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation4/SNP6_allsamplefile.txt"
do_impute2X $OUTPUT_HAPS $GENMAP_FILE $HAPS_FILE $LEGEND_FILE "${impfolder}/chunkresult/SNP6_imp_chr${chr}_PAR1_chunk" "${impfolder}/SNP6_imp_chr${chr}_PAR1.txt" "${impfolder}/SNP6_info_chr${chr}_PAR1.txt" $samplefile

GENMAP_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/genetic_map_chr${chr}_PAR2_combined_b37.txt"
HAPS_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr${chr}_PAR2.hap.gz"
LEGEND_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr${chr}_PAR2.legend_filter1.gz"

#OUTPUT_HAPS=${outfolder}/TCGAnormals_chr${chr}_PAR2.haps
OUTPUT_HAPS=${outfolder}/TCGAnormals_chr${chr}_PAR1.haps
OUTPUT_SAMPLE=${outfolder}/TCGAnormals_chr${chr}_PAR2.sample
OUTPUT_LOG=${outfolder}/TCGAnormals_chr${chr}_PAR2.log
#All three phasing results are the same
#do_shapeitX ${f}_flip $GENMAP_FILE $OUTPUT_HAPS $OUTPUT_SAMPLE $OUTPUT_LOG
do_impute2X $OUTPUT_HAPS $GENMAP_FILE $HAPS_FILE $LEGEND_FILE "${impfolder}/chunkresult/SNP6_imp_chr${chr}_PAR2_chunk" "${impfolder}/SNP6_imp_chr${chr}_PAR2.txt" "${impfolder}/SNP6_info_chr${chr}_PAR2.txt" $samplefile

GENMAP_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/genetic_map_chr${chr}_nonPAR_combined_b37.txt"
HAPS_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr${chr}_NONPAR.hap.gz"
LEGEND_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr${chr}_NONPAR.legend_filter1.gz"

#OUTPUT_HAPS=${outfolder}/TCGAnormals_chr${chr}_nonPAR.haps
OUTPUT_HAPS=${outfolder}/TCGAnormals_chr${chr}_PAR1.haps
OUTPUT_SAMPLE=${outfolder}/TCGAnormals_chr${chr}_nonPAR.sample
OUTPUT_LOG=${outfolder}/TCGAnormals_chr${chr}_nonPAR.log
#do_shapeitX ${f}_flip $GENMAP_FILE $OUTPUT_HAPS $OUTPUT_SAMPLE $OUTPUT_LOG
do_impute2X $OUTPUT_HAPS $GENMAP_FILE $HAPS_FILE $LEGEND_FILE "${impfolder}/chunkresult/SNP6_imp_chr${chr}_nonPAR_chunk" "${impfolder}/SNP6_imp_chr${chr}_nonPAR.txt" "${impfolder}/SNP6_info_chr${chr}_nonPAR.txt" $samplefile


