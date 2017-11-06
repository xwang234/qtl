#!/usr/bin/env bash
#set -e
set -u
set -o pipefail

#reference version
imp=2014OCT

plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink
impute=/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/impute2
shapeit=/fh/fast/stanford_j/Xiaoyu/Tools/shapeit/bin/shapeit
#imputation1: 2012 MAR, imputation2:2014 OCT
if [[ $imp == "2012MAR" ]]; then impfolder=/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation1;fi
if [[ $imp == "2014OCT" ]]; then impfolder=/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation2;fi
outfolder=${impfolder}/plink
#ped,map,flip files:
infolder=/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation1/plink

# parameters
THREAT=4
EFFECTIVESIZE=11418

#impute2
NE=20000


  chr=${1?"chrnumber"}
  echo $chr
  # reference data files
  #2012 MAR
#  GENMAP_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr${chr}_combined_b37.txt"
#  HAPS_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/ALL_1000G_phase1integrated_v3_chr${chr}_impute.hap.gz"
#  LEGEND_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr${chr}_impute.legend.gz"
#  refsample="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3.sample"
  #2014 OCT
  GENMAP_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt"
  HAPS_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz"
  LEGEND_FILE="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz"
  refsample="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/1000GP_Phase3.sample"


  # phase GWAS genotypes
  ## this doesn't work,src/modes/phaser/phaser_algorithm.cpp:150: void phaser::phaseSegment(int): Assertion `conditional_index[segment].size() >= 2' failed
  ## it needs to include reference --input-ref for small amount of samples 
#  OUTPUT_HAPS=${outfolder}/TCGAnormals_chr${chr}.haps
#  OUTPUT_SAMPLE=${outfolder}/TCGAnormals_chr${chr}.sample
#  OUTPUT_LOG=${outfolder}/TCGAnormals_chr${chr}.log
  
#  $shapeit \
#  --input-bed $GWASDATA_BED $GWASDATA_BIM $GWASDATA_FAM \
#  --input-map $GENMAP_FILE \
#  --thread $THREAT \
#  --effective-size $EFFECTIVESIZE \
#  --output-max $OUTPUT_HAPS $OUTPUT_SAMPLE \
#  --output-log $OUTPUT_LOG

  #prefix
  inf=$infolder/TCGAnormals_chr$chr
  f=$outfolder/TCGAnormals_chr$chr
  
onlyrunimp=1
if [[ $onlyrunimp -eq 0 ]]
then
  #make plink files, flip some SNPs on -strand
  $plink --noweb --file $inf --flip $inf.fliplist --recode --out ${f}_flip
  $plink --noweb --file ${f}_flip --make-bed --out ${f}_flip
  
  #prephasing use shapeit
  $shapeit -check -B ${f}_flip --input-map $GENMAP_FILE -R ${HAPS_FILE} ${LEGEND_FILE} $refsample --thread 1 --output-log $f.log
  $shapeit -B ${f}_flip --input-map $GENMAP_FILE -R ${HAPS_FILE} ${LEGEND_FILE} $refsample --thread 1 --exclude-snp $f.snp.strand.exclude -O $f.with.ref
fi

  #imputation
  # main output file
  #not used, strand had been flipped
  #strandfile="/fh/fast/stanford_j/Xiaoyu/QTL/result/imputation/SNP6_strand_chr${chr}.txt"
  #output of prephasing  
  GWAS_HAPS_FILE=$f.with.ref.haps
  lastpos=$(cat $inf.map |tail -n 1 |awk '{print $4}')
  ((lastpos=lastpos+100000))
  firstpos=$(cat $inf.map |head -n 2 |tail -n 1|awk '{print $4}')
  ((firstpos=firstpos-100000))

    
  if [[ $firstpos -lt 1 ]];then firstpos=1;fi
  interval=$((lastpos-$firstpos))  
  nchunks=$(( ($interval+5000000-1)/5000000 ))
  ALL_OUTPUT1="${impfolder}/SNP6_imp_chr${chr}.txt"
  ALL_OUTPUT2="${impfolder}/SNP6_info_chr${chr}.txt"
  if [[ -f $ALL_OUTPUT1 ]];then rm $ALL_OUTPUT1;fi
  if [[ -f $ALL_OUTPUT2 ]];then rm $ALL_OUTPUT2;fi
  for ((chunk=1;chunk<=nchunks;chunk++)) 
  do
    (( CHUNK_START = firstpos+(chunk-1)*5000000))
    (( CHUNK_END = firstpos+chunk*5000000-1))
    OUTPUT_FILE="${impfolder}/chunkresult/SNP6_imp_chr${chr}_chunk${chunk}.txt"
    $impute -use_prephased_g\
       -m $GENMAP_FILE \
       -known_haps_g $GWAS_HAPS_FILE \
       -pgs_miss \
       -filt_rules_l "EUR==0" \
       -h $HAPS_FILE \
       -l $LEGEND_FILE \
       -Ne $NE \
       -k_hap 800 -buffer 500 \
       -int $CHUNK_START $CHUNK_END \
       -o $OUTPUT_FILE \
       -allow_large_regions \
       -phase

    cat $OUTPUT_FILE >> $ALL_OUTPUT1
    cat ${OUTPUT_FILE}_info >> $ALL_OUTPUT2
  done




