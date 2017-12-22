#!/usr/bin/env bash

for i in {1..22}
do
	#echo "sbatch /fh/fast/stanford_j/Xiaoyu/QTL/code/impute_shapeit_TCGAnormals.sh $i"	
	#sbatch /fh/fast/stanford_j/Xiaoyu/QTL/code/impute_shapeit_TCGAnormals.sh $i
	#echo "sbatch /fh/fast/stanford_j/Xiaoyu/QTL/code/impute_shapeit_TCGAallnormals.sh $i"	
	#sbatch /fh/fast/stanford_j/Xiaoyu/QTL/code/impute_shapeit_TCGAallnormals.sh $i
	echo "sbatch /fh/fast/stanford_j/Xiaoyu/QTL/code/impute_shapeit_TCGAalltumors.sh $i"	
	sbatch /fh/fast/stanford_j/Xiaoyu/QTL/code/impute_shapeit_TCGAalltumors.sh $i
	sleep 1s
done

