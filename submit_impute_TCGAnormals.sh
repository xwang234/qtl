#!/usr/bin/env bash

for i in {1..22}
do
	echo "sbatch /fh/fast/stanford_j/Xiaoyu/QTL/code/impute_TCGAnormals.py $i"	
	sbatch /fh/fast/stanford_j/Xiaoyu/QTL/code/impute_TCGAnormals.py $i
	sleep 1s
done

