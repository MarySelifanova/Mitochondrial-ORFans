#!/bin/bash

for file in *bed
do
	gawk -F'\t' 'BEGIN{OFS="\t";}{print $1,$2,$3,$4,$5,$6}' $file > ~/dora/annot/$file
done
