#!/usr/local/bin/bash

#script for retrieve ALT/REF information from other seq file

#first: generate a list of sites without the ALT allele (removed in v2 annotation by yasin)

#take as arguments:
#$1: chr file path
#$2: reference annotation file path 

if [ $# -lt 2 ]
then
	echo "ATTENTION!!!Missing arguments!!"
	echo "Usage:"
	echo "vbi_allele_retriever_launcher.sh <chr_file_path> <ref_annotation_file_path>"
	exit 1
fi

for i in {1..22}
do
chr_file=`ls $1/*chr${i}.*.gz`
bsub -J "vbi_allele_retriever_chr${i}" -o "%J_vbi_allele_retriever_chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement vbi_allele_retriever.sh ${i} $chr_file $1 $2
done

