#!/usr/local/bin/bash

#script to launch the the putative novel site list creator

#take as arguments:
#$1: path of files used to filter 
#$2: path of files to be filtered
#$3: output path

if [ $# -lt 3 ]
then
        echo "ATTENTION!!!Missing arguments!!"
        echo "Usage:"
        echo "putative_novel_list_extractor_launcher.sh <path of files used to filter> <path of files to be filtered> <output_path>"
        exit 1
fi

for i in {1..22}
do
filter_chr_file=`ls $1/*chr${i}.*.map`
tb_filtered_chr_file=`ls $2/*chr${i}.*.gz`

bsub -J "putative_novel_chr${i}" -o "%J_putative_novel_chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement putative_novel_list_extractor.sh ${i} $filter_chr_file $tb_filtered_chr_file $3
done


