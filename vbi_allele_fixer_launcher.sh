#!/usr/local/bin/bash

#script for retrieve ALT/REF information from other seq file

#first: generate a list of sites without the ALT allele (removed in v2 annotation by yasin)

#take as arguments:
#$1: chr file with missing alleles path 
#$2: fixed alleles map file path
#$3: output path

if [ $# -lt 2 ]
then
        echo "ATTENTION!!!Missing arguments!!"
        echo "Usage:"
        echo "vbi_allele_fixer_launcher.sh <chr_file_path_with_missing_alleles> <fixed allele map file path> <output_path>"
        exit 1
fi

for i in {1..22}
do
missing_chr_file=`ls $1/*chr${i}.*.map`
fixing_chr_file=`ls $2/fixed_*chr${i}.map`

bsub -J "vbi_allele_fixer_chr${i}" -o "%J_vbi_allele_fixer_chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement vbi_allele_fixer.sh ${i} $missing_chr_file $fixing_chr_file $3
done


