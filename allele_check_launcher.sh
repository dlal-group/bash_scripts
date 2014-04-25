#!/usr/local/bin/bash

#script for allele check in overlapping novel sites
#take as arguments:
#$1: first chr map file path
#$2: second chr map file path

for i in {1..22}
do
first_file=`ls $1/*chr${i}.*map`
second_file=`ls $2/*chr${i}.*map`
bsub -J "allele_check_chr${i}" -o "%J_allele_check_chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement \
R CMD BATCH "--args $first_file $second_file ${i}" /nfs/users/nfs_m/mc14/Work/r_scripts/alleles_compare.r
done
