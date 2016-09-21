#!/usr/bin/env bash
#
# Script to clean and generate the last final vcf file for the ref panel

vcf=$1
cohort=$2
all_list=$3
outdir=$4
mode=$5

all_list_name=`basename ${all_list}`

bcftools view -h ${vcf} > ${vcf}.header
bcftools view -H ${vcf} | python /nfs/users/nfs_m/mc14/Work/bash_scripts/last_filter.py ${cohort} ${all_list} ${outdir} ${mode}

cat ${vcf}.header ${outdir}/${all_list_name}.${cohort}.${mode}.to_keep.vcf | bgzip -c > ${outdir}/${all_list_name}.${cohort}.${mode}.to_keep.vcf.gz
tabix -f -p vcf ${outdir}/${all_list_name}.${cohort}.${mode}.to_keep.vcf.gz

#clean up a little
rm ${vcf}.header
rm ${outdir}/${all_list_name}.${cohort}.${mode}.to_keep.vcf