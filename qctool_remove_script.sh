#!/usr/bin/env bash

#inner command for the script to convert impute2 output files to bimbam comma sepatarated mean genotype files
#args:
#$1=chr number
#$2=chromosome files geno path
#$3=transformed file output
#Modified to be run as a job array
pop=$1
chr=$2
basefolder=$3

echo "Processing gen file for ${pop} ${chr}.."
tail -n+2 /home/cocca/imputation/31012017_MERGED_TEST/${pop}/chr${chr}_snps_to_filter.list|cut -f 3  >/home/cocca/imputation/31012017_MERGED_TEST/${pop}/chr${chr}_snpsIDS_to_filter.list
qctool -g ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen.gz -excl-rsids /home/cocca/imputation/31012017_MERGED_TEST/${pop}/chr${chr}_snpsIDS_to_filter.list -omit-chromosome -og ${basefolder}/${pop}/MERGED/ALL/FILTERED/chr${chr}.gen.gz 

echo "Processing gen_info file for ${pop} ${chr}.."
fgrep -v -w -f /home/cocca/imputation/31012017_MERGED_TEST/${pop}/chr${chr}_snpsIDS_to_filter.list ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen_info > ${basefolder}/${pop}/MERGED/ALL/FILTERED/chr${chr}.gen_info

echo "moving unfiltered data for ${pop} ${chr}.."
mv ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen.gz ${basefolder}/${pop}/MERGED/ALL/UNFILTERED/chr${chr}.gen.gz
mv ${basefolder}/${pop}/MERGED/ALL/chr${chr}.gen_info ${basefolder}/${pop}/MERGED/ALL/UNFILTERED/chr${chr}.gen_info
