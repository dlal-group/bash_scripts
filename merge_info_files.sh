#!/usr/local/bin/bash

#script to merge info files from separate imputation

#Args:
#$1=first imp path
#$2=second imp path
#$3= out path
first_imp_path=$1
sec_imp_path=$2
outpath=$3

#merge each info file from each chunk in order to create a merged chunk file
#we need to use the .gen file to know how many individuals we have for each imputation
#we assume that each chunk is named like: chr#.##.gen.gz
#use a s first population the one with more sites!
imp1_ind=`zcat ${first_imp_path}/chr${chr}.${chunk}.gen.gz | awk '{print (NF-5)/3}'|uniq`
imp2_ind=`zcat ${sec_imp_path}/chr${chr}.${chunk}.gen.gz | awk '{print (NF-5)/3}'|uniq`

#extract sites with the same position
awk 'NR==FNR{a[$2,$3,$4]=$0;next;}(a[$2,$3,$4])'

