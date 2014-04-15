#!/usr/local/bin/bash

#script for antitumoral drugs sites extraction
#First arg is the file we're working on
file=`sed -n "${LSB_JOBINDEX}p" $1`
# chr=$1
outpath=$2

mkdir -p ${outpath}/GENO
mkdir -p ${outpath}/FREQ/males
mkdir -p ${outpath}/FREQ/females

echo ${file} 
#extract sites by region using antitumoral list
chr=`echo ${file} | cut -f 1 -d " "`
start=`echo ${file} | cut -f 2 -d " "`
end=`echo ${file} | cut -f 3 -d " "`
gene=`echo ${file} | cut -f 4 -d " "`

#extract genotypes
plink --noweb --bfile /lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/ARRAY/fvg_merged/fvg_merged_flipped_clean --chr ${chr} --from-bp ${start} --to-bp ${end} --make-bed --out ${outpath}/GENO/chr${chr}_${gene}

#now calculate also the frequencies in FVG
plink --noweb --bfile /lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/ARRAY/fvg_merged/fvg_merged_flipped_clean --chr ${chr} --from-bp ${start} --to-bp ${end} --freq --out ${outpath}/FREQ/chr${chr}_${gene}_fvgfrq

#now calculate also the frequencies in FVG for males only
plink --noweb --bfile /lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/ARRAY/fvg_merged/fvg_merged_flipped_clean --chr ${chr} --from-bp ${start} --to-bp ${end} --filter-males --freq --out ${outpath}/FREQ/males/chr${chr}_${gene}_fvgfrq

#now calculate also the frequencies in FVG for females only
plink --noweb --bfile /lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/ARRAY/fvg_merged/fvg_merged_flipped_clean --chr ${chr} --from-bp ${start} --to-bp ${end} --filter-females --freq --out ${outpath}/FREQ/females/chr${chr}_${gene}_fvgfrq

