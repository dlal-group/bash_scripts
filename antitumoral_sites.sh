#!/usr/local/bin/bash

#script for antitumoral drugs sites extraction
#First arg is the file we're working on
file=$1
outpath=$2

echo ${file} 
#extract sites by region using antitumoral list
# chr=`echo ${file} | cut -f 1 -d " "`
# start=`echo ${file} | cut -f 2 -d " "`
# end=`echo ${file} | cut -f 3 -d " "`
# gene=`echo ${file} | cut -f 4 -d " "`
# # plink --noweb --bfile /nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/merged/chr${chr}_merged --chr ${chr} --from-bp ${start} --to-bp ${end} --make-bed --out ${outpath}/chr${chr}_${gene}

# #now calculate also the frequencies in FVG
# # plink --noweb --bfile /nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/merged/chr${chr}_merged --chr ${chr} --from-bp ${start} --to-bp ${end} --freq --out ${outpath}/chr${chr}_${gene}_fvgfrq

# #now calculate also the frequencies in FVG for males only
# plink --noweb --bfile /nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/merged/chr${chr}_merged --chr ${chr} --from-bp ${start} --to-bp ${end} --filter-males --freq --out ${outpath}/males/chr${chr}_${gene}_fvgfrq

# #now calculate also the frequencies in FVG for females only
# plink --noweb --bfile /nfs/users/nfs_m/mc14/Work/SANGER/FVG/ANTI_TUMORAL_DRUGS/merged/chr${chr}_merged --chr ${chr} --from-bp ${start} --to-bp ${end} --filter-females --freq --out ${outpath}/females/chr${chr}_${gene}_fvgfrq
