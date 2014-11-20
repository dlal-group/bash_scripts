#!/usr/local/bin/bash

#script to create a info file from imputation with allele information in the last two columns
#Args list:
#$1=chr
#$2=geno path
#$3=output 

chr=$1
geno_path=$2
out_path=$3

(echo "`fgrep certainty ${geno_path}/CHR${chr}/chr${chr}.gen_info` a0 a1"; join -1 1 -2 1 <(zcat ${geno_path}/CHR${chr}/chr${chr}.gen.gz | cut -f 2-5 -d " ") <(fgrep -v certainty ${geno_path}/CHR${chr}/chr${chr}.gen_info | cut -f 2- -d " " ) | awk -v chr=${chr} '{print chr,$1,$2,$7,$8,$9,$10,$11,$12,$13,$14,$15,$3,$4}') > ${out_path}/CHR${chr}/chr${chr}.gen_info_allele
