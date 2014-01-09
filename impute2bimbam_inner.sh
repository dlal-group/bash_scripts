#!/usr/local/bin/bash

#inner command for the script to convert impute2 output files to bimbam comma sepatarated mean genotype files
#args:
#$1=chr number
#$2=chromosome files geno path
#$3=transformed file output
#Modified to be run as a job array
file=$1
chr=$2
outpath=$3
file_name=`basename ${file}`
##################################################################################################################################################
#ATTENTION: to generate the bimbam file we check if the snp id column matches the pattern "rs*" if this is not the case, we set the rsID to "NA"!!!
##################################################################################################################################################
#This generate bimbam with the allele info in the rsID field
# zcat ${file} | sed 's/,rs/|rs/g' | awk -v chr=${chr} '{ snp=(NF-5)/3; if($2 ~/^rs/) s=$2;else s="NA"; printf "chr"chr":"$3"-"s"-"$4"-"$5"," $4 "," $5; for(i=1; i<=snp; i++) printf "," $(i*3+3)*2+$(i*3+4); printf "\n" }' > ${outpath}/${file_name}.bimbam

# this is to generate bimbam without allele information in the rsID field
zcat ${file} | awk -v chr=${chr} '{ snp=(NF-5)/3; printf "chr"chr":"$3","$4","$5; for(i=1; i<=snp; i++) printf "," $(i*3+3)*2+$(i*3+4); printf "\n" }' > ${outpath}/${file_name}.bimbam
# zcat $file | awk '{s=(NF-5)/3;gsub(/,/, "_", $2);printf $2 "," $4 "," $5; for(i=1; i<=s; i++) printf "," $(i*3+3)*2+$(i*3+4); printf "\n"}' | gzip > $outpath/$file_name

#command to create position annotation files for gemma
# chrom=${file_name%%.*}
# chrom_num=${chrom#chr*}

# zcat $file | awk -v chrom=$chr '{gsub(/,/, "_", $2);printf $2 "," $3 "," chrom "\n"}' > $outpath/$file_name.pos
zcat $file | awk -v chrom=$chr '{printf "chr"chr":"$3" "," $3 "," chrom "\n"}' > $outpath/$file_name.pos


