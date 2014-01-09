#!/usr/local/bin/bash

#script to extract sites that have a removed rsID from dbSNP, but that we call in sequences
dbsnp_file_path=$1
current_chr_file_path=$2
output_file_path=$3
chr=$4

#extract sites from seq data
fgrep -w -f <(zcat $dbsnp_file_path | cut -f 1,2 | grep -v CHROM | grep "^${chr}	" ) $current_chr_file_path > $output_file_path.map

#we need to extract also the rs id for our files...so we made also a reverse extraction: extract sites with rsID from dbSNP file
#fgrep -w -f <(cut -f 1,2 $current_chr_file_path | grep -v CHROM) <(zcat $dbsnp_file_path | cut -f 1-5 | grep -v CHROM | grep "^${chr}	") > ${output_file_path}.rsID
fgrep -w -f <(cut -f 1,2 $current_chr_file_path | grep -v CHROM) <(zcat $dbsnp_file_path | grep -v CHROM | grep "^${chr}	") > ${output_file_path}.rsID

