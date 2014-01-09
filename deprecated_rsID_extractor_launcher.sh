#!/usr/local/bin/bash

#script to launch deprecated_rsID_extractor.sh
#Non serve a un chezzo, scritto così!!!
#Args to pass:
#dbsnp_file_path=$1
#current_chr_file_path=$2
#output_file_path=$3
#chr=$4

#for i in {1..22}
#do
#extract sites from seq data
#fgrep -w -f <(zcat $dbsnp_file_path | cut -f 1,2 | grep -v CHROM | grep "^${chr}        " ) $current_chr_file_path > $output_file_path.map

#bsub -J "deprecated_rsID_extraction_chr${i}" -o "%J_deprecated_rsID_extraction_chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
#-q normal deprecated_rsID_extractor.sh $dbsnp_file_path $current_chr_file_path $output_file_path ${i}

#done
