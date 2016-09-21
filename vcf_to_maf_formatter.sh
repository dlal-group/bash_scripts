#!/usr/local/bin/bash

#script we can use to format vcf files in tab spaced file as:
#CHROM\tPOS\tID\tAC\tAN\tAF\tMAF
#Args:
#$1=list of vcf files provided
#$2=out path
#$3= chr

vcf_file=$1
out_path=$2

echo ${vcf_file}
#we need bcftools2 queryto extract info from files, than use akw script to format them
#removed multiallelic sites and formatted the output to avoid missing problem when checking frequencies
(echo "CHROM POS ID REF ALT AC AN TGP_AF TGP_AMR_AF TGP_ASN_AF TGP_AFR_AF TGP_EUR_AF IMP2 VQSLOD AF MAF MINOR";bcftools2 query ${vcf_file} -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/1kg_AF\t%INFO/1kg_AMR_AF\t%INFO/1kg_ASN_AF\t%INFO/1kg_AFR_AF\t%INFO/1kg_EUR_AF\t%INFO/IMP2\t%INFO/VQSLOD\n'| fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " " > ${out_path}/${vcf_file}.maf_table.tab


