#!/usr/local/bin/bash

#Script to create all necessary input files for genotype concordance calculation
#ARGS:
#$1=chr
#$2=vcf file
#$3=out path

if [ $# -lt 3 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "genotype_concordance_input_prep.sh <chr> <vcf_file> <out_path> [<POP/SET name>]\n"
	exit 1
fi

chr=$1
vcf_file=$2
out_path=$3


input_dir=`dirname ${vcf_file}`
out_name=`basename ${vcf_file}`

if [ $# -eq 4 ]
	then
	set_name=$4
	out_file="${out_name}.${set_name}.snps.biallelic.gz"
else
	out_file="${out_name}.snps.biallelic.gz"
fi

if [ -s ${out_path}/${out_file} ]
then
	echo -e "\nATTENTION!!!a file named ${out_file} already exists in ${out_path}!!!!\nThe script is exiting NOW!!"
	exit 1
fi
#We assume we want to check concordance only for BIALLELIC SNPs
#remove the indels and all multiallelic sites
tabix -f -p vcf ${vcf_file}

(zcat ${vcf_file} | grep "^#";tabix ${vcf_file} ${chr} | fgrep -v INDEL | awk '$5!~","' )| bgzip -c > ${out_path}/${out_file}

#index the file
tabix -p vcf ${out_path}/${out_name}.snps.biallelic.gz

#now create the file with the correct genotype format
vcf-query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${out_path}/${out_file} | awk '{printf $1" "$1":"$2" 0 "$2;for(i=5;i<=NF;i++) printf " " $i;printf "\n"}'| tr "./." "0/0 "| tr ".|." "0|0" | tr "|" " " |tr "/" " " > ${out_path}/${out_file}.tped

#and the corrispondent TFAM file (we suppose for now that we are working with the same samples)
zcat ${out_path}/${out_file} | grep "^#CHROM" | cut -f 10- | tr "\t" "\n" | awk '{print $1,$1,0,0,0,0}' > ${out_path}/${out_file}.tfam

#create the reference table for this file
(echo "CHROM	POS	SNP	MAJOR";vcf-query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${out_path}/${out_file} | fgrep -v "CHROM" | awk '{print $1,$2,$1":"$2,$3}' | tr " " "\t" ) > ${out_path}/${out_file}.reference_allele_table.txt

#Next step will be to extract a list of overlapping sites from TPED formatted files...but I need to run this twice to create the two input files I need...
