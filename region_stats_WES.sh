#!/bin/bash

#script to extract stats for all samples for a single region from annotated vcf file
IN_VCF=$1 #input vcf folder
REG_F=$2 #region file
OUT_F=$3 #output folder
FORMAT=$4 #based on format of our file

while read line
do
	case ${FORMAT} in
		UCSC )
			chr=`echo ${line}| cut -f 2 -d " "`
			start=`echo ${line}| cut -f 4 -d " "`
			end=`echo ${line}| cut -f 5 -d " "`
			gene_name=`echo ${line}| cut -f 1 -d " "`
			exon_count=`echo ${line}| cut -f 8 -d " "`
			;;
		GENCODE )
			chr=`echo ${line}| cut -f 3 -d " "`
			start=`echo ${line}| cut -f 5 -d " "`
			end=`echo ${line}| cut -f 6 -d " "`
			gene_name=`echo ${line}| cut -f 2 -d " "`
			exon_count=`echo ${line}| cut -f 9 -d " "`
			;;
	esac
			gene_length=$[end - start]

echo $FORMAT
echo "Analyzing region => ${chr}:${start}-${end} ,"
echo "Gene: ${gene_name},"
echo "NExons: ${exon_count}."

# # we need to read from our region file and save the extracted region, than we're going to extract stats for that region
# bcftools view ${INF}/All.multisampleinitial.allregions.${TYPE}.recalibrated.filtered.vcf -O z -o ${INF}/All.multisampleinitial.allregions.${TYPE}.recalibrated.filtered.vcf.gz
# tabix -p vcf ${INF}/All.multisampleinitial.allregions.${TYPE}.recalibrated.filtered.vcf.gz
# bcftools annotate -a /nfs/users/xe/ggirotto/annotations/All_20150102.vcf.gz -c CHROM,POS,ID,REF,ALT ${INF}/All.multisampleinitial.allregions.${TYPE}.recalibrated.filtered.vcf.gz -O z -o ${INF}/All.multisampleinitial.allregions.${TYPE}.recalibrated.filtered.ann.vcf.gz
# tabix -p vcf ${INF}/All.multisampleinitial.allregions.${TYPE}.recalibrated.filtered.ann.vcf.gz
# bcftools stats -s - ${INF}/All.multisampleinitial.allregions.${TYPE}.recalibrated.filtered.ann.vcf.gz > ${INF}/WES_${TYPE}_stats_ALL.txt
			

done < <(zcat $REG_F)
