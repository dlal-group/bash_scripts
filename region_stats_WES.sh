#!/bin/bash

#script to extract stats for all samples for a single region from annotated vcf file
IN_VCF=$1 #input vcf folder
REG_F=$2 #region file
OUT_F=$3 #output folder
FORMAT=$4 #based on format of our file
TYPE=$5 #variation type : SNP / INDEL

while read line
do
	case ${FORMAT} in
		UCSC )
			chr=`echo ${line}| cut -f 2 -d " "| sed 's/chr//g'`
			start=`echo ${line}| cut -f 4 -d " "`
			end=`echo ${line}| cut -f 5 -d " "`
			gene_name=`echo ${line}| cut -f 1 -d " "`
			exon_count=`echo ${line}| cut -f 8 -d " "`
			;;
		GENCODE )
			chr=`echo ${line}| cut -f 3 -d " "|sed 's/chr//g`
			start=`echo ${line}| cut -f 5 -d " "`
			end=`echo ${line}| cut -f 6 -d " "`
			gene_name=`echo ${line}| cut -f 2 -d " "`
			exon_count=`echo ${line}| cut -f 9 -d " "`
			;;
	esac
			gene_length=$[end - start]

echo "$FORMAT,$TYPE"
echo "Analyzing region => ${chr}:${start}-${end} ,"
echo "Length: ${gene_length}."
echo "Gene: ${gene_name},"
echo "NExons: ${exon_count}."

mkdir -p ${OUT_F}/${gene_name}
OUT_VCF=${OUT_F}/${gene_name}/All.multisampleinitial.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.recalibrated.filtered.vcf.gz


# we need to read from our region file and save the extracted region, than we're going to extract stats for that region
bcftools view ${IN_VCF} -r ${chr}:${start}-${end} -O z -o ${OUT_VCF}
tabix -f -p vcf ${OUT_VCF}
bcftools stats -s - ${OUT_VCF} > ${OUT_F}/${gene_name}/WES.${TYPE}.${FORMAT}.${gene_name}.${chr}.${start}.${end}.stats
			

done < <(zcat $REG_F)
