#!/usr/local/bin/bash
#script to generate a final release after exlusion
#Args:
#$1=vcf file path
#$2=exclusion sample list path
#$3=exclusion sites list path
#$4=output path

vcf=$1
samples=$2
sites=$3
outpath=$4

#In order to exclude sites, we need to tell bcftools what we want to include
mkdir -p ${outpath}

filename=`basename ${vcf}`
#extract a temporary file with chr:pos-pos from the vcf file
bcftools query -f "%CHROM:%POS-%POS\n" ${vcf} -o ${outpath}/${filename}.all_sites.list

#use the exlusion list formatted as CHR:POS-POS to remove sites from the complete file list
fgrep -v -w -f ${sites} ${outpath}/${filename}.all_sites.list | cut -f 1 -d "-" | tr ":" "\t" | awk '{OFS="\t"}{print $1,$2,$2}' > ${outpath}/${filename}.included_sites.list

#now use the sample exclusion list to exclude samples

bcftools view -S ^${samples} -R ${outpath}/${filename}.included_sites.list ${vcf} -O z -o ${outpath}/${filename}

tabix -p vcf ${outpath}/${filename}