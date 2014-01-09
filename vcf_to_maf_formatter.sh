#!/usr/local/bin/bash

#script we can use to format vcf files in tab spaced file as:
#CHROM\tPOS\tID\tAC\tAN\tAF\tMAF
#Args:
#$1=vcf splitted file absolute path
#$2=out path
#$3= chr

chr=$LSB_JOBINDEX
echo ${chr}
#we need vcf-query to extract info from files, than use akw script to format them
(echo "CHROM	POS	ID	AC	AN	AF	MAF";vcf-query $1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${chr}.re_ann.vcf.gz -f '%CHROM\t%POS\t%ID\t%INFO/AC\t%INFO/AN\n'| fgrep -v CHROM | awk '{OFS="\t"}{print $0,$4/$5}' | awk '{OFS="\t"}{if ($6 <= 0.5) print $0,$6;else print $0,(1-$6)}') > $2/chr${chr}_maf_table.tab


