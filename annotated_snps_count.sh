#!/usr/local/bin/bash
#count how many annotated and not snps we have in a indexed vcf file
#the use must set the file suffix
file=$1

if [ $# -lt 1 ]
then
	echo "WRONG argument number!"
	echo "Usage: annotated_snps_count.sh <file>"
	echo "Arguments:"
	echo "file: file you want to count. (eg. re_ann_rsIDs.vcf.gz)"
	echo "----------------------------------------------------------------------------------------"
	exit 1
fi

for i in {1..22} X
do

	#NO_RSID=$(tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${i}.re_ann_rsIDs.vcf.gz chr ${i} | grep "	\.	" | wc -l)
	#W_RS_ID=$(tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${i}.re_ann_rsIDs.vcf.gz chr ${i} | grep -v "	\.	" | wc -l)
	NO_RSID=$(tabix ${file} ${i} | cut -f 1-3 | grep "\.$" | wc -l)
	W_RS_ID=$(tabix ${file} ${i} | cut -f 1-3 | grep -v "\.$" | wc -l)
	echo "$i $W_RS_ID $NO_RSID"
done > snp_count.txt

