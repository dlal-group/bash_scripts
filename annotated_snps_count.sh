#!/usr/local/bin/bash
#count how many annotated and not snps we have in a indexed vcf file
#the use must set the file suffix
suffix=$1

if [ $# -lt 1 ]
then
	echo "WRONG argument number!"
	echo "Usage: annotated_snps_count.sh <file_suffix>"
	echo "Arguments:"
	echo "file_suffix: suffix for the group of files you want to count. (eg. re_ann_rsIDs.vcf.gz)"
	echo "----------------------------------------------------------------------------------------"
	exit 1
fi

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do

	#NO_RSID=$(tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${i}.re_ann_rsIDs.vcf.gz chr ${i} | grep "	\.	" | wc -l)
	#W_RS_ID=$(tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${i}.re_ann_rsIDs.vcf.gz chr ${i} | grep -v "	\.	" | wc -l)
	NO_RSID=$(tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${i}.$suffix chr ${i} | cut -f 1-3 | grep "\.$" | wc -l)
	W_RS_ID=$(tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${i}.$suffix chr ${i} | cut -f 1-3 | grep -v "\.$" | wc -l)
	echo "Chr $i : "
	echo "Annotated : $W_RS_ID , Not Annotated: $NO_RSID"
	echo ""
done > snp_count.txt

