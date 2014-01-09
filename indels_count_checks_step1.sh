#!/usr/local/bin/bash

#FIRST STEP of the pipeline used to extract info about indels

#FIRST:
#Count number of indels for each chr and separate knonwn from novel
#ARGS:
vcf_file=$1

if [ $LSB_JOBINDEX -eq 23 ]
then
	chr="X"
else
	chr=$LSB_JOBINDEX
fi

out_path=$2

mkdir -p $out_path/CHR${chr}

filename=`basename $vcf_file`
suffix=${filename%.vcf*}

#extract sites for each chr
(zcat $vcf_file | egrep "^#";tabix $vcf_file ${chr})| bgzip -c > $out_path/CHR${chr}/${suffix}.${chr}.vcf.gz
tabix $out_path/CHR${chr}/${suffix}.${chr}.vcf.gz

#for each chr count known and novel sites
tot_num=`tabix $out_path/CHR${chr}/${suffix}.${chr}.vcf.gz ${chr} | wc -l`
known_num=`tabix $out_path/CHR${chr}/${suffix}.${chr}.vcf.gz ${chr} | awk '{if ($3 != ".") print $0}' | wc -l`
novel_num=`tabix $out_path/CHR${chr}/${suffix}.${chr}.vcf.gz ${chr} | awk '{if ($3 == ".") print $0}' | wc -l`
echo -e "INDELS count report for CHR ${chr}:\n\nTotal INDELS number = $tot_num\nKnown INDELS number = $known_num\nNovel INDELS number = $novel_num" > $out_path/CHR${chr}/INDELS.${chr}.report
~
