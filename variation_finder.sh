#!/usr/local/bin/bash

#script to find how maf variate before and after imputation
#Args:
#$1=filter file path (assume it is only one file)
#$2=path to files to filter (assume it's a group of files you want to filter out)
#$3=out path

for i in {1..22}
do
	fgrep -w -f <(cut -f 2 -d " " <(fgrep -w -f <(egrep "^${i}	" $1 | cut -f 2) <(cut -f 1-3 -d " " $2/CHR${i}/chr${i}.gen_info_rs_only_maf))) $2/CHR${i}/chr${i}.gen_info_rs_only_maf > $3/chr${i}_low_freq.txt

done


