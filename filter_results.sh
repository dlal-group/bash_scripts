#!/usr/local/bin/bash

#filter results according to define thresholds

#merge chunk files and harmonize snp rs id to be the same as bimbam format:
if [[ ${LSB_JOBINDEX} == "23" ]]
then
	c="X"
else
	c=${LSB_JOBINDEX}
fi

add_stats_path=$1
res_path=$2
trait=$3
out=$4

#filter out sites based on predefined tresholds
maf=0.005
info=0.3
hwe=1e-6
# pval1=1e-5
# pval2=5e-6
# pval3=1e-7

awk -v tm=${maf} -v ti=${info} -v th=${hwe} '$6>=tm && $7>=ti && $8>=th' ${add_stats_path}/${trait}.gemma.${c}.stats.reordered.final > ${out}/${trait}.${c}.maf_info_hwe_filtered_list

#than we want to use this list of filtered sites to extract the relevant results without pvalues treshold
(head -1 ${res_path}/${trait}.chr${c}.tab.assoc.txt;fgrep -w -f <( cut -f 3 -d " " ${out}/${trait}.${c}.maf_info_hwe_filtered_list ) ${res_path}/${trait}.chr${c}.tab.assoc.txt | sort -g -k3,3 ) | cut -f -6,9,11 > ${out}/${trait}.${c}.result.maf_info_hwe_filtered

#create a joint resume table
(echo "CHR RS POS BETA SE P_wald P_score MAF INFO_snptest HWE all_0_freq all_1_freq eff_all_freq";awk 'FNR==NR { a[$3]=$0; next } $2 in a { print $0,a[$2] }' ${out}/${trait}.${c}.maf_info_hwe_filtered_list ${out}/${trait}.${c}.result.maf_info_hwe_filtered |  tr "\t" " " | cut -f -3,5-8,14-16,18- -d " ") > ${out}/${trait}.${c}.result.maf_info_hwe_filtered.join

#create a joint resume with positive control pval info
(awk 'FNR==NR { a[$2]=$0; next }{if($2 in a) {split(a[$2],x," ");print $0,x[8],x[9]}else{print $0,"NA","NA"} }' ${out}/${trait}.${c}.poscon.join ${out}/${trait}.${c}.result.maf_info_hwe_filtered.join |  tr "\t" " " ) > ${out}/${trait}.${c}.result.poscon.join.maf_info_hwe_filtered

#fgrep -v "BETA" ${out}/${trait}.${c}.result.maf_info_hwe_filtered.join | sed 's/^0/23/g' | sort -k3,3n -k1,1n | awk '$6<=1e-1'>> ${out}/${trait}.all.result.maf_info_hwe_filtered.join.plot


#filter out from the "all" file using pvalue
#fgrep -v "BETA" ${out}/${trait}.${c}.result.maf_info_hwe_filtered.join | sed 's/^0/23/g' | awk -v tp=${pval} '$6<=tp' | sort -k1,1n -k3,3n >> ${out}/${trait}.all.result.maf_info_hwe_pval_filtered.join
