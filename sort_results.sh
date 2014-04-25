#!/usr/local/bin/bash

#sort results extracted

add_stats_path=$1
res_path=$2
trait=$3
out=$4

#filter out sites based on predefined tresholds
maf=0.005
info=0.3
hwe=1e-6
pval1=1e-5
pval2=5e-6
pval3=1e-7

#join together all poscon extracted by chr
(head -1 ${out}/${trait}.1.poscon.join ;cat ${out}/${trait}.*.poscon.join|fgrep -v "BETA" | sed 's/^0/X/g' | sort -k1,1g -k3,3g) > ${out}/${trait}.all.poscon.join

#create a filtered file to plot the data
(head -1 ${out}/${trait}.1.result.maf_info_hwe_filtered.join;(cat ${out}/${trait}.*.result.maf_info_hwe_filtered.join | fgrep -v "BETA" | sed 's/^0/23/g' | awk '$6<=1e-1'; cat ${out}/${trait}.*.result.maf_info_hwe_filtered.join | fgrep -v "BETA" | sed 's/^0/23/g' | awk '$6>1e-1' | shuf -n 400000 )| sort -k3,3n -k1,1n)  > ${out}/${trait}.all.result.maf_info_hwe_filtered.join.plot

#filter out from the result files using 3 different pvalues
(head -1 ${out}/${trait}.1.result.maf_info_hwe_filtered.join;cat ${out}/${trait}.*.result.maf_info_hwe_filtered.join | fgrep -v "BETA" | sed 's/^0/23/g' | awk -v tp=${pval1} '$6<=tp' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.maf_info_hwe_pval_${pval1}_filtered.join
(head -1 ${out}/${trait}.1.result.maf_info_hwe_filtered.join;cat ${out}/${trait}.*.result.maf_info_hwe_filtered.join | fgrep -v "BETA" | sed 's/^0/23/g' | awk -v tp=${pval2} '$6<=tp' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.maf_info_hwe_pval_${pval2}_filtered.join
(head -1 ${out}/${trait}.1.result.maf_info_hwe_filtered.join;cat ${out}/${trait}.*.result.maf_info_hwe_filtered.join | fgrep -v "BETA" | sed 's/^0/23/g' | awk -v tp=${pval3} '$6<=tp' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.maf_info_hwe_pval_${pval3}_filtered.join

#filter out from the "all" file using pvalue
(head -1 ${out}/${trait}.1.result.poscon.join.maf_info_hwe_filtered;cat ${out}/${trait}.*.result.poscon.join.maf_info_hwe_filtered | fgrep -v "BETA" | sed 's/^0/23/g' | awk -v tp=${pval1} '$6<=tp' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.poscon.join.maf_info_hwe_pval_${pval1}_filtered
(head -1 ${out}/${trait}.1.result.poscon.join.maf_info_hwe_filtered;cat ${out}/${trait}.*.result.poscon.join.maf_info_hwe_filtered | fgrep -v "BETA" | sed 's/^0/23/g' | awk -v tp=${pval2} '$6<=tp' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.poscon.join.maf_info_hwe_pval_${pval2}_filtered
(head -1 ${out}/${trait}.1.result.poscon.join.maf_info_hwe_filtered;cat ${out}/${trait}.*.result.poscon.join.maf_info_hwe_filtered | fgrep -v "BETA" | sed 's/^0/23/g' | awk -v tp=${pval3} '$6<=tp' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.poscon.join.maf_info_hwe_pval_${pval3}_filtered

#filter out from the "all" file using the poscon columns: extract results that overlap with positive controls
(head -1 ${out}/${trait}.1.result.poscon.join.maf_info_hwe_filtered;cat ${out}/${trait}.*.result.poscon.join.maf_info_hwe_filtered | fgrep -v "BETA" | sed 's/^0/23/g' | awk '$14!="NA" && $15!="NA"' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.poscon.join.maf_info_hwe_filtered
# (head -1 ${out}/${trait}.1.result.poscon.join.maf_info_hwe_filtered;cat ${out}/${trait}.*.result.poscon.join.maf_info_hwe_filtered | fgrep -v "BETA" | sed 's/^0/23/g' | awk -v tp=${pval1} '$6<=tp && $14!="NA" && $15!="NA"' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.poscon.join.maf_info_hwe_pval_${pval1}_filtered
# (head -1 ${out}/${trait}.1.result.poscon.join.maf_info_hwe_filtered;cat ${out}/${trait}.*.result.poscon.join.maf_info_hwe_filtered | fgrep -v "BETA" | sed 's/^0/23/g' | awk -v tp=${pval2} '$6<=tp && $14!="NA" && $15!="NA"' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.poscon.join.maf_info_hwe_pval_${pval2}_filtered
# (head -1 ${out}/${trait}.1.result.poscon.join.maf_info_hwe_filtered;cat ${out}/${trait}.*.result.poscon.join.maf_info_hwe_filtered | fgrep -v "BETA" | sed 's/^0/23/g' | awk -v tp=${pval3} '$6<=tp && $14!="NA" && $15!="NA"' | sort -k1,1n -k3,3n) > ${out}/${trait}.all.result.poscon.join.maf_info_hwe_pval_${pval3}_filtered
