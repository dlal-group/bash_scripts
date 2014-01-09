#!/usr/local/bin/bash

#script to add data to gemma results, instead of using the additional_gemma_stats scripts

#merge chunk files and harmonize snp rs id to be the same as bimbam format:
if [[ ${LSB_JOBINDEX} == "23" ]]
then
	c="X"
else
	c=${LSB_JOBINDEX}
fi
header=`zcat chr1.01.gen.gz.stats.gz | head -1 | cut -f 2,4-11,16- -d " "`
(echo ${header} all_A_freq all_B_freq minor_ALL;(zgrep -v "chromosome" chr${c}.*.gen.gz.stats.gz | sed 's/,rs/|rs/g' | awk -v chr=${c} '{ inf=(NF-3); if($2 ~/^rs/) s=$2;else s="NA"; printf "chr"chr":"$4"-"s"-"$5"-"$6; for(i=1; i<=inf; i++) printf " " $(i+3); printf "\n" }' | cut -f -11,16- -d " " | awk '{print $0,(($8+($9/2))/($8+$9+$10+$11)),(($10+($9/2))/($8+$9+$10+$11))}' | awk '{if($(NF) >= 0.5) print $0,$3;else print $0,$4}'))| gzip -c > chr${c}.gen.gz.stats.gz
