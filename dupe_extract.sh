#!/usr/bin/env bash
#
# Script to extract data formatted in a certain way with information on duplicate data

tab_file=$1
dupe_pos_file=$2
dupe_sites_file=$3

cut -f 2 ${tab_file}| uniq -c| awk 'BEGIN{OFS="\t"}{print $3,$1-1}' > ${dupe_pos_file}
cut -f -2,5- ${tab_file} | uniq -c |awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$1-1}' > ${dupe_sites_file}
