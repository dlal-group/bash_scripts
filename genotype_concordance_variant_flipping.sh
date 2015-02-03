#!/bin/bash

#compare 2 freq tables and write  a list of sites to flip
freq1=$1
freq2=$2
seq_set=$3
#get working directory

#extract alleles mismatch between freq tables
awk 'FNR==NR { a[$2]=$0; next } $2 in a {split(a[$2],b," ");if ($2$3$4!=b[2]b[3]b[4]) print $2}' $freq1 $freq2 > discordance_sites_to_flip.list

#this list is ment to be used on the second genotype set (the sequence one) to adjust strands
#so now we flip the second set
plink2 --file ${seq_set} --flip discordance_sites_to_flip.list --recode --out ${seq_set}.flipped

