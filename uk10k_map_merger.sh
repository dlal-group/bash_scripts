#!/usr/local/bin/bash


#script to merge uk10k map file and remove duplicates for each chr
#ARGS:
#$1 = chr
#$2 = input_path_1
#$3 = input_path_2
#$4 = output_path

cat <(grep -v CHROM $2/UK10K.chr$1.snps.not_overlap.passed.vcf.map) <(grep -v CHROM $3/UK10K.chr$1.snps.not_overlap.passed.vcf.map) > $4/UK10K.chr$1.snps.not_overlap.map

sort -g -k1,1 -k2,2 $4/UK10K.chr$1.snps.not_overlap.map > $4/UK10K.chr$1.snps.not_overlap.sorted.map

rm $4/UK10K.chr$1.snps.not_overlap.map

#Updated for merging not overlapping sites
#cat <(grep -v CHROM $2/UK10K.chr$1.snps.overlap.map) <(grep -v CHROM $3/UK10K.chr$1.snps.overlap.map) > $4/UK10K.chr$1.snps.overlap.map

#sort -g -k1,1 -k2,2 $4/UK10K.chr$1.snps.overlap.map > $4/UK10K.chr$1.snps.overlap.sorted.map

#rm $4/UK10K.chr$1.snps.overlap.map
