#!/usr/local/bin/bash

#script used to check if the overlap with 370K data has all rsID in seq files
file_path=$2
out_path=$3

#this is for the check in biallelic overlap
#grep -w -f <(grep "^$1	" biallelic_overl.SNP.unfilt.geno.seq.VB.map | cut -f 4 ) <(tabix $file_path/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann_rsIDs.vcf.gz chr $1 | cut -f 1-5) > $out_path/chr$1.370k_overlap_ann_check.tab

#this is for the chek in the whole overlap set
fgrep -w -f <(grep "^$1 " all.chr.geno.overlapped_POS.sorted.txt | cut -f 4 -d " ") <(tabix $file_path/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz chr $1 | cut -f 1-5) > $out_path/chr$1.370k_overlap_ann_check.tab
