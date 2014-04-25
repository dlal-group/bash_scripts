#!/usr/local/bin/bash

#script to find not overlapping sites
#for each chr we creare a file with not overlapping entries
tabix ~/lustre110_home/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.vcf.gz chr $1 | fgrep -v -w -f <(cut -f 1,2 ~/lustre110_home/GENOTIPI/COMPARISON/NOT_OVERLAPPING/OVERLAPPING_LISTS/VBI_1KG_all_overlap_chr$1.csv | grep -v "POS") > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf
