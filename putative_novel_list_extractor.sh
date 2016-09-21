#!/usr/local/bin/bash

#script for creation of putative novel site list

#take as arguments:
#$1: chr
#$2: filter file path
#$3: filtered file path
#$4: output_path

#first: remove from not overlap with 1KG all sites overlapping with UK10K

fgrep -v -w -f <(cut -f 1,2 $2) <(tabix $3 chr $1) | bgzip -c > $4/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.vcf.gz 

tabix $4/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.vcf.gz


