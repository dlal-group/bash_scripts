#!/usr/local/bin/bash

#script for fix missing ALT alleles

#take as arguments:
#$1: chr
#$2: chr map file with missing alleles 
#$3: chr map file with fixed alleles
#$4: output file path

#Separate the allele match from the allele mismatch
#First find the allele match
fgrep -v -w -f <(cut -f 1,2 $3) $2 > $4/allele_match.chr$1

#then find the allele mismatch
fgrep -w -f <(cut -f 1,2 $3) $2 > $4/allele_mismatch.chr$1

#now FIX the allele mismatch
fgrep -w -f <(cut -f 1,2 $4/allele_mismatch.chr$1) $3 > $4/fixed_allele_mismatch.chr$1

#then put together again
cat $4/allele_match.chr$1 $4/fixed_allele_mismatch.chr$1 | sort -k1g -k2g > $4/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.vcf.fixed_alleles.map

#now clean up the temp files
rm $4/allele_match.chr$1 $4/allele_mismatch.chr$1 $4/fixed_allele_mismatch.chr$1
