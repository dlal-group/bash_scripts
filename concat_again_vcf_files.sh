#!/usr/local/bin/bash

#script made for concatenate again splitted vcf files

cat <(zcat ../RS_ann_splitted/RE_ANN_no_allele_match/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann_rsIDs.vcf.gz | grep ^# ) <(cat <(tabix ../RS_ann_splitted/RE_ANN_no_allele_match/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann_rsIDs.vcf.gz chr $1) <(tabix ../RS_ann_splitted/W_RSID/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.w_rsIDs.vcf.gz chr $1) | sort -g -k1,1 -k2,2 ) | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.vcf.gz

tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.vcf.gz
