#!/usr/local/bin/bash
#check differencies after reannotation with and without allele match option

diff <(zcat RS_ann_splitted/RE_ANN/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann_rsIDs.vcf.gz | grep -v ^# | cut -f 1-5) <(zcat RS_ann_splitted/RE_ANN_no_allele_match/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann_rsIDs.vcf.gz | grep -v ^# | cut -f 1-5) > ann_check/chr$1_annotation_diff_file.diff

