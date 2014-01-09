#!/usr/local/bin/bash

#A general script to manage TSI files

#cat <(zcat TSI.chr$1.annotate.recode.vcf.gz | grep ^#) <(tabix TSI.chr$1.annotate.recode.vcf.gz chr $1 | grep VT=SNP ) | bgzip -c > ../../TAB/VCF/EUR_ONLY/TSI.chr$1.annotate.recode.snps.vcf.gz

#cat <(zcat TSI.chr$1.annotate.recode.vcf.gz | grep ^#) <(tabix TSI.chr$1.annotate.recode.vcf.gz chr $1 | fgrep EUR_AF ) | bgzip -c > EUR_ONLY/TSI.chr$1.annotate.recode.vcf.gz

#tabix TSI.chr$1.annotate.recode.vcf.gz
#tabix TSI.chr$1.annotate.recode.snps.vcf.gz

#grep -v "rs74733400" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr10.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr10.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs855274" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr12.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr12.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs188690587" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr1.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr1.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs61573637" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr2.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr2.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs55926138" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr6.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr6.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs78198903" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr9.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr9.re_ann.NOT_OVERLAP.vcf.gz

