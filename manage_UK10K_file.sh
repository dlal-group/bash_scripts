#!/usr/local/bin/bash

#A general script to manage UK10K files
#we need ABSOLUTE path to work

cat <(zcat $2/chr$1.vcf.gz | grep ^#) <(tabix $2/chr$1.vcf.gz chr $1 | fgrep -v INDEL ) | bgzip -c > $3/UK10K.chr$1.snps.vcf.gz
tabix $3/UK10K.chr$1.snps.vcf.gz

#cat <(zcat TSI.chr$1.annotate.recode.vcf.gz | grep ^#) <(tabix TSI.chr$1.annotate.recode.vcf.gz chr $1 | fgrep EUR_AF ) | bgzip -c > EUR_ONLY/TSI.chr$1.annotate.recode.vcf.gz

#tabix TSI.chr$1.annotate.recode.vcf.gz
#tabix TSI.chr$1.annotate.recode.snps.vcf.gz

#grep -v "rs74733400" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr10.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr10.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs855274" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr12.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr12.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs188690587" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr1.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr1.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs61573637" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr2.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr2.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs55926138" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr6.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr6.re_ann.NOT_OVERLAP.vcf.gz
#grep -v "rs78198903" esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr9.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr9.re_ann.NOT_OVERLAP.vcf.gz

#fgrep -w -f <(grep -v "CHROM" ~/lustre109_home/GENOTIPI/COMPARISON/VBSEQ_QC/VBI_vs_1KG/QC_OUT/OVERLAP/COMPLETE_OVERLAP/VBI_1KG_all_overlap_chr$1.csv | cut -f 1-2 ) <(tabix ~/lustre109_home/GENOTIPI/COMPARISON/VBSEQ_QC/TSI/ANNOTATED/TAB/VCF/TSI.all.annotate.recode.snps.vcf.gz chr $1) | bgzip -c > ~/lustre109_home/GENOTIPI/COMPARISON/VBSEQ_QC/VBI_vs_1KG/QC_OUT/OVERLAP/COMPLETE_OVERLAP/VCF/VBI_1KG_all_overlap_chr$1.vcf.gz

#find how many overlapping sites have AF_EUR and not
#echo "=========================="
#echo "CHR $1"
#echo "Total VBI/1KG overlapping sites:"
#tabix VBI_1KG_all_overlap_chr$1.vcf.gz chr $1 | wc -l
#echo "Overlapping sites with AF_EUR:"
#tabix VBI_1KG_all_overlap_chr$1.vcf.gz chr $1 | fgrep "EUR_AF" | wc -l
#echo "Overlapping sites without AF_EUR"
#tabix VBI_1KG_all_overlap_chr$1.vcf.gz chr $1 | fgrep -v "EUR_AF" | wc -l
