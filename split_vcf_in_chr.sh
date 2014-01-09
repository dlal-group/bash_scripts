#/!usr/bin/bash
#script for vcf splitting in 22 chr

cat <(zcat esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz | head -200 | grep ^#) <(tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz chr $1) | bgzip -c > SPLITTED_VCF/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.vcf.gz
