#!/usr/local/bin/bash
#split new seq file in SNPS and INDELS
###ATTENTION:the vcftool --keep-only-indels and --remove-indels seems not to work properly!!!So use grep instead
#vcftools --gzvcf esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz --keep-only-indels --recode-INFO-all --out indels_split --recode-to-stream | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.INDELS.vcf.gz
cat <(zcat esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz | head -200 | grep ^#) <(zcat esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz | grep -v ^# | grep INDEL) | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.INDELS.vcf.gz
 
tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.INDELS.vcf.gz

#vcftools --gzvcf esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz --remove-indels --recode-INFO-all --out snps_split --recode-to-stream | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz

cat <(zcat esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz | head -200 | grep ^#) <(zcat esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz | grep -v ^# | grep -v INDEL) | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz

tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz
