#/!usr/bin/bash
#script for splitting annotation in chr

#zcat annots-rsIDs-dbSNPv137-SNV_only.2012-06-16.tab.gz | grep "^$1	" | bgzip -c > dbSNP_splitted/chr$1_dbSNP_b137.tab.gz
zcat annots-rsIDs-dbSNPv137-no_SNV.2012-06-16.tab.gz | grep "^$1	" | bgzip -c > dbSNP_splitted/chr$1_dbSNP_b137_no_SNV.tab.gz
cd dbSNP_splitted
tabix -s 1 -b 2 -e 2 chr$1_dbSNP_b137.tab.gz
tabix -s 1 -b 2 -e 2 chr$1_dbSNP_b137_no_SNV.tab.gz
cd ..
#cat <(zcat esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz | head -200 | grep ^#) <(tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.vcf.gz chr $1) | bgzip -c > SPLITTED_VCF/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.vcf.gz
