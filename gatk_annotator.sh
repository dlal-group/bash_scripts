#!/usr/local/bin/bash

#script to annotate files by chromosome
#Args: 
#Use LSB_JOBINDEX as chr number!
# chr=${LSB_JOBINDEX}

chr=$1
infolder=$2
outfolder=$3

#first fix the header
(tabix -H ${infolder}/${chr}.vcf.gz ${chr}|grep "^##";cat /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140711_ANNOTATED/to_add_to_header2.txt;tabix -H ${infolder}/${chr}.vcf.gz ${chr}| tail -n1;tabix ${infolder}/${chr}.vcf.gz ${chr})| bgzip -c > ${infolder}/${chr}.fixed.vcf.gz;
# (tabix -H 1.vcf.gz 1|grep "^##";cat /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140711_ANNOTATED/to_add_to_header2.txt;tabix -H 1.vcf.gz 1| tail -n1)|
tabix -p vcf ${infolder}/${chr}.fixed.vcf.gz;

bcftools annotate -x INFO/1kg_AA,INFO/1kg_AF,INFO/1kg_AMR_AF,INFO/1kg_ASN_AF,INFO/1kg_AFR_AF,INFO/1kg_EUR_AF ${infolder}/${chr}.fixed.vcf.gz -O z -o ${infolder}/${chr}.fixed.cleaned.vcf.gz
tabix -p vcf ${infolder}/${chr}.fixed.cleaned.vcf.gz


/software/jre1.7.0_25/bin/java -Xmx1000m -Xms1000m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T VariantAnnotator \
--dbsnp /lustre/scratch113/projects/fvg_seq/variant_refinemet/annotations/dbSNP-b138/${chr}.dbsnp_138.vcf.gz \
--variant ${infolder}/${chr}.fixed.cleaned.vcf.gz \
--out ${outfolder}/${chr}.vcf.gz \
-R /lustre/scratch111/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa \
--resource:1kg /lustre/scratch113/projects/fvg_seq/variant_refinemet/annotations/1TGP/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz \
-E 1kg.AA \
-E 1kg.AF \
-E 1kg.AMR_AF \
-E 1kg.ASN_AF \
-E 1kg.AFR_AF \
-E 1kg.EUR_AF

#index with tabix
tabix -p vcf ${outfolder}/${chr}.vcf.gz

zcat ${outfolder}/${chr}.vcf.gz | sed "s/1kg\.//g" | bgzip -c > ${outfolder}/${chr}.clean_annotated.vcf.gz

tabix -p vcf ${outfolder}/${chr}.clean_annotated.vcf.gz

#clean stuff
rm ${infolder}/${chr}.fixed.vcf.gz
rm ${infolder}/${chr}.fixed.cleaned.vcf.gz
rm ${infolder}/${chr}.fixed.vcf.gz.tbi
rm ${infolder}/${chr}.fixed.cleaned.vcf.gz.tbi