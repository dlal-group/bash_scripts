#!/usr/local/bin/bash

#script to annotate files by chromosome
#Args: 
#Use LSB_JOBINDEX as chr number!
chr=${LSB_JOBINDEX}

/software/jre1.7.0_25/bin/java -Xmx1000m -Xms1000m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar -T VariantAnnotator \
--dbsnp /lustre/scratch113/projects/fvg_seq/variant_refinemet/annotations/dbSNP-b138/${chr}.dbsnp_138.vcf.gz \
--variant /lustre/scratch113/projects/fvg_seq/variant_refinemet/20131213_IMPUTE/${chr}.vcf.gz \
--out /lustre/scratch113/projects/fvg_seq/variant_refinemet/20140108_ANNOTATION/${chr}.dbsnp.vcf.gz \
-R /lustre/scratch111/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa \
--resource:1kg /lustre/scratch113/projects/fvg_seq/variant_refinemet/annotations/1TGP/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.${chr}.vcf.gz \
-E 1kg.AF \
-E 1kg.AMR_AF \
-E 1kg.ASN_AF \
-E 1kg.AFR_AF \
-E 1kg.EUR_AF
# --dbsnp /lustre/scratch111/resources/variation/Homo_sapiens/grch37/dbsnp_138.vcf.gz \
# --resource:1kg /lustre/scratch111/resources/variation/Homo_sapiens/grch37/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz \
# --list