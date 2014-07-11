#!/usr/local/bin/bash

#script to annotate files by chromosome
#Args: 
#Use LSB_JOBINDEX as chr number!
# chr=${LSB_JOBINDEX}
in_file=$1
out_file=$2
dbsnp_ann=$3


/software/jre1.7.0_25/bin/java -Xmx1000m -Xms1000m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T VariantAnnotator -nt 2\
--dbsnp ${dbsnp_ann} \
--variant ${in_file} \
--out ${out_file} \
-R /lustre/scratch111/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa \
--resource:1kg /lustre/scratch113/projects/fvg_seq/variant_refinemet/annotations/1TGP/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz \
-E 1kg.AA \
-E 1kg.AF \
-E 1kg.AMR_AF \
-E 1kg.ASN_AF \
-E 1kg.AFR_AF \
-E 1kg.EUR_AF
# --dbsnp /lustre/scratch111/resources/variation/Homo_sapiens/grch37/dbsnp_138.vcf.gz \
# --resource:1kg /lustre/scratch111/resources/variation/Homo_sapiens/grch37/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz \
# --list
tabix -p vcf ${out_file}

zcat ${out_file} | sed "s/1kg\.//g" | bgzip -c > ${out_file}.clean_annotated.vcf.gz

