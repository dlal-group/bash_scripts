#!/bin/bash

#modified script to concat together all the called regions and perform the vqsr filtering

chr=$1
#OUTF=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS
OUTF=$2
mkdir -p ${OUTF}/${chr}

#REF=/users/GD/resource/human/hg19/databases/GATK_resources/bundle/2.8/hg19/ucsc.hg19.fasta <- this file generate errors during the contig header check: mismatch of contig names
REF=/nfs/users/GD/resource/human/hg19/hg19.fasta
DBSNP=/users/GD/resource/human/hg19/databases/dbSNP/dbsnp_138.hg19.vcf
GATK=/users/GD/tools/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar
#PROBE=/users/GD/resource/human/probesets/nimblegene/v3/Target_Regions/SeqCap_EZ_Exome_v3_capture.bed_plus_50
#PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/nimblegen_plus50.bed
#PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/nimblegen_plus50_chr${chr}_r${reg}.bed
GATKRS=/users/GD/resource/human/hg19/databases/GATK_resources/bundle/2.8/hg19
CPU=8

## file containing all the bams
#BAMS=all_pooled.list

##Concat the vcf files toghether for the whole chromosome
#create the list of files
size=`ls $OUTF/${chr}/2.multisampleinitial.allregions.snps.r*.vcf| wc -l`;
i=0;
while [ $i -lt ${size} ]
do
 i=$[i+1];
 ls $OUTF/${chr}/2.multisampleinitial.allregions.snps.r${i}.vcf;
done > $OUTF/${chr}/chr_${chr}_vcf.list

#concat the files back together
vcf-concat -f $OUTF/${chr}/chr_${chr}_vcf.list > $OUTF/${chr}.multisampleinitial.allregions.snps.vcf

