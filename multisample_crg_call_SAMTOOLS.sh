#!/bin/bash

#modified script to work using region for calling the chromosome,so we can work in parallel:
#SAMTOOLS version

chr=$1
reg=$2
## file containing all the bams
#BAMS=/nfs/users/xe/ggirotto/multisample/all_pooled.list
BAMS=$3

#out folder
OUTF=$4

#variant type
VARTYPE=$5

#output mode
OMODE=$6

echo $reg
echo $chr
echo ${VARTYPE}
echo ${OUTF}

mkdir -p ${OUTF}/${chr}

#REF=/users/GD/resource/human/hg19/databases/GATK_resources/bundle/2.8/hg19/ucsc.hg19.fasta <- this file generate errors during the contig header check: mismatch of contig names
REF=/nfs/users/GD/resource/human/hg19/hg19.fasta
DBSNP=/users/GD/resource/human/hg19/databases/dbSNP/dbsnp_138.hg19.vcf
#GATK=/users/GD/tools/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar
GATK=/users/GD/tools/GATK/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar
#PROBE=/users/GD/resource/human/probesets/nimblegene/v3/Target_Regions/SeqCap_EZ_Exome_v3_capture.bed_plus_50
#PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/nimblegen_plus50.bed
#PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/REGIONS/nimblegen_plus50_chr${chr}_r${reg}.bed
PROBE=${reg}
GATKRS=/users/GD/resource/human/hg19/databases/GATK_resources/bundle/2.8/hg19
CPU=4

filename=`basename ${reg}`
reg_name=`echo ${filename%.*}|awk 'BEGIN{FS="_"};{print $(NF)}'`

#we need to work in a different way with respect to GATK: we need to tell what class of variants to skip rather than to keep
case ${VARTYPE} in
	INDEL )
		SKIPTYPE='snps'
		;;
	SNP)
		SKIPTYPE='indels'
		;;
esac

## SAMTOOLS multisample call
samtools2 mpileup -b ${BAMS} -l ${PROBE} -t DP,SP,DV -C50 -pm1 -F0.2 -d 100000 | bcftools call -vmO v -f GQ,GP -o $OUTF/${chr}/${chr}.${reg_name}.multisampleinitial.allregions.${VARTYPE}.vcf


if [ -s $OUTF/${chr}/${chr}.${reg_name}.multisampleinitial.allregions.${VARTYPE}.vcf ]
then
	touch $OUTF/${chr}/${chr}.${reg_name}.multisampleinitial.allregions.${VARTYPE}.done
fi
