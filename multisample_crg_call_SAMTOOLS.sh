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
REF=/nfs/users/GD/resource/human/hg19/hg19.fa
DBSNP=/users/GD/resource/human/hg19/databases/dbSNP/dbsnp_138.hg19.vcf
#GATK=/users/GD/tools/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar
GATK=/users/GD/tools/GATK/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar
#PROBE=/users/GD/resource/human/probesets/nimblegene/v3/Target_Regions/SeqCap_EZ_Exome_v3_capture.bed_plus_50
#PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/nimblegen_plus50.bed
#PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/REGIONS/nimblegen_plus50_chr${chr}_r${reg}.bed
PROBE=${reg}
GATKRS=/users/GD/resource/human/hg19/databases/GATK_resources/bundle/2.8/hg19
CPU=4
#hard code ped file location
# pedfile=/nfs/users/xe/ggirotto/multisample/FVG_WES_samples_no_fam_20102014_IDS_sorted_uniq.ped
pedfile=/nfs/users/xe/ggirotto/multisample/FVG_WES_samples_no_fam_20102014_sex.samples


filename=`basename ${reg}`
reg_name=`echo ${filename%.*}|awk 'BEGIN{FS="_"};{print $(NF)}'`

#we need to work in a different way with respect to GATK: we need to tell what class of variants to skip rather than to keep
case ${VARTYPE} in
	INDEL )
		SKIPTYPE='snps'
		## SAMTOOLS multisample call
		;;
	SNP)
		SKIPTYPE='indels'
		;;
esac
if [ $chr = "X" ]
then
	## SAMTOOLS multisample call
	samtools mpileup -b ${BAMS} -l ${PROBE} -r X:1-60000 -f ${REF} -t DP,SP,DV -s -C50 -pm1 -F0.2 -d 100000 -g -u | bcftools call --skip-variants ${SKIPTYPE} -X -t X:1-60000 -S ${pedfile} -vmO v -f GQ,GP -o $OUTF/${chr}/${chr}.NON_PAR_1.multisampleinitial.allregions.${VARTYPE}.vcf
	samtools mpileup -b ${BAMS} -l ${PROBE} -r X:60001-2699520 -f ${REF} -t DP,SP,DV -s -C50 -pm1 -F0.2 -d 100000 -g -u | bcftools call --skip-variants ${SKIPTYPE} -t X:60001-2699520 -S ${pedfile} -vmO v -f GQ,GP -o $OUTF/${chr}/${chr}.PAR_1.multisampleinitial.allregions.${VARTYPE}.vcf
	samtools mpileup -b ${BAMS} -l ${PROBE} -r X:2699521-154931043 -f ${REF} -t DP,SP,DV -s -C50 -pm1 -F0.2 -d 100000 -g -u | bcftools call --skip-variants ${SKIPTYPE} -X -t X:2699521-154931043 -S ${pedfile} -vmO v -f GQ,GP -o $OUTF/${chr}/${chr}.NON_PAR_2.multisampleinitial.allregions.${VARTYPE}.vcf
	samtools mpileup -b ${BAMS} -l ${PROBE} -r X:154931044-155270560 -f ${REF} -t DP,SP,DV -s -C50 -pm1 -F0.2 -d 100000 -g -u | bcftools call --skip-variants ${SKIPTYPE} -t X:154931044-155270560 -S ${pedfile} -vmO v -f GQ,GP -o $OUTF/${chr}/${chr}.PAR_2.multisampleinitial.allregions.${VARTYPE}.vcf
	bcftools concat $OUTF/${chr}/${chr}.NON_PAR_1.multisampleinitial.allregions.${VARTYPE}.vcf $OUTF/${chr}/${chr}.PAR_1.multisampleinitial.allregions.${VARTYPE}.vcf $OUTF/${chr}/${chr}.NON_PAR_2.multisampleinitial.allregions.${VARTYPE}.vcf $OUTF/${chr}/${chr}.PAR_2.multisampleinitial.allregions.${VARTYPE}.vcf -O v -o $OUTF/${chr}/${chr}.${reg_name}.multisampleinitial.allregions.${VARTYPE}.vcf
else
	## SAMTOOLS multisample call
	samtools mpileup -b ${BAMS} -l ${PROBE} -f ${REF} -t DP,SP,DV -s -C50 -pm1 -F0.2 -d 100000 -g -u | bcftools call --skip-variants ${SKIPTYPE} -vmO v -f GQ,GP -o $OUTF/${chr}/${chr}.${reg_name}.multisampleinitial.allregions.${VARTYPE}.vcf
fi
		

echo "${BAMS}"


if [ -s $OUTF/${chr}/${chr}.${reg_name}.multisampleinitial.allregions.${VARTYPE}.vcf ]
then
	touch $OUTF/${chr}/${chr}.${reg_name}.multisampleinitial.allregions.${VARTYPE}.done
fi
