#!/bin/bash

#modified script to work using region for calling the chromosome,so we can work in parallel

chr=$1
reg=$2
## file containing all the bams
#BAMS=/nfs/users/xe/ggirotto/multisample/all_pooled.list
BAMS=$3

#out folder
OUTF=$4

#variant type
VARTYPE=$5

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
CPU=8

## GATK initial multisample call

#java -jar $GATK -U LENIENT_VCF_PROCESSING -l INFO -R $REF -T UnifiedGenotyper -I $OUTF/$BAMS -nt $CPU -o $OUTF/${chr}/1.mutisampleInitialCall_r${reg}.vcf -A DepthPerAlleleBySample -A QualByDepth -A HaplotypeScore -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A InbreedingCoeff -A Coverage --intervals $PROBE -glm BOTH
#use old UnifiedGenotyper method
java -jar $GATK -U LENIENT_VCF_PROCESSING -l INFO -R $REF -T UnifiedGenotyper \
-I $BAMS \
-nt $CPU \
-o $OUTF/${chr}/1.mutisampleInitialCall_all.vcf \
-A DepthPerAlleleBySample -A QualByDepth -A HaplotypeScore \
-A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand \
-A InbreedingCoeff -A Coverage \
--intervals $PROBE -glm BOTH

#use new Haplotype caller in gVCF mode
#java -jar $GATK -U LENIENT_VCF_PROCESSING -l INFO -R $REF -T UnifiedGenotyper -I $BAMS -nt $CPU -o $OUTF/${chr}/1.mutisampleInitialCall_all.vcf -A DepthPerAlleleBySample -A QualByDepth -A HaplotypeScore -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A InbreedingCoeff -A Coverage --intervals $PROBE -glm BOTH

## Seprate SNPs and INDELs from main multisample call vcf file

#java -jar $GATK -R $REF -T SelectVariants --variant $OUTF/${chr}/1.mutisampleInitialCall_r${reg}.vcf -o $OUTF/${chr}/2.multisampleinitial.allregions.snps.r${reg}.vcf -selectType SNP
# java -jar $GATK -R $REF -T SelectVariants --variant $OUTF/${chr}/1.mutisampleInitialCall_all.vcf -o $OUTF/${chr}.multisampleinitial.allregions.indels.vcf -selectType INDEL
java -jar $GATK -R $REF -T SelectVariants --variant $OUTF/${chr}/1.mutisampleInitialCall_all.vcf -o $OUTF/${chr}.multisampleinitial.allregions.${VARTYPE}.vcf -selectType ${VARTYPE}

if [ -s $OUTF/${chr}.multisampleinitial.allregions.${VARTYPE}.vcf ]
then
	touch $OUTF/${chr}.multisampleinitial.allregions.${VARTYPE}.done
fi
