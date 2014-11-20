#!/bin/bash

#modified script to work using region for calling the chromosome,so we can work in parallel

reg=$1
chr=$2
OUTF=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS
mkdir -p ${OUTF}/${chr}

#REF=/users/GD/resource/human/hg19/databases/GATK_resources/bundle/2.8/hg19/ucsc.hg19.fasta <- this file generate errors during the contig header check: mismatch of contig names
REF=/nfs/users/GD/resource/human/hg19/hg19.fasta
DBSNP=/users/GD/resource/human/hg19/databases/dbSNP/dbsnp_138.hg19.vcf
GATK=/users/GD/tools/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar
#PROBE=/users/GD/resource/human/probesets/nimblegene/v3/Target_Regions/SeqCap_EZ_Exome_v3_capture.bed_plus_50
#PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/nimblegen_plus50.bed
PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/nimblegen_plus50_chr${chr}_r${reg}.bed
GATKRS=/users/GD/resource/human/hg19/databases/GATK_resources/bundle/2.8/hg19
CPU=8

## file containing all the bams
BAMS=all_pooled.list

## GATK initial multisample call

java -jar $GATK -R $REF -T UnifiedGenotyper -I $OUTF/$BAMS -nt $CPU -o $OUTF/${chr}/mutisampleInitialCall_r${reg}.vcf -A DepthPerAlleleBySample -A QualByDepth -A HaplotypeScore -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A InbreedingCoeff -A Coverage --intervals $PROBE -glm BOTH

## Seprate SNPs and INDELs from main multisample call vcf file

java -jar $GATK -R $REF -T SelectVariants --variant $OUTF/${chr}/mutisampleInitialCall_r${reg}.vcf -o $OUTF/${chr}/multisampleinitial.allregions.snps.r${reg}.vcf -selectType SNP
#java -jar $GATK -R $REF -T SelectVariants --variant $OUTF/multisampleinitial.allregions.vcf -o $OUTF/multisampleinitial.allregions.indels.vcf -selectType INDEL


## Variant Recalibration : SNP only
java -jar $GATK -T VariantRecalibrator -R $REF -input $OUTF/multisampleinitial.allregions.snps.vcf -recalFile $OUTF/multisampleinitial.allregions.snps.vcf.recal -tranchesFile $OUTF/multisampleinitial.allregions.snps.vcf.tranches --rscript_file $OUTF/test_VQSR.r --maxGaussians 6 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATKRS/hapmap_3.3.hg19.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 $GATKRS/1000G_omni2.5.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATKRS/dbsnp_138.hg19.vcf -resource:1000g,known=false,training=true,truth=false,prior=10.0 $GATKRS/1000G_phase1.snps.high_confidence.hg19.vcf -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an InbreedingCoeff -mode SNP --target_titv 3.0 --ts_filter_level 92.0 


## Apply Recalibration : SNP only
java -jar $GATK -T ApplyRecalibration -R $REF -input $OUTF/multisampleinitial.allregions.snps.vcf -tranchesFile $OUTF/multisampleinitial.allregions.snps.vcf.tranches -recalFile $OUTF/multisampleinitial.allregions.snps.vcf.recal -o $OUTF/multisampleinitial.allregions.snps.recalibrated.filtered.vcf --ts_filter_level 99.0 -mode SNP


## grep PASS snps from the recalibration
egrep 'PASS|^#' $OUTF/multisampleinitial.allregions.snps.recalibrated.filtered.vcf > $OUTF/multisampleinitial.allregions.snps.recalibrated.filtered.clean.vcf

## variant Evaluation
java -jar $GATK -T VariantEval -R $REF --dbsnp $DBSNP -o $OUTF/report.multisampleinitial.allregions.snps.recalibrated.filtered.clean.vcf --eval $OUTF/multisampleinitial.allregions.snps.recalibrated.filtered.clean.vcf -l INFO --downsampling_type NONE
