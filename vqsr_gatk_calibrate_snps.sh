#!/usr/local/bin/bash
#script to train a new vqsr filter model to test different solution and expand thresholds:
#########################
#     SNPs ONLY!!!		#
#########################
if [[ $# -lt 3 ]]; then
	echo "##############################################"
	echo "Attention!!! Not enough arguments provided!!!"
	echo "##############################################"
	echo "Usage:"
	echo "vqsr_gatk_calibrate_snps.sh <input_file_path> <output_path> <mode>."
	echo -e "The 'mode' option is used to select the GATK version to use:\n1 - GATK v.2.5\n2 - GATK v.2.7\n3 - GATK v.2.8\n5 - Latest GATK v.3.3"
	exit 1
fi


input=$1
OUTF=$2
MODE=$3
REF=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa

mkdir -p $OUTF/LOGS

#those commands only apply if you need to apply the filter from scratch, without the mpileup pipeline


case $MODE in
	1)
	echo "2.5"
	/software/jre1.7.0_25/bin/java -Xmx12000m -Xms12000m -Xss280m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-2.5/GenomeAnalysisTK.jar \
	-T VariantRecalibrator -U LENIENT_VCF_PROCESSING --mode SNP -l INFO -R $REF -input $input \
	--ignore_filter LowQual \
	-recalFile $OUTF/vqsr.sites.snps.vcf \
	-tranchesFile $OUTF/vqsr.sites.snps.tranches \
	--rscript_file $OUTF/snp_VQSR.r \
	--maxGaussians 6 \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/hapmap_3.3.b37.vcf \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/1000G_omni2.5.b37.vcf \
	-resource:1000g,known=false,training=true,truth=false,prior=10.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/1000G_phase1.snps.high_confidence.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.8/dbsnp_138.b37.excluding_sites_after_129.vcf \
	-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff -an DP \
	--target_titv 2.15 --ts_filter_level 96.00 \
	-minNumBad 5000 \
	-percentBad 0.05 \
	-tranche 70.00 -tranche 71.00 -tranche 72.00 -tranche 73.00 -tranche 74.00 -tranche 75.00 -tranche 76.00 -tranche 77.00 -tranche 78.00 -tranche 79.00 -tranche 80.00 -tranche 81.00 -tranche 82.00 -tranche 83.00 -tranche 84.00 -tranche 85.00 -tranche 86.00 -tranche 87.00 -tranche 88.00 -tranche 89.00 -tranche 90.00 -tranche 91.00 -tranche 92.00 -tranche 93.00 -tranche 94.00 -tranche 95.00 -tranche 96.00 -tranche 96.02 -tranche 96.04 -tranche 96.06 -tranche 96.08 -tranche 96.10 -tranche 96.12 -tranche 96.14 -tranche 96.16 -tranche 96.18 -tranche 96.20 -tranche 96.22 -tranche 96.24 -tranche 96.26 -tranche 96.28 -tranche 96.30 -tranche 96.32 -tranche 96.34 -tranche 96.36 -tranche 96.38 -tranche 96.40 -tranche 96.42 -tranche 96.44 -tranche 96.46 -tranche 96.48 -tranche 96.50 -tranche 96.52 -tranche 96.54 -tranche 96.56 -tranche 96.58 -tranche 96.60 -tranche 96.62 -tranche 96.64 -tranche 96.66 -tranche 96.68 -tranche 96.70 -tranche 96.72 -tranche 96.74 -tranche 96.76 -tranche 96.78 -tranche 96.80 -tranche 96.82 -tranche 96.84 -tranche 96.86 -tranche 96.88 -tranche 96.90 -tranche 96.92 -tranche 96.94 -tranche 96.96 -tranche 96.98 -tranche 97.00 -tranche 97.02 -tranche 97.04 -tranche 97.06 -tranche 97.08 -tranche 97.10 -tranche 97.12 -tranche 97.14 -tranche 97.16 -tranche 97.18 -tranche 97.20 -tranche 97.22 -tranche 97.24 -tranche 97.26 -tranche 97.28 -tranche 97.30 -tranche 97.32 -tranche 97.34 -tranche 97.36 -tranche 97.38 -tranche 97.40 -tranche 97.42 -tranche 97.44 -tranche 97.46 -tranche 97.48 -tranche 97.50 -tranche 97.52 -tranche 97.54 -tranche 97.56 -tranche 97.58 -tranche 97.60 -tranche 97.62 -tranche 97.64 -tranche 97.66 -tranche 97.68 -tranche 97.70 -tranche 97.72 -tranche 97.74 -tranche 97.76 -tranche 97.78 -tranche 97.80 -tranche 97.82 -tranche 97.84 -tranche 97.86 -tranche 97.88 -tranche 97.90 -tranche 97.92 -tranche 97.94 -tranche 97.96 -tranche 97.98 -tranche 98.00 -tranche 98.02 -tranche 98.04 -tranche 98.06 -tranche 98.08 -tranche 98.10 -tranche 98.12 -tranche 98.14 -tranche 98.16 -tranche 98.18 -tranche 98.20 -tranche 98.22 -tranche 98.24 -tranche 98.26 -tranche 98.28 -tranche 98.30 -tranche 98.32 -tranche 98.34 -tranche 98.36 -tranche 98.38 -tranche 98.40 -tranche 98.42 -tranche 98.44 -tranche 98.46 -tranche 98.48 -tranche 98.50 -tranche 98.52 -tranche 98.54 -tranche 98.56 -tranche 98.58 -tranche 98.60 -tranche 98.62 -tranche 98.64 -tranche 98.66 -tranche 98.68 -tranche 98.70 -tranche 98.72 -tranche 98.74 -tranche 98.76 -tranche 98.78 -tranche 98.80 -tranche 98.82 -tranche 98.84 -tranche 98.86 -tranche 98.88 -tranche 98.90 -tranche 98.92 -tranche 98.94 -tranche 98.96 -tranche 98.98 -tranche 99.00 -tranche 100.00
	# -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff -an DP \
	# -tranche 70.00 -tranche 71.00 -tranche 72.00 -tranche 73.00 -tranche 74.00 -tranche 75.00 -tranche 76.00 -tranche 77.00 -tranche 78.00 -tranche 79.00 -tranche 80.00 -tranche 81.00 -tranche 82.00 -tranche 83.00 -tranche 84.00 -tranche 85.00 -tranche 86.00 -tranche 87.00 -tranche 88.00 -tranche 89.00 -tranche 90.00 -tranche 91.00 -tranche 92.00 -tranche 93.00 -tranche 94.00 -tranche 95.00 -tranche 96.00 -tranche 96.02 -tranche 96.04 -tranche 96.06 -tranche 96.08 -tranche 96.10 -tranche 96.12 -tranche 96.14 -tranche 96.16 -tranche 96.18 -tranche 96.20 -tranche 96.22 -tranche 96.24 -tranche 96.26 -tranche 96.28 -tranche 96.30 -tranche 96.32 -tranche 96.34 -tranche 96.36 -tranche 96.38 -tranche 96.40 -tranche 96.42 -tranche 96.44 -tranche 96.46 -tranche 96.48 -tranche 96.50 -tranche 96.52 -tranche 96.54 -tranche 96.56 -tranche 96.58 -tranche 96.60 -tranche 96.62 -tranche 96.64 -tranche 96.66 -tranche 96.68 -tranche 96.70 -tranche 96.72 -tranche 96.74 -tranche 96.76 -tranche 96.78 -tranche 96.80 -tranche 96.82 -tranche 96.84 -tranche 96.86 -tranche 96.88 -tranche 96.90 -tranche 96.92 -tranche 96.94 -tranche 96.96 -tranche 96.98 -tranche 97.00 -tranche 97.02 -tranche 97.04 -tranche 97.06 -tranche 97.08 -tranche 97.10 -tranche 97.12 -tranche 97.14 -tranche 97.16 -tranche 97.18 -tranche 97.20 -tranche 97.22 -tranche 97.24 -tranche 97.26 -tranche 97.28 -tranche 97.30 -tranche 97.32 -tranche 97.34 -tranche 97.36 -tranche 97.38 -tranche 97.40 -tranche 97.42 -tranche 97.44 -tranche 97.46 -tranche 97.48 -tranche 97.50 -tranche 97.52 -tranche 97.54 -tranche 97.56 -tranche 97.58 -tranche 97.60 -tranche 97.62 -tranche 97.64 -tranche 97.66 -tranche 97.68 -tranche 97.70 -tranche 97.72 -tranche 97.74 -tranche 97.76 -tranche 97.78 -tranche 97.80 -tranche 97.82 -tranche 97.84 -tranche 97.86 -tranche 97.88 -tranche 97.90 -tranche 97.92 -tranche 97.94 -tranche 97.96 -tranche 97.98 -tranche 98.00 -tranche 98.02 -tranche 98.04 -tranche 98.06 -tranche 98.08 -tranche 98.10 -tranche 98.12 -tranche 98.14 -tranche 98.16 -tranche 98.18 -tranche 98.20 -tranche 98.22 -tranche 98.24 -tranche 98.26 -tranche 98.28 -tranche 98.30 -tranche 98.32 -tranche 98.34 -tranche 98.36 -tranche 98.38 -tranche 98.40 -tranche 98.42 -tranche 98.44 -tranche 98.46 -tranche 98.48 -tranche 98.50 -tranche 98.52 -tranche 98.54 -tranche 98.56 -tranche 98.58 -tranche 98.60 -tranche 98.62 -tranche 98.64 -tranche 98.66 -tranche 98.68 -tranche 98.70 -tranche 98.72 -tranche 98.74 -tranche 98.76 -tranche 98.78 -tranche 98.80 -tranche 98.82 -tranche 98.84 -tranche 98.86 -tranche 98.88 -tranche 98.90 -tranche 98.92 -tranche 98.94 -tranche 98.96 -tranche 98.98 -tranche 99.00 -tranche 99.02 -tranche 99.04 -tranche 99.06 -tranche 99.08 -tranche 99.10 -tranche 99.12 -tranche 99.14 -tranche 99.16 -tranche 99.18 -tranche 99.20 -tranche 99.22 -tranche 99.24 -tranche 99.26 -tranche 99.28 -tranche 99.30 -tranche 99.32 -tranche 99.34 -tranche 99.36 -tranche 99.38 -tranche 99.40 -tranche 99.42 -tranche 99.44 -tranche 99.46 -tranche 99.48 -tranche 99.50 -tranche 99.52 -tranche 99.54 -tranche 99.56 -tranche 99.58 -tranche 99.60 -tranche 99.62 -tranche 99.64 -tranche 99.66 -tranche 99.68 -tranche 99.70 -tranche 99.72 -tranche 99.74 -tranche 99.76 -tranche 99.78 -tranche 99.80 -tranche 99.82 -tranche 99.84 -tranche 99.86 -tranche 99.88 -tranche 99.90 -tranche 99.92 -tranche 99.94 -tranche 99.96 -tranche 99.98 -tranche 100.00
	# -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.8/dbsnp_138.b37.vcf \
	;;

	2)
	echo "2.7"
	/software/jre1.7.0_25/bin/java -Xmx12000m -Xms12000m -Xss280m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar \
	-T VariantRecalibrator -U LENIENT_VCF_PROCESSING --mode SNP -l INFO -R $REF -input $input \
	--ignore_filter LowQual \
	-recalFile $OUTF/vqsr.sites.snps.vcf \
	-tranchesFile $OUTF/vqsr.sites.snps.tranches \
	--rscript_file $OUTF/snp_VQSR.r \
	--maxGaussians 6 \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/hapmap_3.3.b37.vcf \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/1000G_omni2.5.b37.vcf \
	-resource:1000g,known=false,training=true,truth=false,prior=10.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/1000G_phase1.snps.high_confidence.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.8/dbsnp_138.b37.excluding_sites_after_129.vcf \
	-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff -an DP \
	--target_titv 2.15 --ts_filter_level 96.00 \
	--numBadVariants 5000 \
	-tranche 70.00 -tranche 71.00 -tranche 72.00 -tranche 73.00 -tranche 74.00 -tranche 75.00 -tranche 76.00 -tranche 77.00 -tranche 78.00 -tranche 79.00 -tranche 80.00 -tranche 81.00 -tranche 82.00 -tranche 83.00 -tranche 84.00 -tranche 85.00 -tranche 86.00 -tranche 87.00 -tranche 88.00 -tranche 89.00 -tranche 90.00 -tranche 91.00 -tranche 92.00 -tranche 93.00 -tranche 94.00 -tranche 95.00 -tranche 96.00 -tranche 96.02 -tranche 96.04 -tranche 96.06 -tranche 96.08 -tranche 96.10 -tranche 96.12 -tranche 96.14 -tranche 96.16 -tranche 96.18 -tranche 96.20 -tranche 96.22 -tranche 96.24 -tranche 96.26 -tranche 96.28 -tranche 96.30 -tranche 96.32 -tranche 96.34 -tranche 96.36 -tranche 96.38 -tranche 96.40 -tranche 96.42 -tranche 96.44 -tranche 96.46 -tranche 96.48 -tranche 96.50 -tranche 96.52 -tranche 96.54 -tranche 96.56 -tranche 96.58 -tranche 96.60 -tranche 96.62 -tranche 96.64 -tranche 96.66 -tranche 96.68 -tranche 96.70 -tranche 96.72 -tranche 96.74 -tranche 96.76 -tranche 96.78 -tranche 96.80 -tranche 96.82 -tranche 96.84 -tranche 96.86 -tranche 96.88 -tranche 96.90 -tranche 96.92 -tranche 96.94 -tranche 96.96 -tranche 96.98 -tranche 97.00 -tranche 97.02 -tranche 97.04 -tranche 97.06 -tranche 97.08 -tranche 97.10 -tranche 97.12 -tranche 97.14 -tranche 97.16 -tranche 97.18 -tranche 97.20 -tranche 97.22 -tranche 97.24 -tranche 97.26 -tranche 97.28 -tranche 97.30 -tranche 97.32 -tranche 97.34 -tranche 97.36 -tranche 97.38 -tranche 97.40 -tranche 97.42 -tranche 97.44 -tranche 97.46 -tranche 97.48 -tranche 97.50 -tranche 97.52 -tranche 97.54 -tranche 97.56 -tranche 97.58 -tranche 97.60 -tranche 97.62 -tranche 97.64 -tranche 97.66 -tranche 97.68 -tranche 97.70 -tranche 97.72 -tranche 97.74 -tranche 97.76 -tranche 97.78 -tranche 97.80 -tranche 97.82 -tranche 97.84 -tranche 97.86 -tranche 97.88 -tranche 97.90 -tranche 97.92 -tranche 97.94 -tranche 97.96 -tranche 97.98 -tranche 98.00 -tranche 98.02 -tranche 98.04 -tranche 98.06 -tranche 98.08 -tranche 98.10 -tranche 98.12 -tranche 98.14 -tranche 98.16 -tranche 98.18 -tranche 98.20 -tranche 98.22 -tranche 98.24 -tranche 98.26 -tranche 98.28 -tranche 98.30 -tranche 98.32 -tranche 98.34 -tranche 98.36 -tranche 98.38 -tranche 98.40 -tranche 98.42 -tranche 98.44 -tranche 98.46 -tranche 98.48 -tranche 98.50 -tranche 98.52 -tranche 98.54 -tranche 98.56 -tranche 98.58 -tranche 98.60 -tranche 98.62 -tranche 98.64 -tranche 98.66 -tranche 98.68 -tranche 98.70 -tranche 98.72 -tranche 98.74 -tranche 98.76 -tranche 98.78 -tranche 98.80 -tranche 98.82 -tranche 98.84 -tranche 98.86 -tranche 98.88 -tranche 98.90 -tranche 98.92 -tranche 98.94 -tranche 98.96 -tranche 98.98 -tranche 99.00 -tranche 100.00
	# -tranche 70.00 -tranche 71.00 -tranche 72.00 -tranche 73.00 -tranche 74.00 -tranche 75.00 -tranche 76.00 -tranche 77.00 -tranche 78.00 -tranche 79.00 -tranche 80.00 -tranche 81.00 -tranche 82.00 -tranche 83.00 -tranche 84.00 -tranche 85.00 -tranche 86.00 -tranche 87.00 -tranche 88.00 -tranche 89.00 -tranche 90.00 -tranche 91.00 -tranche 92.00 -tranche 93.00 -tranche 94.00 -tranche 95.00 -tranche 96.00 -tranche 96.02 -tranche 96.04 -tranche 96.06 -tranche 96.08 -tranche 96.10 -tranche 96.12 -tranche 96.14 -tranche 96.16 -tranche 96.18 -tranche 96.20 -tranche 96.22 -tranche 96.24 -tranche 96.26 -tranche 96.28 -tranche 96.30 -tranche 96.32 -tranche 96.34 -tranche 96.36 -tranche 96.38 -tranche 96.40 -tranche 96.42 -tranche 96.44 -tranche 96.46 -tranche 96.48 -tranche 96.50 -tranche 96.52 -tranche 96.54 -tranche 96.56 -tranche 96.58 -tranche 96.60 -tranche 96.62 -tranche 96.64 -tranche 96.66 -tranche 96.68 -tranche 96.70 -tranche 96.72 -tranche 96.74 -tranche 96.76 -tranche 96.78 -tranche 96.80 -tranche 96.82 -tranche 96.84 -tranche 96.86 -tranche 96.88 -tranche 96.90 -tranche 96.92 -tranche 96.94 -tranche 96.96 -tranche 96.98 -tranche 97.00 -tranche 97.02 -tranche 97.04 -tranche 97.06 -tranche 97.08 -tranche 97.10 -tranche 97.12 -tranche 97.14 -tranche 97.16 -tranche 97.18 -tranche 97.20 -tranche 97.22 -tranche 97.24 -tranche 97.26 -tranche 97.28 -tranche 97.30 -tranche 97.32 -tranche 97.34 -tranche 97.36 -tranche 97.38 -tranche 97.40 -tranche 97.42 -tranche 97.44 -tranche 97.46 -tranche 97.48 -tranche 97.50 -tranche 97.52 -tranche 97.54 -tranche 97.56 -tranche 97.58 -tranche 97.60 -tranche 97.62 -tranche 97.64 -tranche 97.66 -tranche 97.68 -tranche 97.70 -tranche 97.72 -tranche 97.74 -tranche 97.76 -tranche 97.78 -tranche 97.80 -tranche 97.82 -tranche 97.84 -tranche 97.86 -tranche 97.88 -tranche 97.90 -tranche 97.92 -tranche 97.94 -tranche 97.96 -tranche 97.98 -tranche 98.00 -tranche 98.02 -tranche 98.04 -tranche 98.06 -tranche 98.08 -tranche 98.10 -tranche 98.12 -tranche 98.14 -tranche 98.16 -tranche 98.18 -tranche 98.20 -tranche 98.22 -tranche 98.24 -tranche 98.26 -tranche 98.28 -tranche 98.30 -tranche 98.32 -tranche 98.34 -tranche 98.36 -tranche 98.38 -tranche 98.40 -tranche 98.42 -tranche 98.44 -tranche 98.46 -tranche 98.48 -tranche 98.50 -tranche 98.52 -tranche 98.54 -tranche 98.56 -tranche 98.58 -tranche 98.60 -tranche 98.62 -tranche 98.64 -tranche 98.66 -tranche 98.68 -tranche 98.70 -tranche 98.72 -tranche 98.74 -tranche 98.76 -tranche 98.78 -tranche 98.80 -tranche 98.82 -tranche 98.84 -tranche 98.86 -tranche 98.88 -tranche 98.90 -tranche 98.92 -tranche 98.94 -tranche 98.96 -tranche 98.98 -tranche 99.00 -tranche 99.02 -tranche 99.04 -tranche 99.06 -tranche 99.08 -tranche 99.10 -tranche 99.12 -tranche 99.14 -tranche 99.16 -tranche 99.18 -tranche 99.20 -tranche 99.22 -tranche 99.24 -tranche 99.26 -tranche 99.28 -tranche 99.30 -tranche 99.32 -tranche 99.34 -tranche 99.36 -tranche 99.38 -tranche 99.40 -tranche 99.42 -tranche 99.44 -tranche 99.46 -tranche 99.48 -tranche 99.50 -tranche 99.52 -tranche 99.54 -tranche 99.56 -tranche 99.58 -tranche 99.60 -tranche 99.62 -tranche 99.64 -tranche 99.66 -tranche 99.68 -tranche 99.70 -tranche 99.72 -tranche 99.74 -tranche 99.76 -tranche 99.78 -tranche 99.80 -tranche 99.82 -tranche 99.84 -tranche 99.86 -tranche 99.88 -tranche 99.90 -tranche 99.92 -tranche 99.94 -tranche 99.96 -tranche 99.98 -tranche 100.00

	;;

	3)
	echo "2.8"
	/software/jre1.7.0_25/bin/java -Xmx12000m -Xms12000m -Xss280m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-2.8/GenomeAnalysisTK.jar \
	-T VariantRecalibrator -U LENIENT_VCF_PROCESSING --mode SNP -l INFO -R $REF -input $input \
	--ignore_filter LowQual \
	-recalFile $OUTF/vqsr.sites.snps.vcf \
	-tranchesFile $OUTF/vqsr.sites.snps.tranches \
	--rscript_file $OUTF/snp_VQSR.r \
	--maxGaussians 6 \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/hapmap_3.3.b37.vcf \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/1000G_omni2.5.b37.vcf \
	-resource:1000g,known=false,training=true,truth=false,prior=10.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/1000G_phase1.snps.high_confidence.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.8/dbsnp_138.b37.excluding_sites_after_129.vcf \
	-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff -an DP \
	--target_titv 2.15 \
	-tranche 70.00 -tranche 71.00 -tranche 72.00 -tranche 73.00 -tranche 74.00 -tranche 75.00 -tranche 76.00 -tranche 77.00 -tranche 78.00 -tranche 79.00 -tranche 80.00 -tranche 81.00 -tranche 82.00 -tranche 83.00 -tranche 84.00 -tranche 85.00 -tranche 86.00 -tranche 87.00 -tranche 88.00 -tranche 89.00 -tranche 90.00 -tranche 91.00 -tranche 92.00 -tranche 93.00 -tranche 94.00 -tranche 95.00 -tranche 96.00 -tranche 96.02 -tranche 96.04 -tranche 96.06 -tranche 96.08 -tranche 96.10 -tranche 96.12 -tranche 96.14 -tranche 96.16 -tranche 96.18 -tranche 96.20 -tranche 96.22 -tranche 96.24 -tranche 96.26 -tranche 96.28 -tranche 96.30 -tranche 96.32 -tranche 96.34 -tranche 96.36 -tranche 96.38 -tranche 96.40 -tranche 96.42 -tranche 96.44 -tranche 96.46 -tranche 96.48 -tranche 96.50 -tranche 96.52 -tranche 96.54 -tranche 96.56 -tranche 96.58 -tranche 96.60 -tranche 96.62 -tranche 96.64 -tranche 96.66 -tranche 96.68 -tranche 96.70 -tranche 96.72 -tranche 96.74 -tranche 96.76 -tranche 96.78 -tranche 96.80 -tranche 96.82 -tranche 96.84 -tranche 96.86 -tranche 96.88 -tranche 96.90 -tranche 96.92 -tranche 96.94 -tranche 96.96 -tranche 96.98 -tranche 97.00 -tranche 97.02 -tranche 97.04 -tranche 97.06 -tranche 97.08 -tranche 97.10 -tranche 97.12 -tranche 97.14 -tranche 97.16 -tranche 97.18 -tranche 97.20 -tranche 97.22 -tranche 97.24 -tranche 97.26 -tranche 97.28 -tranche 97.30 -tranche 97.32 -tranche 97.34 -tranche 97.36 -tranche 97.38 -tranche 97.40 -tranche 97.42 -tranche 97.44 -tranche 97.46 -tranche 97.48 -tranche 97.50 -tranche 97.52 -tranche 97.54 -tranche 97.56 -tranche 97.58 -tranche 97.60 -tranche 97.62 -tranche 97.64 -tranche 97.66 -tranche 97.68 -tranche 97.70 -tranche 97.72 -tranche 97.74 -tranche 97.76 -tranche 97.78 -tranche 97.80 -tranche 97.82 -tranche 97.84 -tranche 97.86 -tranche 97.88 -tranche 97.90 -tranche 97.92 -tranche 97.94 -tranche 97.96 -tranche 97.98 -tranche 98.00 -tranche 98.02 -tranche 98.04 -tranche 98.06 -tranche 98.08 -tranche 98.10 -tranche 98.12 -tranche 98.14 -tranche 98.16 -tranche 98.18 -tranche 98.20 -tranche 98.22 -tranche 98.24 -tranche 98.26 -tranche 98.28 -tranche 98.30 -tranche 98.32 -tranche 98.34 -tranche 98.36 -tranche 98.38 -tranche 98.40 -tranche 98.42 -tranche 98.44 -tranche 98.46 -tranche 98.48 -tranche 98.50 -tranche 98.52 -tranche 98.54 -tranche 98.56 -tranche 98.58 -tranche 98.60 -tranche 98.62 -tranche 98.64 -tranche 98.66 -tranche 98.68 -tranche 98.70 -tranche 98.72 -tranche 98.74 -tranche 98.76 -tranche 98.78 -tranche 98.80 -tranche 98.82 -tranche 98.84 -tranche 98.86 -tranche 98.88 -tranche 98.90 -tranche 98.92 -tranche 98.94 -tranche 98.96 -tranche 98.98 -tranche 99.00 -tranche 100.00
	# -tranche 70.00 -tranche 71.00 -tranche 72.00 -tranche 73.00 -tranche 74.00 -tranche 75.00 -tranche 76.00 -tranche 77.00 -tranche 78.00 -tranche 79.00 -tranche 80.00 -tranche 81.00 -tranche 82.00 -tranche 83.00 -tranche 84.00 -tranche 85.00 -tranche 86.00 -tranche 87.00 -tranche 88.00 -tranche 89.00 -tranche 90.00 -tranche 91.00 -tranche 92.00 -tranche 93.00 -tranche 94.00 -tranche 95.00 -tranche 96.00 -tranche 96.02 -tranche 96.04 -tranche 96.06 -tranche 96.08 -tranche 96.10 -tranche 96.12 -tranche 96.14 -tranche 96.16 -tranche 96.18 -tranche 96.20 -tranche 96.22 -tranche 96.24 -tranche 96.26 -tranche 96.28 -tranche 96.30 -tranche 96.32 -tranche 96.34 -tranche 96.36 -tranche 96.38 -tranche 96.40 -tranche 96.42 -tranche 96.44 -tranche 96.46 -tranche 96.48 -tranche 96.50 -tranche 96.52 -tranche 96.54 -tranche 96.56 -tranche 96.58 -tranche 96.60 -tranche 96.62 -tranche 96.64 -tranche 96.66 -tranche 96.68 -tranche 96.70 -tranche 96.72 -tranche 96.74 -tranche 96.76 -tranche 96.78 -tranche 96.80 -tranche 96.82 -tranche 96.84 -tranche 96.86 -tranche 96.88 -tranche 96.90 -tranche 96.92 -tranche 96.94 -tranche 96.96 -tranche 96.98 -tranche 97.00 -tranche 97.02 -tranche 97.04 -tranche 97.06 -tranche 97.08 -tranche 97.10 -tranche 97.12 -tranche 97.14 -tranche 97.16 -tranche 97.18 -tranche 97.20 -tranche 97.22 -tranche 97.24 -tranche 97.26 -tranche 97.28 -tranche 97.30 -tranche 97.32 -tranche 97.34 -tranche 97.36 -tranche 97.38 -tranche 97.40 -tranche 97.42 -tranche 97.44 -tranche 97.46 -tranche 97.48 -tranche 97.50 -tranche 97.52 -tranche 97.54 -tranche 97.56 -tranche 97.58 -tranche 97.60 -tranche 97.62 -tranche 97.64 -tranche 97.66 -tranche 97.68 -tranche 97.70 -tranche 97.72 -tranche 97.74 -tranche 97.76 -tranche 97.78 -tranche 97.80 -tranche 97.82 -tranche 97.84 -tranche 97.86 -tranche 97.88 -tranche 97.90 -tranche 97.92 -tranche 97.94 -tranche 97.96 -tranche 97.98 -tranche 98.00 -tranche 98.02 -tranche 98.04 -tranche 98.06 -tranche 98.08 -tranche 98.10 -tranche 98.12 -tranche 98.14 -tranche 98.16 -tranche 98.18 -tranche 98.20 -tranche 98.22 -tranche 98.24 -tranche 98.26 -tranche 98.28 -tranche 98.30 -tranche 98.32 -tranche 98.34 -tranche 98.36 -tranche 98.38 -tranche 98.40 -tranche 98.42 -tranche 98.44 -tranche 98.46 -tranche 98.48 -tranche 98.50 -tranche 98.52 -tranche 98.54 -tranche 98.56 -tranche 98.58 -tranche 98.60 -tranche 98.62 -tranche 98.64 -tranche 98.66 -tranche 98.68 -tranche 98.70 -tranche 98.72 -tranche 98.74 -tranche 98.76 -tranche 98.78 -tranche 98.80 -tranche 98.82 -tranche 98.84 -tranche 98.86 -tranche 98.88 -tranche 98.90 -tranche 98.92 -tranche 98.94 -tranche 98.96 -tranche 98.98 -tranche 99.00 -tranche 99.02 -tranche 99.04 -tranche 99.06 -tranche 99.08 -tranche 99.10 -tranche 99.12 -tranche 99.14 -tranche 99.16 -tranche 99.18 -tranche 99.20 -tranche 99.22 -tranche 99.24 -tranche 99.26 -tranche 99.28 -tranche 99.30 -tranche 99.32 -tranche 99.34 -tranche 99.36 -tranche 99.38 -tranche 99.40 -tranche 99.42 -tranche 99.44 -tranche 99.46 -tranche 99.48 -tranche 99.50 -tranche 99.52 -tranche 99.54 -tranche 99.56 -tranche 99.58 -tranche 99.60 -tranche 99.62 -tranche 99.64 -tranche 99.66 -tranche 99.68 -tranche 99.70 -tranche 99.72 -tranche 99.74 -tranche 99.76 -tranche 99.78 -tranche 99.80 -tranche 99.82 -tranche 99.84 -tranche 99.86 -tranche 99.88 -tranche 99.90 -tranche 99.92 -tranche 99.94 -tranche 99.96 -tranche 99.98 -tranche 100.00
	
	;;

	4)
	echo "2.5 dbSnp138_excl; No InbreedingCoeff"
	/software/jre1.7.0_25/bin/java -Xmx12000m -Xms12000m -Xss280m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-2.5/GenomeAnalysisTK.jar \
	-T VariantRecalibrator -U LENIENT_VCF_PROCESSING --mode SNP -l INFO -R $REF -input $input \
	--ignore_filter LowQual \
	-recalFile $OUTF/vqsr.sites.snps.vcf \
	-tranchesFile $OUTF/vqsr.sites.snps.tranches \
	--rscript_file $OUTF/snp_VQSR.r \
	--maxGaussians 6 \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/hapmap_3.3.b37.vcf \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/1000G_omni2.5.b37.vcf \
	-resource:1000g,known=false,training=true,truth=false,prior=10.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.5/1000G_phase1.snps.high_confidence.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.8/dbsnp_138.b37.excluding_sites_after_129.vcf \
	-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an DP \
	--target_titv 2.15 --ts_filter_level 96.00 \
	-minNumBad 5000 \
	-percentBad 0.05 \
	-tranche 70.00 -tranche 71.00 -tranche 72.00 -tranche 73.00 -tranche 74.00 -tranche 75.00 -tranche 76.00 -tranche 77.00 -tranche 78.00 -tranche 79.00 -tranche 80.00 -tranche 81.00 -tranche 82.00 -tranche 83.00 -tranche 84.00 -tranche 85.00 -tranche 86.00 -tranche 87.00 -tranche 88.00 -tranche 89.00 -tranche 90.00 -tranche 91.00 -tranche 92.00 -tranche 93.00 -tranche 94.00 -tranche 95.00 -tranche 96.00 -tranche 96.50 -tranche 96.52 -tranche 96.54 -tranche 96.56 -tranche 96.58 -tranche 96.60 -tranche 96.62 -tranche 96.64 -tranche 96.66 -tranche 96.68 -tranche 96.70 -tranche 96.72 -tranche 96.74 -tranche 96.76 -tranche 96.78 -tranche 96.80 -tranche 96.82 -tranche 96.84 -tranche 96.86 -tranche 96.88 -tranche 96.90 -tranche 96.92 -tranche 96.94 -tranche 96.96 -tranche 96.98 -tranche 97.00 -tranche 97.02 -tranche 97.04 -tranche 97.06 -tranche 97.08 -tranche 97.10 -tranche 97.12 -tranche 97.14 -tranche 97.16 -tranche 97.18 -tranche 97.20 -tranche 97.22 -tranche 97.24 -tranche 97.26 -tranche 97.28 -tranche 97.30 -tranche 97.32 -tranche 97.34 -tranche 97.36 -tranche 97.38 -tranche 97.40 -tranche 97.42 -tranche 97.44 -tranche 97.46 -tranche 97.48 -tranche 97.50 -tranche 97.52 -tranche 97.54 -tranche 97.56 -tranche 97.58 -tranche 97.60 -tranche 97.62 -tranche 97.64 -tranche 97.66 -tranche 97.68 -tranche 97.70 -tranche 97.72 -tranche 97.74 -tranche 97.76 -tranche 97.78 -tranche 97.80 -tranche 97.82 -tranche 97.84 -tranche 97.86 -tranche 97.88 -tranche 97.90 -tranche 97.92 -tranche 97.94 -tranche 97.96 -tranche 97.98 -tranche 98.00 -tranche 98.02 -tranche 98.04 -tranche 98.06 -tranche 98.08 -tranche 98.10 -tranche 98.12 -tranche 98.14 -tranche 98.16 -tranche 98.18 -tranche 98.20 -tranche 98.22 -tranche 98.24 -tranche 98.26 -tranche 98.28 -tranche 98.30 -tranche 98.32 -tranche 98.34 -tranche 98.36 -tranche 98.38 -tranche 98.40 -tranche 98.42 -tranche 98.44 -tranche 98.46 -tranche 98.48 -tranche 98.50 -tranche 98.52 -tranche 98.54 -tranche 98.56 -tranche 98.58 -tranche 98.60 -tranche 98.62 -tranche 98.64 -tranche 98.66 -tranche 98.68 -tranche 98.70 -tranche 98.72 -tranche 98.74 -tranche 98.76 -tranche 98.78 -tranche 98.80 -tranche 98.82 -tranche 98.84 -tranche 98.86 -tranche 98.88 -tranche 98.90 -tranche 98.92 -tranche 98.94 -tranche 98.96 -tranche 98.98 -tranche 99.00 -tranche 99.50 -tranche 100.00
	# -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an InbreedingCoeff -an DP \
	# -tranche 70.00 -tranche 71.00 -tranche 72.00 -tranche 73.00 -tranche 74.00 -tranche 75.00 -tranche 76.00 -tranche 77.00 -tranche 78.00 -tranche 79.00 -tranche 80.00 -tranche 81.00 -tranche 82.00 -tranche 83.00 -tranche 84.00 -tranche 85.00 -tranche 86.00 -tranche 87.00 -tranche 88.00 -tranche 89.00 -tranche 90.00 -tranche 91.00 -tranche 92.00 -tranche 93.00 -tranche 94.00 -tranche 95.00 -tranche 96.00 -tranche 96.02 -tranche 96.04 -tranche 96.06 -tranche 96.08 -tranche 96.10 -tranche 96.12 -tranche 96.14 -tranche 96.16 -tranche 96.18 -tranche 96.20 -tranche 96.22 -tranche 96.24 -tranche 96.26 -tranche 96.28 -tranche 96.30 -tranche 96.32 -tranche 96.34 -tranche 96.36 -tranche 96.38 -tranche 96.40 -tranche 96.42 -tranche 96.44 -tranche 96.46 -tranche 96.48 -tranche 96.50 -tranche 96.52 -tranche 96.54 -tranche 96.56 -tranche 96.58 -tranche 96.60 -tranche 96.62 -tranche 96.64 -tranche 96.66 -tranche 96.68 -tranche 96.70 -tranche 96.72 -tranche 96.74 -tranche 96.76 -tranche 96.78 -tranche 96.80 -tranche 96.82 -tranche 96.84 -tranche 96.86 -tranche 96.88 -tranche 96.90 -tranche 96.92 -tranche 96.94 -tranche 96.96 -tranche 96.98 -tranche 97.00 -tranche 97.02 -tranche 97.04 -tranche 97.06 -tranche 97.08 -tranche 97.10 -tranche 97.12 -tranche 97.14 -tranche 97.16 -tranche 97.18 -tranche 97.20 -tranche 97.22 -tranche 97.24 -tranche 97.26 -tranche 97.28 -tranche 97.30 -tranche 97.32 -tranche 97.34 -tranche 97.36 -tranche 97.38 -tranche 97.40 -tranche 97.42 -tranche 97.44 -tranche 97.46 -tranche 97.48 -tranche 97.50 -tranche 97.52 -tranche 97.54 -tranche 97.56 -tranche 97.58 -tranche 97.60 -tranche 97.62 -tranche 97.64 -tranche 97.66 -tranche 97.68 -tranche 97.70 -tranche 97.72 -tranche 97.74 -tranche 97.76 -tranche 97.78 -tranche 97.80 -tranche 97.82 -tranche 97.84 -tranche 97.86 -tranche 97.88 -tranche 97.90 -tranche 97.92 -tranche 97.94 -tranche 97.96 -tranche 97.98 -tranche 98.00 -tranche 98.02 -tranche 98.04 -tranche 98.06 -tranche 98.08 -tranche 98.10 -tranche 98.12 -tranche 98.14 -tranche 98.16 -tranche 98.18 -tranche 98.20 -tranche 98.22 -tranche 98.24 -tranche 98.26 -tranche 98.28 -tranche 98.30 -tranche 98.32 -tranche 98.34 -tranche 98.36 -tranche 98.38 -tranche 98.40 -tranche 98.42 -tranche 98.44 -tranche 98.46 -tranche 98.48 -tranche 98.50 -tranche 98.52 -tranche 98.54 -tranche 98.56 -tranche 98.58 -tranche 98.60 -tranche 98.62 -tranche 98.64 -tranche 98.66 -tranche 98.68 -tranche 98.70 -tranche 98.72 -tranche 98.74 -tranche 98.76 -tranche 98.78 -tranche 98.80 -tranche 98.82 -tranche 98.84 -tranche 98.86 -tranche 98.88 -tranche 98.90 -tranche 98.92 -tranche 98.94 -tranche 98.96 -tranche 98.98 -tranche 99.00 -tranche 99.02 -tranche 99.04 -tranche 99.06 -tranche 99.08 -tranche 99.10 -tranche 99.12 -tranche 99.14 -tranche 99.16 -tranche 99.18 -tranche 99.20 -tranche 99.22 -tranche 99.24 -tranche 99.26 -tranche 99.28 -tranche 99.30 -tranche 99.32 -tranche 99.34 -tranche 99.36 -tranche 99.38 -tranche 99.40 -tranche 99.42 -tranche 99.44 -tranche 99.46 -tranche 99.48 -tranche 99.50 -tranche 99.52 -tranche 99.54 -tranche 99.56 -tranche 99.58 -tranche 99.60 -tranche 99.62 -tranche 99.64 -tranche 99.66 -tranche 99.68 -tranche 99.70 -tranche 99.72 -tranche 99.74 -tranche 99.76 -tranche 99.78 -tranche 99.80 -tranche 99.82 -tranche 99.84 -tranche 99.86 -tranche 99.88 -tranche 99.90 -tranche 99.92 -tranche 99.94 -tranche 99.96 -tranche 99.98 -tranche 100.00
	# -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.8/dbsnp_138.b37.vcf \
	;;
	5)
	echo "Latest GATK"
	java -Xmx4800m -Xms4800m -Xss280m -server -XX:+UseSerialGC -jar /software/hgi/pkglocal/gatk-protected-3.3/GenomeAnalysisTK.jar \
	-T VariantRecalibrator --mode SNP -l INFO -U LENIENT_VCF_PROCESSING -R $REF -input $input \
	--ignore_filter LowQual \
	-recalFile $OUTF/vqsr.sites.snps.vcf \
	-tranchesFile $OUTF/vqsr.sites.snps.tranches \
	--rscript_file $OUTF/snp_VQSR.r \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /lustre/scratch114/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.5/hapmap_3.3.b37.vcf \
	-resource:omni,known=false,training=true,truth=false,prior=12.0 /lustre/scratch114/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.5/1000G_omni2.5.b37.vcf \
	-resource:1000g,known=false,training=true,truth=false,prior=10.0 /lustre/scratch114/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.5/1000G_phase1.snps.high_confidence.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch114/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
	--target_titv 2.3 \
	--maxGaussians 4 \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
	-tranche 99.00 -tranche 99.02 -tranche 99.04 -tranche 99.06 -tranche 99.08 -tranche 99.10 -tranche 99.12 -tranche 99.14 -tranche 99.16 -tranche 99.18 -tranche 99.20 -tranche 99.22 -tranche 99.24 -tranche 99.26 -tranche 99.28 -tranche 99.30 -tranche 99.32 -tranche 99.34 -tranche 99.36 -tranche 99.38 -tranche 99.40 -tranche 99.42 -tranche 99.44 -tranche 99.46 -tranche 99.48 -tranche 99.50 -tranche 99.52 -tranche 99.54 -tranche 99.56 -tranche 99.58 -tranche 99.60 -tranche 99.62 -tranche 99.64 -tranche 99.66 -tranche 99.68 -tranche 99.70 -tranche 99.72 -tranche 99.74 -tranche 99.76 -tranche 99.78 -tranche 99.80 -tranche 99.82 -tranche 99.84 -tranche 99.86 -tranche 99.88 -tranche 99.90 -tranche 99.92 -tranche 99.94 -tranche 99.96 -tranche 99.98 -tranche 99.99 -tranche 100.00
	;;
	*)
	echo "Invalid mode selected"
	;;
esac

#compress and index
bgzip $OUTF/vqsr.sites.snps.vcf
tabix -p vcf $OUTF/vqsr.sites.snps.vcf.gz

