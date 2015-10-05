#!/usr/local/bin/bash
#script to train a new vqsr filter model to test different solution and expand thresholds:
#########################
#     INDELs ONLY!!!		#
#########################
if [[ $# -lt 3 ]]; then
	echo "##############################################"
	echo "Attention!!! Not enough arguments provided!!!"
	echo "##############################################"
	echo "Usage:"
	echo "vqsr_gatk_recalibrate_indels.sh <input_file_path> <output_path> <mode>."
	echo -e "The 'mode' option is used to select the GATK version to use:\n1 - GATK v.2.5\n2 - GATK v.2.7\n3 - GATK v.2.8"
	exit 1
fi

input=$1
OUTF=$2
MODE=$3
REF=/lustre/scratch111/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa

mkdir -p $OUTF/LOGS

case $MODE in
	1)
	echo "2.5"
	/software/jre1.7.0_25/bin/java -Xmx12000m -Xms12000m -Xss280m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-2.5/GenomeAnalysisTK.jar \
	-T VariantRecalibrator -U LENIENT_VCF_PROCESSING --mode INDEL -l INFO -R $REF -input $input \
	--ignore_filter LowQual \
	-recalFile $OUTF/vqsr.sites.indels.vcf \
	-tranchesFile $OUTF/vqsr.sites.indels.tranches \
	--rscript_file $OUTF/indels_VQSR.r \
	-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.2/Mills_and_1000G_gold_standard.indels.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.8/dbsnp_138.b37.vcf \
	-an DP -an FS -an ReadPosRankSum -an MQRankSum \
	-minNumBad 5000 \
	-percentBad 0.05 \
	--target_titv 2.15 --ts_filter_level 90.00 \
	--maxGaussians 6 \
	-tranche 10 -tranche 15 -tranche 20 -tranche 25 -tranche 30 -tranche 35 -tranche 40 -tranche 45 -tranche 50 -tranche 55 -tranche 60 -tranche 65 -tranche 70 -tranche 75 -tranche 80 -tranche 85 -tranche 90 -tranche 95
	
	;;

	2)
	echo "2.7"
	/software/jre1.7.0_25/bin/java -Xmx12000m -Xms12000m -Xss280m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar \
	-T VariantRecalibrator -U LENIENT_VCF_PROCESSING --mode INDEL -l INFO -R $REF -input $input \
	--ignore_filter LowQual \
	-recalFile $OUTF/vqsr.sites.indels.vcf \
	-tranchesFile $OUTF/vqsr.sites.indels.tranches \
	--rscript_file $OUTF/indels_VQSR.r \
	-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.2/Mills_and_1000G_gold_standard.indels.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.8/dbsnp_138.b37.vcf \
	-an DP -an FS -an ReadPosRankSum -an MQRankSum \
	--target_titv 2.15 --ts_filter_level 90.00 \
	--numBadVariants 5000 \
	-tranche 10 -tranche 15 -tranche 20 -tranche 25 -tranche 30 -tranche 35 -tranche 40 -tranche 45 -tranche 50 -tranche 55 -tranche 60 -tranche 65 -tranche 70 -tranche 75 -tranche 80 -tranche 85 -tranche 90 -tranche 95

	;;

	3)
	echo "2.8"
	/software/jre1.7.0_25/bin/java -Xmx12000m -Xms12000m -Xss280m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-2.8/GenomeAnalysisTK.jar \
	-T VariantRecalibrator -U LENIENT_VCF_PROCESSING --mode INDEL -l INFO -R $REF -input $input \
	--ignore_filter LowQual \
	-recalFile $OUTF/vqsr.sites.indels.vcf \
	-tranchesFile $OUTF/vqsr.sites.indels.tranches \
	--rscript_file $OUTF/indels_VQSR.r \
	-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.2/Mills_and_1000G_gold_standard.indels.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch111/resources/variation/grch37/gatk-bundle/2.8/dbsnp_138.b37.vcf \
	-an DP -an FS -an ReadPosRankSum -an MQRankSum \
	--target_titv 2.15 \
	-tranche 10 -tranche 15 -tranche 20 -tranche 25 -tranche 30 -tranche 35 -tranche 40 -tranche 45 -tranche 50 -tranche 55 -tranche 60 -tranche 65 -tranche 70 -tranche 75 -tranche 80 -tranche 85 -tranche 90 -tranche 95
	
	;;
	4)
	echo "Latest GATK"
	java -Xmx4800m -Xms4800m -Xss280m -server -XX:+UseSerialGC -jar /software/hgi/pkglocal/gatk-protected-3.3/GenomeAnalysisTK.jar \
	-T VariantRecalibrator --mode INDEL -l INFO -U LENIENT_VCF_PROCESSING -R $REF -input $input \
	-recalFile $OUTF/vqsr.sites.indels.vcf \
	-tranchesFile $OUTF/vqsr.sites.indels.tranches \
	--rscript_file $OUTF/indels_VQSR.r \
	-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 /lustre/scratch114/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.5/Mills_and_1000G_gold_standard.indels.b37.vcf
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /lustre/scratch114/resources/variation/Homo_sapiens/grch37/gatk-bundle/2.8/b37/dbsnp_138.b37.vcf
	--target_titv 2.3
	--maxGaussians 4
	-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff
	-tranche 10 -tranche 15 -tranche 20 -tranche 25 -tranche 30 -tranche 35 -tranche 40 -tranche 45 -tranche 50 -tranche 55 -tranche 60 -tranche 65 -tranche 70 -tranche 75 -tranche 80 -tranche 85 -tranche 90 -tranche 95 -tranche 98 -tranche 100
	
	;;

	*)
	echo "Invalid mode selected"
	exit 1

	;;
esac

#compress and index
bgzip $OUTF/vqsr.sites.indels.vcf
tabix -p vcf $OUTF/vqsr.sites.indels.vcf.gz

