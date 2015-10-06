#!/usr/local/bin/bash
#script to train a new vqsr filter model:
#############################
#     SNPs and INDELSs!!	#
#############################

#The script has to be run using the ja_runner_par.sh wrapper, to create a job array starting from a file list, where files are listed with the absolute path
input=$1
input_name=`basename ${input}`
OUTF=$2
thr=$3
recal_trance_path=$4

REF=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
echo "Applying recalibration filter with parameters:"
echo -e "Reference file:$REF\nInput file: ${input}\nOutput folder: $OUTF\nTreshold for SNPs call: ${thr}\nTreshold for INDELs call: 90.00 (fixed)\nRecall and tranches path:${recal_trance_path}"

#apply the recalibration for SNPS
echo "SNP mode ..."

# java -Xmx4800m -Xms4800m -Xss280m -server -XX:+UseSerialGC -jar /software/hgi/pkglocal/gatk-protected-3.3/GenomeAnalysisTK.jar \
# -T ApplyRecalibration -U LENIENT_VCF_PROCESSING \
# --ts_filter_level ${thr} \
# --mode SNP \
# -R $REF \
# --input ${input} \
# -recalFile ${recal_trance_path}/vqsr.sites.snps.vcf.gz \
# -tranchesFile ${recal_trance_path}/vqsr.sites.snps.tranches \
# -o $OUTF/${input_name}.vqsr_snp.vcf.gz

# tabix -p vcf -f $OUTF/${input_name}.vqsr_snp.vcf.gz


# echo "INDELs mode ..."
# #apply the recalibration for INDELs
# java -Xmx4800m -Xms4800m -Xss280m -server -XX:+UseSerialGC -jar /software/hgi/pkglocal/gatk-protected-3.3/GenomeAnalysisTK.jar \
# -T ApplyRecalibration -U LENIENT_VCF_PROCESSING \
# --ts_filter_level 95.00 \
# --mode INDEL \
# -R $REF \
# --input $OUTF/${input_name}.vqsr_snp.vcf.gz \
# -recalFile ${recal_trance_path}/vqsr.sites.indels.vcf.gz \
# -tranchesFile ${recal_trance_path}/vqsr.sites.indels.tranches \
# -o $OUTF/${input_name}.vqsr.vcf.gz

# tabix -p vcf -f $OUTF/${input_name}.vqsr.vcf.gz

# bcftools concat -a $OUTF/${input_name}.vqsr_snp.vcf.gz $OUTF/${input_name}.vqsr_indel.vcf.gz | vcf-sort | bgzip -c > $OUTF/${input_name}.vqsr_appl.vcf.gz
bcftools filter -i'FILTER=PASS' -O z -o $OUTF/${input_name}.vqsr.filt.vcf.gz $OUTF/${input_name}.vqsr.vcf.gz

tabix -p vcf -f $OUTF/${input_name}.vqsr.filt.vcf.gz
bcftools stats -d 0,5000,1 -s - $OUTF/${input_name}.vqsr.vcf.gz > $OUTF/${input_name}.vqsr.vcf.gz.vchk
bcftools stats -d 0,5000,1 -s - $OUTF/${input_name}.vqsr.filt.vcf.gz > $OUTF/${input_name}.vqsr.filt.vcf.gz.vchk
