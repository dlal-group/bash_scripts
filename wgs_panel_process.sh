#!/usr/bin/env bash
#
# Script to process the selected files to create a reference panle
# We'll use this to:
# Simply MERGE VCF file -> there will be missing genotypes (phase 1)
# Phase and Impute the merged data -> to remove missing and to correct genotype errors (phase 1)
# To recall all data together from bam files using all sites selected from the merging step (phase 2)
# To refine genotypes calls after recalling (phase 3)
# To create merged reference panel in IMPUTE format for each phase

# vcf=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/22.vcf.gz
# cohort="VBI"
# outdir=/lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL/VBI
#ARGS:
# $1= vcf file for a single chromosome, better if named as "[chr].vcf.gz"
# $2=cohort
# $3=output folder
# $4=mode (snps/indels)

# /nfs/users/nfs_m/mc14/Work/bash_scripts/wgs_panel_process.sh /lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/22.vcf.gz VBI /lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL
set -e

vcf=$1
cohort=$2
outdir=$3/$2
stage=$4

mkdir -p ${outdir}/LOG_${stage}

filename=`basename ${vcf}`
first_suffix="${filename%%.*}"
chr=${first_suffix}

mkdir -p ${outdir}/${chr}
#we're going to split snps and indels, than put them back together again

case ${stage} in
VCF_MERGE )
# define path
;;
IMPUTE_FORMAT )
#first we need to create the legend file
# we need rsID, position, REF,ALT: for consistency we'll assign rsID=chr:pos:ALT to all sites with unknown rsID
echo "bcftools query -f\"%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n\" ${vcf} | awk '{if(\$3!=\".\") print \$3,\$2,\$4,\$5;else if (\$3==\".\") print \"chr\"\$1\":\"\$2\":\"\$5,\$2,\$4,\$5}'| gzip -c > ${outdir}/${chr}.INGI_REF.${cohort}.legend.gz" | bsub -J"extract_${mode}_${cohort}_legend_${first_suffix}" -o"${outdir}/LOG_${stage}/1.%J_extract_${mode}_${cohort}_legend_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal
#than we need to create the hap file we need to create the legend file
echo "bcftools query -f\"[%GT ]\\n\" ${vcf} | tr \"|\" \" \" | gzip -c > ${outdir}/${chr}.INGI_REF.${cohort}.hap.gz" | bsub -J"extract_${mode}_${cohort}_hap_${first_suffix}" -o"${outdir}/LOG_${stage}/2.%J_extract_${mode}_${cohort}_hap_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal


;;
PANEL_MERGE )
#this has to be an iterative process
# We need to tell which and how many cohorts we want to merge and merge them 2 by 2 using the output of the previous step
# to merge the subsequent cohort we can use two files:
	# -one with the cohort order
	# -one with the path for each cohort chr file
	cohorts=$1
	cohorts_files=$2

# ./impute2 -merge_ref_panels \
#  -m ./Example/example.chr22.map \
#  -h ./Example/example.chr22.1kG.haps \
#     ./Example/example.chr22.hm3.haps \
#  -l ./Example/example.chr22.1kG.legend \
#     ./Example/example.chr22.hm3.legend \
#  -merge_ref_panels_output_ref
#  -strand_g ./Example/example.chr22.study.strand \
#  -int 20.4e6 20.5e6 \
#  -Ne 20000 \
#  -o ./Example/example.chr22.two.phased.impute2



;;
esac
