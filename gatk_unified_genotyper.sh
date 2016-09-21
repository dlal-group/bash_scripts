#!/usr/local/bin/bash
#script to create files to use to train the VQSR filter:
#########################
#     SNPs and INDELs	#
#########################
if [[ $# -lt 3 ]]; then
	echo "##############################################"
	echo "Attention!!! Not enough arguments provided!!!"
	echo "##############################################"
	echo "Usage:"
	echo "vqsr_gatk_recalibrate.sh <input_file_path> <output_path> <mode>."
	echo -e "The 'mode' option is used to select the GATK version to use:\n1 - GATK v.2.5\n2 - GATK v.2.7\n3 - GATK v.2.8"
	exit 1
fi

input=$1
OUTF=$2
MODE=$3
REF=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa

mkdir -p $OUTF/LOGS

