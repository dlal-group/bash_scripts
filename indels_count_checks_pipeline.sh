#!/usr/local/bin/bash

#Pipeline used to extract info about indels

#ARGS:
vcf_file=$1
out_path=$2

if [ $# -lt 2 ]
then
	echo "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
	echo "Usage:"
	echo "All paths must been ABSOLUTE!!"
	echo "indels_count_checks_pipeline.sh <vcf_file> <output_path>"
	echo "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"

fi

#FIRST STEP:
#Count number of indels for each chr and separate knonwn from novel

bsub -J "indels_check[1-23]" -o "%J_indels_check.log" -e "%J_indels_check.err" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement /nfs/users/nfs_m/mc14/Work/bash_scripts/indels_count_checks_step1.sh $vcf_file $out_path
