#!/usr/local/bin/bash
#
#Script to split vcf files in rsId annotated vcf and missing rsID annotation vcf
#Written by Max_cocca

#Manually fixed paths
#seq_file_path='/lustre/scratch109/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel'
#splitted_seq_file_path='/lustre/scratch109/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/SPLITTED_VCF'

#User defined from command line path
seq_file_path=$1
splitted_seq_file_path=$2

if [ $# -lt 2 ]
then
	echo "ATTENTION!!Wrong arguments number!!!"
	echo " Usage: split_vcf_by_rsID_annotation_launch.sh <seq_file_path> <splitted_seq_file_path> [W]"
	echo "Options:"
	echo "seq_file_path: Path of sequences file (here will be created all output folders) [REQUIRED]"
	echo "splitted_file_path: Path for chr splitted vcf files [REQUIRED]"
	echo "W : Option to extract sites with rsID [OPTIONAL]"
	echo "--------------------------------------------------------------------------------------------"
	exit 1
fi

mkdir -p RS_ann_splitted

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do

	if [ $# -lt 3 ]
	then
		bsub -J "split_files_by_rs.chr${i}" -o "%J_split_files_by_rs.chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q normal split_vcf_by_rsID_annotation.sh ${i} $seq_file_path $splitted_seq_file_path
	else
		bsub -J "split_files_by_rs_w_rs.chr${i}" -o "%J_split_files_by_rs_w_rs.chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q normal split_vcf_by_rsID_annotation_w_rsID.sh ${i} $seq_file_path $splitted_seq_file_path
	fi
done
