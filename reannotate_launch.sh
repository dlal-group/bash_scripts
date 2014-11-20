#/!usr/local/bin/bash
#script for reannotation of vcf for all chr
#manually set the ABSOLUTE output path
if [ $# -lt 2 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "reannotate_launch.sh <input_path> <output-path> [no_all_match]\n"
	echo -e "<input_path>: Path for the files to annotate."
	echo -e "<output-path>:Desired output path\n"
	echo -e "no_all_match:Add this option if you want to annotate without performing allele match between vcf and annotation file.\n"
	exit 1
fi

input_path=$1
out_path=$2

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	if [ $# -eq 3 ]
	then
		bsub -J "re_ann_file_no_all_chk.chr${i}" -o "LOGS/%J_re_ann_file_no_all_chk.chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q normal reannotate.sh ${i} $out_path $input_path $2

	else
		bsub -J "re_ann_file.chr${i}" -o "LOGS/%J_re_ann_file.chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q normal reannotate.sh ${i} $out_path $input_path

	fi
done
