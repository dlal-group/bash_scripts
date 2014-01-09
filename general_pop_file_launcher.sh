#!/usr/local/bin/bash
#Generic script to launch CUSTOMIZED population file manipulator's scripts: works by chromosome

if [ $# -lt 3 ]
then
	echo "----------------------------------------------------------------------------------------"
	echo "WRONG argument number!"
	echo "Usage: general_pop_file_launcher.sh <pop_script_file_name> <input_path> <output_path>"
	echo "----------------------------------------------------------------------------------------"
	exit 1
fi

filename=$1
input_path=$2
output_path=$3

#for uk10k there is not chrX
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
bsub -J "chr${i}_manipulation" -o "%J_chr${i}_manipulation.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q normal ${filename} ${i} ${input_path} ${output_path}
done

