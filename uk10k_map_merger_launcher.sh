#!/usr/local/bin/bash
#launcher for uk10k_map_merger.sh script: works by chromosome
#Args:
#$1=input_1_path
#$2=input_2_path
#$3=output_path


if [ $# -lt 3 ]
then
        echo "----------------------------------------------------------------------------------------"
        echo "WRONG argument number!"
        echo "Usage: uk10k_map_merger.sh <input_1_path> <input_2_path> <output_path>"
        echo "----------------------------------------------------------------------------------------"
        exit 1
fi

input_1_path=$1
input_2_path=$2
output_path=$3

#for uk10k there is not chrX
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
bsub -J "uk10k_merger_chr${i}" -o "%J_uk10k_merger_chr${i}log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q normal uk10k_map_merger.sh ${i} ${input_1_path} ${input_2_path} ${output_path}
done

