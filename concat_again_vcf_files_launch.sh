#!/usr/local/bin/bash
#launcher for concat_again_vcf_files.sh

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	bsub -J "chr${i}_concat" -o "%J_chr${i}_concat.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement concat_again_vcf_files.sh ${i}
done
