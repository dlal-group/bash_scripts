#/!usr/local/bin/bash
#script for vcf splitting in 22 chr

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
bsub -J "split_ann_file.chr${i}" -o "%J_split_ann_file.chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q normal split_dbSNP_ann_in_chr.sh ${i}
done
