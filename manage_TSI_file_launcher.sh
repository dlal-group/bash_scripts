#!/usr/local/bin/bash

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
#bsub -J "chr${i}_manipulation" -o "%J_chr${i}_manipulation.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
#-q basement manage_TSI_file.sh ${i}
echo "CHR ${i}"
tabix TSI.chr${i}.annotate.recode.snps.vcf.gz chr ${i} | wc -l
done

