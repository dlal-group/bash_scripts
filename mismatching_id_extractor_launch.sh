#!/usr/local/bin/bash

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	while read line
		do echo "RIGA:"
		POS=`echo $line | cut -f 2 -d " "`
		bsub -J "chr${i}_extractor" -o "%J_chr${i}_extractor.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
		-q normal mismatching_id_extractor.sh ${i} $POS
	done < no_rsID_chr${i}.tab
done

