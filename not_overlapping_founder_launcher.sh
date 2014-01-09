 #!/usr/local/bin/bash
 
#script to find not overlapping sites
#for each chr we creare a file with not overlapping entries
for i in 1 2 6 9 10 12
do
	bsub -J "not_overlapping_chr${i}" -o "%J_not_overlapping_chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q basement not_overlapping_founder.sh ${i}
done
