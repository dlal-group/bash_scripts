#/!usr/local/bin/bash
#wrapper for reannotation check for all chr
LOGS=/nfs/users/nfs_m/mc14/lustre110_home/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/LOGS

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	bsub -J "re_ann_diff_check.chr${i}" -o "$LOGS/%J_re_ann_diff_check.chr${i}.log" -e "$LOGS/%J_re_ann_diff_check.chr${i}.err" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q normal reannotate_check.sh ${i}

done
