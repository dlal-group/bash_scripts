#!/usr/local/bin/bash
#launcher script for reannotation check on overlap between seq data and 370k

logs='/nfs/users/nfs_m/mc14/lustre110_home/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/LOGS'
file_path=$1
out_path=$2

if [ $# -lt 2 ]
then
	echo "Wrong argument number!!"
	echo "Usage: annotation_check_launch.sh <seq_file_path> <out_path>"
	echo "(the file for the overlap with gwas data is assumed to be in the current folder and named \'biallelic_overl.SNP.unfilt.geno.seq.VB.map\' or \'all.chr.geno.overlapped_POS.sorted.txt\' if you are checking the complete overlap (not only biallelic sites)"
	echo ""
	exit 1
fi


for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
	bsub -J "overlap_ann_check.chr${i}" -o "$logs/%J_overlap_ann_check.chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q normal annotation_check.sh $i $file_path $out_path
done 
