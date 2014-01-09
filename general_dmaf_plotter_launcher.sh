#!/usr/local/bin/bash

#wrapper for launch general_dmaf_plotter.R script

POP1=$1
POP2=$2
FILEPATH=$3
sign=$4

#first concat files for each range and save them in another folder
mkdir -p $3/ALL
for class in lt02 lte05 lte10 lte20 lte30 lte40 lte50
do
	for i in {1..22}
	do
		cat VBI_vs_EUR_CHR_${i}_$sign.txt.${class} >> $3/ALL/VBI_vs_EUR_all_$sign.txt.${class}
	done
done

#use these files to calculate freq bins
files=`ls $3/ALL/*`

#change working directory
cd $3/ALL

for file in $files
do
filename=`basename $file`

bsub -J "general_dmaf_$file" -o "%J_general_dmaf_$file.log" -M8000000 -R "select[mem>8000] rusage[mem=8000]" \
-q normal R CMD BATCH '--args '$1' '$2' '$filename' '$sign'' /nfs/users/nfs_m/mc14/Work/r_scripts/general_dmaf_plotter.R
done
