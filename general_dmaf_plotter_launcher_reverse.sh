#!/usr/local/bin/bash

#wrapper for launch general_dmaf_plotter.R script
#this works in a different way: this take all chr increase files, than concat them, split them in freq bins and plot dmaf.
POP1=$1
POP2=$2
FILEPATH=$3
sign=$4

#first concat files for each range and save them in another folder
mkdir -p $3/ALL_REVERSE
for i in {1..22}
do
	fgrep -v "CHROM" VBI_vs_EUR_CHR_${i}_$sign.txt >> $3/ALL_REVERSE/VBI_vs_EUR_all_$sign.txt
done

#use these files to calculate freq bins
files=`ls $3/ALL_REVERSE/*`

#change working directory
cd $3/ALL_REVERSE/

for file in $files
do
filename=`basename $file`

bsub -J "general_dmaf_$file" -o "%J_general_dmaf_$file.log" -M8000000 -R "select[mem>8000] rusage[mem=8000]" \
-q normal R CMD BATCH '--args '$1' '$2' '$filename' '$sign'' /nfs/users/nfs_m/mc14/Work/r_scripts/general_dmaf_plotter.R
done
