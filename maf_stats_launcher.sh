#!/usr/local/bin/bash

#launcher for maf stats

if [ $# -lt 1 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "maf_stats_launcher.sh <maf_filename> <pop name> \n"

	exit 1
fi

# submit the job
bsub -J "maf_stats" -o "%J_maf_stats.log" -M8000 -R "select[mem>8000] rusage[mem=8000]" \
-q yesterday R CMD BATCH '--args '$1' '$2'' /nfs/users/nfs_m/mc14/Work/r_scripts/maf_stats.r
