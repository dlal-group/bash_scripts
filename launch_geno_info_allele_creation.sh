#!/usr/local/bin/bash
#Launcher script for geno_info_allele_creation.sh
#
#Args:
#$1=geno path
#$2=output path

geno_path=$1
outpath=$2

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X_PAR1 X_PAR2 X_nonPAR
#for i in 22
do

        mkdir -p ${outpath}/LOGS
        mkdir -p ${outpath}/CHR${i}
	
	bsub -J "chr${chr}_gia" -o "$out_path/LOGS/%J_chr${chr}_gia.log" -e "$out_path/LOGS/%J_chr${chr}_gia.err" -G "team151" -M8000000 -R"select[mem>8000] rusage[mem=8000]" -q normal /nfs/users/nfs_m/mc14/Work/bash_scripts/geno_info_allele_creation.sh ${i} ${geno_path} ${outpath}

done

