#!/usr/local/bin/bash
#Launch vqsrlod plotter

echo "Starting plotter script..."
#-J job_name_spec[index |  start_index-end_index:step,]

#command to extract only vqslod column for plotting..need to be done BEFORE running the script...
#first done on SNPS_FILES content files
#for i in {1..22};do cut -f 1,2,3,11 SNPS_FILES/VBI_chr${i}.snps.csv > SNPS_FILES/VBI_chr${i}.snps.vqslod.csv;done

#then also for MULTIALLELIC folder files..to speed up R upload!!
#for j in {1..22};do cut -f 1,2,3,9 MULTIALLELIC/multiallelic_chr"${j}"_table.csv > MULTIALLELIC/multiallelic_chr"${j}"_table.vqslod.csv;done

bsub -J "vqsrlod_plotter" -o "%J_vqsrlod_out.log"  -M32000000 -R"select[mem>32000] rusage[mem=32000]" \
-q hugemem R CMD BATCH /nfs/users/nfs_m/mc14/Work/r_scripts/vqsrlod_distributions.r

