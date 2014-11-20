#!/usr/local/bin/bash

#script for manhattan plot generation
#take as arguments:
#$1: complete chr map file path
#$2: partial chr map file path

for i in {1..22}
do
all_chr=`ls $1/*chr${i}.*map`
partial_chr=`ls $2/*chr${i}.*map`
bsub -J "manhattan_plot_chr${i}" -o "%J_manhattan_plot_chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement \
R CMD BATCH "--args $all_chr $partial_chr ${i}" /nfs/users/nfs_m/mc14/Work/r_scripts/plot_chr_position_densities.r
done
