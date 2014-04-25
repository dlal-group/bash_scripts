#!/usr/local/bin/bash

#Pipeline to create summary files for analyses for different porpouses:
# 1) Extract positive controls match for each trait
	# - using a master files with positive controls information for each traits we are going to:
	# 	- filter the master file based on phenotype
	# 	- filter result files based on master-file results
	# 	- create a summary table of positive controls matching with:
	# 		- CHR POS rsID P_repl beta se P_pubbl

# 2) Create files for plotting based on MAF and info filtering
	# - extract from the additional calculated stats sites that satisfy our filtering by maf and info
		# - MAF >=0.5%
		# - info >= 0.3
	# - generate qqplots and manhattan plot and calculate lambda for those variants (evaluate if a downsampling is needed to speed up computational time)

# 3) Create resume files based on MAF, info, HWE and Pvalue filtering
	# - filtering based on MAF,info and HWE from the additional stats
	# - extract those sites from the analyses results
	# - filter by Pvalue
###
###STEP1
#positive control files path = /lustre/scratch113/projects/uk10k/users/jh21/scratch3/positive_ctrls/uk10k.poscon

#example
#gwas_summary.sh FVG LDL /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/LIPIDS/NEW_ORDER/output /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/LIPIDS/NEW_ORDER/SUMMARIES

pop=$1
trait=$2
result_path=$3
out_path=$4
# poscon="/lustre/scratch113/projects/uk10k/users/jh21/scratch3/positive_ctrls/uk10k.poscon"
#customized file with chrX coded as X
poscon="/lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/uk10k.poscon"
# rs1976403
mkdir -p ${out_path}/${trait}/LOGS

#extract from trait results, the matching positive controls
#we'll work by chr
#Print command used:
# echo "extract_poscon.sh ${poscon} ${result_path} ${trait} ${out_path}/${trait}"
# echo "filter_results.sh ${add_stats_path} ${result_path} ${trait} ${out_path}/${trait}"
# echo "sort_results.sh ${add_stats_path} ${result_path} ${trait} ${out_path}/${trait}"
# echo "plotting_results.sh ${trait} ${out_path}/${trait} ${out_path}/${trait}/${trait}.all.poscon.join"

bsub -J "poscon_check_${pop}_${trait}[1-23]" -o "${out_path}/${trait}/LOGS/%J_poscon_check.%I.o" -M3000 -R"select[mem>3000] rusage[mem=3000]" -q small -- extract_poscon.sh ${poscon} ${result_path} ${trait} ${out_path}/${trait}

#now create files for plotting and plot them
#we also want a list of the pos controls to higligth on the plots
#we can do this in the meanwhile...
add_stats_path=${result_path}/added_stats/${trait}

bsub -J "filter_results_${pop}_${trait}[1-23]" -o "${out_path}/${trait}/LOGS/%J_filter_results.%I.o" -w"ended(poscon_check_${pop}_${trait})" -M3000 -R"select[mem>3000] rusage[mem=3000]" -q small -- filter_results.sh ${add_stats_path} ${result_path} ${trait} ${out_path}/${trait}

#now sort and merge chr files and filter by pval
bsub -J "sort_results_${pop}_${trait}" -o "${out_path}/${trait}/LOGS/%J_sort_results.o" -w"ended(filter_results_${pop}_${trait}) && ended(poscon_check_${pop}_${trait})" -M3000 -R"select[mem>3000] rusage[mem=3000]" -q small -- sort_results.sh ${add_stats_path} ${result_path} ${trait} ${out_path}/${trait}

#now generate plots
bsub -J "plot_results_${pop}_${trait}" -o "${out_path}/${trait}/LOGS/%J_plot_results.o" -w"ended(sort_results_${pop}_${trait})" -M6000 -R"select[mem>6000] rusage[mem=6000]" -q yesterday -- plotting_results.sh ${trait} ${out_path}/${trait} ${out_path}/${trait}/${trait}.all.poscon.join
