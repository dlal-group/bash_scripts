#!/usr/local/bin/bash
#Calculate Delta maf between population
#It uses a preformatted set of files, splitted by chr containing the overlap between populations
#file preprocessing...

if [ $# -lt 1 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "delta_maf_calc.sh <maf_filename> <pop_1 name> <pop_2 name> \n"

	exit 1
fi

MAF_FILE=$1
mkdir -p CHR_MAF
mkdir -p CHR_MAF/INCREASE
mkdir -p CHR_MAF/DECREASE
mkdir -p CHR_MAF/SAME
mkdir -p CHR_MAF/SUMMARIES
mkdir -p CHR_MAF/RESUMES
mkdir -p CHR_MAF/JPG
mkdir -p ALL_MAF

for i in {1..22}
do
	echo "Extracting chr ${i}.."
	grep "^$i	" ${MAF_FILE} | cut -f 1-3,5,9,10 > CHR_MAF/overlap_maf_chr${i}.tmp
	#we'll have a file with:
	#CHROM	POS	VB_ID	<pop_name>_ID	VBI_MAF <pop_name>_MAF	DELTA_MAF	SIGN
	#calculate delta maf for VBI and other population
	awk '
	{OFS="\t"}
	{
		if ($6 > $5)
			print $0, $6-$5,"-";
		else if ($6 == $5)
			print $0, $6-$5,"=";
		else if ($6 < $5)
			print $0, $5-$6,"+";
	}' CHR_MAF/overlap_maf_chr${i}.tmp > CHR_MAF/overlap_dmaf_chr${i}.csv
	rm CHR_MAF/overlap_maf_chr${i}.tmp

	echo "Delta MAF file created for chr ${i}!"
	echo "Plot job launch..."

	#then submit the plotting job, for each chr
	bsub -J "delta_maf_plot_chr${i}" -o "%J_delta_maf_plot_chr${i}.log" -M8000000 -R "select[mem>8000] rusage[mem=8000]" \
	-q basement R CMD BATCH '--args '$2' '$3' CHR_MAF/overlap_dmaf_chr'${i}'.csv '${i}'' /nfs/users/nfs_m/mc14/lustre109_home/GENOTIPI/COMPARISON/VBSEQ_QC/scripts/DELTA_MAF_plotter.R
	#move output files to chr specific dir
	bsub -J "mv_out_chr${i}" -w "ended(delta_maf_plot_chr${i})" -o "%J_mv_out_chr${i}.log" -M8000000 -R "select[mem>8000] rusage[mem=8000]" \
	-q normal 'mv *CHR_${i}_increase.txt CHR_MAF/INCREASE/;mv *CHR_${i}_decrease.txt CHR_MAF/DECREASE/;mv *.jpg CHR_MAF/JPG/;mv *${i}*resume.txt CHR_MAF/RESUMES/;mv *${i}*same.txt CHR_MAF/SAME/;mv summary*${i}* CHR_MAF/SUMMARIES/'

	#now concat each chr file in a big dmaf file
	cat CHR_MAF/overlap_dmaf_chr${i}.csv >> ALL_MAF/overlap_dmaf_all_chr.csv
done

#then submit the plotting job, for all chr
bsub -J "delta_maf_plot_all_chr" -o "%J_delta_maf_plot_all_chr.log" -M8000000 -R "select[mem>8000] rusage[mem=8000]" \
-q basement R CMD BATCH '--args '$2' '$3' ALL_MAF/overlap_dmaf_all_chr.csv' /nfs/users/nfs_m/mc14/Work/r_scripts/dmaf_all_chr_plotter.r
