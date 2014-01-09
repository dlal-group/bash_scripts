#!/usr/local/bin/bash
#Launch r script vb_vs_other.r via LSF
#file preprocessing...

if [ $# -lt 1 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "vbi_vs_other_launcher.sh <splitted-file prefix> \n"

	exit 1
fi

SPLITTED_FILES_prefix=$1

for i in {1..22};do grep -v CHROM ${SPLITTED_FILES_prefix}${i}.csv >> all_overlap.csv;done
#calculate maf for VBI and other population
awk '
{OFS="\t"}
{
	if ($6 > 0.5 && $8 > 0.5)
		print $0, 1-$6,1-$8;
	else if ($6 > 0.5 && $8 <= 0.5)
		print $0, 1-$6,$8;
	else if ($6 <= 0.5 && $8 > 0.5)
		print $0, $6,1-$8;
	else if ($6 <= 0.5 && $8 <= 0.5)
		print $0,$6,$8;
}' all_overlap.csv > all_overlap_maf.tmp.csv
rm all_overlap.csv

#fix also the MAF for EUR annotated in VBI file
awk '
{OFS="\t"}
{
	if ($7 != "NA") {
		if ($7 > 0.5)
			print $0,1-$7;
		else if ($7 <= 0.5)
			print $0,$7;
		}
	else if ($7 == "NA") {
		print $0,$7;
		}
}' all_overlap_maf.tmp.csv > all_overlap_maf.csv
rm all_overlap_maf.tmp.csv


#the job submission must be manually done after the calculation of all overlap maf for each population
#bsub -J "vb_vs_other" -o "%J_vb_vs_other_out.log" -M8000000 -R "select[mem>8000] rusage[mem=8000]" \
#-q basement 'R CMD BATCH /nfs/users/nfs_m/mc14/lustre109_home/GENOTIPI/COMPARISON/VBSEQ_QC/scripts/vb_vs_other.r'
