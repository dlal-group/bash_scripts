#!/bin/bash

#split files by chr and analyze one chr a time
#ARGS:
#$1 : file for other data
#$2 : file for SEQ
#$3 : indiv list

if [ $# -lt 3 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "genotype_concordance_pipeline.sh <plink file name (no extension) for 1st genotype set (the GWAS dataset)> <plink file name (no extension) for 2nd genotype set (the WGS dataset)> <indiv list file path> [<reference_table_file>]\n"
	exit 1
fi

#26/11/2013
#input file preparation step
# bsub -J "input_prep_$1" -o "%J_input_prep_$1.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q yesterday -- genotype_concordance_input_prep.sh 

outname_f1=`basename $1`.s1
outname_f2=`basename $2`.s2

#recode files in ped format and standardized names
bsub -J "recode_f1" -o "%J_recode_f1.log" -M1000 -R"select[mem>1000] rusage[mem=1000]" \
-q normal "plink2 --bfile $1 --recode --allow-no-sex --out ${outname_f1}"

bsub -J "recode_f2" -o "%J_recode_f2.log" -M1000 -R"select[mem>1000] rusage[mem=1000]" \
-q normal "plink2 --bfile $2 --recode --allow-no-sex --out ${outname_f2}"

#create the freq file for first dataset:this is useful to set the reference only if the first dataset came from gwas data
bsub -J "create_freq_table_d1" -o "%J_create_freq_table_d1.log" -w "ended(recode_f1)" -M1000 -R"select[mem>1000] rusage[mem=1000]" \
-q normal "plink2 --file ${outname_f1} --freq --allow-no-sex --out ${outname_f1}"

#correctly format the frq file
bsub -J "format_freq_table_f1" -o "%J_format_freq_table_f1.log" -w "ended(create_freq_table_d1)" -M4000 -R"select[mem>4000] rusage[mem=4000]" \
-q normal "sed -i 's/^[ 	]*//;s/ \+/	/g' ${outname_f1}.frq"

#create the freq file for second dataset
bsub -J "create_freq_table_d2" -o "%J_create_freq_table_d2.log" -w "ended(recode_f2)" -M1000 -R"select[mem>1000] rusage[mem=1000]" \
-q normal "plink2 --file ${outname_f2} --freq --allow-no-sex --out ${outname_f2}"

#correctly format the frq file
bsub -J "format_freq_table_f2" -o "%J_format_freq_table_f2.log" -w "ended(create_freq_table_d2)" -M4000 -R"select[mem>4000] rusage[mem=4000]" \
-q normal "sed -i 's/^[ 	]*//;s/ \+/	/g' ${outname_f2}.frq"


#sort the same file1 and file2
bsub -J "sort_genotype_table" -o "%J_sort_genotype_table.log" -w "ended(format_freq_table_f*)" -M4000 -R"select[mem>4000] rusage[mem=4000]" \
-q normal "sort -k1 ${outname_f1}.ped > ${outname_f1}.sorted;mv ${outname_f1}.ped ${outname_f1}.unsorted;mv ${outname_f1}.sorted ${outname_f1}.ped;sort -k1 ${outname_f2}.ped > ${outname_f2}.sorted;mv ${outname_f2}.ped ${outname_f2}.unsorted;mv ${outname_f2}.sorted ${outname_f2}.ped"


#create the ref table for  REF/ALT allele:basically this contains the frequencies and the major allele because usually this is the REF allele in GWAS data
if [ $# -lt 4 ]
then
	bsub -J "ref_discordance_table_creator" -w "done(sort_genotype_table)" -o "%J_ref_discordance_table_creator.log" -M8000 -R"select[mem>8000] rusage[mem=8000]" \
-q normal R CMD BATCH "--args ${outname_f1}.frq ${outname_f1}.map" /nfs/users/nfs_m/mc14/Work/r_scripts/gt_discordance.r
fi

# for chr in {1..22}
for chr in 11
do
	if [ $# -eq 4 ]
		#if we provide the 4th argument we are going to use our own REF table:this apply for data as WES or WGS
	then
		bsub -J "gt_calculator_chr${chr}" -o "%J_gt_calculator_chr${chr}.log" -M8000 -R"select[mem>8000] rusage[mem=8000]" \
		-q normal discordance_inner_cicle.sh ${outname_f1} ${outname_f2} ${chr} $3 $4
	else
		#use the automatically created ref discordance table
		bsub -J "gt_calculator_chr${chr}" -w "ended(ref_discordance_table_creator)" -o "%J_gt_calculator_chr${chr}.log" -M8000 -R"select[mem>8000] rusage[mem=8000]" \
		-q normal discordance_inner_cicle.sh ${outname_f1} ${outname_f2} ${chr} $3 ref_discordance_table.txt
	fi
done

bsub -J "last_step" -w "ended(gt_calculator_chr*)" -o "%J_last_step.log" -M8000 -R"select[mem>8000] rusage[mem=8000]" \
-q normal last_discordance_pipeline_step.sh $3
