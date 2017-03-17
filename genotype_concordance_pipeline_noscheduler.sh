#!/bin/bash

#split files by chr and analyze one chr a time
#ARGS:
#$1 : file for other data
#$2 : file for SEQ
#$3 : indiv list

if [ $# -lt 3 ]
then
	echo -e "*******ATTENTION!!!THE SAMPLE LIST AHS TO BE SORTED BY SAMPLE NAME!!!******\n"
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "genotype_concordance_pipeline.sh <plink file name (no extension) for 1st genotype set (the GWAS dataset)> <plink file name (no extension) for 2nd genotype set (the WGS dataset)> <SORTED indiv list file path> [<reference_table_file>]\n"
	exit 1
fi

#26/11/2013
#input file preparation step
# bsub -J "input_prep_$1" -o "%J_input_prep_$1.o" -M2000 -R"select[mem>2000] rusage[mem=2000]" -q yesterday -- genotype_concordance_input_prep.sh 
# module add hgi/plink/1.90b3w

outname_f1=`basename $1`.s1_gwas
outname_f2=`basename $2`.s2_seq

#recode files in ped format and standardized names
plink --bfile $1 --keep $3 --recode --allow-no-sex --out ${outname_f1}

plink --bfile $2 --keep $3 --recode --allow-no-sex --out ${outname_f2}

#create the freq file for first dataset:this is useful to set the reference only if the first dataset came from gwas data
plink --file ${outname_f1} --freq --allow-no-sex --out ${outname_f1}

#correctly format the frq file
sed -i 's/^[ 	]*//;s/ \+/	/g' ${outname_f1}.frq

#create the freq file for second dataset
plink --file ${outname_f2} --freq --allow-no-sex --out ${outname_f2}

#correctly format the frq file
sed -i 's/^[ 	]*//;s/ \+/	/g' ${outname_f2}.frq

#sort the same file1 and file2
sort -k1 ${outname_f1}.ped > ${outname_f1}.sorted;mv ${outname_f1}.ped ${outname_f1}.unsorted;mv ${outname_f1}.sorted ${outname_f1}.ped;sort -k1 ${outname_f2}.ped > ${outname_f2}.sorted;mv ${outname_f2}.ped ${outname_f2}.unsorted;mv ${outname_f2}.sorted ${outname_f2}.ped

#we need to flip data after check on freq tables
$MY_BASH_SCRIPTS/genotype_concordance_variant_flipping.sh ${outname_f1}.frq ${outname_f2}.frq ${outname_f2}

#create the ref table for  REF/ALT allele:basically this contains the frequencies and the major allele because usually this is the REF allele in GWAS data
if [ $# -lt 4 ]
then
	# bsub -J "ref_discordance_table_creator" -w "done(sort_genotype_table)" -o "%J_ref_discordance_table_creator.log" -M8000 -R"select[mem>8000] rusage[mem=8000]" \
	R CMD BATCH "--args ${outname_f1}.frq ${outname_f1}.map" $MY_R_SCRIPTS/gt_discordance.r
fi

chr_list=`cut -f 1 $1.bim| uniq`

# for chr in {1..22}
for chr in ${chr_list}
do
	if [ $# -eq 4 ]
		#if we provide the 4th argument we are going to use our own REF table:this apply for data as WES or WGS
	then
		$MY_BASH_SCRIPTS/discordance_inner_cicle.sh ${outname_f1} ${outname_f2}.flipped ${chr} $3 $4
	else
		#use the automatically created ref discordance table
		$MY_BASH_SCRIPTS/discordance_inner_cicle.sh ${outname_f1} ${outname_f2}.flipped ${chr} $3 ref_discordance_table.txt
	fi
done

$MY_BASH_SCRIPTS/last_discordance_pipeline_step.sh $3
