#!/bin/bash

#Pipeline script for imputation
#this is the case of updated positions with liftover!

#Args (toBe substitute)
#$1 = genotype file path
#$2= reference file path
#$3= output file
#$4= geno mode TRUE/FALSE


function my_trap_handler()
{
        MYSELF="$0"               # equals to my script name
        LASTLINE="$1"            # argument 1: last line of error occurence
        LASTERR="$2"             # argument 2: error code of last command
        echo "${MYSELF}: line ${LASTLINE}: exit status of last command: ${LASTERR}"

        # do additional processing: send email or SNMP trap, write result to database, etc.
	exit 1
}


if [ $# -lt 3 ]
then
	echo -e "**********************\nWRONG ARGUMENT NUMBER!!!\n**********************"
	echo -e "USAGE:\n impute_imputation_pipeline.sh [-g | --geno <genotype files path> ] [-r | --ref <reference files path> ] [-o | --out <output file path>] [-m | --mode <genotypes | prephased ] [-ufia]\n"
	echo -e "-g | --geno <genotype files path> : path for genotypes files"
	echo -e "-r | --ref <reference files path> : path for reference files"
	echo -e "-o | --out <output files path> : output path"
	echo -e "-m | --mode <geno | prephased> : if specified a 'geno' argument, the imputation will be performed using the genotypes otherwise the pre-phased haplotypes\n"
	echo -e "-u : only perform position update\n-f : only perform flipping step\n-i : only perform imputation\n-a : perform all the steps"
exit 1
fi

# trap commands with non-zero exit code
#
trap 'my_trap_handler ${LINENO} $?' ERR


#command to launch subsequent jobs for imputation

#for chr in 16
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
	echo -e "\n*************\nPROCESSING JOBS FOR CHR${chr}!!\n*******************"
	#launch step 1
	#step1_id=`qsub -V -N "step1_chr${chr}" -o "step1_chr${chr}.log" -q workq -- /home/mcocca/bash_scripts/impute_imputation_step1.sh $1 ${chr} $2 $3 no_update`
	step1_id=`qsub -V -N "step1_chr${chr}" -o "step1_chr${chr}.log" -q workq -- /home/mcocca/bash_scripts/impute_imputation_step1.sh $1 ${chr} $2 $3`
	
	echo -e "STEP 1 launched!"

	if [ $# -eq 4 ]
	then
		if [ $4 == "geno" ]
		then
			echo -e "\n*********************\nGENO MODE SELECTED FOR IMPUTATION!!!\n*********************"
			#launch step 2
		        step2_id=`qsub -V -N "step2_chr${chr}" -o "step2_chr${chr}.log" -W depend=afterok:${step1_id} -q workq -- /home/mcocca/bash_scripts/impute_imputation_step2.sh $1 ${chr} $2 $3 $4`
			echo -e "STEP 2 launched!"

			#launch step 3
			qsub -V -N "step3_chr${chr}" -o "step3_chr${chr}.log" -W depend=afterok:${step2_id} -q workq -- /home/mcocca/bash_scripts/impute_imputation_step3.sh $1 ${chr} $2 $3 $4
			echo -e "STEP 3 launched!"
		fi
	else
		#launch step 2
	        step2_id=`qsub -V -N "step2_chr${chr}" -o "step2_chr${chr}.log" -W depend=afterok:${step1_id} -q workq -- /home/mcocca/bash_scripts/impute_imputation_step2.sh $1 ${chr} $2 $3`
		echo -e "STEP 2 launched!"

		#launch step 3
		qsub -V -N "step3_chr${chr}" -o "step3_chr${chr}.log" -W depend=afterok:${step2_id} -q workq -- /home/mcocca/bash_scripts/impute_imputation_step3.sh $1 ${chr} $2 $3
		echo -e "STEP 3 launched!"
	fi

done
