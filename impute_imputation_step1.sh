#!/bin/bash

#Pipeline script for imputation
#this is the case of updated positions with liftover!

#Args (toBe substitute)
#$1= genotype file path
#$2= chr
#$3= reference file path
#$4= output file
#$5= geno mode TRUE/FALSE

#GENO_PATH=$g_path
#CHR=$

function my_trap_handler()
{
        MYSELF="$0"               # equals to my script name
        LASTLINE="$1"            # argument 1: last line of error occurence
        LASTERR="$2"             # argument 2: error code of last command
        echo "${MYSELF}: line ${LASTLINE}: exit status of last command: ${LASTERR}"

        # do additional processing: send email or SNMP trap, write result to database, etc.
	exit 1
}


if [ $# -lt 4 ]
then
	echo -e "**********************\nWRONG ARGUMENT NUMBER!!!\n**********************"
	echo -e "USAGE:\n impute_imputation_pipeline.sh <genotype files path> <chr> <reference files path> <output file path> [imputation mode]\n"
	echo -e "- <genotype files path> : path for genotypes files"
	echo -e "- <chr> : chromosome number"
	echo -e "- <reference files path> : path for reference files"
	echo -e "- <output files path> : output path"
	echo -e "- [imputation mode | update mode ] : if specified a 'geno' argument, the imputation will be performed using the genotypes otherwise the pre-phased haplotypes\n"

exit 1
fi

# trap commands with non-zero exit code
#
trap 'my_trap_handler ${LINENO} $?' ERR


mkdir -p $1/CHR$2

if [[ $2 =~ "^[0-9]+$" ]]
then
	if [ $2 -lt 10 ]
	then
		chr_num="0$2"
	else
		chr_num=$2
	fi
else
	chr_num=$2
fi


tar -xzvf $1/chr${chr_num}.tgz

mv chr${chr_num}.map $1/CHR$2/chr$2.map
mv chr${chr_num}.ped $1/CHR$2/chr$2.ped

#Filter files
#If we use chrX, don't use the mind filter!! (because we don't have internal imputation for this chr)
if [[ $2 =~ "^[0-9]+$" ]]
then
	plink --noweb --file $1/CHR$2/chr$2 --hwe 1e-6 --mind 0.1 --geno 0.1 --maf 0.01 --make-bed --out $1/CHR$2/chr$2.filtered
else
#	and also,change the chr23/25 name to chrX
	sed -i 's/^23/X/g;s/^25/X/g' $1/CHR$2/chr$2.map
	#sed -i 's/^25/23/g' $1/CHR$2/chr$2.map
	plink --noweb --file $1/CHR$2/chr$2 --hwe 1e-6 --geno 0.1 --maf 0.01 --make-bed --out $1/CHR$2/chr$2.filtered
fi

##PART 1: POSITION UPDATE
#create the TO_FLIP folder
mkdir -p $1/CHR$2/TO_FLIP

if [ $# -eq 5 ]
then

	case $5 in
	no_update)
	#if not needed, just convert filtered files to to_flip files
	plink --noweb --bfile $1/CHR$2/chr$2.filtered --make-bed --out $1/CHR$2/TO_FLIP/chr$2.to_flip
	;;
	geno)
	bash /home/mcocca/bash_scripts/update_pos_module.sh $1 $2
	;;
	esac

else
	#launch the module to update positions if needed
	bash /home/mcocca/bash_scripts/update_pos_module.sh $1 $2 
fi

echo -e "\n\nSTARTING FLIPPING PROCEDURE!!!!!\n********************************"
sleep 10

