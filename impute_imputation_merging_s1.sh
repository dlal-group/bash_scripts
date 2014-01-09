#!/usr/local/bin/bash

#Pipeline script for imputation
#this is the case of reference panel merging!

#Args (toBe substitute)
#$1 = genotype file path
#$2= chr
#$3= reference 1 file path
#$4= reference 2 file path
#$5= output file
#$6= geno mode TRUE/FALSE

function my_trap_handler()
{
        MYSELF="$0"               # equals to my script name
        LASTLINE="$1"            # argument 1: last line of error occurence
        LASTERR="$2"             # argument 2: error code of last command
        echo "${MYSELF}: line ${LASTLINE}: exit status of last command: ${LASTERR}"

        # do additional processing: send email or SNMP trap, write result to database, etc.
	exit 1
}


if [ $# -lt 6 ]
then
	echo -e "**********************\nWRONG ARGUMENT NUMBER!!!\n**********************"
	echo -e "USAGE:\n impute_imputation_merging_s1.sh <genotype files path> <chr> <reference 1 files path> <reference 2 files path> <output file path> <chunk file ref>\n"
	echo -e "- <genotype files path> : path for genotypes files"
	echo -e "- <chr> : chromosome number"
	echo -e "- <reference 1 files path> : path for first reference files"
	echo -e "- <reference 2 files path> : path for second reference files"
	echo -e "- <output files path> : output path"
	echo -e "- <chunk file ref> : path to the reference file to use to create the chunks\n"

exit 1
fi

# trap commands with non-zero exit code
#
trap 'my_trap_handler ${LINENO} $?' ERR


##PART 3: IMPUTATION

#function to create the script to be launched for each chunk
function build_chunk_template(){
cat << EOF
#BSUB -J "i_c$1_c$8"
#BSUB -o "%J_i_c$1_c$8.log"
#BSUB -e "%J_i_c$1_c$8.err"


chr=$1
start=$2
end=$3
ref_path1=$4
ref_path2=$5
geno_path=$6
out_path=$7
chunk_n=$8

#ARGS passed
#\$1=chr ($1)
#\$2=start ($2)
#\$3=end ($3)
#\$4=ref_path 1 ($4)
#\$5=ref_path 2 ($5)
#\$6=geno_path ($6)
#\$7=out_path ($7)
#\$8=chunk_n ($8)

/usr/local/bin/bash /nfs/users/nfs_m/mc14/Work/bash_scripts/impute2_merging_launcher_script.sh \${chr} \${start} \${end} \${ref_path1} \${ref_path2} \${geno_path} \${out_path} \${chunk_n}

EOF
}

#now we can launch impute2
GENOTYPES=$1
CHR=$2
REF1_PATH=$3
REF2_PATH=$4
OUTPUT_PATH=$5
#CHUNK_FILE_REF=$6

#create the output folder for current chr
mkdir -p $OUTPUT_PATH/CHR${CHR}

current_out_dir=$OUTPUT_PATH/CHR${CHR}
#we need to generate chromosome chunks

#use same approach used by Cinzia Sala:
#chose what file to use to create the chunks:
#directly passed start end from command line
#chrom_start=`zcat $CHUNK_FILE_REF | grep -v position | awk '{print $2}' | head -n 1`
chrom_start=$6
#chrom_end=`zcat $CHUNK_FILE_REF | grep -v position | awk '{print $2}' | tail -n 1`
chrom_end=$7

#find_start_end.sh $3/ALL_1000G_phase1integrated_v3_chr$2_impute.legend.gz $2 $3

#once done this, we can generate our chunks and launch impute
#initilize the chunk counter
chunk_start=$chrom_start
chunk_count=0
chunk_end=0
 
while [ $chunk_end -ne $chrom_end ]
do
	let chunk_count=$[chunk_count + 1]
	let chunk_end=$[chunk_start + 3000000 - 1]
	
	if [ $chunk_end -gt $chrom_end -o $chunk_end -eq $chrom_end ]
	then
		let chunk_end=$chrom_end
	fi
	
	echo -e "Processing CHROMOSOME ${CHR} ....\nCreated chunk $chunk_count \nStart: $chunk_start \nEnd: $chunk_end"

	start=$chunk_start
	end=$chunk_end
	out_path=$current_out_dir


	if [ $chunk_count -lt 10 ]
	then
		chunk_n="0$chunk_count"
	else
		chunk_n=$chunk_count
	fi
	
	build_chunk_template ${CHR} ${start} ${end} ${REF1_PATH} ${REF2_PATH} ${GENOTYPES} ${out_path} ${chunk_n} > ${out_path}/chr${CHR}.${chunk_n}.${start}_${end}.lsf
	
	#change permission
	chmod ug+x ${out_path}/chr${CHR}.${chunk_n}.${start}_${end}.lsf
	
	#qusub it!
	#bsub ${out_path}/chr${CHR}.${chunk_n}.${start}_${end}.lsf
	
	let chunk_start=$[chunk_end + 1]

done
