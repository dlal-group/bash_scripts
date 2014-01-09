#Script to extract from gemma results common sites with other population
#!/usr/local/bin/bash

#Args:
#$1:gemma results path
#$2:population vcf files path
#$3:output path
#$4:population name

if [ $# -lt 4 ]
then
        echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
        echo -e "extract_gemma_results_from_other.sh <result path> <vcf files path> <output path> <population name> \n"

        exit 1
fi


res_path=$1
geno_path=$2
out_path=$3
pop=$4

for res in `ls $res_path/*.out`
do
	
	#parse res to find the trait name:
	#ATTENTION!! the name of the result MUST be in the format : N.TRAIT.suffix
	filename=${res##*/}
	first=${filename#*.}
	trait=${first%%.*}
	
	mkdir -p $out_path/${trait}
	
	echo "Processing trait ${trait}..."

	for chr in {1..22}
	do
		echo "CHR ${chr}.."
		#create chr result file
		#egrep "^${chr}	" ${res} | cut -f 3 > ${res}_${chr}
		echo "processing file ${res}_${chr}.."
		tabix $geno_path/chr${chr}.vcf.gz chr ${chr} | cut -f 1-7 | fgrep -w -f ${res}_${chr} > $out_path/${trait}/${trait}_${chr}_${pop}.vcf
	done
	
	echo "Done trait ${trait}!"


done
