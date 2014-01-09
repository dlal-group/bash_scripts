#!/usr/local/bin/bash
#Launch VCFtools and qctool through LSF...
# $1: project/study/whateveryouwant name
# $2: vcf-file path
# $3: output path
# $4: min_maf
# $5: max_maf
#TODO: use getopt to read specified options...

if [ $# -lt 5 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "qc_vcftools_wrapper.sh <job name/project name> <vcf-file path> <output-path> <min-maf threshold> <max-maf threshold>\n"

	exit 1
fi

args=("$@")
prjname=${args[0]}
vcfpath=${args[1]}
outpath=${args[2]}
minmaf=${args[3]}
maxmaf=${args[4]}

echo -e "You selected:\n - Project name: $prjname\n - Vcf file path: $vcfpath\n - Output path: $outpath\n"

#FIXME:make this smart:we need to create the directory tree only if there are / in the string
mkdir -p ${outpath%/*}
mkdir -p ${outpath%/*}/LOGS
mkdir -p ${outpath%/*}/FREQ

#The following take too long!!We need to split in chromosome!
#bsub -J "$prjname" -o $prjname.log  -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
#-q basement vcf-stats $vcfpath -p $outpath
for i in {1..22}
do

if [ $# -lt 6 ]
then
	#Unfiltered
	echo "HWE unfiltered selection...."
	#-J job_name_spec[index |  start_index-end_index:step,]
	bsub -J "${prjname}_freq[${i}]" -o "%J_chr%I.log"  -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
	-q basement \
	vcftools --gzvcf $vcfpath --chr ${i} --remove-indels --min-alleles 2 --max-alleles 2 --freq --maf $minmaf --max-maf $maxmaf --out $outpath.chr${i}
else
	hwemax=${args[5]}
	#Filtered with hwe
	echo "HWE filtered selection...."
	bsub -J "${prjname}_freq_chr${i}" -o "${prjname}_freq_chr${i}.log"  -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
	-q basement \
	vcftools --gzvcf $vcfpath --chr ${i} --hwe $hwemax --remove-indels --min-alleles 2 --max-alleles 2 --freq --maf $minmaf --max-maf $maxmaf --out $outpath.chr${i}
fi
done

