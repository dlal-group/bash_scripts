#!/usr/local/bin/bash
#Launch vcf-subset through LSF...
#TODO: use getopt to read specified options...
#The command I have to bsub :
#( zcat /lustre/scratch107/projects/esgi-vbseq/REL-2012-06-07/v1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz | head -200 | grep "^#" ; \
# tabix /lustre/scratch107/projects/esgi-vbseq/REL-2012-06-07/v1/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.vcf.gz 2:125320115-125320115 ) | vcf-subset -a | vcf-annotate --fill-AC-AN
#We need to do it only for multiallelic sites. so we have to split the file in biallelic and multiallelic, fix the multiallelic, merge again the file.
#(e sticazzi non ce lo mettiamo?)
# $1: vcf file path
# $2: output_path
# $3: 
# $4: 
# $5: 

if [ $# -lt 2 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "multiallelic_resolver.sh <vcf-file path> <output-path> [concatenate]\n"
	echo -e "<vcf-file path>: Vcf file path\n"
	echo -e "<output-path>:Desired output path\n"
	echo -e "[concatenate]: if you want ONLY to concatenate existing files!\n"
	exit 1
fi

args=("$@")
vcfpath=${args[0]}
outpath=${args[1]}

echo -e "You selected:\n - Vcf file path: $vcfpath\n - Output path: $outpath\n"

mkdir -p $outpath
if [ $# -eq  2 ]
then
#Now we need to split in chromosome!
#and foreach chr work only on multiallelic sites
for i in {1..22}
do
bsub -J "chr_splitter_${i}" -o "%J_chr${i}.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement "(zcat $vcfpath | head -200 | grep '^#' ; tabix $vcfpath ${i} | awk '{if (\$8 ~ /AC=.+,/) print \$0}') | vcf-subset -a | vcf-annotate --fill-ICF | bgzip -c > $outpath/chr${i}_fixed.vcf.gz"
done
#the same for chr X
bsub -J "chrX_splitter" -o "%J_chrX.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement "(zcat $vcfpath | head -200 | grep '^#' ; tabix $vcfpath X | awk '{if (\$8 ~ /AC=.+,/) print \$0}') | vcf-subset -a | vcf-annotate --fill-ICF | bgzip -c > $outpath/chrX_fixed.vcf.gz"

#in the same time create a file without multialellic sites..we are working also on indels!!
bsub -J "multi_remover" -o "%J_multiremover.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement "(zcat $vcfpath | awk '{if (\$8 !~ /AC=.+,/) print \$0}') | bgzip -c > $outpath/no_multiallelic.vcf.gz"
fi

#now whe need to merge everyfiles back again...
if [ $# -gt 2 ]
then
	echo -e "\nCONCATENATION OPTION SELECTED!!"
	echo -e "this only concatenate already existent files!!"
	cd $outpath
	VCF_FILES=`ls *.gz`
	#this is to fix cos we want to use a variable for the list of files
	bsub -J "chr_concat" -o "%J_chrconcat.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
	-q basement 'vcf-concat chr10_fixed.vcf.gz chr11_fixed.vcf.gz chr12_fixed.vcf.gz chr13_fixed.vcf.gz chr14_fixed.vcf.gz chr15_fixed.vcf.gz chr16_fixed.vcf.gz chr17_fixed.vcf.gz chr18_fixed.vcf.gz chr19_fixed.vcf.gz chr1_fixed.vcf.gz chr20_fixed.vcf.gz chr21_fixed.vcf.gz chr22_fixed.vcf.gz chr2_fixed.vcf.gz chr3_fixed.vcf.gz chr4_fixed.vcf.gz chr5_fixed.vcf.gz chr6_fixed.vcf.gz chr7_fixed.vcf.gz chr8_fixed.vcf.gz chr9_fixed.vcf.gz chrX_fixed.vcf.gz no_multiallelic.vcf.gz | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.MULTIFIXED.vcf.gz'
	cd ..
	
fi



