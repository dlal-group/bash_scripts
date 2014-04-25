#!/usr/local/bin/bash
#Download files from 1KG and extract population individuals, then remove files useless to save space

if [ $# -lt 2 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "population_splitter.sh <population_name> <keeplist_path> \n"
	echo -e "<population_name> : Population name\n"
	echo -e "<keeplist_path>: keeplist path\n"
	exit 1
fi

args=("$@")
popname=${args[0]}
keeplist=${args[1]}

echo -e "You selected:\n - Population: $popname\n - Keeplist path: $keeplist\n"


if [ $# -eq  2 ]
then
#Now we need to work by chromosome!
#all jobs must depend on the previous
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
for i in X
do
#if [ $i -lt 8 ]
#then
	bsub -J "chr${i}_downloader" -o "%J_chr${i}_downloader.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
	-q basement wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
#else
#	h=$[i-2]
#	bsub -J "chr${i}_downloader" -w "ended(chr${h}_extractor)" -o "%J_chr${i}_downloader.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
 #      -q basement wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
#fi
#now we index the file
bsub -J "chr${i}_indexer" -w "ended(chr${i}_downloader)" -o "%J_chr${i}_indexer.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement tabix ALL.chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

#now extract the wanted individuals (and annotate again the file)
bsub -J "chr${i}_extractor" -w "ended(chr${i}_indexer)" -o "%J_chr${i}_extractor.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement /nfs/team151/software/vcftools_0.1.9/cpp/vcftools --gzvcf ALL.chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz --keep $keeplist --recode --recode-INFO-all --out $popname.chr${i}

#now gzip
bsub -J "chr${i}_bgzipper" -w "ended(chr${i}_extractor)" -o "%J_chr${i}_bgzipper.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement bgzip $popname.chr${i}.recode.vcf

#now remove useless files before starting the new download
bsub -J "chr${i}_remover" -w "ended(chr${i}_extractor)" -o "%J_chr${i}_remover.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement rm ALL.chr${i}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

done
	
fi



