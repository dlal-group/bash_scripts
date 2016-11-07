#!/usr/bin/env bash

#script to submit a bunch of gemma jobs
# input COMPRESSED genotype file = /lustre/scratch113/projects/uk10k/users/jh21/imputed/fvg/uk10k1kg.shapeit
# input PLAIN genotype file = /lustre/scratch113/projects/uk10k/users/mc14/imputed/fvg/uk10k1kg.shapeit
# input phenotype file = /nfs/users/nfs_m/mc14/Work/SANGER/FVG/PHENO/ANTROP/new/others/gemma_pheno.txt
# input annotation file = /lustre/scratch113/projects/uk10k/users/jh21/imputed/fvg/uk10k1kg.shapeit/chr9.bimbam.pos
# input relatedness matrix file = /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/kinship/fvg_kin_c.cXX.txt

bimbam_path=$1 #same for genotypes and pos files
pheno_path=$2 #phenotype file path
#we'll read a phenotype list file: each row containd the phenotype name and the phenotype column number in the gemma phenotype file
kinship=$3 #kinship matrix path
platform=$4 #could be SANGER or GENEMONSTER
cov=$5 #covariate file path


if [ $# -lt 4 ]
then
	echo "MISSING ARGUMENTS!!!"
	echo -e "USAGE:\ngemma_launcher.sh <bimbam_path> <pheno_path> <kinship_file_path> <platform[SANGER|GENEMONSTER]> [<cov file path>]"
	exit 1
fi

if [ $# -eq 5 ]
then
	echo "ADDED COVARIATE FILE!!!"
fi

echo "Checking phenotype files existence...."

if [ ! -f ${pheno_path}/gemma_pheno_list.txt ]
then
	echo -e "ATTENTION!!\nMissing ${pheno_path}/gemma_pheno_list.txt file!\nexit from script!"
	exit 1
fi

if [ ! -f ${pheno_path}/gemma_pheno.txt ]
then
	echo -e "ATTENTION!!\nMissing ${pheno_path}/gemma_pheno.txt file!\nexit from script!"
	exit 1
fi

echo "Correct names for phenotype files detected!!"
echo "Cracking!!"
mkdir -p LOGS

while read line
do
	#line="7 TLM"
	trait=`echo $line | cut -f 2 -d " "`
	trait_n=`echo $line | cut -f 1 -d " "`

	if [[ $trait != "ID" ]]
	then
	#here goes the for loop for the chr
		echo "${trait} -> ${trait_n}"
		for chr in {22..22}
		do
			#command for farm 3
			if [ $# -eq 5 ]
			then
				#we want to provide also a covariate file for conditional analyses
				#so we check for a new argument: the covariate file path. If provided we run the analyses with this covariate
				#cova=`echo ${cov#*conditional_}| cut -f 1 -d "."`
				cova=`basename ${cov#*_}| cut -f 1 -d "."`
				echo "Covariate: ${cova}"
				echo "Covariate File: ${cov}"

				outfile=$trait.chr${chr}_${cova}.tab
				##gemma 0.93
				# bsub -J "gemma_${trait}_${chr}" -o "LOGS/%J_gemma_${trait}_${chr}.log" -e "LOGS/%J_gemma_${trait}_${chr}.err" -M3500 -R"select[mem>=3500] rusage[mem=3500]" -q normal -- /nfs/users/nfs_j/jh21/programs/gemma.0.93 -g ${bimbam_path}/chr${chr}.bimbam -p ${pheno_path}/gemma_pheno.txt -a ${bimbam_path}/chr${chr}.bimbam.pos -k ${kinship} -c ${cov} -maf 0 -miss 0 -lmm 4 -n ${trait_n} -o ${outfile}
				##gemma 0.94
				case $platform in
					SANGER )
						bsub -J "gemma_${trait}_${chr}" -o "LOGS/%J_gemma_${trait}_${chr}.log" -e "LOGS/%J_gemma_${trait}_${chr}.err" -M3500 -R"select[mem>=3500] rusage[mem=3500]" -q normal -- gemma -g ${bimbam_path}/chr${chr}.bimbam.gz -p ${pheno_path}/gemma_pheno.txt -a ${bimbam_path}/chr${chr}.bimbam.pos -k ${kinship} -c ${cov} -maf 0 -miss 0 -fa 4 -n ${trait_n} -o ${outfile}
					;;
					GENEMONSTER )
						gemma_cmd=/netapp/nfs/softwares/gemma095alpha/gemma
						echo "${gemma_cmd} -g ${bimbam_path}/chr${chr}.bimbam.gz -p ${pheno_path}/gemma_pheno.txt -a ${bimbam_path}/chr${chr}.bimbam.pos -k ${kinship} -c ${cov} -maf 0 -miss 0 -lm 4 -n ${trait_n} -o ${outfile}" | qsub -cwd -o LOGS/\$JOB_ID_gemma_${trait}_${chr}.log -e LOGS/\$JOB_ID_gemma_${trait}_${chr}.err -V -N gemma_${trait}_${chr} -l h_vmem=12G 
					;;
				esac
			else
				outfile=$trait.chr$chr.tab
				#with the first command we launch gemma on the whole set of snps, regardless of missing genotypes due to panel merging for fvg cohort
				
				## gemma 0.93
				# bsub -J "gemma_${trait}_${chr}" -o "LOGS/%J_gemma_${trait}_${chr}.log" -e "LOGS/%J_gemma_${trait}_${chr}.err" -M3500 -R"select[mem>=3500] rusage[mem=3500]" -q normal -- /nfs/users/nfs_j/jh21/programs/gemma.0.93 -g ${bimbam_path}/chr${chr}.bimbam -p ${pheno_path}/gemma_pheno.txt -a ${bimbam_path}/chr${chr}.bimbam.pos -k ${kinship} -maf 0 -miss 0 -lmm 4 -n ${trait_n} -o ${outfile}
				
				##gemma 0.94
				echo ${platform}
				case $platform in
					SANGER )
						bsub -J "gemma_${trait}_${chr}" -o "LOGS/%J_gemma_${trait}_${chr}.log" -e "LOGS/%J_gemma_${trait}_${chr}.err" -M3500 -R"select[mem>=3500] rusage[mem=3500]" -q normal -- gemma -g ${bimbam_path}/chr${chr}.bimbam.gz -p ${pheno_path}/gemma_pheno.txt -a ${bimbam_path}/chr${chr}.bimbam.pos -k ${kinship} -maf 0 -miss 0 -fa 4 -n ${trait_n} -o ${outfile}
						;;
					GENEMONSTER )
						echo "We're about to submit the analysis..."
						gemma_cmd=/netapp/nfs/softwares/gemma095alpha/gemma
						echo "${gemma_cmd} -g ${bimbam_path}/chr${chr}.bimbam.gz -p ${pheno_path}/gemma_pheno.txt -a ${bimbam_path}/chr${chr}.bimbam.pos -k ${kinship} -maf 0 -miss 0 -lm 4 -n ${trait_n} -o ${outfile}" | qsub -cwd -N "gemma_${trait}_${chr}" -o "LOGS/\$JOB_ID_gemma_${trait}_${chr}.log" -e "LOGS/\$JOB_ID_gemma_${trait}_${chr}.err" -V -l h_vmem=12G 
						;;
				esac
				
				#if we want to remove snps with missing data!!(but the miss option removes the sample or the snp??)
				#bsub -J "gemma_${trait}_${chr}" -o "%J_gemma_${trait}_${chr}.log" -e "%J_gemma_${trait}_${chr}.err" -M7000 -R"select[mem>=7000] rusage[mem=7000]" -q normal -- /nfs/users/nfs_y/ym3/bin/gemma -g ${bimbam_path}/chr${chr}.bimbam -p ${pheno_path}/gemma_pheno.txt -a ${bimbam_path}/chr${chr}.bimbam.pos -k ${kinship} -maf 0 -fa 4 -n ${trait_n} -o $trait.chr$chr.tab
			fi

			#now we want to compress all our outputs to save space
			# but first we want to check that the X chr is coded as X or 23 in the result files, at least
			# if not, we're going to code it as X
			if [[ $chr == "X" ]]
			then
				# chr_code=`cut -f 1 output/${outfile}.assoc.txt | tail -n+2 | sort | uniq`

				# if [[ $chr_code != "X" ]]
				# then
				case $platform in
					SANGER )
						echo "(head -1 output/${outfile}.assoc.txt;awk '{if(\$1 != \"X\") print \"X\",\$0;else print \$1,\$0}' output/${outfile}.assoc.txt | tail -n+2 | tr \" \" \"\t\" | cut -f 1,3-) | gzip -c > output/${outfile}.assoc.txt.gz" | bsub -J "gemma_${trait}_${chr}_shrink" -o "LOGS/%J_gemma_${trait}_${chr}_shrink.log" -e "LOGS/%J_gemma_${trait}_${chr}_shrink.err" -M2000 -R"select[mem>=2000] rusage[mem=2000]" -q normal
						;;
					GENEMONSTER )
						;;
				esac
					# echo "(head -1 output/${outfile}.assoc.txt;awk '{if(\$1 != \"X\") print \"X\",\$0;else print \$1,\$0}' output/${outfile}.assoc.txt | tail -n+2 | tr \" \" \"\t\" | cut -f 1,3-) | gzip -c > output/${outfile}.assoc.txt.gz" | bsub -J "gemma_${trait}_${chr}_shrink" -w "ended(gemma_${trait}_${chr})" -o "LOGS/%J_gemma_${trait}_${chr}_shrink.log" -e "LOGS/%J_gemma_${trait}_${chr}_shrink.err" -M2000 -R"select[mem>=2000] rusage[mem=2000]" -q normal
				# fi
			else
				echo "Submitting compression job..."
				# if [[ ! -f output/${outfile}.assoc.txt ]]; then
					#statements
					# echo "skipping output/${outfile}.assoc.txt, already compressed or not existing!!!"
				# else
				case $platform in
					SANGER )
						bsub -J "gemma_${trait}_${chr}_shrink" -w "ended(gemma_${trait}_${chr})" -o "LOGS/%J_gemma_${trait}_${chr}_shrink.log" -e "LOGS/%J_gemma_${trait}_${chr}_shrink.err" -M1000 -R"select[mem>=1000] rusage[mem=1000]" -q normal -- gzip output/${outfile}.assoc.txt
						;;
					GENEMONSTER )
						echo "gzip output/${outfile}.assoc.txt" | qsub -cwd -N gemma_${trait}_${chr}_shrink -o LOGS/\$JOB_ID_gemma_${trait}_${chr}_shrink.log -e LOGS/\$JOB_ID_gemma_${trait}_${chr}_shrink.err -V -l h_vmem=1G -hold_jid gemma_${trait}_${chr}
						;;
				esac
					# bsub -J "gemma_${trait}_${chr}_shrink" -o "LOGS/%J_gemma_${trait}_${chr}_shrink.log" -e "LOGS/%J_gemma_${trait}_${chr}_shrink.err" -M1000 -R"select[mem>=1000] rusage[mem=1000]" -q normal --  gzip output/${outfile}.assoc.txt
				# fi
			fi
		done
	fi

done < <(cat ${pheno_path}/gemma_pheno_list.txt)

