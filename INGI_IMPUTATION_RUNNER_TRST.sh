#!/usr/bin/env bash

#we should run this like:
# for chr in {3..19}
# do
# echo "$chr $panel"
# ~/Work/bash_scripts/INGI_IMPUTATION_RUNNER_TRST.sh CARL CARL IMPUTE long 12000
# done
postfix=".shapeit"
impute2=/home/cocca/softwares/bin/impute2.3.2
shapeit2=/home/cocca/softwares/bin/shapeit
plink2=plink
chunk_size=3000000
buffer_size=250
window_size=2
thread=8
pop=$1
PANEL=$2
# chr=$3
MODE=$3 #set this to PHASE, if you want to phase and impute; set this to IMPUTE, if you're providing already phased genotypes
q=$4 #selected queue
m=$5 #select memory amount
genmap_dir=/netapp/nfs/resources/1000GP_phase3/impute
# base_out=/netapp/dati/WGS_REF_PANEL/08062016/IMPUTED
base_out=$6
chr=$7
exclude_base=/netapp/dati/INGI_WGS/18112015/${pop}/12112015_FILTERED_REL/LISTS
genotype_base=/netapp/dati/WGS_REF_PANEL/genotypes
refdir=/netapp/dati/WGS_REF_PANEL/REFERENCES/${PANEL}
#run IMPUTE script generator
# for chr in 2 6 11 21
# for chr in 20 21 22
# do
	imputedir=${base_out}/${pop}/${PANEL}$postfix/${chr}
	echo "Generating scripts for ${pop},${chr} (PANEL: ${PANEL}) "
	#added a parameter to discriminate if we are doing a test imputation with samples excluded or not
	# ~/scripts/bash_scripts/impute_INGI_ref_array_TRST.sh ${impute2} ${shapeit2} ${plink2} ${chunk_size} ${buffer_size} ${window_size} ${thread} ${pop} ${PANEL} ${chr} ${MODE} ${q} ${m} ${genmap_dir} ${base_out} ${exclude_base} ${genotype_base} ${refdir}
	~/scripts/bash_scripts/impute_INGI_ref_array_TRST.sh ${impute2} ${shapeit2} ${plink2} ${chunk_size} ${buffer_size} ${window_size} ${thread} ${pop} ${PANEL} ${chr} ${MODE} ${q} ${m} ${genmap_dir} ${base_out} ${exclude_base} ${genotype_base} ${refdir} IMPUTE
	
# done

#job submission
# for chr in 2 6 11 21
# for chr in 20 21 22
# do
	imputedir=${base_out}/${pop}/${PANEL}$postfix/${chr}
	echo "Submitting jobs for ${pop},${chr} (PANEL: ${PANEL}) "
	a_size=`wc -l $imputedir/chr${chr}_command.list| cut -f 1 -d " "`;echo "~/scripts/bash_scripts/ja_runner_TRST.sh -s $imputedir/chr${chr}_command.list"|qsub -t 1-${a_size} -o ${imputedir}/chr${chr}_\$JOB_ID_\$TASK_ID.log -e ${imputedir}/chr${chr}_\$JOB_ID_\$TASK_ID.e -V -N ${pop}_chr${chr} -l h_vmem=${m}
# done
