#!/usr/bin/env bash

#we should run this like:
# for chr in {3..19}
# do
# echo "$chr $panel"
# ~/Work/bash_scripts/INGI_IMPUTATION_RUNNER_TRST.sh CARL CARL ${chr} IMPUTE long 12000
# ~/Work/bash_scripts/INGI_IMPUTATION_RUNNER_TRST.sh FVG FVG ${chr} IMPUTE long 12000
# ~/Work/bash_scripts/INGI_IMPUTATION_RUNNER_TRST.sh VBI VBI ${chr} IMPUTE long 12000
# done

impute2=/home/cocca/softwares/bin/impute2.3.2
shapeit2=/home/cocca/softwares/bin/shapeit
plink2=plink
chunk_size=3000000
buffer_size=250
window_size=2
thread=8
pop=$1
PANEL=$2
chr=$3
MODE=$4 #set this to PHASE, if you want to phase and impute; set this to IMPUTE, if you're providing already phased genotypes
q=$5 #selected queue
m=$6 #select memory amount
genmap_dir=/netapp/nfs/resources/1000GP_phase3/impute
base_out=/netapp/dati/WGS_REF_PANEL/08062016/IMPUTED
exclude_base=/netapp/dati/INGI_WGS/18112015/
genotype_base=/netapp/dati/WGS_REF_PANEL/genotypes
refdir=/netapp/dati/WGS_REF_PANEL/REFERENCES/${PANEL}

#run IMPUTE script generator
~/scripts/bash_scripts/impute_INGI_ref_array.sh ${impute2} ${shapeit2} ${plink2} ${chunk_size} ${buffer_size} ${window_size} ${thread} ${pop} ${PANEL} ${chr} ${MODE} ${q} ${m} ${genmap_dir} ${base_out} ${exclude_base} ${genotype_base} ${refdir}

#job submission
qsub -o /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${chr}/chr${i}_shapeit.log -e /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${i}/chr${chr}_$chunkStr.e -V -N ${pop}_chr${chr}_$chunkStr -pe  $imputedir/chr$chr.$chunkStr.cmd
