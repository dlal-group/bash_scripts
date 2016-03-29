#!/usr/bin/env bash
#
# Script to process the selected files to create a reference panle
# We'll use this to:
# Simply MERGE VCF file -> there will be missing genotypes (phase 1)
# Phase and Impute the merged data -> to remove missing and to correct genotype errors (phase 1)
# To recall all data together from bam files using all sites selected from the merging step (phase 2)
# To refine genotypes calls after recalling (phase 3)
# To create merged reference panel in IMPUTE format for each phase

# vcf=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/22.vcf.gz
# cohort="VBI"
# outdir=/lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL/VBI
#ARGS:
# $1= vcf file for a single chromosome, better if named as "[chr].vcf.gz"
# $2=cohort
# $3=output folder
# $4=mode (snps/indels)


# /nfs/users/nfs_m/mc14/Work/bash_scripts/wgs_panel_process.sh /lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/22.vcf.gz VBI /lustre/scratch113/projects/esgi-vbseq/27112015_INGI_REF_PANEL
set -e

stage=$1
# we need to have stage as $4 when we run the panel conversion
# stage=$4
#we're going to split snps and indels, than put them back together again

case ${stage} in
VCF_MERGE )
# define path
vcf=$2
cohort=$3
outdir=$4/$3

mkdir -p ${outdir}/LOG_${stage}

filename=`basename ${vcf}`
first_suffix="${filename%%.*}"
chr=${first_suffix}

mkdir -p ${outdir}/${chr}
;;
IMPUTE_FORMAT )
vcf=$1
cohort=$2
outdir=$3/$2

mkdir -p ${outdir}/LOG_${stage}

filename=`basename ${vcf}`
first_suffix="${filename%%.*}"
chr=${first_suffix}

mkdir -p ${outdir}/${chr}
echo -e "${vcf}\n${cohort}\n${outdir}\n${filename}\n${first_suffix}\n${chr}"
#first we need to create the legend file
# we need rsID, position, REF,ALT: for consistency we'll assign rsID=chr:pos:ALT to all sites with unknown rsID
# echo "bcftools query -f\"%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n\" ${vcf} | awk '{if(\$3!=\".\") print \$3,\$2,\$4,\$5;else if (\$3==\".\") print \"chr\"\$1\":\"\$2\":\"\$5,\$2,\$4,\$5}'|sed '1 i\rsID position a0 a1'| gzip -c > ${outdir}/${chr}/${chr}.INGI_REF.${cohort}.legend.gz" | bsub -J"extract_${stage}_${cohort}_legend_${first_suffix}" -o"${outdir}/LOG_${stage}/1.%J_extract_${mode}_${cohort}_legend_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal
#than we need to create the hap file we need to create the legend file
# echo "bcftools query -f\"[%GT ]\\n\" ${vcf} | tr \"|\" \" \" | gzip -c > ${outdir}/${chr}/${chr}.INGI_REF.${cohort}.hap.gz" | bsub -J"extract_${mode}_${cohort}_hap_${first_suffix}" -o"${outdir}/LOG_${stage}/2.%J_extract_${stage}_${cohort}_hap_${first_suffix}.o" -M 1000 -R "select[mem>=1000] rusage[mem=1000]" -q normal
bcftools convert -h ${outdir}/${chr}/${chr}.INGI_REF.${cohort}.hap.gz,${outdir}/${chr}/${chr}.INGI_REF.${cohort}.legend.gz,${outdir}/${chr}/INGI_REF.${cohort}.samples ${vcf}
;;
PANEL_MERGE )
#this has to be an iterative process
# We need to tell which and how many cohorts we want to merge and merge them 2 by 2 using the output of the previous step
# to merge the subsequent cohort we can use two files:
# -one with the cohort order
# -one with the path for each cohort chr file

#first:write a chunk file
chr=$2
cohorts=(${@:3})

# cohorts=(TSI CARL FVG VBI)

#for the first round, pop1 is ${cohorts[0]}
# and pop2 is ${cohorts[1]}
# but from the second round, pop1 has to be pop1_pop2
# and pop2 has to be ${cohorts[2]} and so on until there are nop more cohorts in ${cohorts[@]}

# first round
gen_map=/lustre/scratch114/resources/imputation/impute2/2015-05-08/ALL_1000G_phase1interim_jun2011_impute/genetic_map_chr${chr}_combined_b37.txt
buffer=500
chunk_size=3000000

pop1=${cohorts[0]}
pop2=${cohorts[1]}
# pop2="CARL"
iter=0

for p in ${cohorts[@]}
do
if [[ $p == "${pop1}" ]];then
echo "First round"
echo "${pop1},${pop2}"
echo "${iter}"
outdir=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}
mkdir -p ${outdir}
mkdir -p ${outdir}/LOG_${stage}

hap_1=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}/${chr}/${chr}.INGI_REF.${pop1}.hap.gz
hap_2=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop2}/${chr}/${chr}.INGI_REF.${pop2}.hap.gz
leg_1=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}/${chr}/${chr}.INGI_REF.${pop1}.legend.gz
leg_2=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop2}/${chr}/${chr}.INGI_REF.${pop2}.legend.gz

#chunk generation step
/nfs/users/nfs_m/mc14/Work/bash_scripts/chunk_generator.sh ${leg_1} ${leg_2} ${chunk_size} ${pop1} ${pop2} ${chr}

echo -e "Populations: ${pop1} - ${pop2}\nChr: ${chr}\n 
Impute parameters:\n Genetic map: ${gen_map}\n Hap files \n -${hap_1}\n -${hap_2}\n
Legend files: \n -${leg_1}\n -${leg_2}\n
Chunk size: ${chunk_size}\n
Buffer: ${buffer}\n
Output folder: ${outdir}"

# run panel merging
size=`wc -l /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.chunks.txt|cut -f 1 -d " "`;bsub -J "merge_ref_${pop1}_${pop2}_${chr}[1-${size}]" -o "${outdir}/LOG_${stage}/%J_merge_ref_${pop1}_${pop2}_${chr}.%I.o" -M 3000 -R"select[mem>=3000] rusage[mem=3000]" -q normal -- /nfs/users/nfs_m/mc14/Work/bash_scripts/ja_merge_panels.sh /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.chunks.txt ${pop1} ${pop2} ${chr} ${gen_map} ${hap_1} ${hap_2} ${leg_1} ${leg_2} ${buffer}

#merge back panel files
bsub -J"merge_${chr}_${pop1}_${pop2}" -o"${outdir}/LOG_${stage}/%J_merge_${chr}_${pop1}_${pop2}.o" -w"ended(merge_ref_${pop1}_${pop2}_${chr}*)" -M1000 -R"select[mem>=1000] rusage[mem=1000]" -q normal -- /nfs/users/nfs_m/mc14/Work/bash_scripts/chunk_merger.sh ${pop1} ${pop2} ${chr}


elif [[ $p == "${pop2}" ]]; then
echo "${p}, second round"
iter=$[iter +1]
echo "${iter}"

else
echo ${p}

c_pop1=${pop1}
c_pop2=${pop2}
pop1=`echo "${pop1}_${pop2}"`
pop2=${p}
echo "${pop1},${pop2}"
iter=$[iter +1]
echo "${iter}"
outdir=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}
mkdir -p ${outdir}
mkdir -p ${outdir}/LOG_${stage}


hap_1=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}/${chr}/${chr}.INGI_REF.${pop1}.hap.gz
hap_2=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop2}/${chr}/${chr}.INGI_REF.${pop2}.hap.gz
leg_1=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}/${chr}/${chr}.INGI_REF.${pop1}.legend.gz
leg_2=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop2}/${chr}/${chr}.INGI_REF.${pop2}.legend.gz

#chunk generation step
bsub -J"chunks_${chr}_${pop1}_${pop2}" -o"${outdir}/LOG_${stage}/%J_chunks_${chr}_${pop1}_${pop2}.o" -w"ended(merge_${chr}_${c_pop1}_${c_pop2})" -M100 -R"select[mem>=100] rusage[mem=100]" -q normal -- /nfs/users/nfs_m/mc14/Work/bash_scripts/chunk_generator.sh ${leg_1} ${leg_2} ${chunk_size} ${pop1} ${pop2} ${chr}

echo -e "Populations: ${pop1} - ${pop2}\nChr: ${chr}\n 
Impute parameters:\n Genetic map: ${gen_map}\n Hap files \n -${hap_1}\n -${hap_2}\n
Legend files: \n -${leg_1}\n -${leg_2}\n
Chunk size: ${chunk_size}\n
Buffer: ${buffer}\n
Output folder: ${outdir}"

outdir=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}
mkdir -p ${outdir}
mkdir -p ${outdir}/LOG_${stage}

# run panel merging
size=`wc -l /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.chunks.txt|cut -f 1 -d " "`;bsub -J "merge_ref_${pop1}_${pop2}_${chr}[1-${size}]" -o "${outdir}/LOG_${stage}/%J_merge_ref_${pop1}_${pop2}_${chr}.%I.o" -w"ended(chunks_${chr}_${pop1}_${pop2})" -M 3000 -R"select[mem>3000] rusage[mem=3000]" -q normal -- /nfs/users/nfs_m/mc14/Work/bash_scripts/ja_merge_panels.sh /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.chunks.txt ${pop1} ${pop2} ${chr} ${gen_map} ${hap_1} ${hap_2} ${leg_1} ${leg_2} ${buffer}

#merge back panel files
bsub -J"merge_${chr}_${pop1}_${pop2}" -o"${outdir}/LOG_${stage}/%J_merge_${chr}_${pop1}_${pop2}.o" -w"ended(merge_ref_${pop1}_${pop2}_${chr}*)" -M 1000 -R"select[mem>=1000] rusage[mem=1000]" -q normal -- /nfs/users/nfs_m/mc14/Work/bash_scripts/chunk_merger.sh ${pop1} ${pop2} ${chr}


fi
done

;;
CHECK )
# We need to perform some checks on the resultind files
# 1) No presence of duplicated sites with alleles mismatch from different panels 
# 2) ....

# Get single panel legend files and check that overlapping position have the same REF ALT allelels 
# (REMEMBER we're working in MULTIALLELIC SPLITTED MODE)

;;
CONVERT_VCF)
#cpnvert files to VCF format
cohorts=$2
chr=${LSB_JOBINDEX}

bcftools convert -H /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${cohorts}/${chr}/${chr}.INGI_REF.${cohorts}.hap.gz,/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${cohorts}/${chr}/${chr}.INGI_REF.${cohorts}.legend.gz,/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${cohorts}/${cohorts}.REF.samples -O z -o /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${cohorts}/${chr}/${chr}.INGI_REF.${cohorts}.vcf.gz
tabix -f -p vcf /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${cohorts}/${chr}/${chr}.INGI_REF.${cohorts}.vcf.gz
;;
esac
