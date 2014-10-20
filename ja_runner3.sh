#!/bin/bash
# This is the runner file run by qsub on SGE

# Arguments: runner.sh filelist
# Environment variables: SGE_TASK_ID
#mkdir -p LOGS;size=`wc -l result.list|cut -f 1 -d " "`;qsub -t 1-${size} -N "check_${SGE_TASK_ID}" -o "LOGS/check_${SGE_TASK_ID}.o" -e "LOGS/check_${SGE_TASK_ID}.e" -l h_rt=200:00:00 -l virtual_free=5Gb -cwd -q xe-el6 -- ~/Work/bash_scripts/ja_runner.sh result.list


#script to run job arrays in a parametric way
file=`sed -n "${SGE_TASK_ID}p" $1`

filename=`basename ${file}`
#samtools mpileup -D -BQ0 -d10000000 -l /nfs/users/xe/ggirotto/multisample/exome_region_for_coverage_enriched_samtools.txt ${file} | awk '{x+=$4;next}END{print x/NR}' > ${filename}.mean_coverage
samtools mpileup -D -BQ0 -q20 -d5000 -l /nfs/users/xe/ggirotto/multisample/exome_region_for_coverage_enriched_samtools.txt ${file} | awk '{x+=$4;next}END{print x/NR}' > COV_CHK/${filename}.mean_coverage

#java -jar /users/GD/tools/GATK/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T DepthOfCoverage -I ${file} -L /users/GD/resource/human/probesets/nimblegene/v3/Target_Regions/SeqCap_EZ_Exome_v3_capture.bed_plus_150_shore.gatk.bed -R /users/GD/resource/human/hg19/hg19.fasta -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --includeRefNSites -o /nfs/users/xe/ggirotto/multisample/coverage_check/NEW_DATA/${filename}.GATK_DepCov.DATA

