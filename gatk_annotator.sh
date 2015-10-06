#!/usr/local/bin/bash

#script to annotate files by chromosome
#Args: 
#Use LSB_JOBINDEX as chr number!
# chr=${LSB_JOBINDEX}
# 
# Command line example:
# mkdir -p LOGS;size=`wc -l <vcf file list> |cut -f 1 -d " "`; bsub -J "annotate_gatk[1-${size}]" -o "LOGS/%J_annotate_gatk.%I.o" -M5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ja_runner_par.sh -s gatk_annotator.sh <vcf file list> <output folder>


infile=$1
# infolder=$2
outfolder=$2
mkdir -p ${outfolder}
#set those to the last available
dbsnp=/lustre/scratch114/resources/variation/Homo_sapiens/grch37/dbsnp_142_gatk.vcf
ref=/lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa
# TGP=/lustre/scratch114/resources/1000g/release/20130502/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz
TGP=/lustre/scratch114/resources/1000g/release/20130502/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz
infile_name=`basename ${infile}`

#first fix the header AND filter out all sites with a deletion AT THE BEGINNIG of the ALT allele field
# (tabix -H ${infolder}/${chr}.vcf.gz ${chr}|grep "^##";cat /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140711_ANNOTATED/to_add_to_header2.txt;tabix -H ${infolder}/${chr}.vcf.gz ${chr}| tail -n1;tabix ${infolder}/${chr}.vcf.gz ${chr})| awk '$5 !~ "^<DEL>.+"'|dos2unix | bgzip -c > ${infolder}/${chr}.fixed.vcf.gz;
# (tabix -H 1.vcf.gz 1|grep "^##";cat /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140711_ANNOTATED/to_add_to_header2.txt;tabix -H 1.vcf.gz 1| tail -n1;tabix 1.vcf.gz 1)| bgzip -c > ${infolder}/${chr}.fixed.vcf.gz;
# tabix -f -p vcf ${infolder}/${chr}.fixed.vcf.gz;

# bcftools annotate -x INFO/1kg_AA,INFO/1kg_AF,INFO/1kg_AMR_AF,INFO/1kg_ASN_AF,INFO/1kg_AFR_AF,INFO/1kg_EUR_AF,FORMAT/DP,FORMAT/DV,FORMAT/GQ,FORMAT/PL,FORMAT/PQ,FORMAT/PS,FORMAT/SP,FORMAT/DS,FORMAT/GL ${infolder}/${chr}.fixed.vcf.gz -O z -o ${infolder}/${chr}.fixed.cleaned.vcf.gz
# bcftools annotate -x INFO/1kg_AA,INFO/1kg_AF,INFO/1kg_AMR_AF,INFO/1kg_ASN_AF,INFO/1kg_AFR_AF,INFO/1kg_EUR_AF,FORMAT/DP,FORMAT/DV,FORMAT/GQ,FORMAT/PL,FORMAT/PQ,FORMAT/PS,FORMAT/SP,FORMAT/DS,FORMAT/GL ${infile} -O z -o ${infolder}/${infile_name}.fixed.cleaned.vcf.gz
# tabix -f -p vcf ${infolder}/${chr}.fixed.cleaned.vcf.gz


# /software/jre1.7.0_25/bin/java -Xmx1000m -Xms1000m -server -XX:+UseSerialGC -jar /nfs/users/nfs_m/mercury/src/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T VariantAnnotator \
java -Xmx2500m -Xms2500m -server -XX:+UseSerialGC -jar /software/hgi/pkglocal/gatk-protected-3.3/GenomeAnalysisTK.jar -T VariantAnnotator \
--variant ${infile} \
--out ${outfolder}/${infile_name}.ann.vcf.gz \
-R ${ref} \
--dbsnp ${dbsnp} \
--resource:1kg ${TGP} \
-E 1kg.AA \
-E 1kg.EAS_AF \
-E 1kg.AMR_AF \
-E 1kg.AFR_AF \
-E 1kg.SAS_AF \
-E 1kg.EUR_AF
# -E 1kg.AF \

#index with tabix
tabix -p vcf ${outfolder}/${infile_name}.ann.vcf.gz

zcat ${outfolder}/${infile_name}.ann.vcf.gz | sed "s/1kg\.//g" | bgzip -c > ${outfolder}/${infile_name}.clean_ann.vcf.gz

tabix -p vcf ${outfolder}/${infile_name}.clean_ann.vcf.gz

#clean stuff
# rm ${infolder}/${chr}.fixed.vcf.gz
# rm ${infolder}/${chr}.fixed.cleaned.vcf.gz
# rm ${infolder}/${chr}.fixed.vcf.gz.tbi
# rm ${infolder}/${chr}.fixed.cleaned.vcf.gz.tbi

#bit to generate a report....
# PID=$!
# wait $!
# status=$?
# wdir=`pwd -P`
# cmd=`history | tail -n2| head -1| cut -f 2- -d " "`
# info=${infile}
# email=mc14@sanger.ac.uk
# /nfs/users/nfs_m/mc14/Work/bash_scripts/send_report.sh ${status} ${email} ${wdir} ${cmd} ${info}