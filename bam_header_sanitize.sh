#!/usr/local/bin/bash
#
#
#Script to rehead bam files and remove empty CL: field from header
#
bam_path=$1
out_path=$2
old_id=$3
new_id=$4

mkdir -p ${out_path}

#I need to use a trick to rename bams from an old yasin naming convention
# bam_name=`echo ${bam_path} | cut -f 7 -d "/"`
#standard name-retrival procedure
bam_name=`basename ${bam_path}`

#extract header
samtools view -H ${bam_path} > ${out_path}/${bam_name}.header

#sanitize header
sed "s/@RG	ID:${old_id}	SM:${old_id}/@RG	ID:${new_id}	SM:${new_id}/g" ${out_path}/${bam_name}.header > ${out_path}/${bam_name}.header.sanitized

#rehead the bam files
samtools reheader ${out_path}/${bam_name}.header.sanitized ${bam_path} > ${out_path}/${bam_name}.reheaded.bam
#index file
samtools index ${out_path}/${bam_name}.reheaded.bam