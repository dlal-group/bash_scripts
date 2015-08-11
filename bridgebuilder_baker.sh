#!/bin/bash

###############
#BAKER STAGE  #
###############
#ARGS:
#$1=bam file
#Build baker single-ended way to avoid resorting the original BAM according to read names
# module unload hgi/bwa/latest
# module add hgi/bwa/0.5.9
module unload hgi/bwa/0.5.9
module add hgi/bwa/latest
bam_name=$(basename ${1})
if [ ! -f ${2}/${bam_name}\_bridge.bam ]; then
		# /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED/makebridge_bridge.fa ${1} > ${1}_bridge.sai
		bwa aln -q 15 -b /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED/makebridge_bridge.fa ${1} > ${2}/${bam_name}_bridge.sai
		# /software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa samse /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED/makebridge_bridge.fa ${1}_bridge.sai ${1} |samtools view -h -F 4 -Sb - > ${1}_bridge.bam
		bwa samse /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED/makebridge_bridge.fa ${2}/${bam_name}_bridge.sai ${1} |samtools view -h -F 4 -Sb - > ${2}/${bam_name}_bridge.bam
        rm ${2}/${bam_name}\_bridge.sai
fi

