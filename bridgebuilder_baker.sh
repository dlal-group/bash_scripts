#!/bin/bash

###############
#BAKER STAGE  #
###############
#ARGS:
#$1=bam file
#Build baker single-ended way to avoid resorting the original BAM according to read names
if [ ! -f ${1}\_bridge.bam ]; then
		/software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa aln -q 15 -b /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED/makebridge_bridge.fa ${1} > ${1}_bridge.sai
		/software/solexa/bin/aligners/bwa/bwa-0.5.9/bwa samse /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED/makebridge_bridge.fa ${1}_bridge.sai ${1} |samtools view -h -F 4 -Sb - > ${1}_bridge.bam
        rm ${1}\_bridge.sai
fi

