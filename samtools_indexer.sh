#!/usr/local/bin/bash
#indexer script

#ARGS:
#$1=bam filepath

files=(`ls $1/*.bam`)
index=$[LSB_JOBINDEX - 1]

samtools index ${files[$index]}
