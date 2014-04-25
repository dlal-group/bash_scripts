#!/usr/local/bin/bash

#script just to reverse bam files to fastq files

infile=$1

filename=`basename ${infile}`

outfile=`echo ${filename%.*} | tr "#" "_"`

bam2fastq -o ${outfile}_R#.fastq ${infile}
#gzip the output
gzip ${outfile}_R_1.fastq
gzip ${outfile}_R_2.fastq

#check that everything is ok for reads length
zcat ${outfile}_R_1.fastq.gz | awk '{if(NR%4==2) print NR"\t"$0"\t"length($0)}' > ${outfile}_R_1.fastq.readLength
zcat ${outfile}_R_2.fastq.gz | awk '{if(NR%4==2) print NR"\t"$0"\t"length($0)}' > ${outfile}_R_2.fastq.readLength
#check that everything is ok for reads quality
zcat ${outfile}_R_1.fastq.gz | awk '{if(NR%4==0) print NR"\t"$0"\t"length($0)}'  > ${outfile}_R_1.fastq.qualityLength
zcat ${outfile}_R_2.fastq.gz | awk '{if(NR%4==0) print NR"\t"$0"\t"length($0)}'  > ${outfile}_R_2.fastq.qualityLength
#compare both
awk 'NR==FNR{a[$3]++;next}!a[$3]' ${outfile}_R_1.fastq.readLength ${outfile}_R_1.fastq.qualityLength > ${outfile}_R_1.fastq.difflength
awk 'NR==FNR{a[$3]++;next}!a[$3]' ${outfile}_R_2.fastq.readLength ${outfile}_R_2.fastq.qualityLength > ${outfile}_R_2.fastq.difflength


