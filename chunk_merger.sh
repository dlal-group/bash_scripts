#!/usr/bin/env bash
#
#merge back chunks of haplotypes and legend files to generate the gz format
pop1=$1
pop2=$2
chr=$3

(echo "id position a0 a1";(zgrep -h -v position /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.*.legend.gz | sort -g -k2,2 ) ) | gzip -c > /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.legend.gz
(zcat /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.*.hap.gz ) | gzip -c > /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.hap.gz
# (echo "id position a0 a1";(zgrep -v position -h /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.legend.gz | awk -v chrom=${chr} '{print chrom":"$2"_"$3"_"$4,$2,$3,$4}')) | gzip -c > /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.reformatted.legend.gz
