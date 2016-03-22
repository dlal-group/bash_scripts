#!/usr/bin/env bash
#
#merge back chunks of haplotypes and legend files to generate the gz format
pop1=$1
pop2=$2
chr=$3

(echo "id position a0 a1";(fgrep -h -v position /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.*.legend | sort -g -k2,2 ) ) | gzip -c > /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.legend.gz
(cat /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.*.hap ) | gzip -c > /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${pop1}_${pop2}/${chr}/${chr}.INGI_REF.${pop1}_${pop2}.hap.gz
