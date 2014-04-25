#!/usr/local/bin/bash

bam=${1}
stats=$(samtools view $bam\_mkdup.bam|wc -l)
echo $bam $stats > ${bam}_new.stats

#to check total number of reads are the same before and after fix.
old=$(grep $bam $bam.stats|awk '{print $2}')
if [ $stats -eq $old ]; then 
 echo true;
 rm $bam\_mrg.bam
 rm $bam\_mrg.txt
 rm $bam\_rg.txt
 rm $bam\_rgfix.bam
 rm $bam\_rg.bam
 rm $bam\_new.bam
 rm $bam.temp
 touch $bam\_done
else
 echo $bam" total reads "$stats" does not match the record!"
 echo $bam>>rerun.txt
fi
