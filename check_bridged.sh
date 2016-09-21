#!/usr/local/bin/bash


files=`ls */Alignment_result/*new*`

for file in $files
do

(echo $file;samtools flagstat $file;echo ;echo) >> bridged_flagstat_check.txt

done                                                                                                          

