#!/usr/local/bin/bash

#Script to extract some useful info:
#- How many increase above 2%
#- How many increase above 4%
#- How many increase above 10%
#- How many increase above 20%
#- How many decrease below 2%
#- How many decrease below 4%
#- How many decrease below 10%
#- How many decrease below 20%

file_path=$1
sign=$2
cd $file_path
#sign must be one between [increase,decrease,same]
files=`ls *$sign.txt`
mkdir DMAF_BINS
for file in $files
do
   awk '{OFS="\t"}{if ($7 < 0.02 ) print $0}' $file > DMAF_BINS/$file.lt02
   awk '{OFS="\t"}{if ($7 >= 0.02 && $7 <= 0.05) print $0}' $file > DMAF_BINS/$file.lte05
   awk '{OFS="\t"}{if ($7 > 0.05 && $7 <= 0.10) print $0}' $file > DMAF_BINS/$file.lte10
   awk '{OFS="\t"}{if ($7 > 0.10 && $7 <= 0.20) print $0}' $file > DMAF_BINS/$file.lte20
   awk '{OFS="\t"}{if ($7 > 0.20 && $7 <= 0.30) print $0}' $file > DMAF_BINS/$file.lte30
   awk '{OFS="\t"}{if ($7 > 0.30 && $7 <= 0.40) print $0}' $file > DMAF_BINS/$file.lte40
   awk '{OFS="\t"}{if ($7 > 0.40 && $7 <= 0.50) print $0}' $file > DMAF_BINS/$file.lte50
   awk '{OFS="\t"}{if ($7 == 0.50) print $0}' $file > DMAF_BINS/$file.eq50

done


