#!/usr/local/bin/bash

#extract a snp list for different thresholds of a columns value
#Args:
#$1=imput file (complete path)
#$2=output path
#$3=column
#$4=threshold

filename=`basename $1`
mkdir -p $2/$4

#awk -v thr=$4 -v col=$3 '{if($col >= thr) print $0}' $1 > $2/$4/${filename}_$4

qctool -g $1 -omit-chromosome -og - -info $4 1 | gzip -c > $2/$4/${filename}