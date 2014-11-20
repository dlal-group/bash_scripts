#!/usr/local/bin/bash

#script for retrieve ALT/REF information from other seq file

#first: generate a list of sites without the ALT allele (removed in v2 annotation by yasin)

#take as arguments:
#$1: chr
#$2: chr file path
#$3: map output path
#$4: reference annotated file path

#extract chr and position from my files, and save it in a bed-like format
tabix $2 chr $1 | cut -f 1-5 | grep "\.$" | cut -f 1,2 | awk '{printf("%s:%s-%s\n",$1,$2,$2)}'> $3/missing_allele_chr$1.map

#now we need to extract the allele information from the file we want
while read line
do
	tabix $4 chr $line | fgrep -v INDEL | cut -f 1-5 >> $3/fixed_missing_allele_chr$1.map
done < $3/missing_allele_chr$1.map 
