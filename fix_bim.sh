#!/usr/bin/env bash

#function to fix the missing rsID problem on PLINK files
# 19/8/2016 
function fix_bim() {
	orig=$1
	new=${orig}.NEW
	
	awk 'BEGIN{OFS="\t"}{if($2==".") print $1,"chr"$1":"$4,$3,$4,$5,$6;else print $0}' ${orig} > ${new}
	mv ${orig} ${orig}.OLD
	mv ${new} ${orig}
	
}
