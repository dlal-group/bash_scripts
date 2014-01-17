#!/usr/local/bin/bash
#Script to create a table to perform annotation with ANCESTRAL alleles on VCF data
#Args:
#	$1=Allele.bcp.gz
#	$2=SNPAncestralAllele.bcp.gz
#	$3=outpath
allele=$1
ancestral=$2
outpath=$3

#need to use awk to match second column of SNPAncestralAllele.bcp.gz with first column of Allele.bcp.gz
#Extract info from SNPAncestralAllele: get first 2 columns and extract a uniq list of markers and ancestral alleles code
zcat ${ancestral} | cut -f 1,2 | uniq | gzip -c > ${outpath}/${ancestral}.filtered.gz

#now take the allele table to retrieve the ancestral allele info..
awk 'NR==FNR{a[$1]=$2;next;}(a[$2])'