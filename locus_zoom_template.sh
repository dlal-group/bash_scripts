#!/usr/local/bin/bash
#
#Generate locuszoom plots using a file of regions

#Args
#$1=region file (CHR START END)
#$2=result file
#$3=markercol
#$4=pvalcol
#$5=refsnp

# Environment variables: LSB_JOBINDEX
#How to launch the script
#Use this sample command:
# - mkdir -p LOGS;size=`wc -l test_loczoom_tc_fvg_region.list|cut -f 1 -d " "`;bsub -J "lzplot[1-${size}]" -o "LOGS/%J_lzplot.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/locus_zoom_template.sh test_loczoom_tc_fvg_region.list test_loczoom_tc_fvg_2.txt SNPid p_wald FVG chr1:33472299
# - Region file has to be formatted as CHR START END, regardless of file separator
# - result file has to be UNCOMPRESSED and TAB separated (or you can change the separator in the script manually)
# - In the region file you can specify a reference SNP in the forth column, so you dont have to enter it manually each time
# - If you don't provide the ref snp it'll use the one with min pval.


line=`sed -n "${LSB_JOBINDEX}p" $1`
resultfile=$2
markercol=$3
pvalcol=$4
prefix=$5
refsnp=$6

reg_length=`echo ${line} | awk '{print NF}'`
chrom=`echo ${line} | awk '{print $1}'`
start=`echo ${line} | awk '{print $2}'`
end=`echo ${line} | awk '{print $3}'`
echo ${chrom}
echo ${start}
echo ${end}
echo ${reg_length}

if [ $# -eq 6 ]
then
echo "Using a reference snp: ${refsnp}"
locuszoom --no-date --cache None --plotonly \
--metal ${resultfile} \
--chr ${chrom} \
--start ${start} \
--end ${end} \
--refsnp ${refsnp} \
--delim=tab \
--prefix ${prefix} \
--markercol ${markercol} \
--pvalcol ${pvalcol} \
--flank 500kb --build uk10k --pop EUR \
--source UK10K_v4 recombColor=grey showAnnot=F ldColors="gray,lightgray,lightskyblue,green,orange,red,purple" \
--gwas-cat whole-cat_significant-only

else

if [ ${reg_length} -eq 4 ]
then
echo "Using a reference snp: ${refsnp}"
refsnp=`echo ${line} | awk '{print $4}'`
locuszoom --no-date --cache None --plotonly \
--metal ${resultfile} \
--chr ${chrom} \
--start ${start} \
--end ${end} \
--refsnp ${refsnp} \
--delim=tab \
--prefix ${prefix} \
--markercol ${markercol} \
--pvalcol ${pvalcol} \
--flank 500kb --build uk10k --pop EUR \
--source UK10K_v4 recombColor=grey showAnnot=F ldColors="gray,lightgray,lightskyblue,green,orange,red,purple" \
--gwas-cat whole-cat_significant-only

else

locuszoom --no-date --cache None --plotonly \
--metal ${resultfile} \
--chr ${chrom} \
--start ${start} \
--end ${end} \
--delim=tab \
--prefix ${prefix} \
--markercol ${markercol} \
--pvalcol ${pvalcol} \
--flank 500kb --build uk10k --pop EUR \
--source UK10K_v4 recombColor=grey showAnnot=F ldColors="gray,lightgray,lightskyblue,green,orange,red,purple" \
--gwas-cat whole-cat_significant-only
        
fi
fi