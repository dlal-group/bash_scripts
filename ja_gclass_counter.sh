#!/usr/local/bin/bash

# This is the runner file run by bsub
# Arguments: runner.sh filelist
# Environment variables: LSB_JOBINDEX
# mkdir -p LOGS;size=`wc -l result.list|cut -f 1 -d " "`;bsub -J "p_check[1-${size}]" -o "LOGS/%J_p_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner.sh result.list
r_cat=`sed -n "${LSB_JOBINDEX}p" $1`
file=`sed -n "${LSB_JOBINDEX}p" $2`
g_cat=$2

echo ${g_cat} ${r_cat} ${file}


# bsub -o /lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza_2016/gerpinroh/new_counts/outerr/${g_cat}/ooo.${file}.${r_cat}.${g_cat}.out -e /lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza_2016/gerpinroh/new_counts/outerr/${g_cat}/ooo.${file}.${r_cat}.${g_cat}.err -G team151
# python /lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza_2016/gerp/GerpScr/GenotypeSummaryByRegionAndGerp_mergechr.py ${file} /lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza_2016/gerp/BedIndBySize/${r_cat}/${file}.${r_cat}.bed -2 2 > /lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza_2016/gerpinroh/new_counts/inds/${file}/${file}.${r_cat}.${g_cat}.count