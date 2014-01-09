#!/usr/local/bin/bash
#send PLINK merge command through lsf
#Usage:
echo "Send PLINK merge command through lsf.."

#if [$# -lt ]

bsub -J "merge_files.plink" -o merge_files.plink.log  -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement /software/varinf/bin/plink --noweb --file chr01 --merge-list mergelist.txt --recode --out geno_unfiltered

