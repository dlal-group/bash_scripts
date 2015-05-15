#!/usr/local/bin/bash

# This is the runner file run by bsub
# Arguments: runner.sh filelist
# Environment variables: LSB_JOBINDEX
#Example command:
#mkdir -p LOGS;size=`wc -l result.list|cut -f 1 -d " "`;bsub -J "p_check[1-${size}]" -o "LOGS/%J_p_check.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh ~/Work/bash_scripts/check_pvals_test.sh result.list
#mkdir -p LOGS;size=`wc -l match_file.list|cut -f 1 -d " "`;bsub -J "sh_dens[1-${size}]" -o "LOGS/%J_sh_dens.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh ~/Work/bash_scripts/sharingDensity.sh match_file.list map.list
echo "${@}"
while getopts ":ds" opt; do
  case $opt in
    d)
      echo $opt
      echo "Double column list mode triggered!" >&2
      file=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $1}'`
      file2=`sed -n "${LSB_JOBINDEX}p" $3 | awk '{print $2}'`
      echo ${file}
      echo ${file2}
      script=$2
      $script ${file} ${file2} "${@:4}"
      ;;
    s)
      echo $opt
      echo "Single list mode triggered!!" >&2
      file=`sed -n "${LSB_JOBINDEX}p" $3`
      # echo ${file}
      script=$2
      # $script ${file} $4 $5 $6 $7 $8
      $script ${file} "${@:4}"
      ;;
    *)
      echo $opt
    ;;
  esac
done
#bash $script ${file} $3 $4 $5 $6 $7 $8

#06/06/2014
#order inverted ONLY FOR consequence list file splitting
# $script $3 $4 ${file} $5 $6 $7 $8
# /software/vertres/codebase/scripts/bamcheck -c 1,50,1 -d /lustre/scratch113/projects/fvg_seq/F12HPCEUHK0358_HUMpngR/BRIDGED_BAMS/${file}
# zcat ${file} | sed 's/,rs/|rs/g' | awk -v chr=${LSB_JOBINDEX} '{ snp=(NF-5)/3; if($2 ~/^rs/) s=$2;else s="NA"; printf "chr"chr":"$3"-"s"-"$4"-"$5"," $4 "," $5; for(i=1; i<=snp; i++) printf "," $(i*3+3)*2+$(i*3+4); printf "\n" }' > chr${LSB_JOBINDEX}.bimbam


