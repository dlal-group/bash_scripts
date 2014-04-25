#!/usr/local/bin/bash
#Script to create a pipeline for the purge through inbreeding work
#
#
#define the population we are using
pops="FVG VBI TSI CEU"

#retrieve the MODE parameter to select the correct operation
MODE=$1
# Merge different popuplation together
# TODO: add code!!!

#create a case statement to select the operation to do:
case $MODE in
  IBD )
    #Extract IBD information for each population. Here we are using plink2 (1.9)
    echo "Calculate IBD...."
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/fvg_seq/20140410/INGI/FVG
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/fvg_seq/20140410/INGI/VBI
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/fvg_seq/20140410/TGP/TSI
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/fvg_seq/20140410/TGP/CEU
            ;;
      esac
      bsub -J"ibd_${pop}" -o"%J_ibd_${pop}.o" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --vcf ${pop_path}/22.vcf.gz --double-id --biallelic-only --genome gz --parallel 1 2 --out ibd_${pop}
    done
  ;;
  HET )
    #Extract Inbreeding coeff information for each population. Here we are using plink2 (1.9)
    #we'll use the --het option as long as the --ibc option
    echo "Calculate Inbreeding...."
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/fvg_seq/20140410/INGI/FVG/PLINK
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/fvg_seq/20140410/INGI/VBI/PLINK
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/fvg_seq/20140410/TGP/TSI/PLINK
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/fvg_seq/20140410/TGP/CEU/PLINK
            ;;
      esac
        #calculate frequencies, before:
        bsub -J"freq_${pop}" -o"%J_freq_${pop}.o" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --bfile ${pop_path}/22.${pop,,} --freq --nonfounders --out freq_${pop}
        
        #use freq data
        bsub -J"inb_${pop}" -o"%J_inb_${pop}.o" -w "ended(freq_${pop})" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --bfile ${pop_path}/22.${pop,,} --het --read-freq freq_${pop}.frq --out inb_${pop}
        bsub -J"ibc_${pop}" -o"%J_ibc_${pop}.o" -w "ended(freq_${pop})" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --bfile ${pop_path}/22.${pop,,} --ibc --read-freq freq_${pop}.frq --out ibc_${pop}
      done
    ;;
  esac
