#!/usr/local/bin/bash
#Script to create a pipeline for the purge through inbreeding work
#
#
#define the population we are using
pops="FVG VBI TSI CEU"

#retrieve the MODE parameter to select the correct operation
MODE=$1
CHR=$2
#set parameters for beagle:
window=$3
overlap=$4
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
        bsub -J"inb_${pop}" -o"%J_inb_${pop}.o" -w "ended(freq_${pop})" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --bfile ${pop_path}/22.${pop,,} --het --out inb_${pop}
        bsub -J"ibc_${pop}" -o"%J_ibc_${pop}.o" -w "ended(freq_${pop})" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --bfile ${pop_path}/22.${pop,,} --ibc --out ibc_${pop}
      done
    ;;
    ROH )
    echo "Calculate ROH....with BEAGLE and separate population files (no maf filtering)"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/INGI/FVG
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/INGI/VBI
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/TGP/TSI
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/TGP/CEU
            ;;
      esac
        #use freq data
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1230.jar gl=${pop_path}/${CHR}.vcf.gz ibd=true nthreads=2 window=${window} overlap=${overlap} out=${pop}.roh
      
      done

    ;;
    ROH2 )
    echo "Calculate ROH from a unified vcf file....with BEAGLE...we need a file without missing genotypes(NO MAF filter)!!"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
      esac
        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/listpop/all_pop_but_${pop}.txt
        # marker_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140520_ROH/sites_with_missing_genotypes.list
        #use freq data
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1230.jar gtgl=${pop_path}/22.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} excludemarkers=${marker_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1230.jar gt=${pop_path}/22.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1230.jar gt=${pop_path}/${CHR}.nonmissing.maf_gt_05.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${pop}.roh
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1230.jar gt=${pop_path}/${CHR}.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${pop}.roh
      
      done

    ;;
    ROH3 )
    echo "Calculate ROH from a unified vcf file....with BEAGLE...we need a file without missing genotypes(currently is filtered on MAF>5%)!!"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
      esac
        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/listpop/all_pop_but_${pop}.txt
        # marker_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140520_ROH/sites_with_missing_genotypes.list
        #use freq data
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1230.jar gtgl=${pop_path}/22.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} excludemarkers=${marker_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1230.jar gt=${pop_path}/22.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} out=${pop}.roh
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1230.jar gt=${pop_path}/${CHR}.nonmissing.maf_gt_05.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${pop}.roh
      
      done

    ;;
    ROH4 )
    echo "Calculate ROH from different files for each population....using BEAGLE4!(filtered by MAF!!))"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/INGI/FVG/MAF_gt05
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/INGI/VBI/MAF_gt05
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/TGP/TSI/MAF_gt05
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/TGP/CEU/MAF_gt05
            ;;
      esac
        #use freq data
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1230.jar gt=${pop_path}/${CHR}.vcf.gz ibd=true nthreads=2 window=${window} overlap=${overlap} out=${pop}.roh
      
      done

    ;;
    ROH5 )
    echo "Calculate ROH from a unified vcf file....using plink! (NO MAF FILTERING)"
    echo -e "Parameters: \nwindow=${window}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
      esac
        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/listpop/${pop}.keeplist
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --vcf ${pop_path}/${CHR}.nonmissing.vcf.gz --keep ${pop_list} --homozyg --homozyg-window-snp ${window} --out ${pop}.roh
      
      done

    ;;
    ROH5 )
    echo "Calculate ROH from a unified vcf file....using plink! (with MAF FILTERING)"
    echo -e "Parameters: \nwindow=${window}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
      esac
        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/listpop/${pop}.keeplist
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --vcf ${pop_path}/${CHR}.nonmissing.maf_gt_05.vcf.gz --keep ${pop_list} --homozyg --homozyg-window-snp ${window} --out ${pop}.roh
      
      done

    ;;
    ROH6 )
    echo "Calculate ROH from different files for each population....using PLINK!"
    echo -e "Parameters: \nwindow=${window}"
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/INGI/FVG
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/INGI/VBI
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/TGP/TSI
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POPULATIONS/TGP/CEU
            ;;
      esac
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --vcf ${pop_path}/${CHR}.vcf.gz --homozyg --homozyg-window-snp ${window} --out ${pop}.roh
      
      done

    ;;
  esac