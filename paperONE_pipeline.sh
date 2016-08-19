#!/usr/local/bin/bash
#Script to create a pipeline for the purge through inbreeding work
#
#
#define the population we are using
pops="FVG VBI TSI CEU CARL"
# pops="FVG VBI"

#retrieve the MODE parameter to select the correct operation
MODE=$1
CHR=$2
# define the output folder relative to the chromosome, if specified
mkdir -p LOGS

case $MODE in
  HOMDAC)
  sample=$1
  list_path=$4
  cat=$5
  pop=$6
  if [[ -d "${outdir}" ]]; then
  rm -r ${outdir}
  fi
  ;;
  ALTDAC)
  sample=$1
  list_path=$4
  cat=$5
  pop=$6
  if [[ -d "${outdir}" ]]; then
  rm -r ${outdir}
  fi
  ;;
  EXTRFREQ )
  #variant type
  variant=$3

  ;;
  HWECLEAN)
  module add hgi/vcftools/0.1.14
  module add hgi/zlib/1.2.8
  hwe=$3
  q=$4
  mem=$5
  pops="FVG VBI CARL"
  ;;
  INGIROH*)
    window=$3
    overlap=$4
    ld_win=$5
    r2=$6
    mac=$7
    outdir=$8/CHR${CHR}
    mkdir -p ${outdir}
    ;;  
  ROH*)
    #set parameters for beagle:
    window=$3
    overlap=$4
    outdir=$5/CHR${CHR}
    mkdir -p ${outdir}
    ;;
  SPLITCSQ*)
    #set parameters for splitting in consequences
    category=$3
    fixed=$4
    outdir=CHR${CHR}
    mkdir -p ${outdir}
    ;;
  SHARED )
    #set parameters for file input/output
    in_dir=$3
    outdir=CHR${CHR}
    mkdir -p ${outdir}
    ;;
  SHAREDVILL )
    #set parameters for file input/output
    in_file=$3
    outdir=$4/CHR${CHR}
    mkdir -p ${outdir}

    ;;
  MERGEROH )
      #set parameters for file input/output
      in_dir=$3
      outdir=CHR${CHR}
      mkdir -p ${outdir}
      ;;
  IBDCLUST )
      #set parameters for file input/output
      win_size=$3
      density=$4
      in_dir=$5
      outdir=CHR${CHR}
      mkdir -p ${outdir}
      ;;
  MANCLUSTIBD )
      in_dir=$3
      outdir=CHR${CHR}
      ;;
  MAFSPEC )
      #set parameters for file input/output
      input_file=$3
      maf_file_path=$4
      outdir=CHR${CHR}
      mkdir -p ${outdir}
    ;;
  GERMLINE )
    # set up args for germline command line
    MATCH=$3
    HOM=$4
    HET=$5
    BITS=$6
    outdir=CHR${CHR}
    mkdir -p ${outdir}
  ;;
  GERMLINEREG )
    # set up args for germline command line, no need for CHR!
    MATCH=$2
    HOM=$3
    HET=$4
    BITS=$5
    rm -r ${outdir}
  ;;
  PRISHCOUNT )
    # set up args for germline command line, no need for CHR!
    in_path=$3
    outdir=CHR${CHR}
    mkdir -p ${outdir}
  ;;
  EXTRMAF )
    pop_path=$3
    # in_path=$4
    outdir=CHR${CHR}
    mkdir -p ${outdir}
  ;;
  RANDLIST )
  ROUNDS=$2
  SIZE=$3
  outdir=CHR${CHR}
  mkdir -p ${outdir}
  ;;
  MERGEVCF)

  ;;
esac

# Merge different popuplation together
# TODO: add code!!!



#create a case statement to select the operation to do:
case $MODE in
  MERGEVCF)
  #17/08/2016
  #merge files for IBD/ROH calc for SIGU slides
  # chr=${file}
  module purge
  module add hgi/bcftools/1.3.1
  module add hgi/htslib/1.3.1
  module add hgi/samtools/1.3.1


  TGP3_path=/lustre/scratch114/resources/1000g/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
  CARL_path=/lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL/${CHR}.vcf.gz
  FVG_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL/${CHR}.vcf.gz
  VBI_path=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/${CHR}.vcf.gz
  
  mkdir -p /lustre/scratch113/projects/carl_seq/08072016_paperONE/TGP3_INGI_MERGE/${CHR}

  bsub -J"merge_vcf_${CHR}" -o"%J_merge_vcf_${CHR}.o" -q long -M8000 -R"select[mem>=8000] rusage[mem=8000]" -- bcftools merge -m both ${TGP3_path} ${CARL_path} ${FVG_path} ${VBI_path} -O z -o /lustre/scratch113/projects/carl_seq/08072016_paperONE/TGP3_INGI_MERGE/${CHR}/${CHR}.vcf.gz

  ;;
  RANDLIST )
  #create random list for subsampling
  for pop in CARL CEU Erto FVG Illegio Resia Sauris TSI VBI
  do 
    for i in {1..10}
    do 
        mkdir -p /lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES_rounds/${i}/listpop
        sort -R /lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/${pop}_unrelated.list | head -n ${SIZE} > /lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES_rounds/${i}/listpop/${pop}_unrelated.list
    done
  done
  ;;
  ROHENR)
    # check if ROH are enriched in deleterious mutations
    #first select
  ;;
  MERGEROH )
    #merge all ROH to plot a cumulative sites number for each population genome wide
    #we'll work on all populations, included FVG's villages
    pops_updated="FVG VBI TSI CEU CARL Erto Resia Illegio Sauris"
    for pop in $pops_updated
    do
      case $pop in
        FVG )
          pop_path=${in_dir}/CHR${CHR}
          ;;
        VBI )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        TSI )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        CEU )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        CARL )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        Erto )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        Sauris )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        Illegio )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        Resia )
          pop_path=${in_dir}/CHR${CHR}
            ;;
      esac
      # While adding ROH length info to hbd BEAGLE output, filter out with min LOD of 4 and 5 to control for sequence error rate
      awk '{if($1==$3 && $8 >= 4) print $0,$7-$6}' ${pop_path}/${pop}.roh.hbd| tr "\t" " " > ${pop_path}/${pop}.roh.length.4.hbd
      awk '{if($1==$3 && $8 >= 5) print $0,$7-$6}' ${pop_path}/${pop}.roh.hbd| tr "\t" " " > ${pop_path}/${pop}.roh.length.5.hbd

      #While adding length info to IBD data from BEAGLE4, joining col 1 and col 3 to create a unique identifier for the couple and filter out with min LOD of 4 and 5 to control for sequence error rate
      awk '{if($8 >= 4) print $0,$1"_"$3,$7-$6}' ${pop_path}/${pop}.roh.ibd| tr "\t" " " > ${pop_path}/${pop}.roh.length.4.ibd
      awk '{if($8 >= 5) print $0,$1"_"$3,$7-$6}' ${pop_path}/${pop}.roh.ibd| tr "\t" " " > ${pop_path}/${pop}.roh.length.5.ibd

      #now take all the files and join them together to create a whole genome ROH and IBD file for each population/village
      #create a population folder
      mkdir -p ${pop}
      #ROH
      cat ${pop_path}/${pop}.roh.length.4.hbd >> ${pop}/${pop}.WG.roh.length.4.hbd
      cat ${pop_path}/${pop}.roh.length.5.hbd >> ${pop}/${pop}.WG.roh.length.5.hbd
      
      #IBD
      cat ${pop_path}/${pop}.roh.length.4.ibd >> ${pop}/${pop}.WG.roh.length.4.ibd
      cat ${pop_path}/${pop}.roh.length.5.ibd >> ${pop}/${pop}.WG.roh.length.5.ibd

  done
  ;;
  MANCLUSTIBD )
    # manual ibd clustering
    # 1) define window legth: default 1Mb and LOD to filter: 5 (default)
    win=1000000
    LOD=5
    pop_path=${in_dir}/CHR${CHR}
    #remove useless chr folders
    rm -r ${outdir}

    pops_updated="FVG VBI TSI CEU CARL Erto Resia Illegio Sauris"
    for pop in $pops_updated
    do
      mkdir -p ${pop}
      #we need to prepare the input files:
      #first we need to sort by position
      # sort -g -k6,6 -k7,7 ${pop_path}/${pop}.roh.length.${LOD}.ibd > ${outdir}/${pop}.roh.length.${LOD}.sorted.ibd
      sort -g -k6,6 -k7,7 ${pop_path}/${pop}.roh.length.${LOD}.ibd > ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.sorted.ibd

      #we should need also a map file with starting and ending chromosome position...or we can start with the first position
      #in our result ibd file 
      #we'll get our position from the original vcf files used for ibd calculation
      in_vcf=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING/${CHR}.non_missing.vcf.gz
      start=`tabix ${in_vcf} ${CHR}|head -1| cut -f 2`
      end=`tabix ${in_vcf} ${CHR}|tail -1| cut -f 2`

      # define start_w and end_w initial values
      start_w=$start
      end_w=$[$start_w + $win]
      w_n=0
      #now we need to count how many segments we have in a window:
      while [[ $end_w -le $end ]]; do
        #for a segment to be in a window we check only the starting position: if the segment starts in my window,I count it in that window
        # awk -v start=$start_w -v end=$end_w '{if($6 <= end && $6 >= start) print $0}' ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.sorted.ibd > ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.W${w_n}.ibd
        
        # we need to calculate the average recombination rate per window: we'll dot this using the genetic map from TGP
        # because should be enough resolution already
        # select the first and the and last line from the map, anc calculate the mean value from the second column
        avg_rec=`awk -v start=$start_w -v end=$end_w '{if($1 >= start && $1 <= end ) print $0 }' /nfs/team151/reference/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr${CHR}_combined_b37.txt | awk '{sum+=$2;next}END{print sum/NR}'`

        # for a segment to be in a window we check the starting position and how mutch overlap there is with the windows in term of base pairs: 
        # if the segment dont start in the windows but there is still an overlap, I count it in that window
        awk -v start=$start_w -v end=$end_w -v avg_rr=$avg_rec '
        {if(($6 < start) && ($7 >= start && $7 <= end))
          {overlap=$7-start;class=1;}
        else if(($6 >= start && $6 <= end) && ($7 >= start && $7 <= end))
          {overlap=($7-$6);class=2;}
        else if(($6 >= start && $6 <= end) && ($7 > end))
          {overlap=end-$6;class=3;}
        else if(($6 < start && $7 > end))
          {overlap=end-start;class=4;}
        }
        {if(($6 <= end && $6 >= start) || ($7 <= end && $7 >= start)) print $0,overlap,class,avg_rr}' ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.sorted.ibd | tr " " "\t"> ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.W${w_n}.${start_w}_${end_w}.ibd
        
        #lets do a resume of all the pair/samples in each windows: we need to get the last column, sort it, select uniq values and split by "_"
        cut -f 9 -d " " ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.W${w_n}.${start_w}_${end_w}.ibd |sort| uniq| tr "_" "\t" > ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.W${w_n}.pair_file
        (cut -f 1 ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.W${w_n}.pair_file;cut -f 2 ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.W${w_n}.pair_file)| sort | uniq > ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.W${w_n}.sample_file
        
        # create a report for all windows for that population
        lines=`wc -l ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.W${w_n}.${start_w}_${end_w}.ibd| cut -f 1 -d " "`
        echo "${pop}	${w_n}	${start_w}	${end_w} ${lines}	${avg_rec}" >> ${pop}/${pop}.chr${CHR}.roh.length.${LOD}.resume.txt
        w_n=$[$w_n + 1]
        start_w=$[$start_w + $win]
        end_w=$[$start_w + $win]
      done

    done
  ;;
  IBDCLUST )
  #need to work by chromosome!!because it is right...we should have different cluster and sample in cluster, for different chrs!!
  echo "Launch IBD clustering software EMI"
  echo "Using win size: ${win_size}"
  echo "Cluster Density: ${density}"

  pops_updated="FVG VBI TSI CEU CARL Erto Resia Illegio Sauris"
    for pop in $pops_updatedbsub -J"inb_${pop}" -o"%J_inb_${pop}.o" -w "ended(freq_${pop})" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- 
    do
      case $pop in
        FVG )
          pop_path=${in_dir}/CHR${CHR}
          ;;
        VBI )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        TSI )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        CEU )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        CARL )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        Erto )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        Sauris )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        Illegio )
          pop_path=${in_dir}/CHR${CHR}
            ;;
        Resia )
          pop_path=${in_dir}/CHR${CHR}
            ;;
      esac
      #we need to prepare the input files:
      # File format transformation: one for each chromosome 
      awk '{if($8>=4) print $1,$1"."$2-1,$3,$3"."$4-1,$6,$7,$8,0,0}' ${pop_path}/${pop}.roh.length.5.hbd | sort -k 5,1n -k 6,2n > ${outdir}/${pop}.roh.length.5.ibd
      awk '{if($8>=5) print $1,$1"."$2-1,$3,$3"."$4-1,$6,$7,$8,0,0}' ${pop_path}/${pop}.roh.length.5.ibd | sort -k 5,1n -k 6,2n > ${outdir}/${pop}.ibd.length.5.ibd

      #fam file: one for each chromosome
      cat <(cut -f 1 -d  " " ${outdir}/${pop}.roh.length.5.ibd| sort|uniq) <(cut -f 3 -d  " " ${outdir}/${pop}.roh.length.5.ibd| sort|uniq) | sort | uniq | awk '{print $1,$1,0,0,0,-9}' > ${outdir}/${pop}.roh.length.5.fam
      cat <(cut -f 1 -d  " " ${outdir}/${pop}.ibd.length.5.ibd| sort|uniq) <(cut -f 3 -d  " " ${outdir}/${pop}.ibd.length.5.ibd| sort|uniq) | sort | uniq | awk '{print $1,$1,0,0,0,-9}' > ${outdir}/${pop}.ibd.length.5.fam

      # now run EMI for the chromosome with the selected data --> FOR IBD:
      mkdir -p EMI/w_${win_size}k/${density}/${pop}; emi ${outdir}/${pop}.ibd.length.5.ibd -fam ${outdir}/${pop}.ibd.length.5.fam -den ${density} -wgt 7th 5 50 -win $[win_size*1000] bp EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.ibd.out

      # format output file to sort cluster and remove duplicates:
      awk -F "\t" -v OFS="\t" '{$1="clst"; print $0}' EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.ibd.out.clst.tmp|sort -k2n -k3n -k4n | uniq| awk '{$1="clst"NR;print $0}' > EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.ibd.out.clst.uniq.sorted.txt
  
      # Create a file with cluster positions and number of samples:
      awk -v chr=${CHR} '{print $1,chr,$2,$3,(NF-3)/2}' EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.ibd.out.clst.uniq.sorted.txt > EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.ibd.out.clst.uniq.size

      # now run EMI for the chromosome with the selected data --> FOR ROH:
      mkdir -p EMI/w_${win_size}k/${density}/${pop}; emi ${outdir}/${pop}.roh.length.5.ibd -fam ${outdir}/${pop}.roh.length.5.fam -den ${density} -wgt 7th 5 50 -win $[win_size*1000] bp EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.roh.out

      # format output file to sort cluster and remove duplicates:
      awk -F "\t" -v OFS="\t" '{$1="clst"; print $0}' EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.roh.out.clst.tmp|sort -k2n -k3n -k4n | uniq| awk '{$1="clst"NR;print $0}' > EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.roh.out.clst.uniq.sorted.txt
  
      # Create a file with cluster positions and number of samples:
      awk -v chr=${CHR} '{print $1,chr,$2,$3,(NF-3)/2}' EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.roh.out.clst.uniq.sorted.txt > EMI/w_${win_size}k/${density}/${pop}/${pop}.chr${CHR}.roh.out.clst.uniq.size

    done
  ;;
  HET )
    #Extract Inbreeding coeff information for each population. Here we are using plink2 (1.9)
    #we'll use the --het option as long as the --ibc option
    module add hgi/plink/1.90b3w

    echo "Calculate Inbreeding for all pops in a vcf files...."
    pop_path=$2
    pop_count=$3
    pop_base=$4
    outdir=${pop_base}/results
    mkdir -p ${outdir}

    pops_updated=$(cut -f 2 -d " " ${pop_count})
    for pop in $pops_updated
    do
        poplist=${pop_base}/lists/${pop}.ped
        #calculate frequencies, before:
        # bsub -J"freq_${pop}" -o"${outdir}/%J_freq_${pop}.o" -q normal -M4000 -n2 -R"span[hosts=1] select[mem>=4000] rusage[mem=4000]" -- plink --vcf ${pop_path} --biallelic-only --double-id --keep-allele-order --snps-only --keep ${poplist} --freq --nonfounders --out ${outdir}/freq_${pop}
        
        #use freq data
        bsub -J"inb_${pop}" -o"${outdir}/%J_inb_${pop}.o" -q normal -M4000 -n2 -R"span[hosts=1] select[mem>=4000] rusage[mem=4000]" -- plink --vcf ${pop_path} --biallelic-only --double-id --keep-allele-order --snps-only --keep ${poplist} --het --out ${outdir}/inb_${pop}
        bsub -J"ibc_${pop}" -o"${outdir}/%J_ibc_${pop}.o" -q normal -M4000 -n2 -R"span[hosts=1] select[mem>=4000] rusage[mem=4000]" -- plink --vcf ${pop_path} --biallelic-only --double-id --keep-allele-order --snps-only --keep ${poplist} --ibc --out ${outdir}/ibc_${pop}
    done
  ;;

  SHARED )
    #calculate shared and private sites for all populations
    pops_updated="FVG VBI CARL TSI CEU Erto Resia Illegio Sauris"
    for pop in $pops_updated
    do
      case $pop in
        FVG )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          #Private calculation without fixed sites
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$8 !="NA" && $8 > 0 && $8 < 1 && $6 == 0 && $5 == 0'| gzip -c > ${outdir}/FVG_private_chr${CHR}.merged_daf.tab.gz
          # Fixed in Isolate and private
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$8 !="NA" && $8 == 1 && $6 == 0 && $5 == 0'| gzip -c > ${outdir}/FVG_private_chr${CHR}.merged_daf.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$8 !="NA" && $8 > 0 && $8 < 1 && ($6 > 0 || $5 > 0)'| gzip -c > ${outdir}/FVG_shared_chr${CHR}.merged_daf.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$8 !="NA" && $8 == 1 && ($6 > 0 || $5 > 0)'| gzip -c > ${outdir}/FVG_shared_chr${CHR}.merged_daf.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$8 =="NA" && ($6 == "NA" && $5 == "NA")'| gzip -c > ${outdir}/FVG_novel_chr${CHR}.merged_daf.tab.gz
          ;;
        VBI )
          #Private calculation without fixed sites
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$7 !="NA" && $7 > 0 && $7 < 1 && $6 == 0 && $5 == 0'| gzip -c > ${outdir}/VBI_private_chr${CHR}.merged_daf.tab.gz
          # Fixed in Isolate and private
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$7 !="NA" && $7 == 1 && $6 == 0 && $5 == 0'| gzip -c > ${outdir}/VBI_private_chr${CHR}.merged_daf.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$7 !="NA" && $7 > 0 && $7 < 1 && ($6 > 0 || $5 > 0)'| gzip -c > ${outdir}/VBI_shared_chr${CHR}.merged_daf.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$7 !="NA" && $7 == 1 && ($6 > 0 || $5 > 0)'| gzip -c > ${outdir}/VBI_shared_chr${CHR}.merged_daf.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$7 =="NA" && ($6 == "NA" && $5 == "NA")'| gzip -c > ${outdir}/VBI_novel_chr${CHR}.merged_daf.tab.gz
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        CARL )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          #Private calculation without fixed sites
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$9 !="NA" && $9 > 0 && $9 < 1 && $6 == 0 && $5 == 0'| gzip -c > ${outdir}/CARL_private_chr${CHR}.merged_daf.tab.gz
          # Fixed in Isolate and private
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$9 !="NA" && $9 == 1 && $6 == 0 && $5 == 0'| gzip -c > ${outdir}/CARL_private_chr${CHR}.merged_daf.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$9 !="NA" && $9 > 0 && $9 < 1 && ($6 > 0 || $5 > 0)'| gzip -c > ${outdir}/CARL_shared_chr${CHR}.merged_daf.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$9 !="NA" && $9 == 1 && ($6 > 0 || $5 > 0)'| gzip -c > ${outdir}/CARL_shared_chr${CHR}.merged_daf.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_dir}/CHR${CHR}/INGI_chr${CHR}.merged_daf.tab.gz | awk '$9 =="NA" && ($6 == "NA" && $5 == "NA")'| gzip -c > ${outdir}/CARL_novel_chr${CHR}.merged_daf.tab.gz
          ;;
        Erto )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        Resia )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        Illegio )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        Sauris )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
      esac
    done
  ;;
  SHAREDVILL )
    #calculate shared and private sites for all populations with fvg splitted by villages
    # $4 CEU
    # $5 TSI
    # $6 VBI
    # $7 FVG
    # $8 CARL
    # $9 Erto
    # $10 Illegio
    # $11 Resia
    # $12 Sauris

    # pops_updated="FVG VBI CARL TSI CEU Erto Resia Illegio Sauris"
    pops_updated="TSI CEU"
    for pop in $pops_updated
    do
      case $pop in
        TSI )
          #Private calculation without fixed sites
          # zcat ${in_file} | awk '$6 !="NA" && $6 > 0 && $6 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          zcat ${in_file} | cut -f 1,3-| awk '$5 !="NA" && $5 > 0 && ($6 == 0 || $7 == 0 || $8 == 0 || $9 == 0 || $10 == 0 || $11 == 0 || $12 == 0)'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # # Fixed in Isolate and private
          # zcat ${in_file} | awk '$6 !="NA" && $6 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # # Shared calculation without fixed sites
          zcat ${in_file} | cut -f 1,3-| awk '$5 !="NA" && $5 > 0 && $5 < 1 && ($6 > 0 || $7 > 0 || $8 > 0 || $9 > 0 || $10 > 0 || $11 > 0 || $12 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # # Fixed in Isolate and shared
          # zcat ${in_file} | awk '$6 !="NA" && $6 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # # NOVEL in my population
          zcat ${in_file}| cut -f 1,3- | awk '$5 !="NA" && ($6 == "NA" || $7 == "NA" || $8 == "NA" || $9 == "NA" || $10 == "NA" || $11 == "NA" || $12 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        CEU )
          #Private calculation without fixed sites
          # zcat ${in_file} | awk '$4 !="NA" && $4 > 0 && $6 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          zcat ${in_file}| cut -f 1,3- | awk '$4 !="NA" && $4 > 0 && ($6 == 0 || $7 == 0 || $8 == 0 || $9 == 0 || $10 == 0 || $11 == 0 || $12 == 0)'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # # Fixed in Isolate and private
          # zcat ${in_file} | awk '$4 !="NA" && $4 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # # Shared calculation without fixed sites
          zcat ${in_file}| cut -f 1,3-| awk '$4 !="NA" && $4 > 0 && $4 < 1 && ($6 > 0 || $7 > 0 || $8 > 0 || $9 > 0 || $10 > 0 || $11 > 0 || $12 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # # Fixed in Isolate and shared
          # zcat ${in_file} | awk '$4 !="NA" && $4 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in my population
          zcat ${in_file}| cut -f 1,3- | awk '$4 !="NA" && ($6 == "NA" || $7 == "NA" || $8 == "NA" || $9 == "NA" || $10 == "NA" || $11 == "NA" || $12 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        VBI )
          #Private calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$6 !="NA" && $6 > 0 && $6 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file}| cut -f 1,3- | awk '$6 !="NA" && $6 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$6 !="NA" && $6 > 0 && $6 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file}| cut -f 1,3- | awk '$6 !="NA" && $6 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file}| cut -f 1,3- | awk '$6 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        FVG )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          #Private calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$7 !="NA" && $7 > 0 && $7 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file} | cut -f 1,3-| awk '$7 !="NA" && $7 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$7 !="NA" && $7 > 0 && $7 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file}| cut -f 1,3- | awk '$7 !="NA" && $7 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file}| cut -f 1,3- | awk '$7 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          ;;
    
        CARL )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          #Private calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$8 !="NA" && $8 > 0 && $8 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file}| cut -f 1,3- | awk '$8 !="NA" && $8 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$8 !="NA" && $8 > 0 && $8 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file}| cut -f 1,3- | awk '$8 !="NA" && $8 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file}| cut -f 1,3- | awk '$8 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
    
          ;;
        Erto )
          #Private calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$9 !="NA" && $9 > 0 && $9 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file}| cut -f 1,3- | awk '$9 !="NA" && $9 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$9 !="NA" && $9 > 0 && $9 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file}| cut -f 1,3- | awk '$9 !="NA" && $9 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file}| cut -f 1,3- | awk '$9 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          ;;
        Illegio )
         #Private calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$10 !="NA" && $10 > 0 && $10 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file}| cut -f 1,3- | awk '$10 !="NA" && $10 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$10 !="NA" && $10 > 0 && $10 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file}| cut -f 1,3- | awk '$10 !="NA" && $10 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file}| cut -f 1,3- | awk '$10 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz 
          ;;
        Resia )
          #Private calculation without fixed sites
          #Private calculation without fixed sites
          zcat ${in_file} | cut -f 1,3-| awk '$11 !="NA" && $11 > 0 && $11 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file}| cut -f 1,3- | awk '$11 !="NA" && $11 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$11 !="NA" && $11 > 0 && $11 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file}| cut -f 1,3- | awk '$11 !="NA" && $11 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file} | cut -f 1,3-| awk '$11 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          ;;
        Sauris )
         #Private calculation without fixed sites
          zcat ${in_file} | cut -f 1,3-| awk '$12 !="NA" && $12 > 0 && $12 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file}| cut -f 1,3- | awk '$12 !="NA" && $12 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file}| cut -f 1,3- | awk '$12 !="NA" && $12 > 0 && $12 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file}| cut -f 1,3- | awk '$12 !="NA" && $12 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file}| cut -f 1,3- | awk '$12 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz 
          ;;
        
      esac
      echo ${in_file}
      echo ${pop}
      echo ${CHR}
      echo ${outdir}
    done
    
  ;;
  PRISHCOUNT )
    #extract data from different files categories and count to summarize all in Isolates
    echo "Count sites...."
    for pop in $pops
    do
      case $pop in
        FVG )
          # pop_path_pref=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP/CHR${CHR}/${pop}
          ;;
        VBI )
          # pop_path_pref=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP/CHR${CHR}/${pop}
          ;;
        CARL )
          # pop_path_pref=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP/CHR${CHR}/${pop}
          ;;
      esac
      pop_path_pref=${in_path}/CHR${CHR}/${pop}
      #create a folder for results
      mkdir -p ${pop}

      zcat ${pop_path_pref}_private_chr${CHR}.merged_daf.tab.gz >> ${pop}/${pop}_private.merged_daf.tab
      zcat ${pop_path_pref}_private_chr${CHR}.merged_daf.fixed.tab.gz >> ${pop}/${pop}_private.merged_daf.fixed.tab
      zcat ${pop_path_pref}_shared_chr${CHR}.merged_daf.tab.gz >> ${pop}/${pop}_shared.merged_daf.tab
      zcat ${pop_path_pref}_shared_chr${CHR}.merged_daf.fixed.tab.gz >> ${pop}/${pop}_shared.merged_daf.fixed.tab
      zcat ${pop_path_pref}_novel_chr${CHR}.merged_daf.tab.gz >> ${pop}/${pop}_novel.merged_daf.tab

    done
  ;;
  DACMAF )
    #extract data in bed format for different populations in a separate way
    echo "Create bed formatted files..."
    pops_updated="Erto Resia Illegio Sauris"
    # for pop in $pops_updated
    for pop in $pops
    do
      case $pop in
        FVG )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/FVG_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/80_SAMPLES/listpop/FVG_unrelated.list
          ;;
        VBI )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/VBI_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/80_SAMPLES/listpop/VBI_unrelated.list
          ;;
        CARL )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/CARL_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/80_SAMPLES/listpop/CARL_unrelated.list
          ;;
        TSI )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/TSI.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/80_SAMPLES/listpop/TSI_unrelated.list
          ;;
        CEU )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/CEU.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/80_SAMPLES/listpop/CEU_unrelated.list
          ;;
        Erto )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Erto_unrelated.list
            ;;
        Sauris )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Sauris_unrelated.list
            ;;
        Illegio )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Illegio_unrelated.list
            ;;
        Resia )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Resia_unrelated.list
            ;;
      esac
      in_vcf=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140711_ANNOTATED/${CHR}.vcf.gz
      # in_vcf=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POP_MERGED_FILES/FIVE_POPS/20140730_ANNOTATED/${CHR}.clean_annotated.vcf.gz
      out_tab=${pop}.chr${CHR}.tab
      #create files for each population for each chromosome in a separate folder
      echo "ec_dacmacdafmaf2bed.py ${pop_path} ${in_vcf} ${out_tab}; mv ${pop}.chr${CHR}.tab ${outdir}/;gzip ${outdir}/${pop}.chr${CHR}.tab" | bsub -J"dac_exract_${CHR}_${pop}" -o"%J_dac_exract_${CHR}_${pop}.o" -M3000 -R"select[mem>=3000] rusage[mem=3000]" -q normal
      # echo "mv ${pop}.chr${CHR}.tab ${outdir}/;gzip ${outdir}/${pop}.chr${CHR}.tab" | bsub -J"dac_exract_${CHR}_${pop}" -o"%J_dac_exract_${CHR}_${pop}.o" -M3000 -R"select[mem>=3000] rusage[mem=3000]" -q normal

    done
  ;;
  EXTRMAF )
    #extract data in bed format for different populations in a separate way
    echo "Extract from tables for each population, complete info about shared/private sites"
    pops_updated="FVG VBI CARL TSI CEU Erto Resia Illegio Sauris"
    for pop in $pops_updated="FVG VBI CARL TSI CEU Erto Resia Illegio Sauris"
    do
      #those are the files I need to use to extract the info splitted by chr
      # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
      
      private=${pop_path}/${pop}/${pop}_private.merged_daf.tab
      private_fixed=${pop_path}/${pop}/${pop}_private.merged_daf.fixed.tab
      shared=${pop_path}/${pop}/${pop}_shared.merged_daf.tab
      shared_fixed=${pop_path}/${pop}/${pop}_shared.merged_daf.fixed.tab

      # #create the list of variants we need to extract
      # awk '{print $1"O"$2"O"$3}' ${private} | tr " " "\t" | dos2unix > ${private}.list
      # awk '{print $1"O"$2"O"$3}' ${private_fixed} | tr " " "\t" | dos2unix > ${private_fixed}.list
      # awk '{print $1"O"$2"O"$3}' ${shared} | tr " " "\t" | dos2unix > ${shared}.list
      # awk '{print $1"O"$2"O"$3}' ${shared_fixed} | tr " " "\t" | dos2unix > ${shared_fixed}.list

      #create temporary files to do the extraction
      in_file=${pop_path}/CHR${CHR}/INGI_chr${CHR}.merged_maf.tab.gz
      # zcat ${in_file} | awk '{print $1"O"$2"O"$3,$0}' | tr " " "\t" | dos2unix > ${in_file}.tmp

      zcat ${in_file} |awk '{print $1"O"$2"O"$3"O"$4,$0}'|tr " " "\t" | fgrep -w -f <(awk -v cro=${CHR} '{if($1==cro) print $1"O"$2"O"$3"O"$4}' ${private})|cut -f 2- |gzip -c > ${in_file}.${pop}.private.tab.gz
      zcat ${in_file} |awk '{print $1"O"$2"O"$3"O"$4,$0}'|tr " " "\t" | fgrep -w -f <(awk -v cro=${CHR} '{if($1==cro) print $1"O"$2"O"$3"O"$4}' ${private_fixed})|cut -f 2- |gzip -c > ${in_file}.${pop}.private_fixed.tab.gz
      zcat ${in_file} |awk '{print $1"O"$2"O"$3"O"$4,$0}'|tr " " "\t" | fgrep -w -f <(awk -v cro=${CHR} '{if($1==cro) print $1"O"$2"O"$3"O"$4}' ${shared})|cut -f 2- |gzip -c > ${in_file}.${pop}.shared.tab.gz
      zcat ${in_file} |awk '{print $1"O"$2"O"$3"O"$4,$0}'|tr " " "\t" | fgrep -w -f <(awk -v cro=${CHR} '{if($1==cro) print $1"O"$2"O"$3"O"$4}' ${shared_fixed})|cut -f 2- |gzip -c > ${in_file}.${pop}.shared_fixed.tab.gz
      
      # #now grep the file to extract the data we need, using different lists
      # (fgrep -w -f <(grep "^${CHR}" ${private}.list) ${in_file}.tmp)| cut -f 2- | dos2unix |gzip -c > ${in_file}.private.tab.gz
      # (fgrep -w -f <(grep "^${CHR}" ${private_fixed}.list) ${in_file}.tmp)| cut -f 2- | dos2unix |gzip -c > ${in_file}.private_fixed.tab.gz
      # (fgrep -w -f <(grep "^${CHR}" ${shared}.list) ${in_file}.tmp)| cut -f 2- | dos2unix |gzip -c > ${in_file}.shared.tab.gz
      # (fgrep -w -f <(grep "^${CHR}" ${shared_fixed}.list) ${in_file}.tmp)| cut -f 2- | dos2unix |gzip -c > ${in_file}.shared_fixed.tab.gz

      # rm ${in_file}.tmp

    done
  ;;
  MAFSPEC )
    #extract MAF data different populations from a given list
    echo "Extract data for af/maf spectrum plot"
    pops_updated="VBI CARL TSI CEU Erto Resia Illegio Sauris"
   
      maf_file=${maf_file_path}/INGI_villages_chr${CHR}.merged_maf.tab.gz
      out_tmp=`basename ${input_file}`

      (zcat ${maf_file} | cut -f 1,3-|head -1;zcat ${maf_file} | cut -f 1,3- | fgrep -w -f <(awk '{print $3}' ${input_file}) )| tr " " "\t" | gzip -c > ${out_tmp}.maf_file.gz
  ;;
  GERMLINE )
    echo "Calculate IBD using GERMLINE from plink formatted files!!"
    echo -e "Parameters: \nmin_match=${MATCH}\nerr_hom=${HOM}\nerr_het=${HET}\nbits=${BITS}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    # pops_updated="FVG VBI TSI CEU CARL Erto Resia Illegio Sauris"
    pops_updated="ALL"
    for pop in $pops_updated
    do

        # ped_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE
        ped_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/FILTERED
        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/FIVE_POPS/all_pop_but_${pop}.removelist
        
        
        # #we need to create the input file for GERMLINE (WTF!!)
        echo "1" > ${pop}.${CHR}.run
        echo "${ped_path}/${pop}.${CHR}.non_missing.filtered.map" >> ${pop}.${CHR}.run
        echo "${ped_path}/${pop}.${CHR}.non_missing.filtered.ped" >> ${pop}.${CHR}.run
        echo "${outdir}/${pop}.${CHR}.non_missing.filtered" >> ${pop}.${CHR}.run

        # commented to use the ALL population files
        # echo "germline -min_m ${MATCH} -err_hom ${HOM} -err_het ${HET} -bits ${BITS} -h_extend -homoz -from_snp rs62224610 -to_snp rs7410320 < ${pop}.${CHR}.run" | bsub -J"LOGS/ibd_${pop}_${CHR}" -o"LOGS/%J_ibd_${pop}_${CHR}.o" -q basement -M8000 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]"
        echo "germline -min_m ${MATCH} -err_hom ${HOM} -err_het ${HET} -bits ${BITS} -h_extend -homoz < ${pop}.${CHR}.run" | bsub -J"LOGS/ibd_${pop}_${CHR}" -o"LOGS/%J_ibd_${pop}_${CHR}.o" -q basement -M8000 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]"
    done
  ;;
  GERMLINEREG )
    echo "Calculate IBD using GERMLINE from plink formatted files!!"
    echo -e "Parameters: \nmin_match=${MATCH}\nerr_hom=${HOM}\nerr_het=${HET}\nbits=${BITS}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    # pops_updated="FVG VBI TSI CEU CARL Erto Resia Illegio Sauris"
    pops_updated="ALL"
    for pop in $pops_updated
    do

        ped_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/HOM2PED
        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/FIVE_POPS/all_pop_but_${pop}.removelist
        
        for filename in `ls ${ped_path}/*.ped`; do
          #we need to create the input file for GERMLINE (WTF!!)
          file=`basename ${filename}`
          echo "1" > ${file}.run
          echo "${filename%*.ped}.map" >> ${file}.run
          echo "${filename}" >> ${file}.run
          echo "${file}.filtered" >> ${file}.run

          echo "germline -min_m ${MATCH} -err_hom ${HOM} -err_het ${HET} -bits ${BITS} -h_extend -homoz < ${file}.run" | bsub -J"LOGS/ibd_${file}" -o"LOGS/%J_ibd_${file}.o" -q normal -M8000 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]"
        done
    done

  ;;
  HWECLEAN )
    echo "Clean data using VCFTOOLS filtering by HWE p"
    echo -e "Parameters: \nhwe threshold=${hwe}"
    #We'll use PLINK formatted files, which we will reconvert back to VCF after filtering
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL/ALL_FVG_20151113.vcf.gz
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/ALL_VBI_20151113.vcf.gz
            ;;
        CARL )
          pop_path=/lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL/ALL_CARL_20151113.vcf.gz
            ;;
        esac
        #use freq data
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gtgl=${pop_path}/22.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} excludemarkers=${marker_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/22.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} out=${pop}.roh
        log_dir=`dirname ${pop_path}`
        echo "vcftools --gzvcf ${pop_path} --hwe ${hwe} --recode --recode-INFO-all --stdout | bgzip -c > ${pop_path}.hwe_filt.vcf.gz" | bsub -J"filter_hwe_${pop}" -o"${log_dir}/%J_filter_${pop}.o" -q ${q} -M${mem} -R"select[mem>=${mem}] rusage[mem=${mem}]"
        echo "vcftools --gzvcf ${pop_path}.hwe_filt.vcf.gz --het --out ${pop_path}.hwe_filt_HET" | bsub -J"calc_inbreed_${pop}" -o"${log_dir}/%J_calc_inbreed_${pop}.o" -q ${q} -M${mem} -R"select[mem>=${mem}] rusage[mem=${mem}]" -w "ended(filter_hwe_${pop})"
        echo "vcftools --gzvcf ${pop_path}.hwe_filt.vcf.gz --SNPdensity 3000000 --out ${pop_path}.hwe_filt_density" | bsub -J"calc_density_${pop}" -o"${log_dir}/%J_calc_density_${pop}.o" -q ${q} -M${mem} -R"select[mem>=${mem}] rusage[mem=${mem}]" -w "ended(filter_hwe_${pop})"
        # echo "vcftools --gzvcf ${pop_path}.hwe_filt.vcf.gz --SNPdensity 3000000 --out ${pop_path}.hwe_filt_density" | bsub -J"calc_density_${pop}" -o"${log_dir}/%J_calc_density_${pop}.o" -q ${q} -M${mem} -R"select[mem>=${mem}] rusage[mem=${mem}]"
      
    done

  ;;
  ROHVILLAGENEW )
    echo "Calculate ROH from a unified vcf file....with IBDseq...we need a file without missing genotypes(NO MAF filter)!!"
    echo "We'll have also data separate for villages in FVG"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    pops_updated="FVG VBI TSI CEU Erto Resia Illegio Sauris"
    for pop in $pops_updated
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP
            ;;
        CARL )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP
            ;;
        Erto )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP
            ;;
        Sauris )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP
            ;;
        Illegio )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP
            ;;
        Resia )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP
            ;;
      esac

        pop_list=/lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/listpop/HRC_all_samples_but_${pop}.list
        #use freq data
        # bsub -J"roh_${pop}_${CHR}" -o"%J_roh_${pop}_${CHR}.o" -q normal -M8000 -n4 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/beagle.22Feb16.8ef.jar gt=${pop_path}/${chr}.GERP.vcf.gz nthreads=4 excludesamples=${pop_list} out=${outdir}/${pop}.roh
        bsub -J"roh_${pop}_${CHR}" -o"%J_roh_${pop}_${CHR}.o" -q normal -M8000 -n4 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- beagle 5000 gt=${pop_path}/${CHR}.GERP.vcf.gz ibd=true nthreads=4 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
    done

  ;;
  ROHBEAGLE )
    echo "Calculate ROH from different files for each population....using BEAGLE4.1!(filtered by HWE!!))"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    module add hgi/beagle/4.1-22Feb16.8ef
    pops_updated="FVG VBI CARL"
    for pop in $pops_updated
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL/ALL_FVG_20151113.vcf.gz.hwe_filt.vcf.gz
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/ALL_VBI_20151113.vcf.gz.hwe_filt.vcf.gz
            ;;
        CARL )
          pop_path=/lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL/ALL_CARL_20151113.vcf.gz.hwe_filt.vcf.gz
            ;;
      esac
        #use freq data
        bsub -J"roh_${pop}_${CHR}_bg" -o"%J_roh_${pop}_${CHR}_bg.o" -q normal -M8000 -n16 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- beagle 8000 gt=${pop_path} ibd=true ibdlod=1 nthreads=16 window=${window} overlap=${overlap} chrom=${CHR} out=${outdir}/${pop}.roh
    
    done
  ;;
  INGIROH )
    echo "Calculate ROH for INGI populations from splitted vcf files with IBDseq (NO MAF filter)!!"
    # echo "We'll have also data separate for villages in FVG"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}\nr2=${r2}\nMAC=${mac}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    # pops_updated="FVG VBI CARL Erto Resia Illegio Sauris"
    pops_updated="FVG VBI CARL"
    for pop in $pops_updated
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL
            ;;
        CARL )
          pop_path=/lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL
            ;;
        Erto )
          pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL
            ;;
        Sauris )
          pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL
            ;;
        Illegio )
          pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL
            ;;
        Resia )
          pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL
            ;;
      esac

        #we can filter out unrelated after the calculation
        # pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/BEAGLE/all_pop_but_${pop}.txt
        # pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/FIVE_POPS/all_pop_but_${pop}.txt
        # pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/all_pop_but_${pop}.txt
        # pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/listpop/all_pop_but_${pop}.txt
        # marker_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140520_ROH/sites_with_missing_genotypes.list
        #use freq data
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gtgl=${pop_path}/22.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} excludemarkers=${marker_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/22.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.nonmissing.maf_gt_05.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
        # bsub -J"roh_${pop}_${CHR}" -o"%J_roh_${pop}_${CHR}.o" -q normal -M8000 -n4 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz ibd=true nthreads=4 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
        # bsub -J"roh_${pop}_${CHR}" -o"%J_roh_${pop}_${CHR}.o" -q normal -M8000 -n4 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/IBDseq/ibdseq.r1206.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz nthreads=4 excludesamples=${pop_list} out=${outdir}/${pop}.roh
        if [[ $mac -eq 1 ]]; then
          #statements
        bsub -J"roh_${pop}_${CHR}" -o"LOGS/%J_roh_${pop}_${CHR}.o" -q normal -M5000 -n16 -R"span[hosts=1] select[mem>=5000] rusage[mem=5000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/IBDseq/ibdseq.r1206.jar gt=${pop_path}/${CHR}.vcf.gz nthreads=16 ibdlod=1 r2max=${r2} out=${outdir}/${pop}.roh
        else
        bsub -J"roh_${pop}_${CHR}" -o"LOGS/%J_roh_${pop}_${CHR}.o" -q normal -M5000 -n16 -R"span[hosts=1] select[mem>=5000] rusage[mem=5000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/IBDseq/ibdseq.r1206.jar gt=${pop_path}/${CHR}.vcf.gz nthreads=16 ibdlod=1 r2max=${r2} minalleles=${mac} out=${outdir}/${pop}.roh
          
        fi
    done

  ;;
  INGIROHFILT )
    echo "Calculate ROH for INGI populations from splitted vcf files with IBDseq (filtered by HWE)!!"
    # echo "We'll have also data separate for villages in FVG"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}\nr2=${r2}\nMAC=${mac}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    # pops_updated="FVG VBI CARL Erto Resia Illegio Sauris"
    pops_updated="FVG VBI CARL"
    for pop in $pops_updated
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL/ALL_FVG_20151113.vcf.gz.hwe_filt.vcf.gz
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/ALL_VBI_20151113.vcf.gz.hwe_filt.vcf.gz
            ;;
        CARL )
          pop_path=/lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL/ALL_CARL_20151113.vcf.gz.hwe_filt.vcf.gz
            ;;
      esac

        #we can filter out unrelated after the calculation
        # pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/BEAGLE/all_pop_but_${pop}.txt
        # pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/FIVE_POPS/all_pop_but_${pop}.txt
        # pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/all_pop_but_${pop}.txt
        # pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/listpop/all_pop_but_${pop}.txt
        # marker_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140520_ROH/sites_with_missing_genotypes.list
        #use freq data
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gtgl=${pop_path}/22.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} excludemarkers=${marker_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/22.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.nonmissing.maf_gt_05.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
        # bsub -J"roh_${pop}_${CHR}" -o"%J_roh_${pop}_${CHR}.o" -q normal -M8000 -n4 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz ibd=true nthreads=4 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
        # bsub -J"roh_${pop}_${CHR}" -o"%J_roh_${pop}_${CHR}.o" -q normal -M8000 -n4 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/IBDseq/ibdseq.r1206.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz nthreads=4 excludesamples=${pop_list} out=${outdir}/${pop}.roh
        if [[ $mac -eq 1 ]]; then
          #statements
        bsub -J"roh_${pop}_${CHR}" -o"LOGS/%J_roh_${pop}_${CHR}.o" -q normal -M5000 -n16 -R"span[hosts=1] select[mem>=5000] rusage[mem=5000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/IBDseq/ibdseq.r1206.jar gt=${pop_path} chrom=${CHR} nthreads=16 ibdlod=1 r2max=${r2} out=${outdir}/${pop}.roh
        else
        bsub -J"roh_${pop}_${CHR}" -o"LOGS/%J_roh_${pop}_${CHR}.o" -q normal -M5000 -n16 -R"span[hosts=1] select[mem>=5000] rusage[mem=5000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/IBDseq/ibdseq.r1206.jar gt=${pop_path} chrom=${CHR} nthreads=16 ibdlod=1 r2max=${r2} minalleles=${mac} out=${outdir}/${pop}.roh
          
        fi
    done
    ;;
  SPLITCSQ )
    echo "Split a chromosome file using a list of consequences"
    if [ $fixed == "fixed" ]
    then
      outdir=${outdir}_fixed
    else
      outdir=${outdir}_no_fixed
    fi

    for pop in $pops
    do
      if [ $fixed == "fixed" ]
      then
        case $pop in
          FVG )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}.chr${CHR}.tab.gz
            ;;
          VBI )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}.chr${CHR}.tab.gz
              ;;
          TSI )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}.chr${CHR}.tab.gz
              ;;
          CEU )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}.chr${CHR}.tab.gz
              ;;
        esac
      else
        case $pop in
          FVG )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}.chr${CHR}.not_fixed.tab.gz
            ;;
          VBI )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}.chr${CHR}.not_fixed.tab.gz
            ;;
          TSI )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}.chr${CHR}.not_fixed.tab.gz
            ;;
          CEU )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}.chr${CHR}.not_fixed.tab.gz
            ;;
        esac
      fi
      
      mkdir -p ${outdir}/${pop}
      # echo "(zcat ${pop_path}| head -1;zfgrep ${category} ${pop_path} | cut -f 1-6,8-) > ${outdir}/${pop}.${category}." | bsub -J"split_${pop}" -o"%J_slit_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]"
      (zcat ${pop_path}| head -1;zfgrep ${category} ${pop_path} )| cut -f 1-6,8- | gzip -c > ${outdir}/${pop}/${pop}.${category}.${CHR}.tab.gz
    
    done

  ;;
  SPLITCSQPRIV )
    echo "Split a chromosome file using a list of consequences"
    if [ $fixed == "fixed" ]
    then
      outdir=${outdir}_private_fixed
    else
      outdir=${outdir}_private_no_fixed
    fi
    
    for pop in $pops
    do
      if [ $fixed == "fixed" ]
      then
        case $pop in
          FVG )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}_private_chr${CHR}.merged_daf.fixed.tab.gz
            ;;
          VBI )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}_private_chr${CHR}.merged_daf.fixed.tab.gz
            ;;
        esac
      else
        case $pop in
          FVG )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}_private_chr${CHR}.tab.gz
            ;;
          VBI )
            pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/${pop}_private_chr${CHR}.tab.gz
            ;;
        esac
      fi
      
      mkdir -p ${outdir}/${pop}
      echo ${pop_path}
      # echo "(zcat ${pop_path}| head -1;zfgrep ${category} ${pop_path} | cut -f 1-6,8-) > ${outdir}/${pop}.${category}." | bsub -J"split_${pop}" -o"%J_slit_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]"
      (zcat ${pop_path} | head -1;zfgrep ${category} ${pop_path} )| cut -f 1-6,8- | gzip -c > ${outdir}/${pop}/${pop}.private.${category}.${CHR}.tab.gz
    
    done
  ;;

  EXTRFREQ )
    #extract data from a vcf given a MAF threshold, than extract the same variants (if present) in other populations
    #calculate MAF for that population
    # 1) For each population or sub population, extract variants with maf > 1% 
    #   - calculate AF and MAF
    #   - keep ref and ALT alleles
    # 2) Extract those variants for all the others populations
    #   - calculating AF, MAF values
    #   - check alleles concordance
    # 3) Keep track of variants missing from each comparison

    echo "Extract from tables for each population, info about shared sites"
    pops_updated="FVG VBI CARL "
    for pop in $pops_updated
    do
      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL/
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        CARL )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
      esac
      #those are the files I need to use to extract the info splitted by chr
      # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
      sel_var=${outdir}/${pop}/${CHR}_sel_${variant}_${pop}.tab

      #extract the list in bed format only for BIALLELIC sites
      bcftools view -i"MAF<=0.01" -m2 -M2 -v ${variant} ${pop_path}/${CHR}.vcf.gz | bcftools query -f"%CHROM\t%POS\n" -o ${sel_snps}
      
      #extract the list with all informations needed
      echo"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN";bcftools view -i"MAF<=0.01" -m2 -M2 -v ${variant} ${pop_path}/${CHR}.vcf.gz | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n" | awk 'BEGIN{OFS="\t"}{print $0,$7-$6,$6/$7,1-($6/$7)}'


      # #create the list of variants we need to extract
      # awk '{print $1"O"$2"O"$3}' ${private} | tr " " "\t" | dos2unix > ${private}.list
      # awk '{print $1"O"$2"O"$3}' ${private_fixed} | tr " " "\t" | dos2unix > ${private_fixed}.list
      # awk '{print $1"O"$2"O"$3}' ${shared} | tr " " "\t" | dos2unix > ${shared}.list
      # awk '{print $1"O"$2"O"$3}' ${shared_fixed} | tr " " "\t" | dos2unix > ${shared_fixed}.list

      #create temporary files to do the extraction
      in_file=${pop_path}/CHR${CHR}/INGI_chr${CHR}.merged_maf.tab.gz
      # zcat ${in_file} | awk '{print $1"O"$2"O"$3,$0}' | tr " " "\t" | dos2unix > ${in_file}.tmp

      zcat ${in_file} |awk '{print $1"O"$2"O"$3"O"$4,$0}'|tr " " "\t" | fgrep -w -f <(awk -v cro=${CHR} '{if($1==cro) print $1"O"$2"O"$3"O"$4}' ${private})|cut -f 2- |gzip -c > ${in_file}.${pop}.private.tab.gz
      zcat ${in_file} |awk '{print $1"O"$2"O"$3"O"$4,$0}'|tr " " "\t" | fgrep -w -f <(awk -v cro=${CHR} '{if($1==cro) print $1"O"$2"O"$3"O"$4}' ${private_fixed})|cut -f 2- |gzip -c > ${in_file}.${pop}.private_fixed.tab.gz
      zcat ${in_file} |awk '{print $1"O"$2"O"$3"O"$4,$0}'|tr " " "\t" | fgrep -w -f <(awk -v cro=${CHR} '{if($1==cro) print $1"O"$2"O"$3"O"$4}' ${shared})|cut -f 2- |gzip -c > ${in_file}.${pop}.shared.tab.gz
      zcat ${in_file} |awk '{print $1"O"$2"O"$3"O"$4,$0}'|tr " " "\t" | fgrep -w -f <(awk -v cro=${CHR} '{if($1==cro) print $1"O"$2"O"$3"O"$4}' ${shared_fixed})|cut -f 2- |gzip -c > ${in_file}.${pop}.shared_fixed.tab.gz
      
      # #now grep the file to extract the data we need, using different lists
      # (fgrep -w -f <(grep "^${CHR}" ${private}.list) ${in_file}.tmp)| cut -f 2- | dos2unix |gzip -c > ${in_file}.private.tab.gz
      # (fgrep -w -f <(grep "^${CHR}" ${private_fixed}.list) ${in_file}.tmp)| cut -f 2- | dos2unix |gzip -c > ${in_file}.private_fixed.tab.gz
      # (fgrep -w -f <(grep "^${CHR}" ${shared}.list) ${in_file}.tmp)| cut -f 2- | dos2unix |gzip -c > ${in_file}.shared.tab.gz
      # (fgrep -w -f <(grep "^${CHR}" ${shared_fixed}.list) ${in_file}.tmp)| cut -f 2- | dos2unix |gzip -c > ${in_file}.shared_fixed.tab.gz

      # rm ${in_file}.tmp

    done
  ;;
  * )
    echo -e "USAGE:\n
    IBD )\n
      - Extract IBD information for each population. Here we are using plink2 (1.9)\n
    HET )\n
      - Extract Inbreeding coeff information for each population. Here we are using plink2 (1.9)\n
    ROH )\n
      - Calculate ROH....with BEAGLE and separate population files (no maf filtering)\n
    ROH2 )\n
      - Calculate ROH from a unified vcf file....with BEAGLE...we need a file without missing genotypes(NO MAF filter)!!\n
    ROH3 )\n
      - Calculate ROH from a unified vcf file....with BEAGLE...we need a file without missing genotypes(currently is filtered on MAF>5%)!!\n
    ROH4 )\n
      - Calculate ROH from different files for each population....using BEAGLE4!(filtered by MAF!!))\n
    ROH8 )\n
      - Calculate ROH from different files for each population....using IBDseq!(NO MAF filter!!))\n
    
    "
  ;;
  esac