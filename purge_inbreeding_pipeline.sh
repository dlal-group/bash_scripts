#!/usr/local/bin/bash
#Script to create a pipeline for the purge through inbreeding work
#
#
#define the population we are using
pops="FVG VBI TSI CEU CARL"
# pops="FVG VBI"

#retrieve the MODE parameter to select the correct operation
MODE=$2
CHR=$3
# define the output folder relative to the chromosome, if specified
outdir=CHR${CHR}
mkdir -p LOGS
mkdir -p ${outdir}

case $MODE in
  HOMDAC)
  sample=$1
  list_path=$4
  cat=$5
  pop=$6
  rm -r ${outdir}
  ;;
  ROH*)
    #set parameters for beagle:
    window=$3
    overlap=$4
    ;;
  SPLITCSQ*)
    #set parameters for splitting in consequences
    category=$3
    fixed=$4
    ;;
  SHARED )
    #set parameters for file input/output
    in_dir=$3
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
      ;;
  IBDCLUST )
      #set parameters for file input/output
      win_size=$3
      density=$4
      in_dir=$5
      ;;
  MANCLUSTIBD )
      in_dir=$3
      ;;
  MAFSPEC )
      #set parameters for file input/output
      input_file=$3
      maf_file_path=$4
    ;;
  GERMLINE )
    # set up args for germline command line
    MATCH=$3
    HOM=$4
    HET=$5
    BITS=$6
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
  ;;
  EXTRMAF )
    pop_path=$3
    # in_path=$4
  ;;
  RANDLIST )
  ROUNDS=$2
  SIZE=$3
  ;;
esac

# Merge different popuplation together
# TODO: add code!!!



#create a case statement to select the operation to do:
case $MODE in
  HOMDAC)
  #extract counts for alternative alleles count from vcf files fo a subset of variants on a subset of samples
  #out path : /lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/RESULTS/HOMCOUNT
  # pops_updated="FVG VBI CARL TSI CEU Erto Resia Illegio Sauris"
  # for pop in $pops_updated
  #   do
      case $pop in
        FVG )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/FVG_unrelated.list
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/FVG_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/FVG_unrelated.list_1
          ;;
        VBI )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/VBI_unrelated.list
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/VBI_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/VBI_unrelated.list_1
          ;;
        CARL )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/CARL_unrelated.list
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/CARL_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/CARL_unrelated.list_1
          ;;
        TSI )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/TSI.list
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/TSI_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/TSI_unrelated.list_1
          ;;
        CEU )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/CEU.list
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/CEU_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/CEU_unrelated.list_1
          ;;
        Erto )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/Erto_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/Erto_unrelated.list_1
            ;;
        Sauris )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/Sauris_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/Sauris_unrelated.list_1
            ;;
        Illegio )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/Illegio_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/Illegio_unrelated.list_1
            ;;
        Resia )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/Resia_unrelated.list
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/Resia_unrelated.list_1
            ;;
      esac

    # snplist=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza/listsites/syn/syn.17.bed
    snplist=${list_path}/${cat}/${cat}.${CHR}.bed
    mkdir -p /lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/RESULTS/HOMCOUNT/shared/${cat}
    mkdir -p /lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/RESULTS/HOMCOUNT/shared/${cat}/${pop}

    vcf=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POP_MERGED_FILES/FIVE_POPS/20140730_ANNOTATED/${CHR}.clean_annotated.vcf.gz
    out_name=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/RESULTS/HOMCOUNT/shared/${cat}/${pop}_${cat}_${sample}
    shared_bed=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/INPUT_FILES/FIVE_POPS/WG/sharedsites/${pop}_shared_chr${CHR}.bed
    # shared_bed=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/INPUT_FILES/FIVE_POPS/WG/sharedsites/Illegio_shared_chr17.bed
    awk 'FNR==NR { a[$2]=$0; next } $2 in a { print a[$2] }' ${shared_bed} ${snplist} | sort -g -k2,2 |uniq| tr " " "\t" > /lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/RESULTS/HOMCOUNT/shared/${cat}/shared.${pop}.${cat}.${chr}.bed
    # the second file is the one which gives you the items's order
    shared_cat=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/RESULTS/HOMCOUNT/shared/${cat}/shared.${pop}.${cat}.${chr}.bed
    # --derived
    # For use with the previous four frequency and count options only. Re-orders the output file columns so that the ancestral allele appears first.
    # This option relies on the ancestral allele being specified in the VCF file using the AA tag in the INFO field
    # while read line
    # do
      /nfs/team151/software/vcftools/bin/vcftools --gzvcf ${vcf} --bed ${shared_cat} --indv ${sample} --counts --derived --out ${out_name} # --> conte per locus per individuo  
      
      # sample_hom=`tail -n+2 ${out_name}.frq.count | awk '{split($5,ref,":");split($6,alt,":")}{if(ref[2]==2 || alt[2]==2) print ref[2],alt[2]}' | wc -l `
      sample_hom=`tail -n+2 ${out_name}.frq.count | awk '{split($5,aa,":")}{if(aa[2]==2) print aa[2]}' | wc -l `
      tot_snp=`wc -l ${shared_bed}|cut -f 1 -d " "`
      echo "${sample} ${CHR} ${sample_hom} ${pop} ${tot_snp}" > /lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/RESULTS/HOMCOUNT/shared/${cat}/${pop}/${pop}_${cat}_${CHR}_${sample}.tab
    # done < <(cat ${pop_path})
  # done
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
  CSQEXTRACT )
  #extract number of variants per sample belonging to a particular category (syn/missense eventually already stratified by novel/shared)
  /lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza/listsites/miss
  ;;
  ROHENR)
    # check if ROH are enriched in deleterious mutations
    #first select
  ;;
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
      bsub -J"ibd_${pop}" -o"%J_ibd_${pop}.o" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --vcf ${pop_path}/22.vcf.gz --double-id --biallelic-only --genome gz --parallel 1 2 --out ${outdir}/ibd_${pop}
    done
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
        CARL )
          pop_path=/lustre/scratch113/projects/fvg_seq/20140410/INGI/CARL/PLINK
            ;;
      esac
        #calculate frequencies, before:
        bsub -J"freq_${pop}" -o"%J_freq_${pop}.o" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --bfile ${pop_path}/22.${pop,,} --freq --nonfounders --out ${outdir}/freq_${pop}
        
        #use freq data
        bsub -J"inb_${pop}" -o"%J_inb_${pop}.o" -w "ended(freq_${pop})" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --bfile ${pop_path}/22.${pop,,} --het --out ${outdir}/inb_${pop}
        bsub -J"ibc_${pop}" -o"%J_ibc_${pop}.o" -w "ended(freq_${pop})" -q yesterday -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --bfile ${pop_path}/22.${pop,,} --ibc --out ${outdir}/ibc_${pop}
      done
  ;;
  SHARED )
    #calculate shared and private sites for all populations
    pops_updated="FVG VBI CARL TSI CEU Erto Resia Illegio Sauris"
    for pop in $pops_updated="FVG VBI CARL TSI CEU Erto Resia Illegio Sauris"
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

    pops_updated="FVG VBI CARL TSI CEU Erto Resia Illegio Sauris"
    for pop in $pops_updated
    do
      case $pop in
        TSI )
          #Private calculation without fixed sites
          # zcat ${in_file} | awk '$6 !="NA" && $6 > 0 && $6 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # # Fixed in Isolate and private
          # zcat ${in_file} | awk '$6 !="NA" && $6 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # # Shared calculation without fixed sites
          zcat ${in_file}|cut -f 1,3- | awk '$5 !="NA" && $5 > 0 && $5 < 1 && ($6 > 0 || $7 > 0 || $8 > 0 || $9 > 0 || $10 > 0 || $11 > 0 || $12 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # # Fixed in Isolate and shared
          # zcat ${in_file} | awk '$6 !="NA" && $6 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # # NOVEL in my population
          zcat ${in_file} |cut -f 1,3-| awk '$5 !="NA" && ($6 == "NA" || $7 == "NA" || $8 == "NA" || $9 == "NA" || $10 == "NA" || $11 == "NA" || $12 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        CEU )
          #Private calculation without fixed sites
          # zcat ${in_file} | awk '$4 !="NA" && $4 > 0 && $6 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # # Fixed in Isolate and private
          # zcat ${in_file} | awk '$4 !="NA" && $4 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # # Shared calculation without fixed sites
          zcat ${in_file}|cut -f 1,3- | awk '$4 !="NA" && $4 > 0 && $4 < 1 && ($6 > 0 || $7 > 0 || $8 > 0 || $9 > 0 || $10 > 0 || $11 > 0 || $12 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # # Fixed in Isolate and shared
          # zcat ${in_file} | awk '$4 !="NA" && $4 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in my population
          zcat ${in_file} |cut -f 1,3-| awk '$4 !="NA" && ($6 == "NA" || $7 == "NA" || $8 == "NA" || $9 == "NA" || $10 == "NA" || $11 == "NA" || $12 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        VBI )
          #Private calculation without fixed sites
          zcat ${in_file}|cut -f 1,3- | awk '$6 !="NA" && $6 > 0 && $6 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file}|cut -f 1,3- | awk '$6 !="NA" && $6 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file}|cut -f 1,3- | awk '$6 !="NA" && $6 > 0 && $6 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file}|cut -f 1,3- | awk '$6 !="NA" && $6 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file}|cut -f 1,3- | awk '$6 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          ;;
        FVG )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          #Private calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$7 !="NA" && $7 > 0 && $7 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file} |cut -f 1,3- | awk '$7 !="NA" && $7 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$7 !="NA" && $7 > 0 && $7 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file} |cut -f 1,3- | awk '$7 !="NA" && $7 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file} |cut -f 1,3- | awk '$7 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          ;;
    
        CARL )
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
          #Private calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$8 !="NA" && $8 > 0 && $8 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file} |cut -f 1,3- | awk '$8 !="NA" && $8 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$8 !="NA" && $8 > 0 && $8 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file} |cut -f 1,3- | awk '$8 !="NA" && $8 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file} |cut -f 1,3- | awk '$8 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
    
          ;;
        Erto )
          #Private calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$9 !="NA" && $9 > 0 && $9 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file} |cut -f 1,3- | awk '$9 !="NA" && $9 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$9 !="NA" && $9 > 0 && $9 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file} |cut -f 1,3- | awk '$9 !="NA" && $9 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file} |cut -f 1,3- | awk '$9 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          ;;
        Illegio )
         #Private calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$10 !="NA" && $10 > 0 && $10 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file} |cut -f 1,3- | awk '$10 !="NA" && $10 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$10 !="NA" && $10 > 0 && $10 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file} |cut -f 1,3- | awk '$10 !="NA" && $10 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file} |cut -f 1,3- | awk '$10 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz 
          ;;
        Resia )
          #Private calculation without fixed sites
          #Private calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$11 !="NA" && $11 > 0 && $11 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file} |cut -f 1,3- | awk '$11 !="NA" && $11 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$11 !="NA" && $11 > 0 && $11 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file} |cut -f 1,3- | awk '$11 !="NA" && $11 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file} |cut -f 1,3- | awk '$11 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz
          ;;
        Sauris )
         #Private calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$12 !="NA" && $12 > 0 && $12 < 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and private
          zcat ${in_file} |cut -f 1,3- | awk '$12 !="NA" && $12 == 1 && $5 == 0 && $4 == 0'| gzip -c > ${outdir}/${pop}_private_chr${CHR}.merged_frq.fixed.tab.gz
          # Shared calculation without fixed sites
          zcat ${in_file} |cut -f 1,3- | awk '$12 !="NA" && $12 > 0 && $12 < 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.tab.gz
          # Fixed in Isolate and shared
          zcat ${in_file} |cut -f 1,3- | awk '$12 !="NA" && $12 == 1 && ($5 > 0 || $4 > 0)'| gzip -c > ${outdir}/${pop}_shared_chr${CHR}.merged_frq.fixed.tab.gz
          # NOVEL in Isolate
          zcat ${in_file} |cut -f 1,3- | awk '$12 !="NA" && ($5 == "NA" && $4 == "NA")'| gzip -c > ${outdir}/${pop}_novel_chr${CHR}.merged_frq.tab.gz 
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
      # case $pop in
      #   FVG )
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
      #     ;;
      #   VBI )
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
      #     ;;
      #   CARL )
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
      #     ;;
      #   TSI )
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
      #     ;;
      #   CEU )
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF/FIVE_POP
      #     ;;
      # esac
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
    pops_updated="FVG VBI CARL TSI CEU Erto Resia Illegio Sauris"
    # for pop in $pops
    # for pop in $pops_updated
    # do
      # case $pop in
      #   FVG )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/FVG_unrelated.list
      #     ;;
      #   VBI )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/VBI_unrelated.list
      #     ;;
      #   CARL )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/CARL_unrelated.list
      #     ;;
      #   TSI )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/TSI.list
      #     ;;
      #   CEU )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/CEU.list
      #     ;;
      #   Erto )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Erto_unrelated.list
      #       ;;
      #   Sauris )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Sauris_unrelated.list
      #       ;;
      #   Illegio )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Illegio_unrelated.list
      #       ;;
      #   Resia )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Resia_unrelated.list
      #       ;;
      # esac
      # maf_file=${maf_file_path}/INGI_chr${CHR}.merged_maf.tab.gz
      # maf_file=${maf_file_path}/INGI_villages_chr${CHR}.merged_daf.tab.gz
      maf_file=${maf_file_path}/INGI_villages_chr${CHR}.merged_maf.tab.gz
      out_tmp=`basename ${input_file}`

      (zcat ${maf_file} | cut -f 1,3-|head -1;zcat ${maf_file} | cut -f 1,3- | fgrep -w -f <(awk '{print $3}' ${input_file}) )| tr " " "\t" | gzip -c > ${out_tmp}.maf_file.gz

    # done
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
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gl=${pop_path}/${CHR}.vcf.gz ibd=true nthreads=2 window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
      
    done

  ;;
  ROH2 )
    echo "Calculate ROH from a unified vcf file....with BEAGLE...we need a file without missing genotypes(NO MAF filter)!!"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    pops_updated="FVG VBI TSI CEU CARL Erto Resia Illegio Sauris"
    for pop in $pops_updated
    do

      # case $pop in
      #   FVG )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POP_MERGED_FILES/FIVE_POPS/20140818_NONMISSING
      #     ;;
      #   VBI )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
      #       ;;
      #   TSI )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
      #       ;;
      #   CEU )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
      #       ;;
      #   CARL )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES
      #     # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
      #       ;;
      #   Erto )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
      #       ;;
      #   Sauris )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
      #       ;;
      #   Illegio )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
      #       ;;
      #   Resia )
      #     pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
      #       ;;
      # esac
        pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/POP_MERGED_FILES/FIVE_POPS/20140818_NONMISSING
        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/listpop/IBD_lists/all_pop_but_${pop}.txt
        # commented to use the ALL population files
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
        bsub -J"roh_${pop}_${CHR}" -o"%J_roh_${pop}_${CHR}.o" -q normal -M8000 -n4 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz ibd=true nthreads=4 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
      
    done

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
  ROHVILLAGE2 )
    echo "Calculate ROH from a unified vcf file....with BEAGLE...we need a file without missing genotypes(NO MAF filter)!!"
    echo "We'll have also data separate for villages in FVG"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    pops_updated="FVG VBI TSI CEU CARL Erto Resia Illegio Sauris"
    for pop in $pops_updated
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        CARL )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        Erto )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        Sauris )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        Illegio )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        Resia )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
      esac
        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/BEAGLE/all_pop_but_${pop}.txt
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
        bsub -J"roh_${pop}_${CHR}" -o"%J_roh_${pop}_${CHR}.o" -q normal -M8000 -n4 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz ibd=true nthreads=4 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
      
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
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gtgl=${pop_path}/22.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} excludemarkers=${marker_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/22.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} out=${pop}.roh
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.nonmissing.maf_gt_05.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
      
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
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.vcf.gz ibd=true nthreads=2 window=${window} overlap=${overlap} out=${outdir}/${pop}.roh
    
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
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --vcf ${pop_path}/${CHR}.nonmissing.vcf.gz --keep ${pop_list} --homozyg --homozyg-window-snp ${window} --out ${outdir}/${pop}.roh
    
    done

  ;;
  ROH6 )
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
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --vcf ${pop_path}/${CHR}.nonmissing.maf_gt_05.vcf.gz --keep ${pop_list} --homozyg --homozyg-window-snp ${window} --out ${outdir}/${pop}.roh
      
    done

  ;;
  ROH7 )
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
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- plink2 --vcf ${pop_path}/${CHR}.vcf.gz --homozyg --homozyg-window-snp ${window} --out ${outdir}/${pop}.roh
      
    done

  ;;
  ROH8 )
    echo "Calculate ROH from a unified vcf file....with IBDseq...we need a file without missing genotypes(NO MAF filter)!!"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    for pop in $pops
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES
          # pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140419_ANNOTATED
            ;;
      esac
        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/all_pop_but_${pop}.txt
        # pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/ALL/listpop/all_pop_but_${pop}.txt
        # marker_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/POP_MERGED_FILES/20140520_ROH/sites_with_missing_genotypes.list
        #use freq data
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gtgl=${pop_path}/22.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} excludemarkers=${marker_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/22.nonmissing.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} out=${pop}.roh
        # bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/beagle_4/b4.r1274.jar gt=${pop_path}/${CHR}.nonmissing.maf_gt_05.vcf.gz ibd=true nthreads=2 excludesamples=${pop_list} window=${window} overlap=${overlap} out=${pop}.roh
        bsub -J"roh_${pop}" -o"%J_roh_${pop}.o" -q normal -M8000 -n2 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/IBDseq/ibdseq.r1206.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz nthreads=2 excludesamples=${pop_list} out=${outdir}/${pop}.roh
      
    done

  ;;
  ROHVILLAGE8 )
    echo "Calculate ROH from a unified vcf file....with IBDseq...we need a file without missing genotypes(NO MAF filter)!!"
    echo "We'll have also data separate for villages in FVG"
    echo -e "Parameters: \nwindow=${window}\noverlap=${overlap}"
    #use he same vcf file for all the samples but change the sample list of individuals toi exclude from the analysis
    pops_updated="FVG VBI TSI CEU CARL Erto Resia Illegio Sauris"
    for pop in $pops_updated
    do

      case $pop in
        FVG )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
          ;;
        VBI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        TSI )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        CEU )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        CARL )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        Erto )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        Sauris )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        Illegio )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
        Resia )
          pop_path=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/POP_MERGED_FILES/FIVE_POPS/20140801_NONMISSING
            ;;
      esac

        pop_list=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/listpop/BEAGLE/all_pop_but_${pop}.txt
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
        bsub -J"roh_${pop}_${CHR}" -o"%J_roh_${pop}_${CHR}.o" -q normal -M8000 -n4 -R"span[hosts=1] select[mem>=8000] rusage[mem=8000]" -- java -Xms5000m -Xmx5000m -jar /nfs/team151/software/IBDseq/ibdseq.r1206.jar gt=${pop_path}/${CHR}.non_missing.vcf.gz nthreads=4 excludesamples=${pop_list} out=${outdir}/${pop}.roh
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
    ROH5 )\n
      - Calculate ROH from a unified vcf file....using plink! (NO MAF FILTERING)\n
    ROH6 )\n
      - Calculate ROH from a unified vcf file....using plink! (with MAF FILTERING)\n
    ROH7 )\n
      - Calculate ROH from different files for each population....using PLINK!\n
    "
  ;;
  esac