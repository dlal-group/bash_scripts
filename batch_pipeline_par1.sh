#!/usr/local/bin/bash

#Pipeline for batch effect analysis check
#What we would do (are doing):
###############################################################
## 1. find the extreme site
# CC analyses for BGI vs SANGER seq center for TWINS and ALSPAC and for the JOINT set:
#
#Two different maf tresholds:
#1) ALL sites
#2) MAF >= 0.01

#TFor the whole genome (unpruned) set and for the LD pruned set
# cohort="ALSPAC"
# cohort="TWINSUK"
# for cohort in TWINSUK ALSPAC JOINT
for cohort in JOINT
do
    if [ $cohort == "TWINSUK" ]
    then
        cohort_path="twins"
    fi
    if [ $cohort == "ALSPAC" ]
    then
        cohort_path="alspac"
    fi
    if [ $cohort == "JOINT" ]
    then
        cohort_path="joint_NEW"
    fi


    outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic_all
    mkdir -p ${outdir}

    for chr in {1..22}
    do
        if [ $cohort == "JOINT" ]
        then
            # bsub -J "plink_run_2" -o "%J_plink_run_2.o" -P "analysis-rd" -M 7000000 -R"select[mem>7000] rusage[mem=7000]" -q normal -- plink --noweb \
            bsub -J "plink_run_${chr}" -o "${outdir}/%J_plink_run_${chr}.o" -P "analysis-rd" -M 6000000 -R"select[mem>6000] rusage[mem=6000]" -q normal -- plink --noweb \
            --bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
            --chr ${chr} \
            --allow-no-sex \
            --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
            --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.centre \
            --covar /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.pop \
            --logistic \
            --adjust \
            --out ${outdir}/chr${chr}.txt

        else

            bsub -J "plink_run_${chr}" -o "${outdir}/%J_plink_run_${chr}.o" -P "analysis-rd" -M 6000000 -R"select[mem>6000] rusage[mem=6000]" -q normal -- plink --noweb \
            --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed \
            --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY \
            --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam \
            --allow-no-sex \
            --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
            --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.centre \
            --logistic \
            --adjust \
            --out ${outdir}/chr${chr}.txt
        fi
    done
done


###############################################################
##2. remove the extreme site, then prune the data
#format results and extract data
#Extract sites with different threshold for pval:
# 1- < 5e-8
# 2- < 5e-6
# 3- < 1e-5
# 4- < 1e-2
# for cohort in ALSPAC TWINSUK JOINT
for cohort in JOINT
do
    if [ $cohort == "TWINSUK" ]
    then
        cohort_path="twins"
    fi
    if [ $cohort == "ALSPAC" ]
    then
        cohort_path="alspac"
    fi
    if [ $cohort == "JOINT" ]
    then
        cohort_path="joint_NEW"
    fi

    indir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/logistic_all
    outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/logistic_all/REMOVELIST

    mkdir -p ${outdir}

    for thr in  5e-8 5e-6 1e-5 1e-2
    do
        #for chr in {1..22}; do
	for chr in X; do
            if [ $cohort == "JOINT" ]
            then
                awk -v pthr=${thr} '{if($9<= pthr && $5 == "ADD") print $2}' ${indir}/chr${chr}.txt.assoc.logistic  > ${outdir}/chr${chr}.txt.assoc.logistic.${cohort}.${thr}.removelist
            else
                awk -v pthr=${thr} '{if($9<= pthr) print $2}' ${indir}/chr${chr}.txt.assoc.logistic  > ${outdir}/chr${chr}.txt.assoc.logistic.${cohort}.${thr}.removelist
            fi
        done
        
#        for chr in {1..22}; do
 #           cat ${outdir}/chr${chr}.txt.assoc.logistic.${cohort}.${thr}.removelist
  #      done > ${outdir}/All.txt.assoc.logistic.${cohort}.${thr}.removelist
    done
done



############################################################### 
##3. calculate IBS metrics, done separately, in parallel
#first calculate freq:
#By CHR
for cohort in ALSPAC TWINSUK JOINT
do
    if [ $cohort == "TWINSUK" ]
    then
        cohort_path="twins"
    fi
    if [ $cohort == "ALSPAC" ]
    then
        cohort_path="alspac"
    fi
    if [ $cohort == "JOINT" ]
    then
        cohort_path="joint_NEW"
    fi

    # indir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/logistic
    outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/MDS

    mkdir -p ${outdir}
    for chr in {1..22}; do
        bsub -J "plink_freq_${chr}" -o "${outdir}/%J_plink_freq_${chr}.o" -P "analysis-rd" -M 6000000 -R"select[mem>6000] rusage[mem=6000]" -q normal -- plink --noweb \
        --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed \
        --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY \
        --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam \
        --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
        --freq \
        --out ${outdir}/chr${chr}.txt
    done
done

#genome wide
for cohort in ALSPAC TWINSUK JOINT
do
    if [ $cohort == "TWINSUK" ]
    then
        cohort_path="twins"
    fi
    if [ $cohort == "ALSPAC" ]
    then
        cohort_path="alspac"
    fi
    if [ $cohort == "JOINT" ]
    then
        cohort_path="joint_NEW"
    fi

    outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input
    mkdir -p ${outdir}

    bsub -J "plink_freq_${cohort}" -o "${outdir}/%J_plink_freq_${cohort}.o" -P "analysis-rd" -M 80000000 -R"select[mem>80000] rusage[mem=80000]" -q normal -- plink --noweb \
    --bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged_unpruned \
    --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
    --freq \
    --out ${outdir}/${cohort}_all_chr_merged_unpruned

done

##################ALSPAC/TWINSUK by chr
# for chr in 4 6 8 9 12 20
# for chr in 4
for cohort in ALSPAC TWINSUK JOINT
for cohort in ALSPAC
do
    if [ $cohort == "TWINSUK" ]
    then
        cohort_path="twins"
    fi
    if [ $cohort == "ALSPAC" ]
    then
        cohort_path="alspac"
    fi
    if [ $cohort == "JOINT" ]
    then
        cohort_path="joint_NEW"
    fi


    for chr in {1..22}
    do
        exc_thr="1e-2"
        filter="FILTERED_1e2_1"
        indir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/MDS
        outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/MDS/${filter}/CHR${chr}
        mkdir -p ${outdir}
        a=$[`ls ${indir}/tmp.list*| wc -l`-1]
        let i=0;let j=0;

        while [ $i -le $a ]
        do
            while [ $j -le $a ]
            do

            # bsub -J "plink_run_${i}_${j}_${chr}" -o "%J_plink_run_${i}_${j}_${chr}.o" -M 4000 -R"select[mem>4000] rusage[mem=4000]" -q yesterday -- plink --noweb \
            bsub -J "${outdir}/plink_run_${i}_${j}_${chr}_${exc_thr}" -o "${outdir}/%J_plink_run_${i}_${j}_${chr}_${exc_thr}.o" -P "analysis-rd" -M 7000000 -R"select[mem>7000] rusage[mem=7000]" -q normal -- plink --noweb \
            --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed \
            --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY \
            --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam \
            --exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/logistic/REMOVELIST/All.txt.assoc.logistic.${cohort}.${exc_thr}.removelist \
            --maf 0.01 \
            --allow-no-sex \
            --read-freq ${indir}/chr${chr}.freq.frq \
            --genome \
            --genome-lists ${indir}/tmp.list`printf "%03i\n" $i` \
            ${indir}/tmp.list`printf "%03i\n" $j` \
            --out ${outdir}/data.sub.$i.$j

            let j=$j+1
            done
        let i=$i+1
        let j=$i
        done
    done
done

###############################################################
############ ALSPAC / TWINSUK GWA ###
# for cohort in ALSPAC TWINSUK JOINT
for cohort in JOINT
do
    if [ $cohort == "TWINSUK" ]
    then
        cohort_path="twins"
    fi
    if [ $cohort == "ALSPAC" ]
    then
        cohort_path="alspac"
    fi
    if [ $cohort == "JOINT" ]
    then
        cohort_path="joint_NEW"
    fi

    # exc_thr="1e-2"
    # exc_thr="1e-5"
    # exc_thr="5e-6"
    exc_thr="5e-8"

    # filter="UNFILTERED_1"
    # filter="FILTERED_1e2_1"
    # filter="FILTERED_1e5_1"
    # filter="FILTERED_5e6_1"
    # filter="FILTERED_5e8_1"
    # filter="FILTERED_KLAUDIA_1"
    filter="FILTERED_KLAUDIA"
    # filter="FILTERED_1e2"
    # filter="FILTERED_1e5"
    # filter="FILTERED_5e6"
    # filter="FILTERED_5e8"

    indir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/MDS

    outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/MDS/GW/${filter}

    maf=0.01 #only common variants
    #max_maf=0.05
    mkdir -p ${outdir}

    a=$[`ls ${indir}/tmp.list*| wc -l`-1]

    let i=0;let j=0;

    while [ $i -le $a ]
    do
        while [ $j -le $a ]
        do

            bsub -J "plink_run_${i}_${j}" -o "${outdir}/%J_plink_run_${i}_${j}.o" -P "analysis-rd" -M 16000000 -R"select[mem>16000] rusage[mem=16000]" -q normal -- plink --noweb \
            --bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged_unpruned_1 \
            --exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/logistic/REMOVELIST/All.txt.assoc.logistic.${cohort}.${exc_thr}.removelist \
            --maf ${maf} \
            --allow-no-sex \
            --read-freq /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/${cohort}_all_chr_merged_unpruned_1.frq \
            --genome \
            --genome-lists ${indir}/tmp.list`printf "%03i\n" $i` \
            ${indir}/tmp.list`printf "%03i\n" $j` \
            --out ${outdir}/data.sub.$i.$j

            let j=$j+1
        done

        let i=$i+1
        let j=$i

    done
done

cat data.sub*genome > data.genome
rm data.sub.*

##4. calculate mds for 4 groups (sequencing center x cohorts)
#####################################################################################################################
############################################ MDS
########## TWINSUK / ALSPAC JOINT
# for cohort in ALSPAC TWINSUK
# for cohort in ALSPAC TWINSUK JOINT
for cohort in JOINT
do
    if [ $cohort == "TWINSUK" ]
    then
        cohort_path="twins"
    fi
    if [ $cohort == "ALSPAC" ]
    then
        cohort_path="alspac"
    fi
    if [ $cohort == "JOINT" ]
    then
        cohort_path="joint_NEW"
    fi

    exc_thr="1e-2"
    # exc_thr="1e-5"
    # exc_thr="5e-6"
    # exc_thr="5e-8"

    # filter="FILTERED_KLAUDIA_1"
    filter="FILTERED_KLAUDIA"
    # filter="FILTERED_1e2_WR"
    # filter="FILTERED_1e2_1"
    # filter="FILTERED_1e5_1"
    # filter="FILTERED_5e6_1"
    # filter="FILTERED_5e8_1"
    # filter="UNFILTERED_1"

    genoset="PRUNED"
    # genoset="UNPRUNED"

    indir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/${genoset}/MDS
    outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/${genoset}/MDS/GW/${filter}

    maf=0.01 #only common variants
    #max_maf=0.05
    mkdir -p ${outdir}

    bsub -J "plink_mds_all" -o "${outdir}/%J_plink_mds_all.o" -M8000000 -P "analysis-rd" -R"select[mem>=8000] rusage[mem=8000]" -q normal -- plink --noweb \
    --bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
    --allow-no-sex \
    --exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/All.chr.from_klaudia.snplist.removelist \
    --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
    --read-genome ${outdir}/data.genome \
    --cluster --mds-plot 10 \
    --out ${outdir}/all_chr_merged
done


    --maf ${maf} \
    --bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
    --exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/${genoset}/logistic/REMOVELIST/All.txt.assoc.logistic.${cohort}.${exc_thr}.removelist \
