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
# for cohort in JOINT
for cohort in TWINSUK ALSPAC JOINT
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
            --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed \
            --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY \
            --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam \
            --allow-no-sex \
            --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
            --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.centre \
            --covar /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples_pheno.pop \
            --logistic \
            --adjust \
            --out ${outdir}/chr${chr}.txt

        else

            bsub -J "plink_run_${chr}" -o "${outdir}/%J_plink_run_${chr}.o" -P "analysis-rd" -M 6000000 -R"select[mem>6000] rusage[mem=6000]" -q normal -- plink --noweb \
            --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bed \
            --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.bim.COPY \
            --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.fam \
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
# for cohort in JOINT
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

    indir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic
    outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic/REMOVELIST

    mkdir -p ${outdir}

    for thr in  5e-8 5e-6 1e-5 1e-2
    do
        for chr in {1..22}; do
            if [ $cohort == "JOINT" ]
            then
                awk -v pthr=${thr} '{if($9<= pthr && $5 == "ADD") print $2}' ${indir}/chr${chr}.txt.assoc.logistic  > ${outdir}/chr${chr}.txt.assoc.logistic.${cohort}.${thr}.removelist
            else
                awk -v pthr=${thr} '{if($9<= pthr) print $2}' ${indir}/chr${chr}.txt.assoc.logistic  > ${outdir}/chr${chr}.txt.assoc.logistic.${cohort}.${thr}.removelist
            fi
        done
        
        for chr in {1..22}; do
            cat ${outdir}/chr${chr}.txt.assoc.logistic.${cohort}.${thr}.removelist
        done > ${outdir}/All.txt.assoc.logistic.${cohort}.${thr}.removelist
    done
done



############################################################### 
##3. calculate IBS, done separately, in parallel
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
    outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS

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

    bsub -J "plink_freq_${chr}" -o "${outdir}/%J_plink_freq_${chr}.o" -P "analysis-rd" -M 8000000 -R"select[mem>8000] rusage[mem=8000]" -q normal -- plink --noweb \
    --bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged_unpruned \
    --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
    --out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/${cohort}_unpruned_gwa_data_freq \
    
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
for cohort in ALSPAC TWINSUK JOINT
for cohort in TWINSUK ALSPAC
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

    # filter="FILTERED_1e2_1"
    # filter="FILTERED_1e5_1"
    # filter="FILTERED_5e6_1"
    filter="FILTERED_5e8_1"

    indir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/MDS

    outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/UNPRUNED/MDS/GW/${filter}/

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

################JOINT COHORT
for chr in 20
do
    mkdir -p /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint/MDS/FILTERED_1e5/CHR${chr}

let i=0 a=36;let j=0;

while [ $i -le $a ]
do
while [ $j -le $a ]
do

bsub -J "plink_run_${i}_${j}_${chr}" -o "%J_plink_run_${i}_${j}_${chr}.o" -M 6000 -R"select[mem>6000] rusage[mem=6000]" -q yesterday -- plink --noweb \
--bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed \
--bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY \
--fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam \
--allow-no-sex \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint/REMOVELIST/chr${chr}.txt.assoc.JOINT.1e-5.removelist \
--read-freq /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint/MDS/chr${chr}.freq.frq \
--genome \
--genome-lists /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint/MDS/tmp.list`printf "%03i\n" $i` \
/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint/MDS/tmp.list`printf "%03i\n" $j` \
--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint/MDS/FILTERED_1e5/CHR${chr}/data.sub.$i.$j

let j=$j+1
done

let i=$i+1
let j=$i

done
done

cat data.sub*genome > data.genome
rm data.sub.*


##4. calculate mds for 4 groups (sequencing center x cohorts)
###########################################################################
#######################TWINSUK/ALSPAC
# filter="FILTERED_5e6_1"
# filter="FILTERED_1e5_1"
# filter="UNFILTERED_1"
# cohort_path="twins"
# cohort="TWINSUK"
cohort_path="alspac"
cohort="ALSPAC"
filter="FILTERED_5e8_1"

for chr in {1..22}
do
outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/${filter}/CHR${chr}
mkdir -p ${outdir}

bsub -J "plink_mds_${chr}" -o "%J_plink_mds_${chr}.o" -M2000 -R"select[mem>=2000] rusage[mem=2000]" -q yesterday -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--chr ${chr} \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic/REMOVELIST/chr${chr}.txt.assoc.logistic.${cohort}.5e-8.removelist \
--allow-no-sex \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
--read-genome /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/${filter}/CHR${chr}/data.genome \
--cluster --mds-plot 10 \
--out ${outdir}/chr${chr}
done

--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic/REMOVELIST/chr${chr}.txt.assoc.logistic.${cohort}.5e-6.removelist \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic/REMOVELIST/chr${chr}.txt.assoc.logistic.${cohort}.1e-5.removelist \
#######################ALSPAC
for chr in 20
do
bsub -J "plink_mds_${chr}" -o "%J_plink_mds_${chr}.o" -M2000 -R"select[mem>=2000] rusage[mem=2000]" -q yesterday -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--chr ${chr} \
--allow-no-sex \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/3621_ALSPAC_samples.keeplist \
--read-genome /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/UNFILTERED/CHR${chr}/data.genome \
--cluster --mds-plot 10 \
--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/UNFILTERED/CHR${chr}/chr${chr}
done

--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/REMOVELIST/chr${chr}.txt.assoc.ALSPAC.5e-6.removelist \
##############################JOINT
for chr in 1 {3..22}
do
bsub -J "plink_mds_${chr}" -o "%J_plink_mds_${chr}.o" -M3000 -R"select[mem>=3000] rusage[mem=3000]" -q yesterday -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--chr ${chr} \
--allow-no-sex \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_JOINT_samples.keeplist \
--read-genome /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/UNFILTERED/CHR${chr}/data.genome \
--cluster --mds-plot 10 \
--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/UNFILTERED/CHR${chr}/chr${chr}
done


--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint/REMOVELIST/chr${chr}.txt.assoc.JOINT.1e-5.removelist \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint/REMOVELIST/chr${chr}.txt.assoc.JOINT.5e-8.removelist \

##############################GWAS
# cohort="TWINSUK"
# cohort="ALSPAC"

bsub -J "plink_run" -o "%J_plink_run.o" -M 6000 -R"select[mem>6000] rusage[mem=6000]" -q yesterday -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/3621_${cohort}_samples.keeplist \
--freq \
--out gwa_data_freq

cohort="JOINT"
cohort_path="joint_NEW"
filter="FILTERED_1e2_WR"
maf=0.01 #only common variants
# max_maf=0.05
outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/MDS/GWA/${filter}
mkdir -p ${outdir}

let i=0 a=36;let j=0;

while [ $i -le $a ]
do
while [ $j -le $a ]
do

bsub -J "plink_run_${i}_${j}" -o "${outdir}/%J_plink_run_${i}_${j}.o" -P "analysis-rd" -M 4000000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--allow-no-sex \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/MDS/GWA/All.txt.assoc.logistic.${cohort}.1e-2.joint_WR_gwa.removelist \
--read-freq /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/MDS/GWA/gwa_data_freq.frq \
--maf ${maf} \
--genome \
--genome-lists /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/MDS/tmp.list`printf "%03i\n" $i` \
/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/MDS/tmp.list`printf "%03i\n" $j` \
--out ${outdir}/data.sub.$i.$j

let j=$j+1
done

let i=$i+1
let j=$i

done

#add for exclusion
#change for alspac
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/MDS/GWA/All.txt.assoc.logistic.${cohort}.1e-2.joint_gwa.removelist \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/logistic/REMOVELIST/All.txt.assoc.logistic.JOINT.1e-5.removelist \
--read-freq /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/GW/gwa_data_freq.frq \
--genome-lists /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/tmp.list`printf "%03i\n" $i` \
/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/tmp.list`printf "%03i\n" $j` \
--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/GW/UNFILTERED/data.sub.$i.$j
#
#
#####################################################################################################################
############################################ MDS
########## TWINSUK / ALSPAC
# cohort_path="alspac"
# cohort="ALSPAC"
# cohort="TWINSUK"
# cohort_path="twins"
cohort="JOINT"
cohort_path="joint_NEW"
filter="FILTERED_1e2_WR"
# outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/GW/${filter}
outdir=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/MDS/GWA/${filter}
mkdir -p ${outdir}

bsub -J "plink_mds_all" -o "${outdir}/%J_plink_mds_all.o" -M4000000 -P "analysis-rd" -R"select[mem>=4000] rusage[mem=4000]" -q normal -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--allow-no-sex \
--maf 0.01 \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/MDS/GWA/All.txt.assoc.logistic.${cohort}.1e-2.joint_WR_gwa.removelist \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
--read-genome /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/MDS/GWA/${filter}/data.genome \
--cluster --mds-plot 10 \
--out ${outdir}/all_chr_merged

####to exclude sites
--read-genome /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/GW/${filter}/data.genome \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/GW/All.txt.assoc.logistic.${cohort}.1e-2.joint_WR_gwa.removelist \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/assoc/REMOVELIST/All.txt.assoc.TWINSUK.1e-5.removelist \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/assoc/REMOVELIST/All.txt.assoc.TWINSUK.5e-6.removelist \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/assoc/REMOVELIST/All.txt.assoc.TWINSUK.1e-5.removelist \


########## ALSPAC

bsub -J "plink_mds_all" -o "%J_plink_mds_all.o" -M4000000 -P "analysis-rd" -R"select[mem>=4000] rusage[mem=4000]" -q yesterday -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--allow-no-sex \
--maf 0.01 \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/assoc/REMOVELIST/All.txt.assoc.ALSPAC.1e-5.removelist \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_ALSPAC_samples.keeplist \
--read-genome /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/GW/FILTERED_1e5_1/data.genome \
--cluster --mds-plot 10 \
--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/GW/FILTERED_1e5_1/all_chr_merged

###to exclude
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/assoc/REMOVELIST/All.txt.assoc.ALSPAC.5e-6.removelist \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/assoc/REMOVELIST/All.txt.assoc.ALSPAC.5e-8.removelist \

########## JOINT

bsub -J "plink_mds_all" -o "%J_plink_mds_all.o" -M4000 -R"select[mem>=4000] rusage[mem=4000]" -q yesterday -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--allow-no-sex \
--maf 0.01 \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/logistic/REMOVELIST/All.txt.assoc.logistic.JOINT.1e-5.removelist \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_JOINT_samples.keeplist \
--read-genome /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/GWA/FILTERED_1e5_1/data.genome \
--cluster --mds-plot 10 \
--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/GWA/FILTERED_1e5_1/all_chr_merged


####to exclude
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/logistic/REMOVELIST/All.txt.assoc.logistic.JOINT.5e-6.removelist \
--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/logistic/REMOVELIST/All.txt.assoc.logistic.JOINT.5e-8.removelist \




##############################GWAS - splitted by chr
# cohort="TWINSUK"
# cohort="ALSPAC"

# cohort="JOINT"
# plink --noweb \
# --bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
# --allow-no-sex --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_${cohort}_samples.keeplist \
# --chr ${chr}
# --freq \
# --out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/UNFILTERED/CHR${chr}/chr${chr}_freq

for chr in {1..22} X
do
echo "PROCESSING CHR${chr} "
echo "############################"
mkdir -p /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/UNFILTERED/CHR${chr}

let i=0 a=36;let j=0;

while [ $i -le $a ]
do
while [ $j -le $a ]
do

bsub -J "plink_run_${i}_${j}_${chr}" -o "%J_plink_run_${i}_${j}_${chr}.o" -M 4000 -R"select[mem>4000] rusage[mem=4000]" -q yesterday -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--chr ${chr} \
--allow-no-sex \
--read-freq /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/UNFILTERED/CHR${chr}/chr${chr}_freq.frq \
--genome \
--genome-lists /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/tmp.list`printf "%03i\n" $i` \
/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/tmp.list`printf "%03i\n" $j` \
--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/UNFILTERED/CHR${chr}/data.sub.$i.$j

let j=$j+1
done

let i=$i+1
let j=$i

done
done

bsub -J "plink_mds_all" -o "%J_plink_mds_all.o" -M6000 -R"select[mem>=6000] rusage[mem=6000]" -q yesterday -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--allow-no-sex \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_JOINT_samples.keeplist \
--read-genome /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/GWA/FILTERED_5e8_1/data.genome \
--cluster --mds-plot 10 \
--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/GWA/UNFILTERED_1/all_chr


--exclude /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/logistic/REMOVELIST/All.txt.assoc.logistic.JOINT.5e-8.removelist \
--out /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/GWA/FILTERED_5e8_1/all_chr
###############################################################
## 5. plots
plink \
--noweb \
--bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chrX.bed \
--bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chrX.bim.COPY \
--fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chrX.fam \
--allow-no-sex \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/TWINSUK_samples.keeplist \
--pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/TWINSUK_samples_pheno.centre \
--assoc \
--out chrX.txt

sed -i 's/ \+/ /g' chr*.txt.assoc; sed -i 's/^ //g' chr*.txt.assoc; sed -i 's/ $//g' chr*.txt.assoc; sed -i 's/ /\t/g' chr*.txt.assoc
sed -i 's/ \+/ /g' chr*.txt.assoc.adjusted;sed -i 's/^ //g' chr*.txt.assoc.adjusted; sed -i 's/ $//g' chr*.txt.assoc.adjusted; sed -i 's/ /\t/g' chr*.txt.assoc.adjusted


/lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/merged_gwas/MDS/PCA

let i=0 a=15;let j=0;

while [ $i -le $a ]
do
while [ $j -le $a ]
do

bsub -J "plink_run_${i}_${j}" -o "%J_plink_run_${i}_${j}.o" -M 4000 -R"select[mem>4000] rusage[mem=4000]" -q basement -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/merged_gwas/fvg_merged_geno_overlap \
--allow-no-sex \
--read-freq /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/merged_gwas/MDS/fvg_merged_geno_overlap.frq \
--genome \
--genome-lists /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/merged_gwas/MDS/tmp.list`printf "%03i\n" $i` \
/lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/merged_gwas/MDS/tmp.list`printf "%03i\n" $j` \
--out /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/merged_gwas/MDS/PCA/data.sub.$i.$j

let j=$j+1
done

let i=$i+1
let j=$i

done


bsub -J "plink_mds_all" -o "%J_plink_mds_all.o" -M4000 -R"select[mem>=4000] rusage[mem=4000]" -q basement -- plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/merged_gwas/fvg_merged_geno_overlap \
--allow-no-sex \
--read-genome /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/merged_gwas/MDS/PCA/data.genome \
--cluster --mds-plot 10 \
--out /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/merged_gwas/MDS/PCA/fvg_merged_geno_overlap



#####################REMOVE region for batch effect
# filter="FILTERED_5e6_1"
# filter="FILTERED_1e5_1"
# filter="UNFILTERED_1"
# filter="FILTERED_5e8_1"
# chr=4
# cohort_path="alspac"
# cohort="ALSPAC"
cohort_path="twins"
cohort="TWINSUK"
for chr in 4 6 8 9 12 20
do
# cat /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic/REMOVELIST/chr${chr}.txt.assoc.logistic.${cohort}.1e-2.removelist <(tail -n+2 /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/chr${chr}_cluster_center.txt.assoc.logistic.mod | awk '{if($9 <= 0.01) print $2}') | sort | uniq > /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/${cohort}_CHR${chr}_joint_1e-2_gwa.removelist
cat /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic/REMOVELIST/chr${chr}.txt.assoc.logistic.${cohort}.1e-2.removelist /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/${cohort}_${chr}_detail.RR.list | sort | uniq > /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/${cohort}_CHR${chr}_joint_1e-2_WRR_gwa.removelist
done
# cat /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic/REMOVELIST/chr${chr}.txt.assoc.logistic.${cohort}.1e-2.removelist <(tail -n+2 /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/${cohort}_${chr}_detail.region | awk '{if($9 <= 0.01) print $1}') | sort | uniq > /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/${cohort}_CHR${chr}_joint_1e-2_gwa.removelist
# cat /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic/REMOVELIST/chr${chr}.txt.assoc.logistic.${cohort}.1e-2.removelist <(tail -n+2 /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/chr${chr}_cluster_center.txt.assoc.logistic.mod| awk '{if($9 <= 0.00001) print $1}' )| sort | uniq > /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/${cohort}_CHR${chr}_joint.removelist
# cat /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/logistic/REMOVELIST/chr${chr}.txt.assoc.logistic.${cohort}.1e-5.removelist <(tail -n+2 /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/chr${chr}_cluster_center.txt.assoc.logistic.mod| awk '{if($9 <= 0.00001) print $1}' )| sort | uniq > /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/${cohort_path}/PRUNED/MDS/UNFILTERED_1/CHR${chr}/logistic/${cohort}_CHR${chr}_joint.removelist

#create removelists
cat CHR4/logistic/TWINSUK_4_detail.RR.list CHR6/logistic/TWINSUK_6_detail.RR.list CHR8/logistic/TWINSUK_8_detail.RR.list CHR9/logistic/TWINSUK_9_detail.RR.list CHR12/logistic/TWINSUK_12_detail.RR.list CHR20/logistic/TWINSUK_20_detail.RR.list > TWINSUK_all_detail.RR.list

cat CHR4/logistic/ALSPAC_4_detail.RR.list CHR6/logistic/ALSPAC_6_detail.RR.list CHR8/logistic/ALSPAC_8_detail.RR.list CHR9/logistic/ALSPAC_9_detail.RR.list CHR12/logistic/ALSPAC_12_detail.RR.list CHR20/logistic/ALSPAC_20_detail.RR.list > ALSPAC_all_detail.RR.list

#removelist for gwa sets
cat /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/logistic/REMOVELIST/All.txt.assoc.logistic.ALSPAC.1e-2.removelist /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/UNFILTERED_1/ALSPAC_all_detail.RR.list | sort | uniq > /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/GW/All.txt.assoc.logistic.ALSPAC.1e-2.joint_WR_gwa.removelist

cat /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/logistic/REMOVELIST/All.txt.assoc.logistic.TWINSUK.1e-2.removelist /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/UNFILTERED_1/TWINSUK_all_detail.RR.list | sort | uniq > /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/GW/All.txt.assoc.logistic.TWINSUK.1e-2.joint_WR_gwa.removelist

#create the removelist for the JOINT dataset
cat <(sort /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/GW/All.txt.assoc.logistic.ALSPAC.1e-2.joint_WR_gwa.removelist | uniq) <(sort /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/GW/All.txt.assoc.logistic.TWINSUK.1e-2.joint_WR_gwa.removelist | uniq) | sort |uniq > All.txt.assoc.logistic.JOINT.1e-2.joint_WR_gwa.removelist

#extract data from the regions
chr=20
outdir=CHR${chr}/logistic/GENO_REGION


plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_TWINSUK_samples.keeplist \
--chr ${chr} \
--from-bp 1267130 --to-bp 1793951 \
--make-bed  \
--out ${outdir}/chr${chr}_1267130_1793951


chr=20
outdir=CHR${chr}/logistic/GENO_REGION

plink --noweb \
--bfile /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/input/all_chr_merged \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_ALSPAC_samples.keeplist \
--chr ${chr} \
--from-bp 1267130 --to-bp 1762433 \
--make-bed  \
--out ${outdir}/chr${chr}_1267130_1762433


TWINSUK         
CHR     START   END SIZE
8   6949495 13311013    6.361518
9   38795517    70000798    31.205281
12  7564190 8756791 1.192601
20  1267130 1793951 0.526821
            
            
            
ALSPAC          
CHR     START   END SIZE
8   6948647 12703085    5.754438
9   39013118    74031576    35.018458
12  8398794 8726581 0.327787
20  1267130 1762433 0.495303
