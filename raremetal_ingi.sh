#!/usr/local/bin/bash
#
#
#Renamed samples in vcf files for all cohorts
#bsub -J"rename samples" -o"%J_rename.o" -M2000 -R"select[mem>=2000] rusage[mem=2000]" -q yesterday -- bcftools reheader -s conversion_table_EGa_clinic_id.txt CARL_20140710.all.vcf.gz -o CARL_20140710.all.renamed.vcf.gz

#bsub -J"rename samples" -o"%J_rename.o" -M2000 -R"select[mem>=2000] rusage[mem=2000]" -q yesterday -- bcftools reheader -s conversion_table_EGa_clinic_id.txt CARL_20140710.all.vcf.gz -o CARL_20140710.all.renamed.vcf.gz
#

######
# extract phenotype info for WGS samples and create the ped formatted file from gemma's format
#FVG:
# awk 'FNR==NR { a[$1]=$2; next } $1 in a { print a[$1],a[$1],0,0,$4 }' /lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140517_RELEASE/WGS_FVG_id_conversion_complete_sorted.list /lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/fvg_seq.samples > /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/fvg_seq_sex.samples
# lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_gemma_pheno.txt | cut -f -5,7- -d " " > /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno.txt 
# lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_gemma_pheno.txt | cut -f -5,7- -d " " > /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno.txt 
# lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_gemma_pheno.txt | cut -f -5,7- -d " " > /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno.txt 
# via Torinlustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_gemma_pheno.txt | cut -f -5,7- -d " " > /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno.txt 

# #VBI:
# awk 'FNR==NR { a[$1]=$2; next } $1 in a { print a[$1],a[$1],0,0,$4 }' /lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/20140518_RELEASE/VBI_SEQ2HSR_table.txt /lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/vbi_seq_sex.samples > /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/vbi_seq_sex.samples
# # awk 'FNR==NR { a[$1]=$0; next } $1 in a { print a[$1],$0 }' /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/vbi_seq_sex.samples /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_gemma_pheno.txt | cut -f -5,7- -d " " > /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno.txt 

# #CARL:
# awk 'FNR==NR { a[$1]=$2; next } $1 in a { print a[$1],a[$1],0,0,$4 }' /lustre/scratch113/projects/carl_seq/20140410/20140410_VQSR2.5_reapply_138_excl/20140710_RELEASE/conversion_table_EGa_clinic_id.txt /lustre/scratch113/projects/carl_seq/20140410/20140410_VQSR2.5_reapply_138_excl/carl_seq_sex.samples > /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/carl_seq_sex.samples
# awk 'FNR==NR { a[$1]=$0; next } $1 in a { print a[$1],$0 }' /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/carl_seq_sex.samples /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/PHENO/LIPIDS_gemma_pheno.txt | cut -f -5,7- -d " " > /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/PHENO/LIPIDS_raremetalworker_pheno.txt 

#now run raremetalworker for the selected region

# /nfs/team151/software/raremetal_4.13.1/raremetal/bin/raremetalworker
#
#Run with: 
#
#FVG:
#mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt|cut -f 1 -d " "`;bsub -J "raremetalw[1-${size}]" -o "LOGS/%J_raremetal_w.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/raremetal_ingi.sh /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140517_RELEASE/FVG_20140517.all.renamed.vcf.gz /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/RAREMETALWORKER /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/raremetal_regions.list /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/RAREMETALWORKER/GRM/FVG_WGS.Empirical.Kinship.gz
#mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt|cut -f 1 -d " "`;bsub -J "raremetalw[1-${size}]" -o "LOGS/%J_raremetal_w.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/raremetal_ingi.sh /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140517_RELEASE/FVG_20140517.all.renamed.vcf.gz /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/RAREMETALWORKER/04132015 /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/raremetal_regions.list /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/RAREMETALWORKER/GRM/FVG_WGS.Empirical.Kinship.gz
#mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt|cut -f 1 -d " "`;bsub -J "raremetalw[1-${size}]" -o "LOGS/%J_raremetal_w.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/raremetal_ingi.sh /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140517_RELEASE/FVG_20140517.all.renamed.vcf.gz /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/RAREMETALWORKER/04282015/INGI /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/raremetal_regions.list /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/RAREMETALWORKER/GRM/FVG_WGS.Empirical.Kinship.gz
#HRC-FVG genotypes data:
#bsub -J"rename samples" -o"%J_rename.o" -M2000 -R"select[mem>=2000] rusage[mem=2000]" -q yesterday -- bcftools reheader -s /lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140517_RELEASE/WGS_FVG_id_conversion_complete_sorted.list /lustre/scratch113/projects/fvg_seq/hrc-subset/HRC.r1.FVG.GRCh37.ALL.shapeit3.mac5.genotypes.vcf.gz -o /lustre/scratch113/projects/fvg_seq/hrc-subset/HRC.r1.FVG.GRCh37.ALL.shapeit3.mac5.genotypes.renamed.vcf.gz
#mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt|cut -f 1 -d " "`;bsub -J "raremetalw[1-${size}]" -o "LOGS/%J_raremetal_w.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/raremetal_ingi.sh /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/fvg_seq/hrc-subset/HRC.r1.FVG.GRCh37.ALL.shapeit3.mac5.genotypes.renamed.vcf.gz /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/RAREMETALWORKER/04282015/HRC /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/raremetal_regions.list /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/RAREMETALWORKER/GRM/FVG_WGS.Empirical.Kinship.gz
#VBI:
#mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno_list.txt|cut -f 1 -d " "`;bsub -J "raremetalw[1-${size}]" -o "LOGS/%J_raremetal_w.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/raremetal_ingi.sh /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/20140518_RELEASE/VBI_20140518.all.renamed.vcf.gz /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/RAREMETALWORKER /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/raremetal_regions.list /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/RAREMETALWORKER/GRM/VBI_WGS.Empirical.Kinship.gz
#mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno_list.txt|cut -f 1 -d " "`;bsub -J "raremetalw[1-${size}]" -o "LOGS/%J_raremetal_w.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/raremetal_ingi.sh /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/20140518_RELEASE/VBI_20140518.all.renamed.vcf.gz /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/RAREMETALWORKER/04282015/INGI /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/raremetal_regions.list /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/RAREMETALWORKER/GRM/VBI_WGS.Empirical.Kinship.gz
#HRC-VBI genotypes data:
#bsub -J"rename samples" -o"%J_rename.o" -M2000 -R"select[mem>=2000] rusage[mem=2000]" -q yesterday -- bcftools reheader -s /lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/20140518_RELEASE/VBI_SEQ2HSR_table.txt /lustre/scratch113/projects/esgi-vbseq/hrc-subset/HRC.r1.VBSEQ.GRCh37.ALL.shapeit3.mac5.genotypes.vcf.gz -o /lustre/scratch113/projects/esgi-vbseq/hrc-subset/HRC.r1.VBSEQ.GRCh37.ALL.shapeit3.mac5.genotypes.renamed.vcf.gz
#mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno_list.txt|cut -f 1 -d " "`;bsub -J "raremetalw[1-${size}]" -o "LOGS/%J_raremetal_w.%I.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/raremetal_ingi.sh /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/esgi-vbseq/hrc-subset/HRC.r1.VBSEQ.GRCh37.ALL.shapeit3.mac5.genotypes.renamed.vcf.gz /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/RAREMETALWORKER/04282015/HRC /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/raremetal_regions.list /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_VBI/RAREMETALWORKER/GRM/VBI_WGS.Empirical.Kinship.gz

#CARL:
#mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/PHENO/LIPIDS_raremetalworker_pheno_list.txt|cut -f 1 -d " "`;bsub -J "raremetalw[1-${size}]" -o "LOGS/%J_raremetal_w.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner_par.sh -s ~/Work/bash_scripts/raremetal_ingi.sh /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/carl_seq/20140410/20140410_VQSR2.5_reapply_138_excl/20140710_RELEASE/CARL_20140710.all.renamed.vcf.gz /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/PHENO/LIPIDS_raremetalworker_pheno.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/PHENO/LIPIDS_raremetalworker_pheno_list.txt /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/RAREMETALWORKER /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/raremetal_regions.list /lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_CARLANTINO/RAREMETALWORKER/GRM/CARL_WGS.Empirical.Kinship.gz
#################################### Run RAREMETALWORKER for 107 regions #####################################
## Trait

trait=$2

## INPUT files
vcfpath=$3 #/lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140517_RELEASE/FVG_20140517.all.renamed.vcf.gz
ped_file=$4 #/lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno.txt
dat=$5 #/lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/PHENO/LIPIDS_raremetalworker_pheno_list.txt
out_folder=$6 #/lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/RAREMETALWORKER
regions=$7 #/lustre/scratch113/projects/uk10k/users/mc14/UK10K_replica/INGI_FVG/raremetal_regions.list
kinship=$8 #kinship matrix path (calculated from WGS data with RAREMETALWORKER!!)

mkdir -p $out_folder
echo $vcfpath
echo $ped_file
echo $dat
echo $out_folder
echo $regions

#regions=/nfs/users/nfs_k/kw8/uk10k/lipids/data/meta-skat_lipids_independent_regions_for_replication_chrX.txt
cat $regions | while read -a line
do
    echo $trait
    chr=${line[0]}
    start=${line[1]}
    end=${line[2]}
    echo "$chr $start $end"

    #select vcf file
    vcf=${out_folder}/${trait}_chr${chr}"_"${start}"_"${end}.vcf.gz

    # create vcf for the region extracting the right samples
    cut -f 1 -d " " ${ped_file} > ${out_folder}/${trait}_samples.list
    bcftools view ${vcfpath} -r ${chr}:${start}-${end} -S ${out_folder}/${trait}_samples.list -O z -o ${vcf}
    tabix -p vcf -f ${vcf}
    #extact samples' code to sort phenotype file
    # tabix -H $vcf | cut -f 10- | tr "\t" "\n" > $out_folder/WGS_samples.list
    zcat ${vcf}| fgrep "#CHROM" | cut -f 10- | tr "\t" "\n" > ${out_folder}/WGS_${trait}_samples.list
    #sort phenotype file based on genotypes...
    awk 'FNR==NR { a[$1]=$0; next } $1 in a { print a[$1]}' ${ped_file} ${out_folder}/WGS_${trait}_samples.list > ${ped_file}.${trait}.sorted.txt

    outfile=${out_folder}/${trait}"_chr"${chr}"_"${start}"_"$end

    # /nfs/team151/software/raremetal_4.13.1/raremetal/bin/raremetalworker --ped $ped_file.sorted.txt --kinGeno --kinSave --prefix $outfile --vcf $vcf
    #/nfs/team151/software/raremetal_4.13.1/raremetal/bin/raremetalworker --ped $ped_file.sorted.txt --dat $dat --cpu 4 --kinGeno --kinSave --kinOnly --vcf $vcfpath --prefix $outfile
    # /nfs/team151/software/raremetal_4.13.1/raremetal/bin/raremetalworker --ped $ped_file.sorted.txt --dat $dat --kinFile $outfile.Empirical.Kinship.gz --vcf $vcf --traitName $trait --prefix $outfile
    /nfs/team151/software/raremetal_4.13.1/raremetal/bin/raremetalworker --ped ${ped_file}.${trait}.sorted.txt --dat ${dat} --kinFile ${kinship} --vcf ${vcf} --traitName ${trait} --prefix ${outfile}

    bgzip -f ${outfile}.${trait}.singlevar.score.txt
    bgzip -f ${outfile}.${trait}.singlevar.cov.txt
    tabix -f -s 1 -b 2 -e 2 -c "#" ${outfile}.${trait}.singlevar.score.txt.gz
    tabix -f -s 1 -b 2 -e 2 -c "#" ${outfile}.${trait}.singlevar.cov.txt.gz

done

