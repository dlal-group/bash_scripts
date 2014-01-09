#!/usr/local/bin/bash

#general script to launch fast jobs

#cat esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf.gz

#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf.gz
#zcat esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf | wc -l
#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf.gz chr $1 >> esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.vcf
#cat ~/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/VBI/VBI.vcf.header esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.vcf > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.HEADED.vcf

#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.vcf.gz chr $1 | wc -l

#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf.gz chr $1 | cut -f 1-3 | fgrep "	." > NO_RSID_LISTS/VBI_1kg_NOT_OVERLAPP_no_rs_id_chr$1.tab
#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf.gz chr $1 | fgrep -w -f <(cut -f 1,2 NO_RSID_LISTS/VBI_1kg_NOT_OVERLAPP_no_rs_id_chr$1.tab) | bgzip -c > NO_RSID_VCF/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.vcf.gz

#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.vcf.gz
#cat ~/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/VBI/VBI.vcf.header esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.NO_RSID.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.NO_RSID.vcf.gz
#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.NO_RSID.vcf.gz
#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf.gz chr $1 | cut -f 1-3 | fgrep -v "	." > WITH_RSID/chr$1.re_ann.NOT_OVERLAP.rsID.tab
#fgrep -w -f <(cut -f 3 LISTS/chr$1.re_ann.NOT_OVERLAP.rsID.tab) <(tabix ../esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.vcf.gz chr $1) | bgzip -c > chr$1.re_ann.NOT_OVERLAP.rsID.vcf.gz
#tabix chr$1.re_ann.NOT_OVERLAP.rsID.vcf.gz chr $1 | fgrep AF_EUR | bgzip -c > chr$1.re_ann.NOT_OVERLAP.rsID.AF_EUR.vcf.gz
#fgrep -w -f <(tabix ../chr$1.re_ann.NOT_OVERLAP.rsID.vcf.gz chr $1 | cut -f 3) <(tabix ~/lustre110_home/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/dbSNP_splitted/chr$1_dbSNP_b137.tab.gz chr $1) > chr$1.re_ann.NOT_OVERLAP.rsID.dbSNP.list

#command used to create a map file for plot
#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.vcf.gz chr $1 | cut -f 1-3 > MAP/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.map

#command to create files for plot for not overlapping sites with DP
#vcf-query <(zcat ../esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.NO_RSID.vcf.gz | grep ^# ; tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.vcf.gz chr $1) -f "%CHROM\t%POS\t%ID\t%INFO/AN\t%INFO/AC\t%INFO/DP\t%INFO/HWE\n" > MAP/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.map

#command to create files for plot all sites with DP
#vcf-query ../esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.vcf.gz -f "%CHROM\t%POS\t%ID\t%INFO/AN\t%INFO/AC\t%INFO/DP\t%INFO/HWE\n" > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.map

#command to extract sites with rsId
#fgrep -v "	." <(fgrep -w PASS ../UK10K.chr$1.snps.overlap.map | cut -f 1-3) > UK10K.chr$1.snps.overlap.rsID.map

#command to extract all sites with rsID
#fgrep -v "	." <(cut -f 1-3 ../UK10K.chr$1.snps.overlap.map) > UK10K.chr$1.snps.overlap.rsID.map
#fgrep -v "	." <(cut -f 1-3 ../UK10K.chr$1.snps.overlap.passed.map) > UK10K.chr$1.snps.overlap.passed.rsID.map

#command to extract all sites PASSED in UK10K 
#fgrep -w PASS ../UK10K.chr$1.snps.overlap.map > UK10K.chr$1.snps.overlap.passed.map

#command to extract only non SNV sites from dbSNP137 annotation file
#vcf-query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/dbSNPBuildID\t%INFO/VC\n' 00-All.vcf.gz | grep -v SNV |sort -k1,1d -k2,2n -k4,4d -k5,5d | coalesce-duplicates | bgzip -c > annots-rsIDs-dbSNPv137-no_SNV.2012-06-16.tab.gz


#command to create files for plot for not overlapping sites with DP
#vcf-query <(zcat ~/lustre110_home/GENOTIPI/COMPARISON/NOT_OVERLAPPING/NOT_OVERLAPPING_LIST/with_genotypes/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.vcf.gz | grep ^# ; tabix ../esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.vcf.gz chr $1) -f "%CHROM\t%POS\t%ID\t%INFO/AN\t%INFO/AC\t%INFO/DP\t%INFO/HWE\n" > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.former_SNPs_removed.map

#command to extract allelels
#tabix ../esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_former_rsID.vcf.gz chr $1 | cut -f 1-5 > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_former_rsID.map

#command to create tab spaced files for putative novel sites with DP and HWE
#vcf-query <(zcat ~/lustre110_home/GENOTIPI/COMPARISON/NOT_OVERLAPPING/NOT_OVERLAPPING_LIST/with_genotypes/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.vcf.gz | grep ^# ; tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.vcf.gz chr $1) -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/DP\t%INFO/HWE\n" > MAP/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.map
#
#
#command to extract all sites with ac>2
#awk '{OFS="\t"}{if ($7 > 2) print $0}' esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.map > ac_gt2/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.ac_gt2.map

#command to extract sites with ac=2 but only heterozygotic
#vcf-query <(zcat ~/lustre110_home/GENOTIPI/COMPARISON/NOT_OVERLAPPING/NOT_OVERLAPPING_LIST/with_genotypes/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.vcf.gz | grep ^# ; tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.vcf.gz chr $1 | fgrep AC=2 | fgrep -v "1|1") -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\t%INFO/DP\t%INFO/HWE\n" > MAP/ac_eq2/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.ac_eq2.map

#command to extract all sites from a list
#fgrep -w -f <(cut -f 1,2 MAP/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.gte2_indivs.map) <(tabix /nfs/users/nfs_m/mc14/lustre110_home/GENOTIPI/COMPARISON/NOT_OVERLAPPING/PUTATIVE_NOVEL/NEW_RUN/former_rsID_filtered/UK10K_FILTERED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.vcf.gz chr $1) | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr$1.re_ann.NOT_OVERLAP.NO_RSID.no_UK10K.gte2_indivs.vcf.gz

#command toi extract sites from 1000Gp files
#vcf-query /lustre/scratch109/sanger/mc14/GENOTIPI/REF_PANEL/1000G/VCF/ALL.chr${1}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\n" > /lustre/scratch109/sanger/mc14/GENOTIPI/REF_PANEL/1000G/VCF/ALL.chr${1}.phase1_release_v3.20101123.snps_indels_svs.genotypes.tab

#command to cancat vcffiles
#vcf-concat /lustre/scratch110/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr7.re_ann.vcf.gz /lustre/scratch110/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/INDELS/CHR7/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.INDELS.7.vcf.gz | bgzip -c > /lustre/scratch109/sanger/mc14/GENOTIPI/REF_PANEL/VBI/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.chr7.vcf.gz

#command to create haplotype formatted file
#vcf-query /lustre/scratch109/sanger/mc14/GENOTIPI/REF_PANEL/1000G/VCF/ALL.chr7.phase1_release_v3.20101123.snps_indels_svs.genotypes.final.vcf.gz -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GTR]\n' | tr "|" " " | tr "/" " " | tr "\t" " "| cut -f 6- -d " " | gzip -c > /lustre/scratch109/sanger/mc14/GENOTIPI/REF_PANEL/1000G/VCF/HAPS/CHR7/chr7.haps.gz

#(zgrep "^#" /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr7.re_ann.vcf.gz;tabix /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.gz ${1} | awk '($5==".")' ) | bgzip -c > /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.to_fix.gz

#tabix -f -p vcf /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.to_fix.gz

#(zgrep "^#" /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr7.re_ann.vcf.gz;tabix /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.gz ${1} | awk '($5!=".")' ) | bgzip -c > /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.OK.gz

#tabix -f -p vcf /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.OK.gz

#commands to fix the missing alt allele problem on SNPs only!
#mkdir -p /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/SPLITTED/fixed/

#(zgrep "^#" /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/SPLITTED/tofix/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.to_fix.gz;join -t "	" -1 2 -2 2 <(tabix /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/SPLITTED/tofix/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.to_fix.gz ${1}) /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/allele_fixing/fixed_missing_allele_chr${1}.map | awk '{printf "%s\t%s\t%s\t%s\t%s\t",$2,$1,$3,$4,$NF; for(i=6;i<=(NF-4);i++) printf "%s\t",$i; printf "\n"}') | bgzip -c > /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/SPLITTED/fixed/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.fixed.gz

#tabix -f -p vcf /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/SPLITTED/fixed/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.fixed.gz

#command to join back fixed files
#(zgrep "^#" /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/SPLITTED/OK/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.OK.gz;(tabix /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/SPLITTED/OK/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.OK.gz ${1};tabix /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/SPLITTED/fixed/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.vcf.fixed.gz ${1}) | sort -g -k2,2) | bgzip -c > /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/FIXED_ALT/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.alt_fixed.vcf.gz

#tabix -f -p vcf /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/FIXED_ALT/esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${1}.re_ann.alt_fixed.vcf.gz

#download all 1KGP ref vcf files
#cd /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/REF_PANEL/1000G/VCF
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr${1}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr${1}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.tbi

#merge geno and legend files
#for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19_20 21_22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53_54
#do
#fgrep -h -v exp_freq_a1 /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/COMBINED_IMPUTATION/CHR7/chr7.${i}.geno_info >> /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/COMBINED_IMPUTATION/CHR7/chr7.geno_info
#done


#for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53
#do
#zcat /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/COMBINED_IMPUTATION/1000GP_UK10K/CHR7/chr7.${i}.geno.gz | gzip -c >> /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/COMBINED_IMPUTATION/1000GP_UK10K/CHR7/chr7.geno.gz
# done

#gzip -c /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/COMBINED_IMPUTATION/CHR7/chr7.geno.gz > /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/COMBINED_IMPUTATION/CHR7/chr7.geno

# for chr in {1..22} X
# do echo "${chr}"
# fgrep -v beta $1.chr${chr}.tab.assoc.txt >> $1.all.tab.assoc.txt
# done
# (fgrep "beta" $1.chr9.tab.assoc.txt;cat $1.all.tab.assoc.txt)| gzip -c > $1.all.tab.assoc.txt.gz

#to create a file with all results
#we can use a compressed mode or a plain mode
# type=`file -ib WEIGHT.chr19.tab.assoc.txt.gz | cut -f 1 -d ";" | cut -f 2 -d "/"`
# if [ $# -eq 2 ]
# then
# 	(zcat $1.chr22.tab.assoc.txt.gz | head -1 ;for chr in {1..22} X;do zgrep -v beta $1.chr${chr}.tab.assoc.txt;done) > $1.all.tab.assoc.txt
# else
	# (head -1 $1.chr22_*.tab.assoc.txt;for chr in {1..22} X;do fgrep -v beta $1.chr${chr}_*.tab.assoc.txt;done) > $1.all.tab.assoc.txt
	#(head -1 $1.chr22.tab.assoc.txt;for chr in {1..22} X;do fgrep -v beta $1.chr${chr}.tab.assoc.txt;done) > $1.all.tab.assoc.txt
	# thr=5e-2
	# mkdir -p SUMMARY
	# (head -1 $1.chr22.tab.assoc.txt;awk -v thr=${thr} '$9<=thr' $1.all.tab.assoc.txt) > SUMMARY/$1.all.tab.assoc.${thr}.txt
	# (head -1 $1.chr22.tab.assoc.txt;zcat $1.all.tab.assoc.txt.gz | awk -v thr=${thr} '$9<=thr') > SUMMARY/$1.all.tab.assoc.${thr}.txt
	#(zcat chr22.txt.assoc.gz | head -1;for chr in {1..22};do zfgrep -v CHISQ chr${chr}.txt.assoc.gz;done) > All.tab.assoc.txt
	#(head -1 chr5.txt.assoc.gz_chi2.plots ;for chr in {1..22};do fgrep -v pval chr5.txt.assoc.gz_chi2.plots;done) > All.txt.assoc.gz_chi2.plots
# fi

#for loop to create imputation scripts
# echo "Processing chr $1"
# for file in /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/REF_PANEL/MERGED/VBI_1000GP/MERGED/CHR$1/HAPS/CHR$1/$1.*.*.hap.gz
# do
# chunk=`echo ${file} | tr "." "\t"| cut -f 2`
# echo "Processing $chunk"
# impute_imputation_step3_combined.sh /nfs/users/nfs_m/mt11/1785_VBI_pre_phasing_370k/plain $1 /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/REF_PANEL/MERGED/VBI_1000GP/MERGED/CHR$1/HAPS/CHR$1 /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/COMBINED_IMPUTATION/VBI_1000GP ${chunk}
# done

#to extract data for plots
# thr=1e-1
# awk -v thr=${thr} '$9<=thr' $1.all.tab.assoc.txt | cut -f -3,5,6,9 > $1.all.tab.assoc.${thr}.to_plot

#to extract indels and snps count from bimbam files
# for chr in {1..22} X
# do
# indels=`cut -f -3 -d ',' chr${chr}.bimbam | tr "," " " |awk '{if(length($2)>1 || length($3)>1) print $0}' |wc -l`
# snps=`cut -f -3 -d ',' chr${chr}.bimbam | tr "," " " |awk '{if(length($2)==1 && length($3)==1) print $0}' |wc -l`
# tot=$[indels + snps]
# echo "CHR${chr} ${indels} ${snps} ${tot}"
# done

#we want to extract also filtered numbers of sites
#with those filters:
#INFO ($7) >= 0.3
#HWE ($14) >1e-6 
#MAF ($12) >0.05'
# file=$1
# findels=`zcat ${file} | tail -n+2 |awk '$7 >= 0.3 && $14>1e-6 && $12>0.05'| awk 'length($3)>1 || length($4)>1' | wc -l`
# fsnps=`zcat ${file} | tail -n+2 |awk '$7 >= 0.3 && $14>1e-6 && $12>0.05'| awk 'length($3)==1 && length($4)==1' | wc -l`
# ftot=$[findels + fsnps]
# echo "${file} ${ftot} ${fsnps} ${findels}"

#quick command to reformat bimbams
# zcat /lustre/scratch113/projects/uk10k/users/jh21/imputed/vb/uk10k1kg.shapeit/chr$1.bimbam.gz | tr "," " "|awk '{gsub(/-/," ",$1);print}'| cut -f 1,5- -d " "| tr " " "," > /lustre/scratch113/projects/uk10k/users/mc14/imputed/vbi/chr$1.bimbam

#count sitest for different categories:
#SNPS: biallelic and multiallelic missing data
#indels: biallelic multiallelic missing data
# missing_snps_bi=`tabix $1.vcf.gz $1| grep "\./\." | awk '$8!~/INDEL/ && $5!~/,/'|wc -l`
# missing_snps_multi=`tabix $1.vcf.gz $1| grep "\./\." | awk '$8!~/INDEL/ && $5~/,/'|wc -l`
# missing_indels_bi=`tabix $1.vcf.gz $1| grep "\./\." | awk '$8~/INDEL/ && $5!~/,/'|wc -l`
# missing_indels_multi=`tabix $1.vcf.gz $1| grep "\./\." | awk '$8~/INDEL/ && $5~/,/'|wc -l`

# echo "CHR$1 ${missing_snps_bi} ${missing_snps_multi} ${missing_indels_bi} ${missing_indels_multi}"

# #USE THE SAME COMMAND TO EXTRACT THE ACTUAL VCF file to process and check for values in the PL/GL columns
# (zcat $1.vcf.gz | grep "^#";tabix $1.vcf.gz $1| grep "\./\." | awk '$8!~/INDEL/ && $5!~/,/' )| bgzip -c > missing_snps_bi.$1.vcf.gz
# tabix -p vcf missing_snps_bi.$1.vcf.gz
# bcftools2 query -f '%CHROM\t%POS\t%REF\t%ALT[\t%PL\t%GL]\n' missing_snps_bi.$1.vcf.gz > missing_snps_bi.$1.vcf.tab
# (zcat $1.vcf.gz | grep "^#";tabix $1.vcf.gz $1| grep "\./\." | awk '$8!~/INDEL/ && $5~/,/' )| bgzip -c > missing_snps_multi.$1.vcf.gz
# tabix -p vcf missing_snps_multi.$1.vcf.gz
# bcftools2 query -f '%CHROM\t%POS\t%REF\t%ALT[\t%PL\t%GL]\n' missing_snps_multi.$1.vcf.gz > missing_snps_multi.$1.vcf.tab
# (zcat $1.vcf.gz | grep "^#";tabix $1.vcf.gz $1| grep "\./\." | awk '$8~/INDEL/ && $5!~/,/' )| bgzip -c > missing_indels_bi.$1.vcf.gz
# tabix -p vcf missing_indels_bi.$1.vcf.gz
# bcftools2 query -f '%CHROM\t%POS\t%REF\t%ALT[\t%PL\t%GL]\n' missing_indels_bi.$1.vcf.gz > missing_indels_bi.$1.vcf.tab
# (zcat $1.vcf.gz | grep "^#";tabix $1.vcf.gz $1| grep "\./\." | awk '$8~/INDEL/ && $5~/,/' )| bgzip -c > missing_indels_multi.$1.vcf.gz
# tabix -p vcf missing_indels_multi.$1.vcf.gz
# bcftools2 query -f '%CHROM\t%POS\t%REF\t%ALT[\t%PL\t%GL]\n' missing_indels_multi.$1.vcf.gz > missing_indels_multi.$1.vcf.tab

#extract variants to check if the excluded indels overlap with SNPs
# tabix $1.vcf.gz $1 | cut -f -5 |fgrep -w -f <(cut -f 2 missing_indels_bi.$1.vcf.tab ) > missing_indels_bi.$1.dupchk.tab

#create chr files from 1KGP ALL file and index them
# tabix -h /lustre/scratch111/resources/variation/Homo_sapiens/grch37/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz $1| bgzip > ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.$1.vcf.gz

# tabix -p vcf ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.$1.vcf.gz

#do the same for dbSNP-b138
tabix -h /lustre/scratch113/projects/fvg_seq/variant_refinemet/annotations/dbSNP-b138/00-All.vcf.gz $1 | bgzip > $1.dbsnp_138.vcf.gz

tabix -p vcf $1.dbsnp_138.vcf.gz