#!/bin/bash


# echo '~/scripts/uk10k/map_chr_pos_rsID_chrom_tabix.sh ${LSB_JOBINDEX}' | bsub -J "chr[1-22]%10" -o ~/map_chr_pos_rsID_chrom_tabix_chr%I.out


##### Assign dbSNP IDs to UK10K 809 TwinsUK samples #####

## chomosome
chr=$1

## UK10K split by chromosome 
#uk10k=/lustre/scratch107/projects/uk10k/cohorts/Twins/REL-2011-06-08/v3/vcf/chroms/chr$chr".vcf.gz"
uk10k=/lustre/scratch107/projects/uk10k/cohorts/Twins/REL-2011-06-08/v3/vcf/chroms_sites_filtered/chr$chr".vcf.gz"
## dbSNP134
#dbsnp=/lustre/scratch103/sanger/kw8/dbSNP/chroms/chr$chr".vcf.gz"
## dbSNP135
tabix /lustre/scratch103/sanger/kw8/dbSNP/dbsnp135/00-All.vcf.gz $chr > /lustre/scratch103/sanger/kw8/dbSNP/tmp_chr$chr.txt
gzip /lustre/scratch103/sanger/kw8/dbSNP/tmp_chr$chr.txt
wait
dbsnp=/lustre/scratch103/sanger/kw8/dbSNP/tmp_chr$chr.txt.gz

## output files
#outfile=/lustre/scratch107/projects/uk10k/cohorts/Twins/REL-2011-06-08/v3/vcf/chr_pos_rsID_2011-11-28/chr$chr".txt"
#outfile=/lustre/scratch107/projects/uk10k/cohorts/Twins/REL-2011-06-08/v3/vcf/chr_pos_rsID_2011-12-16_dbSNP135/chr$chr".txt"
outfile=/lustre/scratch107/projects/uk10k/cohorts/Twins/REL-2011-06-08/v3/vcf/chr_pos_rsID_2012-01-02_dbSNP135/chr$chr".vcf"

## for each chromosome walk through UK10K and dbSNP in parallel
## change output file in perl script !!!
perl ~/scripts/uk10k/map_chr_pos_rsID_chrom.pl $chr $uk10k $dbsnp > $outfile
wait
rm $dbsnp