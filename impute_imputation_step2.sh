#!/bin/bash

#Pipeline script for imputation
#this is the case of updated positions with liftover!

#Args (toBe substitute)
#$1 = genotype file path
#$2= chr
#$3= reference file path
#$4= output file
#$5= geno mode TRUE/FALSE

#GENO_PATH=$g_path
#CHR=$

function my_trap_handler()
{
        MYSELF="$0"               # equals to my script name
        LASTLINE="$1"            # argument 1: last line of error occurence
        LASTERR="$2"             # argument 2: error code of last command
        echo "${MYSELF}: line ${LASTLINE}: exit status of last command: ${LASTERR}"

        # do additional processing: send email or SNMP trap, write result to database, etc.
	exit 1
}


if [ $# -lt 4 ]
then
	echo -e "**********************\nWRONG ARGUMENT NUMBER!!!\n**********************"
	echo -e "USAGE:\n impute_imputation_pipeline.sh <genotype files path> <chr> <reference files path> <output file path> [imputation mode]\n"
	echo -e "- <genotype files path> : path for genotypes files"
	echo -e "- <chr> : chromosome number"
	echo -e "- <reference files path> : path for reference files"
	echo -e "- <output files path> : output path"
	echo -e "- [imputation mode] : if specified a 'geno' argument, the imputation will be performed using the genotypes otherwise the pre-phased haplotypes\n"

exit 1
fi

# trap commands with non-zero exit code
#
#trap 'my_trap_handler ${LINENO} $?' ERR


##PART 2: SNP FLIPPING

#Extract frequencies from genotypes
plink --noweb --bfile $1/CHR$2/TO_FLIP/chr$2.to_flip --freq --out $1/CHR$2/TO_FLIP/chr$2_freq_to_flip

#Remove duplicated spaces from the frq file (A1=minor A2=major; the freq in the file refers to A1):
sed -i 's/^ \+//g' $1/CHR$2/TO_FLIP/chr$2_freq_to_flip.frq
sed -i 's/ \+/\t/g' $1/CHR$2/TO_FLIP/chr$2_freq_to_flip.frq

#Print a file with snps and alleles
grep -v CHR $1/CHR$2/TO_FLIP/chr$2_freq_to_flip.frq | cut -f 2,3,4,5 > $1/CHR$2/TO_FLIP/chr$2_snp_check.list

#Extract alleles and frequencies from ref files (A0=ref_allele A1:alt allele;extracted using eur freq;first freq is A1 freq, second freq is maf, so could be the same as aaf)
mkdir -p $3/CHR$2

zcat $3/ALL_1000G_phase1integrated_v3_chr$2_impute.legend.gz | fgrep -v position | awk '{OFS="\t"}{print $1,$3,$4,$8,$12}' > $3/CHR$2/$2_ref.list

#Extract overlapping snps and allelels from the reference file:
fgrep -w -f <(cut -f 1 $1/CHR$2/TO_FLIP/chr$2_snp_check.list) $3/CHR$2/$2_ref.list > $3/CHR$2/$2_overlap_allele.map

#we need to work only on overlapping sites.Extract overlapping sites from genotypes also
fgrep -w -f <(cut -f 1 $3/CHR$2/$2_overlap_allele.map) $1/CHR$2/TO_FLIP/chr$2_snp_check.list > $1/CHR$2/TO_FLIP/chr$2_overlap_snp_check.list

#Extract from genotypes all A/T C/G sites
awk '{OFS="\t"}{ if (($2=="T" && $3=="A")||($2=="A" && $3=="T")||($2=="C" && $3=="G")||($2=="G" && $3=="C")) print $0, "ambig" ; else print $0 ;}' $1/CHR$2/TO_FLIP/chr$2_overlap_snp_check.list | grep ambig > $1/CHR$2/TO_FLIP/chr$2_snp_check_ATCG.list
awk '{OFS="\t"}{ if (($2=="T" && $3=="A")||($2=="A" && $3=="T")||($2=="C" && $3=="G")||($2=="G" && $3=="C")) print $0, "ambig" ; else print $0 ;}' $1/CHR$2/TO_FLIP/chr$2_overlap_snp_check.list | grep -v ambig > $1/CHR$2/TO_FLIP/chr$2_snp_check_noATCG.list

#Extract the same snps from ref files:
fgrep -w -f <(cut -f 1 $1/CHR$2/TO_FLIP/chr$2_snp_check_noATCG.list) $3/CHR$2/$2_overlap_allele.map > $3/CHR$2/$2_overlap_allele_noATCG.list
fgrep -v -w -f <(cut -f 1 $1/CHR$2/TO_FLIP/chr$2_snp_check_noATCG.list) $3/CHR$2/$2_overlap_allele.map > $3/CHR$2/$2_overlap_allele_ATCG.list

#Check if all alleles matches, if not extract the snp list of mismatching alleles and use them to flip the snps (first file always the ref file)
#The result from join command will have the following header:
#RS_ID REF_REF_all ALT_REF_all Alt_all_freq_ref Maf_ref Minor_all_geno Major_all_geno MAF_geno


#FIRST CHECK FOR NON AMBIGUOUS SNPS
join <(sort -k1 $3/CHR$2/$2_overlap_allele_noATCG.list) <(sort -k1 $1/CHR$2/TO_FLIP/chr$2_snp_check_noATCG.list) | awk '{OFS="\t"}{if((($2 == $6) && ($3 == $7 )) || (($2 == $7) && ($3 == $6))) print $0; else print $0,"flip";}' | fgrep -w "flip" > $1/CHR$2/TO_FLIP/chr$2_to_flip_noATGC.list

#NOW CHECK FOR ATCG SNPS
join <(sort -k1 $3/CHR$2/$2_overlap_allele_ATCG.list) <(sort -k1 $1/CHR$2/TO_FLIP/chr$2_snp_check_ATCG.list) | awk '{OFS="\t"}{if ($4 > $5) minor=$2;else minor=$3; if ($6 != minor) print $0}' > $1/CHR$2/TO_FLIP/chr$2_to_flip_ATGC.list

#NOW JOIN the 2 lists
cat $1/CHR$2/TO_FLIP/chr$2_to_flip_noATGC.list $1/CHR$2/TO_FLIP/chr$2_to_flip_ATGC.list | cut -f 1 > $1/CHR$2/TO_FLIP/chr$2_to_flip.list 

#now create the strand files to be used in imputation using pre-phased haplotypes
awk '{OFS=" "}{print $1,"-"}' $1/CHR$2/TO_FLIP/chr$2_to_flip.list > $1/CHR$2/TO_FLIP/chr$2_neg_strand.list
fgrep -v -w -f <(cut -f 1 -d " " $1/CHR$2/TO_FLIP/chr$2_neg_strand.list) $1/CHR$2/TO_FLIP/chr$2_overlap_snp_check.list | awk '{OFS=" "}{print $1,"+"}' > $1/CHR$2/TO_FLIP/chr$2_pos_strand.list

#create the impute input directory
mkdir -p $1/CHR$2/IMPUTE_INPUT

#Join the stranded lists
#extract positions for snps to create the real strand file
#POSITIVE
fgrep -w -f <(cut -f 1 -d " " $1/CHR$2/TO_FLIP/chr$2_pos_strand.list) $1/CHR$2/TO_FLIP/chr$2.to_flip.bim | awk '{OFS=" "}{print $4,"+"}' > $1/CHR$2/TO_FLIP/chr$2_pos.strand

#NEGATIVE
fgrep -w -f <(cut -f 1 -d " " $1/CHR$2/TO_FLIP/chr$2_neg_strand.list) $1/CHR$2/TO_FLIP/chr$2.to_flip.bim | awk '{OFS=" "}{print $4,"-"}' > $1/CHR$2/TO_FLIP/chr$2_neg.strand

#Join the stranded lists
cat $1/CHR$2/TO_FLIP/chr$2_pos.strand $1/CHR$2/TO_FLIP/chr$2_neg.strand > $1/CHR$2/IMPUTE_INPUT/chr$2.strand

#this can be skipped,try to let impute adjust for strand issue!
#if [ $# -eq 5 ]
#then
#IF YOU USE THE GENOTYPES (NO PREPHASING) YOU NEED TO DO ALSO THOSE STEPS
#Use the generated flipping list with plink and create the final file for impute

#	plink --noweb --bfile $1/CHR$2/TO_FLIP/chr$2.to_flip --flip $1/CHR$2/TO_FLIP/chr$2_to_flip.list --recode --out $1/CHR$2/IMPUTE_INPUT/chr$2_to_convert

#else
	plink --noweb --bfile $1/CHR$2/TO_FLIP/chr$2.to_flip --recode --out $1/CHR$2/IMPUTE_INPUT/chr$2_to_convert

#fi

echo -e "FLIPPING PROCEDURE FINISHED!!\n\nSTARTING IMPUTATION PROCEDURE!!!!!\n********************************"
sleep 10



