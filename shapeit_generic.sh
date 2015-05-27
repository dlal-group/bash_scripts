#!/usr/local/bin/bash
#
#################################################
#	CUSTOMIZED FOR CHRX imputation by max!!!!!		#
#################################################
# For CGI:
#	1. remove extra_str; remove PRE-PHASING file check
#	2. -g $scratch113/imputed/$geno/wgs/chr$chr.gwas.gen.gz -k 100

## from Shane: /lustre/scratch106/projects/uk10k/RELEASE/UK10K_COHORT/REL-2012-06-02/v3/

k_hap=10000 #numbers of haplotypes used by Impute
geno=$1
genodir=$2
refname=uk10k1kg ## 1kg, uk10k, uk10k1kg
postfix=".shapeit" ## "" or ".shapeit"
by_chunk=Y  ## "Y" or "N"
scratch113=/lustre/scratch113/projects/carl_seq/05272015_MERGED_REF_PANEL/GWAS
# refdir=$scratch113/references_panel
phasedir=$3
# phasedir=/lustre/scratch113/projects/carl_seq/05272015_MERGED_REF_PANEL/GWAS/FVG/shapeit
# imputedir=$scratch113/imputed/$geno/${refname}$postfix
mkdir -p ${phasedir}
# mkdir -p ${imputedir}

impute2=/nfs/team151/software/impute_v2.3.1_x86_64_static/impute
shapeit2=/nfs/team151/software/shapeit.v2.r790/shapeit
chunk_size=3000000
buffer_size=250
window_size=2
thread=8
extra_str="-verbose" #"-verbose" #"-phase"
# extra_str="-exclude_samples_g /nfs/team151/jh21/data/uk10kgwas/uk10kgwas.1000.excluded -sample_h $refdir/$refname/$refname.sample_h.sample -exclude_samples_h /nfs/team151/jh21/data/uk10kgwas/uk10kgwas.1000.included"
# extra_str="-exclude_samples_g $scratch113/references_panel/uk10k/uk10k.sample.ids"

# for chr in {1..22}; do
# for chr in X_PAR1 X_PAR2 X ; do
for chr in {1..22} X_PAR1 X_PAR2 X; do
	if [[ "$chr" == "X_PAR1" ]]; then # (60,001 - 2,699,520)
		plink_str="--chr X --from-bp 60001 --to-bp 2699520"
		chrX_phase_str=""
		chrX_impute_str="-chrX -Xpar"
	elif [[ "$chr" == "X" ]]; then
		plink_str="--chr X --from-bp 2699521 --to-bp 154931043"
		chrX_phase_str="--chrX" 
		chrX_impute_str="-chrX" ## -sample_g already specified elsewhere
	elif [[ "$chr" == "X_PAR2" ]]; then # (154,931,044-155,260,560)
		plink_str="--chr X --from-bp 154931044 --to-bp 999999999"
		chrX_phase_str=""
		chrX_impute_str="-chrX -Xpar"
	else
		plink_str="--chr $chr"
		chrX_phase_str=""
		chrX_impute_str=""
	fi

	### step 1: pre-phase ###
	if [[ ! -e $phasedir/chr$chr.hap.gz ]]; then
		echo phase $geno chr$chr
		echo -e "#!/usr/local/bin/bash
		\necho \"Starting on : \$(date); Running on : \$(hostname); Job ID : \$LSB_JOBID\"
		\nplink --noweb --bfile $genodir/$geno  $plink_str --make-bed --out chr$chr\n\n
		\n$shapeit2 --thread $thread --window $window_size --states 200 --effective-size 11418 -B chr$chr --input-map /lustre/scratch114/resources/imputation/impute2/1000GP_Phase3_v1a/genetic_map_chr${chr}_combined_b37.txt --output-log chr$chr.shapeit --output-max chr$chr.hap.gz chr$chr.sample $chrX_phase_str
		" > $phasedir/chr$chr.cmd
		cd $phasedir
		bsub -J $geno.shapeit.chr$chr -q long -o chr$chr.shapeit.log -e chr$chr.shapeit.err -n$thread -R "span[ptile=$thread] select[mem>7000] rusage[mem=7000]" -M7000 < chr$chr.cmd
		continue
	fi

done
