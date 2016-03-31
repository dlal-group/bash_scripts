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
# refname=uk10k1kg ## 1kg, uk10k, uk10k1kg
postfix=".shapeit" ## "" or ".shapeit"
by_chunk=N  ## "Y" or "N"
# phasedir=$scratch113/imputed/$geno/shapeit

impute2=/nfs/team151/software/impute_v2.3.1_x86_64_static/impute2
shapeit2=/nfs/team151/software/shapeit.v2.r790/shapeit
plink2=/nfs/team151/software/plink2_18_April_2015/plink
chunk_size=3000000
buffer_size=250
window_size=2
thread=8
extra_str="-verbose" #"-verbose" #"-phase"
pop=$1
PANEL=$2
chr=$3
MODE=$4 #set this to PHASE, if you want to phase and impute; set this to IMPUTE, if you're providing already phased genotypes

# imputedir=/lustre/scratch113/projects/carl_seq/05272015_MERGED_REF_PANEL/IMPUTED/${pop}/${PANEL}$postfix
# imputedir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/14102015_MERGED_REF_PANEL/IMPUTED/${pop}/${PANEL}$postfix
imputedir=/lustre/scratch113/projects/esgi-vbseq/31032016_IMPUTATION/${pop}/${PANEL}$postfix
# mkdir -p ${phasedir}
mkdir -p ${imputedir}

case $pop in
	VBI)
	extra_str_excl_samples="-exclude_samples_g /lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/LISTS/VBI_impute_exclude_sample.list"
	# extra_str_excl_snps="-exclude_snps_g /lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/snplist/${pop}_chr${chr}.exclude -impute_excluded"
	genodir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/MERGED_REF_PANEL_Feb2015/VBI_geno
	phasedir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/MERGED_REF_PANEL_Feb2015/VBI_geno
	;;
	FVG)
	extra_str_excl_samples="-exclude_samples_g /lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL/LISTS/FVG_impute_exclude_sample.list"
	# extra_str_excl_snps="-exclude_snps_g /lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/snplist/${pop}_chr${chr}.exclude -impute_excluded"
	genodir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/GWAS/FVG/shapeit
	phasedir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/GWAS/FVG/shapeit
	;;
	CARL)
	extra_str_excl_samples="-exclude_samples_g /lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL/LISTS/CARL_impute_exclude_sample.list"
	# extra_str_excl_snps="-exclude_snps_g /lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/snplist/${pop}_chr${chr}.exclude -impute_excluded"
	genodir=/lustre/scratch113/teams/soranzo/users/jh21/imputed/carl/shapeit
	phasedir=/lustre/scratch113/teams/soranzo/users/jh21/imputed/carl/shapeit
	;;
	INCIPE2 )
	# extra_str_excl_snps="-exclude_snps_g /lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/snplist/${pop}_chr${chr}.exclude -impute_excluded"
	genodir=/lustre/scratch113/teams/soranzo/users/jh21/imputed/incipe2/shapeit
	phasedir=/lustre/scratch113/teams/soranzo/users/jh21/imputed/incipe2/shapeit
	;;
esac

refdir=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/${PANEL}
refhap=$refdir/${chr}/${chr}.INGI_REF.${PANEL}.hap.gz
reflegend=$refdir/${chr}/${chr}.INGI_REF.${PANEL}.legend.gz

# case $PANEL in
# 	FVG )
# 	refdir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/MERGED_REF_PANEL_Feb2015/SHAPEIT/FVG_HAP_LEGEND
# 	refhap=$refdir/${PANEL}.chr$chr.hap.gz
# 	reflegend=$refdir/${PANEL}.chr$chr.legend.gz
# 	;;
# 	CARL )
# 	refdir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/MERGED_REF_PANEL_Feb2015/SHAPEIT/CARL_HAP_LEGEND
# 	refhap=$refdir/${PANEL}.chr$chr.hap.gz
# 	reflegend=$refdir/${PANEL}.chr$chr.legend.gz
# 	;;
# 	VBI )
# 	refdir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/MERGED_REF_PANEL_Feb2015/SHAPEIT/VBI_HAP_LEGEND
# 	refhap=$refdir/${PANEL}.chr$chr.hap.gz
# 	reflegend=$refdir/${PANEL}.chr$chr.legend.gz
# 	;;
# 	INGI )
# 	refdir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/MERGED_REF_PANEL_Feb2015/MERGER_CARL_VBI_FVG
# 	refhap=$refdir/${PANEL}.chr$chr.hap.gz
# 	reflegend=$refdir/${PANEL}.chr$chr.legend.gz
# 	;;
# 	1000Gph1 )
# 	refdir=/lustre/scratch114/resources/imputation/impute2/2015-05-08/ALL_1000G_phase1interim_jun2011_impute
# 	refhap=$refdir/ALL_1000G_phase1interim_jun2011_chr${chr}_impute.hap.gz
# 	reflegend=$refdir/ALL_1000G_phase1interim_jun2011_chr${chr}_impute.legend.gz
# 	;;
# 	INGI_1000GPh3 )
# 	refdir=/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/MERGED_REF_PANEL_Feb2015/MERGER_INGI_1000GPh3
# 	refhap=$refdir/${PANEL}.chr$chr.hap.gz
# 	reflegend=$refdir/${PANEL}.chr$chr.legend.gz
# 	;;
# 	1000GP_Phase3)
# 	refdir=/lustre/scratch114/resources/imputation/impute2/1000GP_Phase3_v1a
# 	refhap=$refdir/${PANEL}_chr$chr.hap.gz
# 	reflegend=$refdir/${PANEL}_chr$chr.legend.gz
# 	;;
# 	INGI_1000GPh3_UK10K )
# 	refdir=/lustre/scratch114/resources/imputation/impute2/1000GP_Phase3_v1a
# 	refhap=$refdir/${PANEL}_chr$chr.hap.gz
# 	reflegend=$refdir/${PANEL}_chr$chr.legend.gz
# 	;;
# esac

extra_str=`echo $extra_str_excl_snps $extra_str_excl_samples`

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
if [ $MODE=="PHASE" ]; then
	echo phase $geno chr$chr
	# echo -e "#!/usr/local/bin/bash
	# \necho \"Starting on : \$(date); Running on : \$(hostname); Job ID : \$LSB_JOBID\"
	# \n$plink2 --bfile $genodir/$geno  $plink_str --make-bed --out chr$chr\n\n
	# \n$shapeit2 --thread $thread --window $window_size --states 200 --effective-size 11418 -B chr$chr --input-map $scratch113/references_panel/1kg/genetic_map_chr${chr}_combined_b37.txt --output-log chr$chr.shapeit --output-max chr$chr.hap.gz chr$chr.sample $chrX_phase_str
	# " > $phasedir/chr$chr.cmd
	# cd $phasedir
	# bsub -J $geno.shapeit.chr$chr -q long -o chr$chr.shapeit.log -e chr$chr.shapeit.err -n$thread -R "span[ptile=$thread] select[mem>18000] rusage[mem=18000]" -M18000 < chr$chr.cmd
	# continue
fi

### step 2: impute ###
chr_begin=`zcat $reflegend | awk 'NR==2 {printf \$2}'`
chr_end=`zcat $reflegend | tail -1 | awk '{printf \$2}'`
let "chunk_num=($chr_end - $chr_begin)/$chunk_size" # bash rounds automatically
if [[ $chunk_num <1 ]]; then
	chunk_num=1
fi
for chunk in `seq 1 $chunk_num`; do
	chunkStr=`printf "%02d" $chunk`
	if [[ -e $imputedir/chr$chr.$chunkStr.log ]]; then
		continue
	fi
	if [[ $chunk -gt 1 ]]; then
		chunk_begin=`echo "$chr_begin+($chunk-1)*$chunk_size+1" | bc`
	else
		chunk_begin=$chr_begin
	fi
	if [[ $chunk -eq $chunk_num ]]; then
		mem=5000
		queue=normal
		chunk_end=$chr_end
	else
		mem=5000
		queue=normal
		chunk_end=`echo "$chr_begin+($chunk*$chunk_size)" | bc`
	fi
	if [[ $chr.$chunkStr == 8.02 ]]; then
		mem=36000
		queue=hugemem
	fi
	if [[ $by_chunk == "Y" ]]; then
		refhap=$refdir/$refname/chr$chr.${chunkStr}$postfix.hap.gz
		reflegend=$refdir/$refname/chr$chr.${chunkStr}$postfix.legend.gz
	fi
	gen_map=/lustre/scratch114/resources/imputation/impute2/2015-05-08/ALL_1000G_phase1interim_jun2011_impute/genetic_map_chr${chr}_combined_b37.txt
	echo -e "#!/usr/local/bin/bash
	\n$impute2 -allow_large_regions -m ${gen_map} -h $refhap -l $reflegend -known_haps_g $phasedir/chr$chr.hap.gz -sample_g $phasedir/chr$chr.sample $extra_str -use_prephased_g -k_hap $k_hap -int $chunk_begin $chunk_end -Ne 20000 -buffer $buffer_size -o chr$chr.$chunkStr.gen $chrX_impute_str
	\ngzip -f chr$chr.$chunkStr.gen
	\nif [[ -e chr$chr.$chunkStr.gen_allele_probs ]]; then
	\ngzip chr$chr.$chunkStr.gen_allele_probs chr$chr.$chunkStr.gen_haps
	\nfi
	\nN_info=\`awk 'NR>1' chr$chr.$chunkStr.gen_info | wc -l | awk '{printf \$1}'\`
	\nN_gen=\`zcat chr$chr.$chunkStr.gen.gz | wc -l | awk '{printf \$1}'\`
            if [[ \$N_info != \$N_gen ]]; then
                    echo \"chr$chr $chunkStr: \$N_info for info, \$N_gen for gen\" > $imputedir/chr$chr.$chunkStr.ERR
            fi
            " > $imputedir/chr$chr.$chunkStr.cmd
	cd $imputedir
	bsub -J ${refname}$postfix.$geno.chr$chr.$chunkStr -q $queue -R "select[mem>$mem] rusage[mem=$mem]" -M${mem} -o $imputedir/chr$chr.$chunkStr.log -e $imputedir/chr$chr.$chunkStr.err < $imputedir/chr$chr.$chunkStr.cmd
done
