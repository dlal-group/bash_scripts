t%INFO/DP4\t%INFO/DP\t%INFO/HOB\t%INFO/ICB\t%INFO/IDV


bcftools annotate -x INFO/IMF,INFO/INDEL,INFO/MQ0F,INFO/MQSB,INFO/SGB,INFO/VQSLOD,INFO/culprit,INFO/CSQ,INFO/AN,INFO/AC,INFO/REF_PANEL
bcftools annotate -x "^INFO/DP4,INFO/DP,INFO/HOB,INFO/ICB,INFO/IDV",FORMAT | bcftools view -G

bcftools annotate -x "^INFO/DP4,INFO/DP,INFO/HOB,INFO/ICB,INFO/IDV",FORMAT -r 4:21261508-21261508 4.vcf.gz.indel_REF.vcf.gz| bcftools view -G -H| cut -f -2,4-5,8


 bcftools annotate -x "^INFO/DP4,INFO/DP,INFO/HOB,INFO/ICB,INFO/IDV",FORMAT -r 4:21261508-21261508 4.vcf.gz.indel_REF.vcf.gz| bcftools view -G -H| cut -f -2,4-5,8 


cut -f -2 test_last_chr4_REF_keep.tab| uniq -c|tr " " "\t"| awk 'BEGIN{OFS="\t"}{print $3,$1-1}' > test_4.vcf.gz.indel_REF_keep_dupe_pos.tab
cut -f -2,5- test_last_chr4_REF_keep.tab | uniq -c |tr " " "\t"| awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1-1}' > test_4.vcf.gz.indel_REF_keep_dupe_sites.tab

4.vcf.gz.indel_REF_keep_dupe_pos.tab
21261508        3
43167342        3
43369293        3
45090421        3
55378598        4
66659244        3
68661342        3
83728802        3
155137787       3
162616576       3
168255473       3
172145410       3


4.vcf.gz.indel_REF_keep_dupe_sites.tab

4.vcf.gz.indel_clean_REF.vcf.gz
4.vcf.gz.snp_REF_keep.tab.CARL.snp.to_keep.tab
nn
4.vcf.gz.indel_REF_keep.tab

${outdir}/PANEL/${filename}.${mode}_REF_keep.tab

/nfs/users/nfs_m/mc14/Work/bash_scripts/dupe_extract.sh ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab ${outdir}/PANEL/${filename}.${mode}_REF_keep_dupe_pos.tab ${outdir}/PANEL/${filename}.${mode}_REF_keep_dupe_sites.tab
python /nfs/users/nfs_m/mc14/Work/bash_scripts/panel_check.py ${cohort} ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab ${outdir}/PANEL/${filename}.${mode}_REF_keep_dupe_pos.tab ${outdir}/PANEL/${filename}.${mode}_REF_keep_dupe_sites.tab ${outdir}/PANEL/ ${mode}

${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab.gz

${outdir}/PANEL/${filename}.${mode}_KEEP_REF.vcf.gz

${outdir}/PANEL/${filename}.${mode}_clean_REF.vcf.gz
/nfs/users/nfs_m/mc14/Work/bash_scripts/last_vcf.sh ${outdir}/PANEL/${filename}.${mode}_clean_REF.vcf.gz ${cohort} ${outdir}/PANEL/${filename}.${mode}_REF_keep.tab.${cohort}.${mode}.to_keep.tab ${outdir}/PANEL ${mode}


clean_snp_vcf=`tabix ${pop}/PANEL/${i}.vcf.gz.snp_clean_REF.vcf.gz ${i} | wc -l`
final_clean_indels_vcf=`tabix ${pop}/PANEL/${i}.vcf.gz.indel_REF_keep.tab.${cohort}.indel.to_keep.tab.${cohort}.indel.to_keep.vcf.gz ${i}| wc -l`


bcftools view -i "INFO/GERP > -2" -s XX773282 -r 22:42580333-42872043 /lustre/scratch113/projects/esgi-vbseq/25082015_purging/26082015_ANNOTATED/GERP/CLASS/22.GERP.CLASS.vcf.gz 

for i in {1..22}
do
echo "/lustre/scratch113/projects/esgi-vbseq/09022016_PANEL_SESOURCES/TGP3/TSI/TRIMMED/29022016_NORMALIZED/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/CARL/PANEL/MERGED/${i}.vcf.gz.ALL_REF.vcf.gz /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/FVG/PANEL/MERGED/${i}.vcf.gz.ALL_REF.vcf.gz /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/VBI/PANEL/MERGED/${i}.vcf.gz.ALL_REF.vcf.gz"
done > pop_vcf_files.list

gen_map=/lustre/scratch114/resources/imputation/impute2/2015-05-08/ALL_1000G_phase1interim_jun2011_impute/genetic_map_chr1_combined_b37.txt
hap_1=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/TSI/1/1.INGI_REF.TSI.hap.gz
hap_2=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL/1/1.INGI_REF.CARL.hap.gz
leg_1=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/TSI/1/1.INGI_REF.TSI.legend.gz
leg_2=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL/1/1.INGI_REF.CARL.legend.gz
start_pos=10177
chunk_size=3000000
end_pos=$[start_pos + chunk_size]
out_ref=/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/TSI_CARL/1/1.INGI_REF.TSI_CARL.${start_pos}-${end_pos}
buffer=500


/nfs/team151/software/impute_v2.3.2_x86_64_static/impute2 -allow-large-regions -m ${gen_map} -h ${hap_1} ${hap_2} -l ${leg_1} ${leg_2} -merge_ref_panels -merge_ref_panel_output_ref ${out_ref} -int ${start_pos} ${end_pos} -Ne 20000 -buffer ${buffer}

chr=2
bcftools convert -H ${chr}/${chr}.INGI_REF.CARL.hap.gz,${chr}/${chr}.INGI_REF.CARL.legend.gz,/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL/CARL.REF.samples -O z -o ${chr}/${chr}.INGI_REF.CARL.vcf.gz
bcftools convert -H ${chr}/${chr}.INGI_REF.TSI_CARL.hap.gz,${chr}/${chr}.INGI_REF.TSI_CARL.reformatted.legend.gz,/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/TSI_CARL/TSI_CARL.REF.samples -O z -o ${chr}/${chr}.INGI_REF.TSI_CARL.vcf.gz
bcftools convert -H ${chr}/${chr}.INGI_REF.TSI_CARL.hap.gz,${chr}/${chr}.INGI_REF.TSI_CARL.reformatted.legend.gz,/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/TSI_CARL/TSI_CARL.REF.samples -O z -o ${chr}/${chr}.INGI_REF.TSI_CARL.vcf.gz
bcftools convert -H ${chr}/${chr}.INGI_REF.TSI_CARL_FVG.hap.gz,${chr}/${chr}.INGI_REF.TSI_CARL_FVG.reformatted.legend.gz,/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/TSI_CARL_FVG/TSI_CARL_FVG.REF.samples -O z -o ${chr}/${chr}.INGI_REF.TSI_CARL_FVG.vcf.gz
bcftools convert -H ${chr}/${chr}.INGI_REF.TSI_CARL_FVG_VBI.hap.gz,${chr}/${chr}.INGI_REF.TSI_CARL_FVG_VBI.legend.gz,/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/TSI_CARL_FVG_VBI/TSI_CARL_FVG_VBI.REF.samples -O z -o ${chr}/${chr}.INGI_REF.TSI_CARL_FVG_VBI.vcf.gz
bcftools convert -H ${chr}/${chr}.INGI_REF.CARL_FVG_VBI_TSI.hap.gz,${chr}/${chr}.INGI_REF.CARL_FVG_VBI_TSI.legend.gz,/lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TSI/CARL_FVG_VBI_TSI.REF.samples -O z -o ${chr}/${chr}.INGI_REF.CARL_FVG_VBI_TSI.vcf.gz

##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele. Format: AA|REF|ALT|IndelType. AA: Ancestral allele, REF:Reference Allele, ALT:Altern
EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF,AA

'3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'


/software/hgi/pkglocal/ensembl-vep-release-83-GRCh37/bin/variant_effect_predictor.pl --quiet --regulatory --sift b --polyphen b --plugin Condel,/software/vertres/bin-external/VEP_plugins/config/Condel/config/,b --symbol --format vcf --force_overwrite --cache --dir /nfs/team151/software/VEP_83_GRCh37 -o /lustre/scratch113/projects/esgi-vbseq/02032016_INGI_REF_PANEL/IMPUTE/CARL_FVG_VBI_TSI/30032016_BEAUTIFY/22/47654436-48083118.vep.part


bcftools annotate -c CHROM,POS,ID,REF,ALT,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF,AA -a /lustre/scratch114/resources/1000g/release/20130502/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz 


/lustre/scratch105/vrpipe/refs/human/ncbi37/imputation-ref-panels/uk10k+1000g-phase3/impute2/

chr=22
chunkStr=11
bsub -J rerun.chr$chr.$chunkStr -q long -R "select[mem>12000] rusage[mem=12000]" -M12000 -o chr$chr.$chunkStr.log -e chr$chr.$chunkStr.err < chr$chr.$chunkStr.cmd

chr20.20.log
chr21.12.log
chr22.11.log

384236	30669	353567	265836	184208

SNPID   RSID    chromosome      position        A_allele        B_allele        minor_allele    major_allele    AA      AB      BB      AA_calls    AB_calls        BB_calls        MAF     HWE     missing missing_calls   information


for i in 1 2 20 21 22
do
ls ${i}/chr${i}*.gz |wc -l
ls ${i}/chr${i}*.cmd |wc -l
done

for pop in CARL FVG INCIPE2 VBI
do
pop=INCIPE2
queue=long
mem=30000
chr=2
PANEL=CARL_FVG_VBI_TGP3_ALL
postfix=".shapeit"
chunkStr=28
imputedir=/lustre/scratch113/projects/esgi-vbseq/31032016_IMPUTATION/${pop}/${PANEL}$postfix/${chr}
cd $imputedir
bsub -J ${PANEL}$postfix.chr$chr.$chunkStr -q $queue -R "select[mem>${mem}] rusage[mem=${mem}]" -M ${mem} -o $imputedir/chr$chr.$chunkStr.log -e $imputedir/chr$chr.$chunkStr.err < $imputedir/chr$chr.$chunkStr.cmd
done


for pop in carl fv vbi incipe2
do
for pan in CARL_FVG_VBI.shapeit CARL_FVG_VBI_TSI.shapeit CARL_FVG_VBI_TGP3_ALL.shapeit uk10k1kg.ref TGP3_ALL.shapeit EUR.shapeit CARL.shapeit FVG.shapeit VBI.shapeit
do
for chr in 2
do
echo $pop $pan $chr

bsub -J "merge_$pop_$pan_$chr" -o "LOGS/%J_merge_$pop_$pan_$chr.o" -M 2000 -R"select[mem>2000] rusage[mem=2000]" -q normal -- ~/Work/bash_scripts/INGI_ref_checks.sh $pop $pan $chr

done
done
done

for pop in CARL FVG VBI
do
# for panel in CARL_FVG_VBI.shapeit CARL_FVG_VBI_TGP3_ALL.shapeit CARL_FVG_VBI_TSI.shapeit EUR.shapeit TGP3_ALL.shapeit
for panel in CARL_FVG_VBI_TGP3_ALL.shapeit
do
for chr in {1..22}
do
genz=`ls ${pop}/${panel}/${chr}/*.gen.gz|wc -l | cut -f 1 -d " "`
cmdz=`ls ${pop}/${panel}/${chr}/*.cmd|wc -l | cut -f 1 -d " "`
if [ $genz -eq $cmdz ];then
echo "$pop $panel $chr $genz $cmdz OK"
else
echo "$pop $panel $chr $genz $cmdz WARNING"
fi
done
done
done > 07072016_impute_check.txt


for panel in CARL_FVG_VBI_UK10K_TGP3_ALL
do
for chr in $(seq 22)
do
# legs=$(ls ${chr}.INGI_REF.${panel}.*.legend.gz)
# last=${leg[-1]}
# IFS="." read -ra ZNAME <<< "$last"
# last=${ZNAME[-3]}
last=`wc -l ${panel}/${chr}/${chr}.chunks.txt| cut -f 1 -d " "`
legendz=`ls ${panel}/${chr}/${chr}.INGI_REF.${panel}.*.legend.gz|wc -l | cut -f 1 -d " "`
hapz=`ls ${panel}/${chr}/${chr}.INGI_REF.${panel}.*.hap.gz|wc -l | cut -f 1 -d " "`
if [ $legendz -eq $last ];then
echo "$pop $panel $chr $legendz $hapz $last OK"
else
echo "$pop $panel $chr $legendz $hapz $last WARNING"
fi
done
done

#convert to VCF files
for pop in CARL FVG VBI
do
for chr in {1..22}
do
mkdir -p /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${chr}/vcf/
plink --bfile /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${chr}/chr${chr} --recode vcf-iid --out /netapp/dati/WGS_REF_PANEL/genotypes/${pop}/merged/cleaned/${chr}/vcf/chr${chr}

done
done

