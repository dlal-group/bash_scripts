#Comparisons and Stats
#Comparison between INGI and TGP and UK10K and within INGI

#bcftools isec 

#echo "bcftools isec ~/carl_seq/variant_refinement/13102015_RELEASE/ALL_CARL_20151013.vcf.gz
# /lustre/scratch113/projects/esgi-vbseq/08092015/13102015_RELEASE/ALL_VBI_20151013.vcf.gz
# /lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/20140518_RELEASE/VBI_20140518.all.vcf.gz
#   -O v -n -3 -p UNION"| bsub -J"isec" -o"%J_isec.o" -M 5000 -R"select[mem>=5000] rusage[mem=5000]" -q yesterday


#Extract samples from OUTBRED pops in TGP ph3
#
#for pop in EUR TSI 
#do
##for list in `ls ${pop}/*.list`
##do
#echo ${list}
#for chr in {1..22}
#do
##bsub -J "%J_pop_extr" -o "${pop}/%J_pop_extr.o" -M4000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- bcftools view -S ${list} -c 1 /lustre/scratch114/resources/1000g/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -O z -o /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz 
#bsub -J "%J_pop_extr" -o "${pop}/%J_pop_extr.o" -M4000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- tabix -f -p vcf /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz 
#done
##done
#done


#Intersection between isolates and TGP pops, single pop comparison

for ipop in CARL VBI FVG
do
case $ipop in
CARL)
pop_path=/lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL
;;
VBI)
pop_path=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL
;;
FVG)
pop_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL
;;
esac

for pop in EUR TSI 
do
echo ${ipop} ${pop}
for chr in {1..22}
do
mkdir -p /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${ipop}/${pop}/UNION
bsub -J "%J_${ipop}_${pop}_isec" -o "${ipop}/%J_${ipop}_${pop}_isec.o" -M4000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- bcftools isec /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz ${pop_path}/${chr}.vcf.gz -O z -n -2 -p /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${ipop}/${pop}/UNION/${chr}
done
done
done


#Intersection between isolates and TGP pops, all isolates comparison
#plus intersection with Uk10K + TGP3
#uk10k_tgp3=/lustre/scratch116/vr/ref/human/GRCh37/imputation/uk10k+1000g-phase3/vcf


for pop in EUR TSI
for pop in UK10K
for pop in TGPph3
for pop in UK10K_TGPph3
do
echo ${pop}
pop_path_carl=/lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL
pop_path_vbi=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL
pop_path_fvg=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL
for chr in {1..22}
do
# chr=10
# pop=EUR
# tgp3_file=/lustre/scratch114/resources/1000g/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
uk10k_tgp3=/lustre/scratch116/vr/ref/human/GRCh37/imputation/uk10k+1000g-phase3/vcf/${chr}.bcf
# tgp3_file=/lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
#uk10k_file=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/vcf_sites_filtered/chr${chr}.sites.vcf.gz
mkdir -p /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION
#bsub -J "%J_multi_${pop}_isec" -o "${pop}/%J_multi_${pop}_isec.o" -M4000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- bcftools isec /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/${pop}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz ${pop_path_carl}/${chr}.vcf.gz ${pop_path_vbi}/${chr}.vcf.gz ${pop_path_fvg}/${chr}.vcf.gz -O z -n -4 -p /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/${pop}/UNION/${chr}
#bsub -J "%J_multi_${pop}_isec" -o "${pop}/%J_multi_${pop}_isec.o" -M4000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- bcftools isec ${uk10k_file} ${pop_path_carl}/${chr}.vcf.gz ${pop_path_vbi}/${chr}.vcf.gz ${pop_path_fvg}/${chr}.vcf.gz -O z -n -4 -p /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}
# bsub -J "%J_multi_${pop}_isec" -o "${pop}/%J_multi_${pop}_isec.o" -M4000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- bcftools isec ${tgp3_file} ${pop_path_carl}/${chr}.vcf.gz ${pop_path_vbi}/${chr}.vcf.gz ${pop_path_fvg}/${chr}.vcf.gz -O z -n =4 -p /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}_test_eq
# bsub -J "%J_multi_${pop}_isec" -o "${pop}/%J_multi_${pop}_isec.o" -M4000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- bcftools isec ${tgp3_file} ${pop_path_carl}/${chr}.vcf.gz ${pop_path_vbi}/${chr}.vcf.gz ${pop_path_fvg}/${chr}.vcf.gz -O z -n +4 -p /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}_test_gt
bsub -J "%J_multi_${pop}_isec_${chr}" -o "${pop}/%J_multi_${pop}_isec_${chr}.o" -M4000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- bcftools isec ${uk10k_tgp3} ${pop_path_carl}/${chr}.vcf.gz ${pop_path_vbi}/${chr}.vcf.gz ${pop_path_fvg}/${chr}.vcf.gz -O z -n -4 -p /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}
done
done


#Comparisons between INGI

echo ${pop}
for chr in {1..22}
do
pop_path_carl=/lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL
pop_path_vbi=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL
pop_path_fvg=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL
mkdir -p /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION
bsub -J "%J_multi_INGI_isec" -o "INGI/%J_multi_INGI_isec.o" -M4000 -R"select[mem>4000] rusage[mem=4000]" -q normal -- bcftools isec ${pop_path_carl}/${chr}.vcf.gz ${pop_path_vbi}/${chr}.vcf.gz ${pop_path_fvg}/${chr}.vcf.gz -O z -n -3 -p /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION/${chr}
done


#######count of overlapping sites with TGP phase 3

#Extract for each chr , how many sites we have in common between all populations and all sites belonging only to the isolates (in shared between 1,2 or 3 pops)
#Population order:
#TGP TSI EUR UK10K
#CARL CARL CARL CARL
#VBI VBI VBI VBI
#FVG FVG FVG FVG

#union_path=/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION
#union_path=/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/UK10K/UNION
#union_path=/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/TGPph3/UNION
#union_path=/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/EUR/UNION
#union_path=/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/TSI/UNION

for pop in EUR TSI UK10K
for pop in TGPph3
for pop in EUR TSI UK10K TGPph3
do
union_path=/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION
echo ${pop}
(echo "chr all_common tgp tgp_vb tgp_fvg tgp_carl tgp_vb_fvg tgp_vb_carl tgp_fvg_carl only_ingi only_vb vb_fvg vb_carl only_fvg fvg_carl only_carl";
(for chr in {1..22}
do

all_common=`awk '$(NF)=="1111"' ${union_path}/${chr}/sites.txt | wc -l`

out=`awk '$(NF)=="1000"' ${union_path}/${chr}/sites.txt | wc -l`
out_vb=`awk '$(NF)=="1010"' ${union_path}/${chr}/sites.txt | wc -l`
out_fvg=`awk '$(NF)=="1001"' ${union_path}/${chr}/sites.txt | wc -l`
out_carl=`awk '$(NF)=="1100"' ${union_path}/${chr}/sites.txt | wc -l`
out_vb_fvg=`awk '$(NF)=="1011"' ${union_path}/${chr}/sites.txt | wc -l`
out_vb_carl=`awk '$(NF)=="1110"' ${union_path}/${chr}/sites.txt | wc -l`
out_fvg_carl=`awk '$(NF)=="1101"' ${union_path}/${chr}/sites.txt | wc -l`


only_ingi=`awk '$(NF)=="0111"' ${union_path}/${chr}/sites.txt | wc -l`
vb=`awk '$(NF)=="0010"' ${union_path}/${chr}/sites.txt | wc -l`
vb_fvg=`awk '$(NF)=="0011"' ${union_path}/${chr}/sites.txt | wc -l`
vb_carl=`awk '$(NF)=="0110"' ${union_path}/${chr}/sites.txt | wc -l`

fvg=`awk '$(NF)=="0001"' ${union_path}/${chr}/sites.txt | wc -l`
fvg_carl=`awk '$(NF)=="0101"' ${union_path}/${chr}/sites.txt | wc -l`
carl=`awk '$(NF)=="0100"' ${union_path}/${chr}/sites.txt | wc -l`

echo "${chr} ${all_common} ${out} ${out_vb} ${out_fvg} ${out_carl} ${out_vb_fvg} ${out_vb_carl} ${out_fvg_carl} ${only_ingi} ${vb} ${vb_fvg} ${vb_carl} ${fvg} ${fvg_carl} ${carl}"
done)) > ${union_path}/${pop}_comparison.tab
done


###################comparison within INGI #############################
#Populations order:
#CARL
#VBI
#FVG

union_path=/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/INGI/UNION
(echo "chr all_common vb vb_fvg vb_carl fvg fvg_carl carl";
(for chr in {1..22}
do

all_common=`awk '$(NF)=="111"' ${union_path}/${chr}/sites.txt | wc -l`
vb=`awk '$(NF)=="010"' ${union_path}/${chr}/sites.txt | wc -l`
vb_fvg=`awk '$(NF)=="011"' ${union_path}/${chr}/sites.txt | wc -l`
vb_carl=`awk '$(NF)=="110"' ${union_path}/${chr}/sites.txt | wc -l`
fvg=`awk '$(NF)=="001"' ${union_path}/${chr}/sites.txt | wc -l`
fvg_carl=`awk '$(NF)=="101"' ${union_path}/${chr}/sites.txt | wc -l`
carl=`awk '$(NF)=="100"' ${union_path}/${chr}/sites.txt | wc -l`


echo "${chr} ${all_common} ${vb} ${vb_fvg} ${vb_carl} ${fvg} ${fvg_carl} ${carl}"
done)) > ${union_path}/INGI_comparison.tab

#####################################
#Merge all pop together to calculate ROH/IBD

	mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/EUR_INGI_MERGE/chr.list|cut -f 1 -d " "`;bsub -J "merge_chr[1-${size}]" -o "LOGS/%J_merge_chr.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q long -- ~/Work/bash_scripts/ja_runner.sh /lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/EUR_INGI_MERGE/chr.list

#######################################
#Extract AF data for overlapping sites
#extract freq for each files with overlap

for pop in EUR TSI TGPph3 UK10K
for pop in TGPph3
do

for chr in {1..22}
do
case ${pop} in 
EUR )
pop_path=/lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/EUR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
;;
TSI )
pop_path=/lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/TGP3/TSI/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
;;
TGPph3 )
pop_path=/lustre/scratch114/resources/1000g/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
;;
UK10K )
pop_path=/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/vcf_sites_filtered/chr${chr}.sites.vcf.gz
;;
esac

carl_path=/lustre/scratch113/projects/carl_seq/variant_refinement/12112015_FILTERED_REL/${chr}.vcf.gz
fvg_path=/lustre/scratch113/projects/fvg_seq/16092015/12112015_FILTERED_REL/${chr}.vcf.gz
vbi_path=/lustre/scratch113/projects/esgi-vbseq/08092015/12112015_FILTERED_REL/${chr}.vcf.gz

echo "(echo \"CHROM POS AC AN ALT_F MAF\";bcftools query -R /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}/sites.txt -f \"%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n\" ${pop_path} | awk '{if(\$3~\",\") ;else print \$0,\$3/\$4}'| awk '{if(\$5<=0.5) print \$0,\$5;else print \$0,1-\$5}')| tr \"\t\" \" \" > /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}/${pop}_freq.txt" | bsub -J "freq_calc_${pop}" -o "%J_freq_calc_${pop}.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal
echo "(echo \"CHROM POS AC AN ALT_F MAF\";bcftools query -R /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}/sites.txt -f \"%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n\" ${carl_path} | awk '{if(\$3~\",\");else print \$0,\$3/\$4}'| awk '{if(\$5<=0.5) print \$0,\$5;else print \$0,1-\$5}')| tr \"\t\" \" \" > /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}/carl_freq.txt" | bsub -J "freq_calc_${pop}_carl" -o "%J_freq_calc_${pop}_carl.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal
echo "(echo \"CHROM POS AC AN ALT_F MAF\";bcftools query -R /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}/sites.txt -f \"%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n\" ${fvg_path} | awk '{if(\$3~\",\");else print \$0,\$3/\$4}'| awk '{if(\$5<=0.5) print \$0,\$5;else print \$0,1-\$5}')| tr \"\t\" \" \" > /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}/fvg_freq.txt" | bsub -J "freq_calc_${pop}_fvg" -o "%J_freq_calc_${pop}_fvg.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal
echo "(echo \"CHROM POS AC AN ALT_F MAF\";bcftools query -R /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}/sites.txt -f \"%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n\" ${vbi_path} | awk '{if(\$3~\",\");else print \$0,\$3/\$4}'| awk '{if(\$5<=0.5) print \$0,\$5;else print \$0,1-\$5}')| tr \"\t\" \" \" > /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}/vbi_freq.txt" | bsub -J "freq_calc_${pop}_vbi" -o "%J_freq_calc_${pop}_vbi.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal

done
done


############################################################################
---->>>> RUN THIS ONE!!!!!
#It's in ja_runner.sh!!
#for pop in EUR TSI UK10K
#do
#for chr in {1..22}
#do
##pop="UK10K"
##chr=21
#pop_path="/lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/${pop}/UNION/${chr}"
#awk '{if(length($3)==length($4)) print $0}' ${pop_path}/sites.txt > ${pop_path}/sites_snp.txt
#awk 'FNR==NR{a[$2]=$6;next}{if(a[$2]) print $0,a[$2];else print $0,"NA"}' ${pop_path}/${pop}_freq.txt ${pop_path}/sites_snp.txt > ${pop_path}/sites_${pop}.txt
#awk 'FNR==NR{a[$2]=$6;next}{if(a[$2]) print $0,a[$2];else print $0,"NA"}' ${pop_path}/carl_freq.txt ${pop_path}/sites_${pop}.txt > ${pop_path}/sites_${pop}_carl.txt
#awk 'FNR==NR{a[$2]=$6;next}{if(a[$2]) print $0,a[$2];else print $0,"NA"}' ${pop_path}/vbi_freq.txt ${pop_path}/sites_${pop}_carl.txt > ${pop_path}/sites_${pop}_carl_vbi.txt
#awk 'FNR==NR{a[$2]=$6;next}{if(a[$2]) print $0,a[$2];else print $0,"NA"}' ${pop_path}/fvg_freq.txt ${pop_path}/sites_${pop}_carl_vbi.txt > ${pop_path}/sites_${pop}_carl_vbi_fvg.txt
#
#done
#done

mkdir -p LOGS;size=`wc -l /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/chr.list|cut -f 1 -d " "`;bsub -J "merge_chr[1-${size}]" -o "LOGS/%J_merge_chr.%I.o" -M 5000 -R"select[mem>5000] rusage[mem=5000]" -q normal -- ~/Work/bash_scripts/ja_runner.sh /lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/chr.list


##################################################################################
#Extract data to use for maf plot

# for pop in EUR TSI UK10K TGPph3
for pop in TGPph3
do
    for chr in {1..22}
    do
    pop_path="/lustre/scratch113/projects/esgi-vbseq/16112015_TRIESTE/${pop}/UNION/${chr}"
        for ipop in fvg carl vbi
        do
        echo ${chr} ${pop} ${ipop}
            case ${ipop} in
            carl)
            awk '$7!="NA"' ${pop_path}/sites_${pop}_carl_vbi_fvg.txt |  tr "\t" " "| cut -f -6,7 -d " " | sort -u -t ' ' -g -k2,2 > ${pop_path}/sites_${pop}_carl_only.txt
            ;;
            vbi)
            awk '$8!="NA"' ${pop_path}/sites_${pop}_carl_vbi_fvg.txt | tr "\t" " "| cut -f -6,8 -d " " | sort -u -t ' ' -g -k2,2 > ${pop_path}/sites_${pop}_vbi_only.txt
            ;;
            fvg)
            awk '$9!="NA"' ${pop_path}/sites_${pop}_carl_vbi_fvg.txt | tr "\t" " "| cut -f -6,9 -d " " | sort -u -t ' ' -g -k2,2 > ${pop_path}/sites_${pop}_fvg_only.txt
            ;;
            esac
        done
    done
done
