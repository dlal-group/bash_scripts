#!/usr/local/bin/bash
 
#general script to launch fast jobs


mkdir -p LOGS/
for i in {1..6} {8..22} X
# for i in $@
do
#bsub -J "chr${i}_generic_manipulation" -o "%J_chr${i}_generic_manipulation.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
#-q normal general_script.sh ${i}
#tabix  esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.chr${i}.re_ann.NOT_OVERLAP.NO_RSID.vcf.gz chr ${i} >> esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.NO_RSID.vcf
#cat ~/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/VBI/VBI.vcf.header esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.vcf | bgzip -c > esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.vcf.gz
#echo "CHR ${i}"
#tabix esgi-vbseq.vqsr.beagle.impute2.anno.20120607.csq.SNPS.re_ann.NOT_OVERLAP.NO_RSID.vcf.gz chr ${i} | wc -l
# bsub -J "generic_manipulation_chr${i}" -o "%J_generic_manipulation_chr${i}.log" -e "%J_generic_manipulation_chr${i}.err" -M3000000 -R"select[mem>3000] rusage[mem=3000]" \
# -q yesterday general_script.sh ${i}
echo "Element: ${i}"
# bsub -J "generic_manipulation_${i}" -o "LOGS/%J_generic_manipulation_${i}.log" -e "LOGS/%J_generic_manipulation_${i}.err" -M1000 -R"select[mem>1000] rusage[mem=1000]" \
bsub -J "generic_manipulation_${i}" -o "LOGS/%J_generic_manipulation_${i}.log" -e "LOGS/%J_generic_manipulation_${i}.err" -n 4,8 -M1000 -R"span[hosts=1] select[mem>1000] rusage[mem=1000]" -q small -- general_script.sh ${i}
	# for file in /lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/COMBINED_IMPUTATION/VBI_1000GP/CHR${chr}/IMPUTE_INPUT/*.lsf
	# do
	# 	echo $file
	# 	bsub -M5000 -R"select[mem>5000] rusage[mem=5000]" -q normal < ${file}
	# done
done
