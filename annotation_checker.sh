#!/usr/local/bin/bash
#args:
#$1=input file
#$2=output path
#$3=population name
chr=$1
in_path=$2
out_path=$3
pop=$4

mkdir -p ${out_path}

(echo "CHR #SNPs_w_ann #SNPs_wo_ann #SNPs_SINGLETONs #SNPs_TGP_overlap #SNPs_UK10K_overlap #SNP_private #SNP_private_no_sing #INDELs_w_ann #INDELs_wo_ann #INDELs_SINGLETONs #INDELs_TGP_overlap #INDELs_UK10K_overlap #INDELs_private #INDELs_private_no_sing";
for i in ${chr}
do

	#extract some numbers for general checks
	#Whole variant set
	#SNPs
	#Novel variants
	snp_wo_anno=`tabix ${in_path} ${i} | fgrep -v "INDEL" | awk '$3=="."' | wc -l`
	snp_w_anno=`tabix ${in_path} ${i} | fgrep -v "INDEL" | awk '$3!="."' | wc -l`
	#singletons
	snp_singleton1=`tabix ${in_path} ${i} | fgrep -v "INDEL" | fgrep "SINGLETON" | wc -l`
	snp_singleton2=`tabix ${in_path} ${i} | fgrep -v "INDEL" | fgrep -v "SINGLETON"| fgrep "AC=499" | wc -l`
	#overlap with TGP and uk10k
	snp_tgp=`tabix ${in_path} ${i} | fgrep -v "INDEL" | fgrep "TGP" | wc -l`
	snp_uk10k=`tabix ${in_path} ${i} | fgrep -v "INDEL" | fgrep "UK10K" | wc -l`
	#Putative novel (no rsID and no in tgp or uk10k)
	snp_novel=`tabix ${in_path} ${i} | fgrep -v "INDEL" | awk '$3=="."' | fgrep -v "UK10K" | fgrep -v "TGP" | wc -l`
	snp_novel_no_singleton=`tabix ${in_path} ${i} | fgrep -v "INDEL" | awk '$3=="."' | fgrep -v "UK10K" | fgrep -v "TGP" | fgrep -v "SINGLETON" | wc -l`


	#INDELs
	#Novel variants
	indels_wo_anno=`tabix ${in_path} ${i} | fgrep "INDEL" | awk '$3=="."' | wc -l`
	indels_w_anno=`tabix ${in_path} ${i} | fgrep "INDEL" | awk '$3!="."' | wc -l`
	#singletons
	indels_singletons1=`tabix ${in_path} ${i} | fgrep "INDEL" | fgrep "SINGLETON" | wc -l`
	indels_singletons2=`tabix ${in_path} ${i} | fgrep "INDEL" | fgrep -v "SINGLETON"| fgrep "AC=499" | wc -l`
	#overlap with TGP and uk10k
	indels_tgp=`tabix ${in_path} ${i} | fgrep "INDEL" | fgrep "TGP" | wc -l`
	indels_uk10k=`tabix ${in_path} ${i} | fgrep "INDEL" | fgrep "UK10K" | wc -l`
	#Putative novel (no rsID and no in tgp or uk10k)
	indels_novel=`tabix ${in_path} ${i} | fgrep "INDEL" | awk '$3=="."' | fgrep -v "UK10K" | fgrep -v "TGP" | wc -l`
	indels_novel_no_singleton=`tabix ${in_path} ${i} | fgrep "INDEL" | awk '$3=="."' | fgrep -v "UK10K" | fgrep -v "TGP" | fgrep -v "SINGLETON"| wc -l`
		
	
	echo "${i} ${snp_w_anno} ${snp_wo_anno} ${snp_singletons} ${snp_tgp} ${snp_uk10k} ${snp_novel} ${snp_novel_no_singleton} ${indels_w_anno} ${indels_wo_anno} ${indels_singletons} ${indels_tgp} ${indels_uk10k} ${indels_novel} ${indels_novel_no_singleton}"
 
done) > ${out_path}/${pop}_${chr}_count_table.txt

