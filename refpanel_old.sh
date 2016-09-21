#! /usr/local/bin/bash

#/nfs/team151/uk10k/shared_files/UK10K_3798-samples_annotated_final_20130205.txt
g1kdir=/lustre/scratch107/projects/uk10k/users/jh21/references_panel/1kg
uk10kvcfdir=/lustre/scratch106/projects/uk10k/RELEASE/UK10K_COHORT/REL-2012-06-02/v2
uk10kdir=/lustre/scratch107/projects/uk10k/users/jh21/references_panel/uk10k
uk10k1kgdir=/lustre/scratch107/projects/uk10k/users/jh21/references_panel/uk10k1kg
vcftools=/nfs/users/nfs_j/jh21/programs/vcftools_0.1.10/bin/vcftools

for chr in {20..20}; do
        chrX_vcf_str=""
        chrX_tped_str=""
        if [[ "$chr" == "X" ]]; then
                chrX_vcf_str=" && \$2 > 2699520 && \$2 < 154931044"
                chrX_tped_str=" | sed 's/\\\t\([A-Z]\+\)\\\t/\\\t\1 \1\\\t/g' |  sed 's/\\\t\([A-Z]\+\)\\\t/\\\t\1 \1\\\t/g' | sed 's/\\\t\([A-Z]\+\)$/\\\t\1 \1/m' "
        fi

	#### create list of variants for exclusion ###
	echo -e "#!/usr/local/bin/bash
	\necho \"SNP CHR POS rsID REF ALT AN\" > $g1kdir/chr$chr.stats
	\nzcat $g1kdir/vcfs/ALL.chr$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz | awk '\$1 !~/^#/ {print \"chr\"\$1\":\"\$2, \$1,\$2,\$3,\$4,\$5 \"\\\t\" \$8}' | sed 's/\\\t.*AC=\\\([0-9]*\\\).*/ \\\1/' >> $g1kdir/chr$chr.stats
	\necho \"SNP CHR POS rsID REF ALT AN\" > $uk10kdir/chr$chr.stats
	\nzcat $uk10kvcfdir/$chr.beagle.anno.csq.20130126.vcf.gz | awk '\$1 !~/^#/ {print \"chr\"\$1\":\"\$2, \$1,\$2,\$3,\$4,\$5 \"\\\t\" \$8}' | sed 's/\\\t.*AC=\\\([0-9]*\\\).*/ \\\1/' >> $uk10kdir/chr$chr.stats
	\nawk '{if (g[\$3] !=\"Y\") print \$0; g[\$3]=\"Y\"}' $g1kdir/chr$chr.stats > $g1kdir/chr$chr.stats.uniq
	\nawk '{if (g[\$3] !=\"Y\") print \$0; g[\$3]=\"Y\"}' $uk10kdir/chr$chr.stats > $uk10kdir/chr$chr.stats.uniq
	\nplink --noweb --id-match $g1kdir/chr$chr.stats.uniq POS $uk10kdir/chr$chr.stats.uniq POS --out chr$chr
	\nsed -i 's/  */ /g' chr$chr.matched
	\necho \"reason SNP CHR POS ID REF ALT AC\" > $g1kdir/chr$chr.exclude.txt
	\nawk '{ if (\$5\$6 ~/,/) print \"multi-allelic\",\$0; else if (\$4 ~/^esv/) print \"esv_variants\",\$0; else if (g[\$3]==\"Y\") print \"dup_positions\",\$0; g[\$3]=\"Y\" }' $g1kdir/chr$chr.stats >> $g1kdir/chr$chr.exclude.txt
	\n## if ((\$7==1 || \$7==2183) && \$8==\"NA\") print \"singletons\",\$1,\$2,\$3,\$4,\$5,\$6,\$7;  KEEP SINGLETONS FOR NOW
	\nawk '{ if (\$5\$6 !=\$12\$13 && \$5\$6 !=\$13\$12 && \$1 !=\"NA\" && \$8 !=\"NA\") print \"mismatched-alleles\",\$1,\$2,\$3,\$4,\$5,\$6,\$7;}' chr$chr.matched >> $g1kdir/chr$chr.exclude.txt
	\nawk 'NR >1 {print \$3,\$4}' $g1kdir/chr$chr.exclude.txt | sort -n -k 2,2 | uniq > $g1kdir/chr$chr.exclude
	\necho \"reason SNP CHR POS ID REF ALT AN AC\" > $uk10kdir/chr$chr.exclude.txt
	\nawk '{ if (\$5\$6 ~/,/) print \"multi-allelic\",\$0; else if (\$4 ~/^esv/) print \"esv_variants\", \$0; else if (g[\$3]==\"Y\") print \"dup_positions\",\$0; g[\$3]=\"Y\" }' $uk10kdir/chr$chr.stats >> $uk10kdir/chr$chr.exclude.txt
	\nawk '{ if ((\$14==1 || \$14==7561) && \$1==\"NA\") print \"singletons\",\$8,\$9,\$10,\$11,\$12,\$13,\$14 }' chr$chr.matched >> $uk10kdir/chr$chr.exclude.txt
	\nawk 'NR >1 {print \$3,\$4}' $uk10kdir/chr$chr.exclude.txt | sort -n -k 2,2 | uniq > $uk10kdir/chr$chr.exclude
	" > $uk10k1kgdir/chr$chr.cmd
	cd $uk10k1kgdir
	bsub -J "uk10k1kg.ref.chr$chr" -q basement -o chr$chr.LOG -e chr$chr.ERR -R "select[mem>12000] rusage[mem=12000]" -M12000000 < chr$chr.cmd

	#### create 1000GP statistics and PLINK ###
	echo -e "#!/usr/local/bin/bash
	\nn_excluded=\`wc -l chr$chr.exclude | awk '{printf \$1}'\`
	\ngzip -f chr$chr.exclude
	\nzcat chr$chr.exclude.gz vcfs/ALL.chr$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz | awk -v n_excl=\$n_excluded '{if(NR<=n_excl) pos[\$2]=\"Y\"; else { if (\$1 ~/^#/ || (pos[\$2] !=\"Y\" $chrX_vcf_str)) print \$0} }' |  bgzip -f > chr$chr.vcf.gz
	\ntabix -p vcf chr$chr.vcf.gz
        \nzcat chr$chr.vcf.gz | awk 'BEGIN{print \"id position a0 a1\"}{ if (\$1 ~/^#/); else if (\$3==\".\") print \"chr\"\$1\":\"\$2, \$2,\$4,\$5; else print \$3,\$2,\$4,\$5}' | gzip -f > chr$chr.legend.gz
	\nvcf-query -f '[\\\t%GTR]\\\n' chr$chr.vcf.gz | sed 's/^\s*//' | sed 's/\// /g' | sed 's/|/ /g' | gzip -f > chr$chr.hap.gz
	\nvcf-query -f '%CHROM %ID 0 %POS [\\\t%GT]\\\n' chr$chr.vcf.gz | sed 's/\// /g' | sed 's/|/ /g' $chrX_tped_str | awk '{if (\$2==\".\") \$2=\"chr\"\$1\":\"\$4; print \$0}' > chr$chr.tped
	\nvcf-query -l chr$chr.vcf.gz >  chr$chr.sample.ids
	\nawk '{print \$1,\$1,0,0,0,0}' chr$chr.sample.ids > chr$chr.tfam
	\nplink --noweb --tfile chr$chr --make-bed --out chr$chr
	\nmv chr$chr.bim chr$chr.bim.COPY
	\nawk '{\$2=\"chr\"\$1\":\"\$4; print \$0}' chr$chr.bim.COPY | sed 's/chr23/chrX/' > chr$chr.bim
	" > $g1kdir/chr$chr.cmd
	if [[ $chr -eq 20 ]]; then
		echo "plink --noweb --bfile chr$chr --recode --tab --out chr$chr" >> $g1kdir/chr$chr.cmd
	fi
	cd $g1kdir
	bsub -J "g1k.ref.chr$chr" -w "ended(uk10k1kg.ref.chr$chr)" -q basement -o chr$chr.LOG -e chr$chr.ERR -R "select[mem>12000] rusage[mem=12000]" -M12000000 < chr$chr.cmd

	#### create UK10K ref panel, statistics and PLINK ###
	echo -e "#!/usr/local/bin/bash
	\nn_excluded=\`wc -l chr$chr.exclude | awk '{printf \$1}'\`
	\ngzip -f chr$chr.exclude
	\nzcat chr$chr.exclude.gz $uk10kvcfdir/$chr.beagle.anno.csq.20130126.vcf.gz | awk -v n_excl=\$n_excluded '{if(NR<=n_excl) pos[\$2]=\"Y\"; else { if (\$1 ~/^#/ || (pos[\$2] !=\"Y\" $chrX_vcf_str)) print \$0} }' |  bgzip -f > chr$chr.vcf.gz
	\ntabix -p vcf chr$chr.vcf.gz
	\nzcat chr$chr.vcf.gz | awk 'BEGIN{print \"id position a0 a1\"}{ if (\$1 ~/^#/); else if (\$3==\".\") print \"chr\"\$1\":\"\$2, \$2,\$4,\$5; else print \$3,\$2,\$4,\$5}' | gzip -f > chr$chr.legend.gz
	\nvcf-query -f '[\\\t%GTR]\\\n' chr$chr.vcf.gz | sed 's/^\s*//' | sed 's/\// /g' | sed 's/|/ /g' | gzip -f > chr$chr.hap.gz
	\nvcf-query -f '%CHROM %ID 0 %POS [\\\t%GT]\\\n' chr$chr.vcf.gz | sed 's/\// /g' | sed 's/|/ /g' $chrX_tped_str | awk '{if (\$2==\".\") \$2=\"chr\"\$1\":\"\$4; print \$0}' > chr$chr.tped
	\nvcf-query -l chr$chr.vcf.gz >  chr$chr.sample.ids
        \nawk '{print \$1,\$1,0,0,0,0}' chr$chr.sample.ids > chr$chr.tfam
        \n/nfs/users/nfs_j/jh21/scripts/library/join_file.pl -i \"chr$chr.sample.ids,SPACE,0 /nfs/team151/jh21/files/master_twins_gwas_seqs.ids,SPACE,1\" -a 1 -o chr$chr.sample.ids.merged
        \nawk '{ if (\$2==\"NA\") print \$1; else print \$2 }' chr$chr.sample.ids.merged | sed 's/[ab] / /' > chr$chr.sample.ids
	\nplink --noweb --tfile chr$chr --make-bed --out chr$chr
	\nplink --noweb --bfile chr$chr --remove /nfs/team151/jh21/files/samples_exclusion_twinsalspac_REL2012_06_02.txt --make-bed --out chr$chr.clean
	\nmv chr$chr.bim chr$chr.bim.COPY
	\nawk '{\$2=\"chr\"\$1\":\"\$4; print \$0}' chr$chr.bim.COPY > chr$chr.bim
	" > $uk10kdir/chr$chr.cmd
	if [[ $chr -eq 20 ]]; then
		echo "plink --noweb --bfile chr$chr --recode --tab --out chr$chr" >> $uk10kdir/chr$chr.cmd
	fi
	cd $uk10kdir
	bsub -J "uk10k.ref.chr$chr" -w "ended(uk10k1kg.ref.chr$chr)" -q basement -o chr$chr.LOG -e chr$chr.ERR -R "select[mem>12000] rusage[mem=12000]" -M12000000 < chr$chr.cmd

	#### create combined reference panel
	echo -e "
	\nplink --noweb --bfile $uk10kdir/chr$chr --bmerge $g1kdir/chr$chr.bed $g1kdir/chr$chr.bim $g1kdir/chr$chr.fam --geno 0.1 --make-bed --out chr$chr
	\nplink --noweb --bfile chr$chr --freq --out chr$chr.frq
        \n# /nfs/users/nfs_j/jh21/scripts/impute2go/refpanel.pl $chr $g1kdir/chr$chr.vcf.gz $uk10kdir/chr$chr.vcf.gz $uk10k1kgdir all ## LIST THE SECONDARY FILE FIRST
	" > $uk10k1kgdir/chr$chr.cmd2
	if [[ $chr -eq 20 ]]; then
		echo "plink --noweb --bfile chr$chr --recode --tab --out chr$chr" >> $uk10k1kgdir/chr$chr.cmd2
	fi
	cd $uk10k1kgdir
#	echo process chr$chr
#	bsub -J "ref.chr$chr" -w "ended(*.ref.chr$chr)" -q basement -o chr$chr.LOG -e chr$chr.ERR -R "select[mem>12000] rusage[mem=12000]" -M12000000 < chr$chr.cmd2
done
