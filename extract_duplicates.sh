#!/usr/local/bin/bash

#extract duplicates for imputation results:
#from info and geno files
#and extract also only snps and save in a file
#$3 is CHR


#first only for position (col 3)
#take in account only variation with rs annotation
awk  '{if (x[$3]) { x_count[$3]++; print $0; if (x_count[$3] == 1) { print x[$3] } } x[$3] = $0}' <(fgrep rs $1.gen_info) > $2.info_only_pos

#then also for rsId (cols 2,3)
awk  '{if (x[$2,$3]) { x_count[$2,$3]++; print $0; if (x_count[$2,$3] == 1) { print x[$2,$3] } } x[$2,$3] = $0}' <(fgrep rs $1.gen_info ) > $2.info_pos_and_rs_id

#take in account only variation with rs annotation
awk  '{if (x[$3]) { x_count[$3]++; print $0; if (x_count[$3] == 1) { print x[$3] } } x[$3] = $0}' <(zcat $1.gen.gz | fgrep rs | cut -f 1-5 -d " ") > $2.gen_only_pos

#then also for rsId (cols 2,3)
awk  '{if (x[$2,$3]) { x_count[$2,$3]++; print $0; if (x_count[$2,$3] == 1) { print x[$2,$3] } } x[$2,$3] = $0}' <(zcat $1.gen.gz | fgrep rs | cut -f 1-5 -d " ") > $2.gen_pos_and_rs_id

#then also for pos and alleles length (cols 3,4,5)
awk  '{if (x[$3,length($4),length($5)]) { x_count[$2,length($4),length($5)]++; print $0; if (x_count[$2,length($4),length($5)] == 1) { print x[$2,length($4),length($5)] } } x[$2,length($4),length($5)] = $0}' <(zcat $1.gen.gz | fgrep rs | cut -f 1-5 -d " ") > $2.gen_pos_and_alleles

#now remove the duplicates by rs ids (col 2)
fgrep -v -w -f <(cut -f 2 -d " " $2.gen_pos_and_rs_id) $1.gen_info_allele > $1.geno_info_no_dup

#extract only snp!!
zcat $1.gen.gz | fgrep rs | cut -f 1-5 -d " "| awk '{if (((length($4) == 1) && ($4 != "-")) && ((length($5) == 1) && ($5 != "-"))) print $0}' > $1.gen_rs_only

#extract from info file only rs snps
(fgrep position $1.gen_info;fgrep -w -f <(cut -f 2 -d " " $1.gen_rs_only) $1.gen_info) > $1.gen_info_rs_only

#calculate and reformat the file adding a MAF column at the end
(echo "CHROM `fgrep position $1.gen_info_rs_only` MAF";fgrep -v position $1.gen_info_rs_only | awk -v chr=$3 '{OFS=" "}{if ($4 <= 0.5) print chr,$0,$4;else print chr,$0,(1-$4)}') | cut -f 1,3- -d " " > $1.gen_info_rs_only_maf

