#!/bin/bash

#grep "SAMPLE" QC_OUT/CHR${i}/sample_concordance_discordance_table_chr1.txt > QC_OUT/sample_header.txt
#We read the files created for samples by chr and merge them toghether
for i in {1..22}
do
if [ ! -s QC_OUT/sample_header.txt ]
then
	grep "SAMPLE" QC_OUT/CHR${i}/sample_concordance_discordance_table_chr${i}.txt > QC_OUT/sample_header.txt
fi

grep -v "SAMPLE" QC_OUT/CHR${i}/sample_concordance_discordance_table_chr${i}.txt >> QC_OUT/all_samples_all_chr_concordance_discordance_table_nohead.txt
done

cat QC_OUT/sample_header.txt QC_OUT/all_samples_all_chr_concordance_discordance_table_nohead.txt > QC_OUT/all_samples_all_chr_concordance_discordance_table.txt

#now call the script to calculate the overall GC/NRD by sample
python2.7 ~/Work/bash_scripts/non_ref_discordance_by_sample.py QC_OUT/all_samples_all_chr_concordance_discordance_table.txt $1

#move files:
mv all_sample_concordance_discordance* QC_OUT/

#create a unified fiele for all sites
(echo "CHR POS OGC NRD"
for i in {1..22}
do
  grep -v CHR CHR${i}/site_concordance_discordance_table_chr${i}.txt
done) | tr " " "\t" > all_sites_all_chr_concordance_discordance_table.txt

