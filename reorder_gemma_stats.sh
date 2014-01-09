#!/usr/local/bin/bash

#reorder the qctool output
out_path=$1
stat_filename=$2
gen_info=$3
chr=$4

#(echo "CHR POS rsID all_A all_B minor major MAF info_qctool info_impute HWE";join <(cut -f 2-8,15,16,19 $out_path/$stat_filename | grep -v RSID ) <(cut -f 2,5 -d " " $gen_info | grep -v rs_id ) | awk '{print $2,$3,$1,$4,$5,$6,$7,$8,$9,$11,$10}') > $out_path/$stat_filename.reordered
#reorder for snptest format
#(echo "CHR POS rsID all_0 all_1 MAF info_snptest info_impute HWE QUAL_snptest QUAL_impute";join -1 1 -2 1 <(cut -f 2-6,8,9,18,20 -d " " $stat_filename | grep -v chromosome | sort -k1,1) <(cut -f 2,5,6 -d " " $gen_info | grep -v rs_id | sort -k1,1) | awk -v chr=${chr} '{print chr,$3,$1,$4,$5,$8,$7,$10,$9,$6,$11}'| sort -g -k2,2 ) > $stat_filename.reordered


#(echo "CHR POS rsID all_A all_B minor major MAF info_qctool info_impute HWE";join <(cut -f 2-8,15,16,19 $out_path/$stat_filename | grep -v RSID ) <(cut -f 2,5 -d " " $gen_info | grep -v rs_id ) | awk '{print $2,$3,$1,$4,$5,$6,$7,$8,$9,$11,$10}') > $out_path/$stat_filename.reordered
#reorder for snptest format
echo "Join and format..."
#(echo "CHR POS rsID all_0 all_1 MAF info_snptest info_impute HWE QUAL_snptest QUAL_impute";/usr/bin/join -1 1 -2 1 <(cut -f 2-6,8,9,18,20 -d " " $stat_filename | grep -v chromosome | sort -k1,1) <(cut -f 2,5,6 -d " " $gen_info | grep -v rs_id | sort -k1,1) | awk -v chr=${chr} '{print chr,$3,$1,$4,$5,$8,$7,$10,$9,$6,$11}'| sort -g -k2,2 ) > $stat_filename.reordered
#(echo "CHR POS rsID all_0 all_1 MAF info_snptest info_impute HWE QUAL_snptest QUAL_impute";awk "{l=\$0;id=\$2;while(getline < \"${stat_filename}\") {if(id == \$2) print \$0,l;break}}" $gen_info | grep -v chromosome | grep -v rs_id | awk -v chr=${chr} '{print $0}'| sort -g -k2,2 ) > $stat_filename.reordered
#(echo "CHR POS rsID all_0 all_1 MAF info_snptest HWE QUAL_snptest";fgrep -v chromosome $stat_filename | cut -f 2-6,8-9,18,20 -d " " | awk -v chr=${chr} '{OFS=" "}{print chr,$3,$1,$4,$5,$8,$7,$9,$6}'| sort -g -k2,2) > $stat_filename.reordered
(echo "CHR POS rsID all_0 all_1 MAF info_snptest HWE QUAL_snptest all_0_freq all_1_freq eff_all_freq";fgrep -v chromosome $stat_filename | cut -f 2-6,8-9,14-18,20 -d " " | awk -v chr=${chr} '{OFS=" "}{print chr,$3,$1,$4,$5,$12,$7,$13,$6,(($8*2)+$9)/(($8+$9+$10+$11)*2),(($10*2)+$9)/(($8+$9+$10+$11)*2),(($8*2)+$9)/(($8+$9+$10+$11)*2)}'| sort -g -k2,2) > $stat_filename.reordered

(echo "CHR POS rsID all_0 all_1 MAF info_snptest HWE QUAL_snptest all_0_freq all_1_freq eff_all_freq";tail -n+2 $stat_filename.reordered | awk '{print "chr"$1":"$2,$3,$0}' | awk '{
if($1==$2) 
  printf $3" "$4" "$1"-NA-"$6"-"$7;
else
  printf $3" "$4" "$1"-"$2"-"$6"-"$7;

{for(i=6;i<=NF;i++) 
  printf " " $i;
  printf "\n";}
 }') > $stat_filename.reordered.final


#in case of merged imputation, there is no $gen_info file, so the next passage will be skipped
if [ -f $gen_info ]; then

	echo -e "CRACKING!!\nInfo file from IMPUTE available!!"

	(echo "CHR POS rsID exp_freq_a1 info_impute QUAL_impute";fgrep -w -f <(cut -f 3 -d " " $stat_filename.reordered | grep -v chromosome ) $gen_info | awk -v chr=${chr} '{print chr,$3,$2,$4,$5,$6}'| sort -g -k2,2 ) > $stat_filename.gen_info.reordered
	header=`join -1 3 -2 3 $stat_filename.reordered $stat_filename.gen_info.reordered | cut -f 1-12,15- -d " " | head -1`
	(echo "CHR POS RSid all_0 all_1 MAF info_snptest HWE QUAL_snptest all_0_freq all_1_freq eff_all_freq exp_freq_a1 info_impute QUAL_impute";join -1 3 -2 3 $stat_filename.reordered $stat_filename.gen_info.reordered | cut -f 1-12,15- -d " " | sort -g -k3,3 | tail -n+2 | awk '{print "chr"$2":"$3,$1,$0}'| awk '{
	 	if($1==$2) 
	 	  printf $4" "$5" "$1"-NA-"$6"-"$7;
	 	else
	 	  printf $4" "$5" "$1"-"$2"-"$6"-"$7;
	 
	 	{for(i=6;i<=NF;i++) 
	 	  printf " " $i;
	 	  printf "\n";}
	 	 }') > $stat_filename.merged.final
else
	echo -e "D'OH!!!!\nInfo file from IMPUTE NOT available!!"
	echo -e "Information extraction terminated with file : ${stat_filename}.reordered.final"
fi   
