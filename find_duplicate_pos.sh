#!/usr/local/bin/bash

#script to find a common set of varints between info file from imputation
#we perform a check on the allele columns added to the info_file

#Args list:
#$1=geno_file_path
#$2=first set path
#$3=second set path
#$4=out_file path



geno_file_path=$1
first_set_file_path=$2
second_set_file_path=$3
out_file_path=$4

        # find common variants in both set based on position and allele match                 
        #awk 'FNR==NR{A[$3,$11]=1;next} A[$3,$11]' 370K/CHR${i}/chr${i}.geno_info_allele 700K/CHR${i}/chr${i}.geno_info_allele > chr${i}.370in700k
        awk 'FNR==NR{A[$3,$4,$5]=1;next} A[$3,$11,$12]' $geno_file_path $first_set_file_path > $out_file_path.first_set
        #awk 'FNR==NR{A[$3,$11]=1;next} A[$3,$12]' 370K/CHR${i}/chr${i}.geno_info_allele 700K/CHR${i}/chr${i}.geno_info_allele > chr${i}.370in700k.2
        awk 'FNR==NR{A[$3,$4,$5]=1;next} A[$3,$12,$11]' $geno_file_path $first_set_file_path > $out_file_path.first_set.flipped

        #awk 'FNR==NR{A[$3,$11]=1;next} A[$3,$11]' 700K/CHR${i}/chr${i}.geno_info_allele 370K/CHR${i}/chr${i}.geno_info_allele > chr${i}.700in370k
        awk 'FNR==NR{A[$3,$4,$5]=1;next} A[$3,$11,$12]' $geno_file_path $second_set_file_path > $out_file_path.second_set
        #awk 'FNR==NR{A[$3,$11]=1;next} A[$3,$12]' 700K/CHR${i}/chr${i}.geno_info_allele 370K/CHR${i}/chr${i}.geno_info_allele > chr${i}.700in370k.2
        awk 'FNR==NR{A[$3,$4,$5]=1;next} A[$3,$12,$11]' $geno_file_path $second_set_file_path > $out_file_path.second_set.flipped

        #fix the flipped alleles in info files:
        #invert expected a1 frequency with 1 - exp_af1
        fsf=`wc -l $out_file_path.first_set.flipped|cut -f 1 -d " "`
        if [ $fsf -gt 0 ]
        then
                awk '{print $1,$2,$3,(1-$4),$5,$6,$7,$8,$9,$10,$12,$11}' $out_file_path.first_set.flipped > $out_file_path.first_set.flipped.fixed
        else
                cp $out_file_path.first_set.flipped $out_file_path.first_set.flipped.fixed
        fi

        ssf=`wc -l $out_file_path.second_set.flipped | cut -f 1 -d " "`
        if [ $ssf -gt 0 ]
        then
                awk '{print $1,$2,$3,(1-$4),$5,$6,$7,$8,$9,$10,$12,$11}' $out_file_path.second_set.flipped > $out_file_path.second_set.flipped.fixed
        else
                cp $out_file_path.second_set.flipped $out_file_path.second_set.flipped.fixed
        fi

        #metto tutto insieme e poi rimuovo le righe completamente identiche
        cat $out_file_path.first_set $out_file_path.first_set.flipped.fixed > $out_file_path.first_set.merged
        cat $out_file_path.second_set $out_file_path.second_set.flipped.fixed > $out_file_path.second_set.merged

        #we need to merge back all togheter and order by postion and alleles!!  
        awk '!($0 in a){a[$0];print}' $out_file_path.first_set.merged | sort -k3,3g -k11,11 -k12,12> $out_file_path.first_set
        awk '!($0 in a){a[$0];print}' $out_file_path.second_set.merged| sort -k3,3g -k11,11 -k12,12> $out_file_path.second_set

        # trovo le posizioni duplicate (indel e rs) nel risultato del confronto
        awk  '{if (x[$3]) { x_count[$3]++; print $0; if (x_count[$3] == 1) { print x[$3] } } x[$3] = $0}' $out_file_path.first_set > $out_file_path.duplicates

        #remove useless files
        rm $out_file_path.*.flipped $out_file_path.*.flipped.fixed $out_file_path.*.merged
