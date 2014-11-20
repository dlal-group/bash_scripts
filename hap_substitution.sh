#! /bin/bash
# Print all unique fields
echo $allele0

COL=`awk '{if (NR==1) print NF}' aplo10snp`
for (( col=1; col<=$COL; col++ ))
do
echo "Type the FIRST allele to substitute (0 allele), followed by [ENTER]:"
read allele0
cut -f $col -d ' ' | sed -i '{s/'${allele0}'/0/}' aplo10snp;
#awk '(all0='$allele0')(col=$'$col'){gsub(/all0/,"0",col);print}' aplo10snp
#echo "Welcome $col imes.."
	#awk '(c=$'$col')!(c in a){a[c];print c}' prova.txt | awk 'ORS=NR%2?" ":"\n"'
#	cut -f $col -d ' ' prova.txt | awk '!($0 in a){a[$0]=$0} END { for (i in a) print(a[i]) }' | awk '{sub(/--/,"\n");print}' | awk 'ORS=NR%2?" ":"\n"'
#echo "Type the SECOND allele to substitute (1 allele), followed by [ENTER]:"
#read allele1
#awk 'BEGIN{(all1='$allele1',col=$'$col')}{gsub(/all1/,"1",col);print}' aplo10snp
echo "Type the SECOND allele to substitute (1 allele), followed by [ENTER]:"
read allele1
cut -f $col -d ' ' | sed -i '{s/'${allele1}'/1/}' aplo10snp
done 

