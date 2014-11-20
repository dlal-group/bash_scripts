#! /bin/bash
# Print all unique fields
echo "Type the year that you want to check (4 digits), followed by [ENTER]:"
read year
echo $year
#COL=`awk '{if (NR==1) print NF}' prova.txt`
#for (( col=1; col<=$COL; col++ ))
#do
##echo "Welcome $col imes.."
#	#awk '(c=$'$col')!(c in a){a[c];print c}' prova.txt | awk 'ORS=NR%2?" ":"\n"'
#	cut -f $col -d ' ' prova.txt | awk '!($0 in a){a[$0]=$0} END { for (i in a) print(a[i]) }' | awk '{sub(/--/,"\n");print}' | awk 'ORS=NR%2?" ":"\n"'
#done 

