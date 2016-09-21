#Script to extract from gemma results common sites with other population
#!/usr/local/bin/bash

#Args:
#$1:first set results path
#$2:second set results path 
#$3:output path
#$4:population name

if [ $# -lt 4 ]
then
        echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
        echo -e "merge_extracted_results.sh <first set result path> <second set result path> <output path> <population name> \n"

        exit 1
fi


first_path=$1
second_path=$2
out_path=$3
pop=$4

for trait in `ls $first_path/`
do
	echo "Processing trait ${trait}..."
	#create the output path
	mkdir -p $out_path/${trait}

	#create a file list for each pathand trait
	list1=( `ls $first_path/${trait}` )
	list2=( `ls $second_path/${trait}` )
	count="0"
	
	#echo ${#list1[@]}

	while [ $count -lt ${#list1[@]} ]
	do
	
		echo ${list1[$count]}
		echo ${list2[$count]}
		res=${list2[$count]}

		#find all common sites
		
grep -w -f <( sort -g -k 2 "$first_path/${trait}/${list1[$count]}" | cut -f 2 ) "$second_path/${trait}/${list2[$count]}" | cut -f 1-7 > $out_path/${trait}/${res%_*}.common
		#join -1 2 -2 2 <(sort -k 2 "$first_path/${trait}/${list1[$count]}" | cut -f 1-7) <(sort -k 2 "$second_path/${trait}/${list2[$count]}" | cut -f 1-7) > $out_path/${trait}/${res%_*}.common
		
		#now find what is not in common
		(fgrep -v -w -f <(cut -f 2 $out_path/${trait}/${res%_*}.common ) $first_path/${trait}/${list1[$count]} | cut -f 1-7;fgrep -v -w -f <(cut -f 1 -d " " $out_path/${trait}/${res%_*}.common ) $second_path/${trait}/${list2[$count]} | cut -f 1-7) > $out_path/${trait}/${res%_*}.diff
		#now merge back all
		#(awk '{OFS="\t"}{print $2,$1,$3,$4,$5,$6,$7}' $out_path/${trait}/${res%_*}.common;cat $out_path/${trait}/${res%_*}.diff) | sort -g -k 2 > $out_path/${trait}/${res%_*}.merged
		(cat $out_path/${trait}/${res%_*}.common $out_path/${trait}/${res%_*}.diff) | sort -g -k 2 > $out_path/${trait}/${res%_*}.merged

		count=$[$count+1]
	
	done	
	echo "Done trait ${trait}!"
	echo $count
	unset count

done
