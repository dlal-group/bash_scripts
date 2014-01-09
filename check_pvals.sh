#!/usr/local/bin/bash
#We can use this same launcher to check data using BASH or R
#We will pass a flag for the R or BASH way
#Args:
#$1 = file to process
#$2 = gz flag
#$3 = trait
if [ $# -lt 1 ]
then
	echo "MISSING ARGUMENTS!!!"
	echo -e "USAGE:"
	echo -e "BASH mode:\ncheck_pvals.sh <file_path>"
	echo -e "R mode:\ncheck_pvals.sh <file_path> <gzip flag> <trait name>"
	echo -e "<file_path>: file with results to check\n<gzip flag>: [TRUE/FALSE] used to tell R if the file is in gz format\n<trait name> : name of the trait"
	exit 1
fi


if [ $# -lt 2 ]
then
	echo "Shell Mode"
	echo -e "Processing file $1"
	#check pval calculating lambda and qqplot
	(echo "chr	rs	pos	beta	se	p_wald	chi2";cut -f -3,5,6,9 $1 | tail -n+2 | awk '{OFS=" "}{print $0,(($4/$5)*($4/$5))}' ) | tr " " "\t" | gzip -c > $1_chi2.gz

	# lets try to calculate the same formula with awk : lambda=median(chi2,na.rm=T)/qchisq(0.5,1)
	# qchisq(0.5,1) = 0.4549364 
	#This calculate the MEDIAN value in a NUMERICALLY SORTED column without NA:
	# awk '{count[NR] = $2;}END{if (NR % 2) {print count[(NR + 1) / 2];} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;}}' 

	md=`zcat $1_chi2.gz | tail -n+2 | cut -f 7 | fgrep -v "NA" | sort -g | awk '{count[NR] = $1}END{if (NR % 2) {print count[(NR + 1) / 2];} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;}}'`
	median=`printf "%f\n" $md`
	l=$( bc -l <<< "${median}/0.4549364")
	lambda=$(echo "$l")
	check=$( bc -l <<< $l'>'1.01)
	echo "#######################################"
	echo "The inflation factor fo this trait is:"
	echo "$l"

	if [ $check -eq 1 ]
	then
		echo "You need to correct your p-values!!"
		#calculate the corrected chi2
		# (echo "chr	rs	pos	beta	se	p_wald	chi2	chi2_gc";zcat $1_chi2.gz | tail -n+2 | awk -v lam=$lambda '{print $0, $7/lam}' )| tr " " "\t" | gzip -c > $1_chi2_gc.gz
	else
		echo "You don't need to correct your p-values!!"
	fi
	echo "#######################################"
else
	echo "R mode!"
	if [ $# -lt 3 ]
	then
		echo "MISSING ARGUMENTS for R mode!!!"
		echo -e "USAGE:"
		echo -e "R mode:\ncheck_pvals.sh <file_path> <gzip flag> <trait name>"
		echo -e "<file_path>: file with results to check\n<gzip flag>: [TRUE/FALSE] used to tell R if the file is in gz format\n<trait name> : name of the trait"
		exit 1
	fi
	echo -e "Processing file $1\non trait $3"

	#check pval calculating lambda and qqplot
	(echo "chr	rs	pos	beta	se	p_wald	chi2";cut -f -3,5,6,9 $1 | tail -n+2 | awk '{OFS=" "}{print $0,(($4/$5)*($4/$5))}' ) | tr " " "\t"| gzip -c > $1_chi2.gz
	#Call the R script
	R CMD BATCH '--args '$1_chi2.gz' 'TRUE' '$3'' ~/Work/r_scripts/result_check.r $3.Rout
fi

