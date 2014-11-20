#!/bin/bash
exec > >(tee marchini_formatter.log)
#
#
#Script to create files for phasing with shapeit v2
#Usage:

# ./marchini_formatter.sh <inputfile> <geno plain file with path>

####################################################################################
if [ $# -lt 2 ]
then
	echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
	echo -e "./marchini_formatter.sh <inputfile> <geno plain file with path> \n"
	exit
fi
#MAP/DAT creation
#Read the geno argument, extract info for chr
#
args=("$@")
inputfile=${args[0]}
genofile=${args[1]}
inputfilename=$(basename $inputfile)
genofilename=$(basename $genofile)

echo "Create output folder..."
mkdir -p $inputfilename

echo "Copy map file in output folder..."
cp ${genofile%.*}.map $inputfilename/

#add coded sex, beacuse in input file the sex is coded M for MALE and F for FEMALE
awk '{if($5 == "F") print $1,$2,$3,$4,2,$6}' $inputfile > $inputfilename.female_coded.tmp
awk '{if($5 == "M") print $1,$2,$3,$4,1,$6}' $inputfile > $inputfilename.male_coded.tmp

#Sort input list by the Clinic id column
echo "Sort input list by the Clinic id column"
cat $inputfilename.female_coded.tmp $inputfilename.male_coded.tmp | sort -k 6,6 > $inputfilename.sorted.tmp

#Now remove useless columns from genotypes file in plink format
echo "Now remove useless columns from genotypes file in plink format.."
cut -f 1,7- -d " " $genofile > $genofilename.partial.tmp
#we now have the XX id and the genotypes foreach individual
echo "Check and extract untyped individuals..."
#now we have to check if we have XXids for non genotyped people
grep -v -w -f <(cut -f 1 -d " " $genofilename.partial.tmp) $inputfilename.sorted.tmp > not_geno.list
echo "Now extract from the input ped file the information for untyped indivs..."
#now extract from the input ped file the information for untyped indivs
awk '{print $6,$1,$2,$3,$4,$5}' not_geno.list > $inputfilename.nogeno

echo "Join the inputfile with genofile..."
#now join the inputfile with genofile
#Add extracted genotypes to the formatted output created before: return ONLY JOINED LINES
#and remove the first column that contains the XX id
#join -1 6 -2 1 $inputfile.sorted.tmp $genofile.partial.tmp | cut -f 2- -d " " > complete_chr${i}.ped.tmp
join -1 6 -2 1 $inputfilename.sorted.tmp $genofilename.partial.tmp > $genofilename.last_step
echo "Adding information about untyped individuals at the end of the genotype file"
#add information about untyped individuals at the end of the genotype file
cat $genofilename.last_step $inputfilename.nogeno | cut -f 2- -d " " | sort -g -k 1,1 > $genofilename.final
echo "Move final file in output folder"
mv $genofilename.final $inputfilename/
cd $inputfilename/
mv $genofilename.final $genofilename
cd ..
echo "DONE!!"

