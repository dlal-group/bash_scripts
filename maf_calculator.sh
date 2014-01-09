#!/usr/local/bin/bash
#Calculate Delta maf between population
#It uses a preformatted set of files, splitted by chr
#file preprocessing...

if [ $# -lt 1 ]
  then
  echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
  echo -e "maf_calculator.sh <AC col number> <AN col number> \n"

  exit 1
fi

FILES=`ls *`

mkdir -p CHR_MAF

#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
for i in $FILES
do
  echo "Calculating af and maf for ${i}.."
  #we have a file with:
  #CHROM  POS ID <some_columns> AN <some_columns> AC <other columns>
  #calculate af and maf
  awk -v AC=$1 -v AN=$2 '
  {OFS="\t"}
  {
    print $0,$(AC)/$(AN)
  }' ${i} \
  | \
  awk '
  {OFS="\t"}
  {
    if ($(NF) > 0.5)
     print $0, 1-$(NF);
   else if ($(NF) <= 0.5)
     print $0, $(NF);
   }' > CHR_MAF/${i}.freq

   #We now have a file with:
   #CHROM  POS ID <some_columns> AN <some_columns> AC <other columns> AF MAF
   echo "MAF file created for chr ${i}!"

done
