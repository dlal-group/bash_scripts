#!/usr/local/bin/bash
#$ -S /bin/sh
#$ -cwd
#$ -l mem=8G,time=03:00:00
#!/bin/sh


HOM=$4

HET=$5

BITS=$6

#MORE="-h_extend -homoz"

#MORE="-h_extend -homoz" #to use in case the phasing is not so accurate

MORE="-haploid"
ARGUMENTS=" -min_m 0.5 -err_hom $HOM -err_het $HET -bits $BITS $MORE"

# --- SOFTWARE SOURCES --- #

PIPE=/SOFTWARE/germline/

GERMLINE="$PIPE/bin/germline"

# --- SOFTWARE SOURCES --- #


if [ $# -lt 3 ]; then

  echo "---"
  echo " GERMLINE IBD detection script"
  echo " $PIPE gline.sh"

  echo "---"

  echo -e "Usage:\t$0 [phased ped file] [map file] [output name] [optional parameters]"

  echo -e "Output:\tGenerates [output] match file with IBD segments"

  exit

fi


PED=$1

MAP=$2

OUT=`echo $3.$HOM.$HET.$BITS.$MORE | tr ' ' '.'`

#OUT=$3.$HOM.$HET.$BITS.$MORE


for f in $GERMLINE $PED $MAP; do

  if [ ! -f $f ]; then

    echo "The file $f does not exist, please check that it is referenced properly"

    exit

  fi

done


echo -e "1\tRunning GERMLINE"

echo '1' > $OUT.run

echo $MAP >> $OUT.run

echo $PED >> $OUT.run

echo $OUT >> $OUT.run


cat $OUT.run


$GERMLINE $ARGUMENTS < $OUT.run


rm $OUT.run 

