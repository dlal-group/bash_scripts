#!/bin/bash

cat tiroide.phl > tiroide_temp
X=0
DONE=false
FILE="tiroide_temp"
until $DONE
do
        X=`expr ${X} + 1`
        #subset -Ftiroide_temp -K0 > tsub"$X".phl  
	head -n2 tiroide_temp > tsub"$X".phl  
        grep -wv -f tsub"$X".phl tiroide_temp > new_tiroide_temp
        mv new_tiroide_temp tiroide_temp
        echo $X
	if [[ -s $FILE ]] ; then
		echo "$FILE has data."		
	else
		echo "$FILE is empty."
		echo "DONE!"
		DONE=true
	fi ;
done
