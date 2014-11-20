#!/usr/local/bin/bash

#args:
#$1=transcode_filename
#$2=transcoded_filename

#this works for a transcoder file which has space as delimiter and useful values in 1 and 5 colmn. The substitution is performed IN PLACE!!!A bkup copy of original file is done.
#make a backup copy of transcoded file
cp $2 $2.old

cut -f 1,5 $1 | while read src rep
do
	sed -i "s/ $src / $rep /g" $2
done


