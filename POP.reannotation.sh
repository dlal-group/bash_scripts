#!/usr/local/bin/bash
#reannotate all 

if [ $# -lt 2 ]
then
        echo -e "\nError!!Missing arguments\n\n****** USAGE *****"
        echo -e "POP.reannotation.sh <file_path> <output_path>\n"
        echo -e "<file_path>: population files path\n"
        echo -e "<output_path>: population files output path\n"
        echo -e "All paths MUST BE ABSOLUTE!!\n"
        exit 1
fi

args=("$@")
filepath=${args[0]}
outfilepath=${args[1]}

mkdir -p $outfilepath
cd $filepath
files=`ls *.gz`

#Now we need to work by chromosome!
echo "Reannotate population chr files files with correct AN,AC count"
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
for file in $files
do
#now we want annotate again the file
echo $file
bsub -J "${file}_annotate" -o "%J_${file}_annotate.log" -M8000000 -R"select[mem>8000] rusage[mem=8000]" \
-q basement "zcat ${file} | vcf-annotate --fill-ICF | bgzip -c > ${outfilepath}/${file}.gz"
done

