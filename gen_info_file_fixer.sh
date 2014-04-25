#!/usr/local/bin/bash
geno_file=$1
info_file=$2
chr=$3

#fix geno and info file
#remove duplicates
#gunzip -c ${geno_file} | awk '!($0 in a){a[$0];print}' > ${geno_file}.gen.tmp
gunzip -c ${geno_file}  > ${geno_file}.gen.tmp

#insert the right values for rsId instead of missing
awk -v chr_n=${chr} 'BEGIN{s=0}{
if($1 == "." && $2 == ".")
    printf "%s %s","chr"chr_n":"$3":"s,"chr"chr_n":"$3":"s;
else 
    printf "%s %s",$1":"s,$2":"s;
s=s + 1
{for(i=3;i<=NF;i++) 
  printf " " $i;
  printf "\n";}
      }' ${geno_file}.gen.tmp | gzip -c > ${geno_file}

#fix info file
#awk '!($0 in a){a[$0];print}' ${info_file} > ${info_file}.tmp
fgrep -v position ${info_file} > ${info_file}.tmp

awk -v chr_n=${chr} 'BEGIN{t=0}{
if($1 == "." && $2 == ".")
    printf "%s %s","chr"chr_n":"$3":"t,"chr"chr_n":"$3":"t;
else
    printf "%s %s",$1":"t,$2":"t;
t=t + 1
{for(i=3;i<=NF;i++)
  printf " " $i;
  printf "\n";}
      }' ${info_file}.tmp > ${info_file}.tmp2

cat <(fgrep -w position ${info_file}) ${info_file}.tmp2 > ${info_file}
#rm ${info_file}.tmp;rm ${info_file}.tmp2


