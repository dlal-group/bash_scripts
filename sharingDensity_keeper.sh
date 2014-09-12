#!/usr/local/bin/bash

MATCH=$1 #germline format
# MATCH=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/IBD/GERMLINE/CHR10/FVG.10.non_missing.match

MAP=$2 #plink format in cM
# MAP=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/FVG.10.non_missing.map


# N=$3 #numero individui
N=`wc -l ${MAP%.map}.ped | cut -f 1 -d " "`

first=`head -n 1 $MAP | awk '{print $3}'`
last=`tail -n 1 $MAP | awk '{print $3}'`
minDens=0;
resolution=0.5 #dimensione bin

echo -e "Arguments and parameters:
\rMatch file-> $MATCH
\rMap file -> $MAP
\rN samples -> $N
\rfirst position -> $first
\rLast position -> $last
\rminDens-> $minDens
\rResolution -> $resolution"

#haplotype extension se il phasing non è buono -h_extend file da usare from germline
# T=`head -n 1 $MATCH | awk '{a=substr($2,length($2)-1,2); if (a==".0" || a==".1") print 1; else print 0; }'`

# if [ "$T" -eq "1" ]; then

# cat $MATCH | sed 's/\.[01]//' | sed 's/\.[01]//'

# else

# cat $MATCH

# # fi | awk '{id1=tolower(substr($2,length($2),1)); id2=tolower(substr($4,length($4),1)); if ((id1=="a"||id1=="b")&&(id2=="a"||id2=="b")) print; }' \
# fi | awk -v first=$first -v last=$last -v resolution=$resolution -v map=$MAP -v N=$N -v minDens=$minDens '
# BEGIN{
#   while (getline < map) {
#     gen[$2]=$3;
#   }
#   begin=int(first/resolution);
#   finish=int(last/resolution);
# }
# #snp usati per calcolare ibd diviso lunghezza segmento
# $10/$11>=minDens { 
#   cnt++;
#   # if (cnt%10000==0) print > "/dev/stderr"; 
#   start=int(gen[$8]/resolution);
#   end=int(gen[$9]/resolution);
#   start_w[gen[$8]]=$6
#   end_w[gen[$9]]=$7
#   # print start,end > "/dev/stderr";
#   for (i=start; i<=end; i++) {
#     dens[i]=dens[i]+1;
#     # print start, i, end;
#   }
#   # print gen[$8],gen[$9],start_w[gen[$8]], end_w[gen[$9]] > "check_f.check";
# }
# END{
#   for (i=begin; i<=finish; i++) {
#     print i*resolution "\t" 0+dens[i]/(N*(N-1)/2 - N/2);
#   } 
# }' > $MATCH.shareDens

# La sharing density probabilità che una coppia in quel bin abbia IBD
# prendo histogrammi tolgo regioni che deviano dalla densità media di un 5sd
# now we have the sharingDensity, let's calculate sd and mean for exclusion
# avg=`awk '{sum+=$2} END { printf "%.8f", sum/NR}' $MATCH.shareDens`
# std=`awk -v mu=${avg} '{sumsq += ($2 - mu)^2} END { printf "%.8f", sqrt((sumsq)/NR)}' $MATCH.shareDens`

# max_d=`echo $( bc -l <<< "scale=8; ${avg}+5*${std}") | awk '{printf "%.8f", $1}'`

# extract regions with excess of sharing and retrieve their coordinates from the map file
# awk -v max_sd=${max_d} '$2>=max_sd' $MATCH.shareDens | awk -v map=${MAP} '{}'
# while read line
# do
#   start_r=`echo $line | awk '{print $1*0.5}'`
#   end_r=`echo $line | awk '{print ($1+0.5)*0.5}'`
#   awk -v rstrt=$start_r -v rnd=$end_r '{if($3 >= rstrt && $3<= rnd ) print $2}' $MAP
# done < <(cat $MATCH.shareDens.to_include ) > $MATCH.keepsnps

#modify to extract start and end snp to extract region from plink
for reg_file in `ls ${MATCH}.shareDens_R*.to_include`
do
    start_r=`head -1 ${reg_file} | awk '{print $1*0.5}'`
    end_r=`tail -1 ${reg_file} | awk '{print ($1)*0.5}'`
    reg_start=`awk -v rstrt=$start_r -v rnd=$end_r '{if($3 >= rstrt && $3<= rnd ) print $2}' $MAP|head -1`
    reg_end=`awk -v rstrt=$start_r -v rnd=$end_r '{if($3 >= rstrt && $3<= rnd ) print $2}' $MAP|tail -1`
    echo "${reg_start} ${reg_end}" > ${reg_file}.keepsnps
done

# fgrep -v -f $MATCH.excluded_regions $MATCH
# done < <(awk -v max_sd=${max_d} '$2>=max_sd' CEU.22.non_missing.match.shareDens) > excluded_regions.txt

# MATCH=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/IBD/GERMLINE/CHR22/CEU.22.non_missing.match
# MAP=/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/IBD_INPUT/GERMLINE/CEU.22.non_missing.map
# start_r=`awk -v max_sd=${max_d} '$2>=max_sd' $MATCH.shareDens | awk '{print $1*0.5}'`
# end_r=`awk -v max_sd=${max_d} '$2>=max_sd' $MATCH.shareDens | awk '{print ($1+0.5)*0.5}'`

# avg=`awk '{sum+=$2} END { printf "%.8f", sum/NR}' CEU.22.non_missing.match.shareDens`
# std=`awk -v mu=${avg} '{sumsq += ($2 - mu)^2} END { printf "%.8f", sqrt((sumsq)/NR)}' CEU.22.non_missing.match.shareDens` 
# start_r=`awk -v max_sd=${max_d} '$2>=max_sd' CEU.22.non_missing.match.shareDens | awk '{print $1*0.5}'`
# end_r=`awk -v max_sd=${max_d} '$2>=max_sd' CEU.22.non_missing.match.shareDens | awk '{print ($1+0.5)*0.5}'`


# numero medio coppie di individui vs size in cM -> lo si fa tra unrelated e tra related e vedere cosa succede
