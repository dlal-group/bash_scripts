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
T=`head -n 1 $MATCH | awk '{a=substr($2,length($2)-1,2); if (a==".0" || a==".1") print 1; else print 0; }'`

if [ "$T" -eq "1" ]; then

cat $MATCH | sed 's/\.[01]//' | sed 's/\.[01]//'

else

cat $MATCH

# fi | awk '{id1=tolower(substr($2,length($2),1)); id2=tolower(substr($4,length($4),1)); if ((id1=="a"||id1=="b")&&(id2=="a"||id2=="b")) print; }' \
fi |
cat $MATCH |  awk -v first=$first -v last=$last -v resolution=$resolution -v map=$MAP -v N=$N -v minDens=$minDens '
BEGIN{
  while (getline < map) {
    gen[$2]=$3;
  }
  begin=int(first/resolution);
  finish=int(last/resolution);
}
#snp usati per calcolare ibd diviso lunghezza segmento
$10/$11>=minDens { 
  cnt++;
  # if (cnt%10000==0) print > "/dev/stderr"; 
  start=int(gen[$8]/resolution);
  end=int(gen[$9]/resolution);
  start_w[cnt]=$6
  end_w[cnt]=$7
  # print start,end > "/dev/stderr";
  for (i=start; i<=end; i++) {
    dens[i]=dens[i]+1;
    # print start, i, end;
  }
}
END{
  for (i=begin; i<=finish; i++) {
    print i*resolution "\t" 0+dens[i]/(N*(N-1)/2 - N/2) "\t" start_w[i] "\t" end_w[i];
  } 
}' > $MATCH.shareDens

# prendo histogrammi tolgo regioni che deviano dalla densità media di un 5sd
# numero medio coppie di individui vs size in cM -> lo si fa tra unrelated e tra related e vedere cosa succede
  


