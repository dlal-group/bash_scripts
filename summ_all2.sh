#!/usr/local/bin/bash

head -1 GLUadjBMI.chr1.tab.assoc.txt >> GLUadjBMI.all_chr.assoc.txt;for i in {1..22} X ; do echo "GLUadjBMI.chr${i}.tab.assoc.txt";fgrep -v "beta" GLUadjBMI.chr${i}.tab.assoc.txt >> GLUadjBMI.all_chr.assoc.txt;done;awk '{print $1,$2,$5,$6,$9}' GLUadjBMI.all_chr.assoc.txt | gzip -c > GLUadjBMI.all_chr.assoc.txt.gz