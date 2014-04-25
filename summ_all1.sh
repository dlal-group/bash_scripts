#!/usr/local/bin/bash

head -1 GLU.chr1.tab.assoc.txt >> GLU.all_chr.assoc.txt;for i in {1..22} X ; do echo "GLU.chr${i}.tab.assoc.txt";fgrep -v "beta" GLU.chr${i}.tab.assoc.txt >> GLU.all_chr.assoc.txt;done;awk '{print $1,$2,$5,$6,$9}' GLU.all_chr.assoc.txt | gzip -c > GLU.all_chr.assoc.txt.gz