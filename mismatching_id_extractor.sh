#!/usr/local/bin/bash

#extract missing rsID sites from dbSNP annotation file

tabix ../../dbSNP_splitted/chr$1_dbSNP_b137.tab.gz chr $1:$2-$2 >> chr$1_dbSNP_mismatch.txt

