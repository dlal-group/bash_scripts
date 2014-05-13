#!/software/bin/python

import numpy
import re
import gzip
import sys
import math 
import os
import scipy.stats
import random 

"""
- take as input file a list of "CHR POSITION " (test-list); mandatory no header; chromosome and position must appear in this order 
- give as output a list of sites on the chromosme passed as argument matched for allele  frequency to the test-list
- require a file with allele frequencies in the minumun format CHROM POS DAF where 
	CHROM <- chromosome 
	POS <- position 
	DAF <- alle efrequency 
	CHR POS and DAF in the header are mandatory, can  be changed in the code 
	
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chromosome=sys.argv[1] # chromosome to match 
missingvaluekeyword=sys.argv[2]
pathtoinfile=sys.argv[3] ##test-list file 
pathtodaffile=sys.argv[4] ## allele frequency file 



#~~~~~~~~~~~~~~read the test lists  file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_site={}
site_freq=[]
forbidden_sites=[]
 
for chr in range (1 , 23):
	if not chr in test_site: test_site[chr]=[]

for line in open( pathtoinfile, 'r'):
	x=line.rstrip().split()
	test_site[int(x[0])].append(x[1])
	if chr==chromosome: forbidden_sites.append(x[1])

#~~~~~~~~~~~~~~~~~find out alle frequency distribution for test sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	


# freqfile=gzip.open( pathtodaffile, 'r') #remove gzip if required 
freqfile=open( pathtodaffile, 'r') #remove gzip if required
for line in freqfile: 
	x=line.strip().split()
	if re.match('CHROM', line): 
		coldaf=x.index('DAF')
		colchr=x.index('CHROM')
		colpos=x.index('POS')
	else: 
		if x[colpos] in test_site[int(x[colchr])]: 
			if not re.search (missingvaluekeyword, line ): 
				site_freq.append(float(x[coldaf]))
				print 'test', x[colchr], x[colpos], x[coldaf]
			else: 
				print 'test', x[colchr], x[colpos], missingvaluekeyword
freqfile.close()

#print len(site_freq)

#~~~~~calculate stats and distribution parameters for test-list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nbtot=len(site_freq)
lower=min(site_freq)
upper=max(site_freq)
nbofquartiles=10 ##subdivisions in 10 quartiles increase/reduce this number to get more fine/coarse resolution 
punti=range(0, 100, nbofquartiles ) 
scores=[]
for pp in punti:  
	valore=scipy.stats.scoreatpercentile(site_freq, pp)
	scores.append(valore)
scores.append(upper)
#print scores 
quartili={}
contaquartili=0
for item in scores[:-1]:
	quartili[contaquartili]=(item, scores[scores.index(item)+1])
	contaquartili+=1
frak=[]
for quar in quartili:
	count=0
	for val in site_freq:
		if (( val>=quartili[quar][0]) and (val<quartili[quar][1]) ): count +=1
	frak.append(count)	

#~~~~~~~keep information for allel frequencies on chromosome "outchromosome" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dic_freqchr={}
for quar in quartili: dic_freqchr[quar]=[]
listaposizioni=[]
diz_frequenze={}

# freqfile=gzip.open( pathtodaffile, 'r') #remove gzip if required 
freqfile=open( pathtodaffile, 'r') 
for line in freqfile: 
	x=line.rstrip().split('\t')
	if re.match('CHROM', line): 
		coldaf=x.index('DAF')
		colchr=x.index('CHROM')
		colpos=x.index('POS')
	else: 
		if not re.search('na', line):
			if x[colchr]==chromosome:  #change here if you want to have results for more than one chromosome but there might be memory issues 
				posizione=x[colpos]; frequenza=x[coldaf]
				if not  (x[colpos] in forbidden_sites): 
					listaposizioni.append(posizione)
					diz_frequenze[posizione]=frequenza

freqfile.close()

#print len(listaposizioni) 

sitiutili=set(listaposizioni)
for posizioneutile in sitiutili: 
	for quar in quartili: 
		if (( float(diz_frequenze[posizioneutile])>=quartili[quar][0]) and (float(diz_frequenze[posizioneutile])<quartili[quar][1]) ):
			dic_freqchr[quar].append((posizioneutile, diz_frequenze[posizioneutile]))
			
#print dic_freqchr
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for quar in dic_freqchr: 
	countmatched=0

	maxnuberof=frak.pop(0)
	#print maxnuberof
	random.shuffle(dic_freqchr[quar])
	for pair in  dic_freqchr[quar]:
		countmatched+=1
		#print countmatched , pair, pair[0], pair[1] 
		if countmatched <= maxnuberof: 
			print 'match', chromosome, pair[0], pair[1] 
		#else: break 
