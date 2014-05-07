#!/usr/bin/python
import re
import gzip
import sys
import numpy
import scipy
import os

"""
*** USAGE ***
hd1_dac_mac_1kg.py vcfinput individuals.list 

"""

inputvcf=sys.argv[1]
individualssamplelist=sys.argv[2]

#~~~~~~~~~~~~~ output files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~ routines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def frequencies_anc_known_confidence(  genlist, ref, alt, anc) :
        R=A=M=0
        D=0; confidence=''
        #swap={'a':'A', 't':'T', 'g':'G', 'c':'C'}
        for item in genlist:
                x=item.split(':')
                r=x[0].count('0'); a=x[0].count('1'); m=x[0].count('.')
                R=R+r; A=A+a; M=M+m
                if ref==anc: D=A
                else: D=R
                    #upperanc=swap[anc]
                    #if ref==upperanc: D=A
                    #else: D=R
	mac=min(R,A )
        return R, A, M, D, mac  #, confidence

def frequencies_mac(  genlist, ref, alt) :
	R=A=M=0
	for item in genlist:
		x=item.split(':')
		r=x[0].count('0'); a=x[0].count('1'); m=x[0].count('.')
		R=R+r; A=A+a; M=M+m
	mac=min(R,A )
	return R, A, M, mac

def remove_dups(seq):
	x = {}
    	for y in seq: x[y] = 1
    	u = x.keys(); u.sort()
    	return u 

def findpopranges(dic_index, pop, dic_indiv_pop): 
	list=[]
	for ind in dic_index: 
		if dic_indiv_pop[ind]==pop:
			list.append(dic_index[ind])
	#return  pop, min(list), max(list)+1
		list.sort()
	return  pop, list

#~~~~~~~~~~~~~ find out which samples are in the vcf file and in which columns  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampleind=[]
for line in open (individualssamplelist , 'r') : 
	x=line.rstrip().split()
	sampleind.append(x[0])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print '#CHROM\tPOZ\tPOS\tID\tREF\tALT\tINFO\tREC\tALC\tDAC\tMAC\tDAF\tMAF'

sampleind_index=[]
for line in gzip.open(inputvcf, 'r'):
	if re.match('#CHROM', line): 
		header=line.rstrip().split() 
		for indiv in sampleind:
			sampleind_index.append(header.index(indiv))
		#break 
        elif re.match('\d+\t', line):
		confidence=''
                z=line.split()
		chr=z[0]; poz=int(z[1])-1; position=z[1]; vid=z[2]; ref=z[header.index('REF')]; alt=z[header.index('ALT')]; infofield=z[header.index('INFO')]
		infosplit=infofield.split(';')
		for ii in infosplit: 
			if re.match('AA=', ii): 
				iisplitted=ii.split('=')
				ancestralallele=iisplitted[1]
		temporary_genotypes=[]		
		for id in sampleind_index: temporary_genotypes.append(z[id])
		#print temporary_genotypes

		#print '##', ref, alt, ancestralallele
		if not  (re.search('\.', ancestralallele)   or re.search('-' , ancestralallele) or re.search('N', ancestralallele)  ):
			alleles_count=frequencies_anc_known_confidence(temporary_genotypes, ref, alt, ancestralallele) 
			rac=alleles_count[0]; alc=alleles_count[1]; dac=alleles_count[3]; mac=alleles_count[4]
		else: 
			alleles_count=frequencies_mac(  temporary_genotypes, ref, alt) 
			rac=alleles_count[0]; alc=alleles_count[1]; dac='na'; mac=alleles_count[3]
		for item in  [chr,poz,position, vid, ref, alt, infofield]:  #z[0:8]: 
                        print '%s\t' %(item),
                print '%s\t%s\t%s\t%s\t' %(rac,alc, dac, mac),
                if not dac=='na': print '%.50f\t' %(int(dac)/float(len(sampleind_index)*2)),
                else: print 'na\t',
                print int(mac)/float(len(sampleind_index)*2)
	
