#!/usr/bin/env python
import re
import gzip
import sys
import numpy
import scipy
import os
# import pdb
"""
*** USAGE ***
ec_dac_mac_1kg.py individuals.list (from vcf header) vcfinput outfile

"""
# individualssamplelist="/netapp/dati/INGI_WGS/18112015/CARL/12112015_FILTERED_REL/30092016_CONV_ID/CARL_wgs.samples"
individualssamplelist=sys.argv[1]
# inputvcf="/netapp/dati/INGI_WGS/18112015/CARL/12112015_FILTERED_REL/30092016_CONV_ID/1.vcf.gz"
inputvcf=sys.argv[2]
# outprefix="FVG.chr22.tab"
outfile=sys.argv[3]

#~~~~~~~~~~~~~ output files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~ routines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def frequencies_anc_known_confidence(genlist, ref, alt, anc) :
	R=A=M=0
	D=0; confidence=''
	swap={'a':'A', 't':'T', 'g':'G', 'c':'C'}
	for item in genlist:
		x=item.split(':')
		r=x[0].count('0'); a=x[0].count('1'); m=x[0].count('.')
		R=R+r; A=A+a; M=M+m #here we're not taking in account that a site can be missing for an entire population
		if ref==anc: D=A
		elif alt==anc: D=R
		else:
			upperanc=swap[anc]
			if ref==upperanc: D=A
			elif alt==upperanc: R=A

	if (R==0 and A==0):
		mac="NA"
	else:
		mac=min(R,A )

	return R, A, M, D, mac  #, confidence

def frequencies_mac(genlist, ref, alt) :
	R=A=M=0
	for item in genlist:
		x=item.split(':')
		r=x[0].count('0'); a=x[0].count('1'); m=x[0].count('.')
		R=R+r; A=A+a; M=M+m #here we're not taking in account that a site can be missing for an entire population
	if (R==0 and A==0):
		mac="NA"
	else:
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
out_file=open(outfile , 'w')
print >> out_file, '#CHROM\tPOZ\tPOS\tID\tREF\tALT\tAA\tREC\tALC\tDAC\tMAC\tDAF\tMAF'

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

		#skip all multiallelic sites
		if !re.match(",", alt):
			for ii in infosplit:
				if re.match('AA=', ii): 
					# pdb.set_trace()
					aasplitted=ii.split('=')
					iisplitted=aasplitted[1].split('|')
					ancestralallele=iisplitted[0]
					info_ancestral=ii
					break
				else:
					ancestralallele='.'
					info_ancestral='.'
			

			temporary_genotypes=[]		

			for id in sampleind_index: temporary_genotypes.append(z[id])
			#print temporary_genotypes

			#print '##', ref, alt, ancestralallele
			if not (re.search('\.', ancestralallele) or re.search('-' , ancestralallele) or re.search('N', ancestralallele)):
				# pdb.set_trace()
				if ((ref == ancestralallele)	or (alt == ancestralallele) or (ref == ancestralallele.upper())	or (alt == ancestralallele.upper())):
					alleles_count=frequencies_anc_known_confidence(temporary_genotypes, ref, alt, ancestralallele.upper())
					rac=alleles_count[0]; alc=alleles_count[1]; dac=alleles_count[3]; mac=alleles_count[4]
				else:
					alleles_count=frequencies_mac(  temporary_genotypes, ref, alt) 
					rac=alleles_count[0]; alc=alleles_count[1]; dac='NA'; mac=alleles_count[3]

			else: 
				alleles_count=frequencies_mac(  temporary_genotypes, ref, alt) 
				rac=alleles_count[0]; alc=alleles_count[1]; dac='NA'; mac=alleles_count[3]

			# for item in  [chr,poz,position, vid, ref, alt, infofield]:  #z[0:8]: 
			for item in  [chr,poz,position, vid, ref, alt, info_ancestral]:  #z[0:8]: 
				print >> out_file, '%s\t' %(item),

			print >> out_file, '%s\t%s\t%s\t%s\t' %(rac,alc, dac, mac),

			if not dac=='NA': print >> out_file, '%.6f\t' %(int(dac)/float(len(sampleind_index)*2)),
			else: print >> out_file, 'NA\t',

			if not mac=='NA': print >> out_file, int(mac)/float(len(sampleind_index)*2)
			else: print >> out_file, 'NA\t'


