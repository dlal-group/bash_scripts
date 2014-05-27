#!/software/bin/python

import gzip 
import re 
import sys 
from sys import stdout
"""
usage python ec_mergeXac.py mychr myvarofinterest mycontofinterest > myfilemerged.bed

"""
chr=sys.argv[1]
xac=sys.argv[2]
contofinterest=sys.argv[3]

dic_xac={}
dic_type={}
sitelist=[]

vartype=['SNP', 'INDEL', 'MULTI', 'NA']
for vt in vartype: 
	dic_type[vt]=[]
	dic_xac[vt]={}

poplist={'EUR':['CEU', 'FIN', 'GBR', 'TSI'], 'ASN': ['CHB', 'CHS', 'JPT'], 'AFR': ['ASW', 'LWK', 'YRI'], 'CON':['AFR', 'ASN', 'EUR'], 'INGI':['CEU','TSI','VBI','FVG']} 

for pop in poplist[contofinterest]:
	for vt in vartype:  dic_xac[vt][pop]={}

	for line in gzip.open('%s.chr%s.not_fixed.tab.gz'%(pop, chr ), 'r'):
		if re.match('CHR', line): 
			y=line.split('\t')
			colofinterest=y.index(xac)
			#print colofinterest

		if re.match('\d+\t', line):
			y=line.split('\t')
			posizione=y[2]; sitelist.append(posizione)
			# if re.search('VT=SNP' ,line ): dic_type['SNP'].append(posizione); dic_xac['SNP'][pop][posizione]=y[colofinterest]
			# elif re.search('VT=INDEL' ,line ): dic_type['INDEL'].append(posizione); dic_xac['INDEL'][pop][posizione]=y[colofinterest]
			# elif re.search ('VT=SV' , line ) : dic_type['SV'].append(posizione); dic_xac['SV'][pop][posizione]=y[colofinterest]
			# else: dic_type['NA'].append(posizione); dic_xac['NA'][pop][posizione]=y[colofinterest]
			if (len(y[4]) == 1 and len(y[5]) == 1): dic_type['SNP'].append(posizione); dic_xac['SNP'][pop][posizione]=y[colofinterest]
			elif (len(y[4]) != len(y[5]) and re.search(",",y[5]) ): dic_type['MULTI'].append(posizione); dic_xac['MULTI'][pop][posizione]=y[colofinterest]
			elif (len(y[4]) != len(y[5]) and (not (re.search(",",y[5])))): dic_type['INDEL'].append(posizione); dic_xac['INDEL'][pop][posizione]=y[colofinterest]
			elif re.search ('VT=SV' , line ) : dic_type['SV'].append(posizione); dic_xac['SV'][pop][posizione]=y[colofinterest]
			else: dic_type['NA'].append(posizione); dic_xac['NA'][pop][posizione]=y[colofinterest]




stdout.write('#CHR\tPOZ\tPOS\tVT')
# for pop in poplist[contofinterest]: print '%s' %(pop),
for pop in poplist[contofinterest]: stdout.write('\t%s' % pop)
print '\r'

for site in set(sitelist):
	for vt in vartype:
		if site in dic_xac[vt][poplist[contofinterest][0]]:  
			# print '%s\t%s\t%s\t%s' %(chr,int(site)-1,site,vt)
			stdout.write('%s\t%s\t%s\t%s' %(chr,int(site)-1,site,vt))
			for pop in poplist[contofinterest]: 
				# print '%s' %(dic_xac[vt][pop][site]),
				if site in dic_xac[vt][pop]
					stdout.write('\t%s' %(dic_xac[vt][pop][site]))
				else
					stdout.write('\tna')
			print '\r'
 
