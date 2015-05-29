#!/software/bin/python

import gzip 
import re 
import sys 
from sys import stdout
"""
usage python ec_mergeXac.py mychr myvarofinterest mycontofinterest > myfilemerged.bed

"""
file1=sys.argv[1]
file2=sys.argv[2]
xac=sys.argv[3]

dic_xac={}

for line in gzip.open(file1, 'r'):
	y=line.split('\t')
	colofinterest=y.index(xac)
	#print colofinterest

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
		for main_pop in poplist[contofinterest]:
			if site in dic_xac[vt][main_pop]:  
				# print '%s\t%s\t%s\t%s' %(chr,int(site)-1,site,vt)
				stdout.write('%s\t%s\t%s\t%s' %(chr,int(site)-1,site,vt))
				for pop in poplist[contofinterest]: 
					# print '%s' %(dic_xac[vt][pop][site]),
					if site in dic_xac[vt][pop]:
						stdout.write('\t%s' %(dic_xac[vt][pop][site]))
					else:
						stdout.write('\tna')
				print '\r'
				break
		# elif site in dic_xac[vt][poplist[contofinterest][1]]:
		# 	stdout.write('%s\t%s\t%s\t%s' %(chr,int(site)-1,site,vt))
		# 	for pop in poplist[contofinterest]: 
		# 		# print '%s' %(dic_xac[vt][pop][site]),
		# 		if site in dic_xac[vt][pop]:
		# 			stdout.write('\t%s' %(dic_xac[vt][pop][site]))
		# 		else:
		# 			stdout.write('\tna')
		# 	print '\r'
		# elif site in dic_xac[vt][poplist[contofinterest][2]]:
		# 	stdout.write('%s\t%s\t%s\t%s' %(chr,int(site)-1,site,vt))
		# 	for pop in poplist[contofinterest]: 
		# 		# print '%s' %(dic_xac[vt][pop][site]),
		# 		if site in dic_xac[vt][pop]:
		# 			stdout.write('\t%s' %(dic_xac[vt][pop][site]))
		# 		else:
		# 			stdout.write('\tna')
		# 	print '\r'
		# elif site in dic_xac[vt][poplist[contofinterest][3]]:
		# 	stdout.write('%s\t%s\t%s\t%s' %(chr,int(site)-1,site,vt))
		# 	for pop in poplist[contofinterest]: 
		# 		# print '%s' %(dic_xac[vt][pop][site]),
		# 		if site in dic_xac[vt][pop]:
		# 			stdout.write('\t%s' %(dic_xac[vt][pop][site]))
		# 		else:
		# 			stdout.write('\tna')
		# 	print '\r'
