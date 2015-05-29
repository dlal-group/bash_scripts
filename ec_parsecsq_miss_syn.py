import sys 
import gzip
import re 

def derived_or_ancestral( ref, alt, anc) :         
	D=0; state=0 
	if ref==anc:D=alt; state="alt"
	else: D=ref; state="ref"
	return D, state  	

myvcf=sys.argv[2]
mykeyword=sys.argv[1]  # synonymus_variant  missense_variant
for line in gzip.open (myvcf, 'r'): 
	if re.search('#CHROM', line):
		header=line.split()
		
	if not re.match('#', line) :
		z=line.split()
		myref=z[header.index('REF')]; myalt=z[header.index('ALT')]; 
		myinfo=z[header.index('INFO')].split(";"); mychr=z[0];  mypos=int(z[1]) ; #print myinfo  
		w=[y for y in myinfo if re.match("AA", y) if re.search("AA=", z[header.index('INFO')] ) ] ; #print w 
		if len(w)>0 : 
			myancestral=w[0].split("=")[1].upper() ; #print myancestral
			if myancestral not in ['A', 'T', 'G', 'C']: mystate=['NA', 'NA']
			else: mystate=derived_or_ancestral( myref, myalt,  myancestral) 
		else: mystate=['NA', 'NA'] 
		mycsq=[  x for  x in myinfo if re.match("CSQ", x)] 
		if len(mycsq)>0 :

			#print mypos , mycsq 
			mytranscripts=mycsq[0].split("+") 
			for trans in mytranscripts: 
				transname=trans.split(":")[0].lstrip("CSQ=")  
				mymissense=[ x for  x in trans.split(":")  if re.match(mykeyword, x) ]
				#print mymissense
				if len(mymissense) >0 : 
					myoutline=[ mychr, mypos-1 ,mypos,  mymissense[0].split(",")[0] , mystate[0], mystate[1], transname]
					print "\t".join(map(str, myoutline) ) 
