#Script to create a table to perform annotation with ANCESTRAL alleles on VCF data
import sys 
import gzip
import re 
import os 
import resource

usage=resource.getrusage(resource.RUSAGE_SELF)
print "mem=%s mb" %((usage[2]*resource.getpagesize())/(1024**2))

allele_name=sys.argv[1]
ancestral_name=sys.argv[2]
dbsnp_data=sys.argv[3]
outpath=sys.argv[4]

usage=resource.getrusage(resource.RUSAGE_SELF)
print "mem=%s mb" %((usage[2]*resource.getpagesize())/(1024**2))

#read the tables
allele=gzip.open(allele_name,'r')
dbsnp_table=gzip.open(dbsnp_data,'r')
#This has the form= rsID Allele someDATA

#need to use awk to match second column of SNPAncestralAllele.bcp.gz with first column of Allele.bcp.gz
#Extract info from SNPAncestralAllele: get first 2 columns and extract a uniq list of markers and ancestral alleles code
os.system("zcat %s | cut -f 1,2 | uniq | gzip -c > %s/%s.filtered.gz" %(ancestral_name,outpath,ancestral_name))

usage=resource.getrusage(resource.RUSAGE_SELF)
print "mem=%s mb" %((usage[2]*resource.getpagesize())/(1024**2))

ancestral=gzip.open("%s/%s.filtered.gz" %(outpath,ancestral_name),'r')

print "allele=%s mb" %(sys.getsizeof(allele)/(1024**2))
print "dbsnp_table=%s mb" %(sys.getsizeof(dbsnp_table)/(1024**2))
print "ancestral=%s mb" %(sys.getsizeof(ancestral)/(1024**2))

#this ha the form = rsID ANC_allele_id
usage=resource.getrusage(resource.RUSAGE_SELF)
print "mem=%s mb" %((usage[2]*resource.getpagesize())/(1024**2))


allels = {}
for allline in allele:
	y = allline.split()
	allels[y[0]] = y[1]
	if len(allels)%100000 == 0:
		usage=resource.getrusage(resource.RUSAGE_SELF)
		print "mem=%s mb" %((usage[2]*resource.getpagesize())/(1024**2))


usage2=resource.getrusage(resource.RUSAGE_SELF)
print "mem=%s mb" %((usage2[2]*resource.getpagesize())/(1024**2) )
print "allels=%s mb" %(sys.getsizeof(allels)/(1024**2))

table = {}
for ancline in ancestral:
	x = ancline.split()
	#we need to format all those lines with a double ancestor in a way like rsID anc_rs1,anc_rs2 anc_all1,anc_all2

	if len(table)%30000 == 0:
		usage=resource.getrusage(resource.RUSAGE_SELF)
		print "mem=%s mb" %((usage[2]*resource.getpagesize())/(1024**2) )
		print "table=%s mb (length %s)" %(sys.getsizeof(table)/(1024**2),len(table))

	if x[0] in table:
		table[x[0]] = [";rs".join([table[x[0]][0],x[1]]),";".join([table[x[0]][1],allels[x[1]]])]
	else:
		table[x[0]] = ["rs"+x[1],allels[x[1]]]
	
# if i%5000000 == 0:
dbsnp = {}
for snpline in dbsnp_table:
	snp = snpline.split()
	dbsnp[snp[2]] = [snp[0],snp[1],snp[3],snp[4]]
	if len(dbsnp)%3000000 == 0:
		usage=resource.getrusage(resource.RUSAGE_SELF)
		print "mem=%s mb" %((usage[2]*resource.getpagesize())/(1024**2) )
		print "dbsnp=%s mb (length %s)" %(sys.getsizeof(dbsnp)/(1024**2),len(dbsnp))

out=open('%s/Ancestral_ann_table.txt' %(outpath), 'w')
sys.stdout=out

# print "CHROM POS ID REF ALT anc_rsID ANC_ALL"
# for k in sorted(table.iterkeys(),key=int):
# 	#this will have the form: rsID => CHROM,POS,REF,ALT

# 		# if 'rs'+k in dbsnp:
# 		if 'rs'+k == snp[2]:
# 			# print "%s %s rs%s %s %s %s %s"  %(dbsnp['rs'+k][0],dbsnp['rs'+k][1],k,dbsnp['rs'+k][2],dbsnp['rs'+k][3],table[k][0],table[k][1])
# 			print "%s %s rs%s %s %s %s %s"  %(snp[0],snp[1],k,snp[3],snp[4],table[k][0],table[k][1])

print "CHROM POS ID REF ALT anc_rsID ANC_ALL"
for k in sorted(table.iterkeys(),key=int):
	#this will have the form: rsID => CHROM,POS,REF,ALT
	# for snpline in dbsnp_table:
	# 	snp = snpline.split()
	if 'rs'+k in dbsnp:
		# if 'rs'+k == snp[2]:
			# dbsnp[snp[2]] = [snp[0],snp[1],snp[3],snp[4]]
		print "%s %s rs%s %s %s %s %s"  %(dbsnp['rs'+k][0],dbsnp['rs'+k][1],k,dbsnp['rs'+k][2],dbsnp['rs'+k][3],table[k][0],table[k][1])
			# print "%s %s rs%s %s %s %s %s"  %(snp[0],snp[1],k,snp[3],snp[4],table[k][0],table[k][1])
			# dbsnp_table.seek(0)
			# break

			

# del ancestral

# print in a table format in a file
