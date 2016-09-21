import sys
import gzip
import re
import os
import pdb

#define a function to calculate the overall concordance and NRD
def g_concordance(table):
	overall_conc_den=sum([table['SRR_ORR'],table['SRR_ORA'],table['SRR_OAA'],table['SRA_ORR'],table['SRA_ORA'],table['SRA_OAA'],table['SAA_ORR'],table['SAA_ORA'],table['SAA_OAA']])
	overall_conc_num= sum([table['SRR_ORR'],table['SRA_ORA'],table['SAA_OAA']])
	if overall_conc_den!= 0:
		overall_concordance=float(overall_conc_num)/float(overall_conc_den)
	else:
		overall_concordance='NA'
		
	return overall_concordance

#define a function to calculate the complementary allele
def comp_ref(allele):
	if allele == 'A':
		complementary='T'
	elif allele == 'T':
		complementary='A'
	elif allele == 'C':
		complementary='G'
	elif allele == 'G':
		complementary='C'
	elif allele == '0':
		complementary='.'
		
	return complementary

#define a function to calculate the complementary genotype
def comp_geno(genotype):
	a1=genotype.strip().split(' ')[0]
	a2=genotype.strip().split(' ')[1]
	c_a1=comp_ref(a1)
	c_a2=comp_ref(a2)
	g_complementary=(' ').join([c_a1,c_a2])
	return g_complementary

def nr_discordance(table):
	non_ref_den=sum([table['SRR_ORA'],table['SRR_OAA'],table['SRA_ORR'],table['SRA_ORA'],table['SRA_OAA'],table['SAA_ORR'],table['SAA_ORA'],table['SAA_OAA']])
	non_ref_num=sum([table['SRA_ORR'],table['SRA_OAA'],table['SAA_ORR'],table['SAA_ORA'],table['SRR_ORA'],table['SRR_OAA']])
	#if the site is monomorphic for the ref allele,we cant calculate the NRD, so set the value at NA
	if non_ref_den!= 0:
		non_ref_disc=float(non_ref_num)/float(non_ref_den)
	else:
		non_ref_disc='NA'
	return non_ref_disc

def print_table(table):
	print ' \t|RR\t|RA\t|AA\t|NA'
	print '---\t---\t---\t---\t---\t'
	print 'RR\t|%s\t|%s\t|%s\t|%s' %(table['SRR_ORR'],table['SRR_ORA'],table['SRR_OAA'],table['SRR_ONA'])
	print '---\t---\t---\t---\t---\t'
	print 'RA\t|%s\t|%s\t|%s\t|%s' %(table['SRA_ORR'],table['SRA_ORA'],table['SRA_OAA'],table['SRA_ONA'])
	print '---\t---\t---\t---\t---\t'
	print 'AA\t|%s\t|%s\t|%s\t|%s' %(table['SAA_ORR'],table['SAA_ORA'],table['SAA_OAA'],table['SAA_ONA'])
	print '---\t---\t---\t---\t---\t'
	print 'NA\t|%s\t|%s\t|%s\t|%s' %(table['SNA_ORR'],table['SNA_ORA'],table['SNA_OAA'],table['SNA_ONA'])
	print '---\t---\t---\t---\t---\t'
	return

def table_counter(first_table_el,second_table_el):
	el_cohord={}
	if first_table_el==0:
		if second_table_el==0:
			el_cohord['SRR_ORR']+=1
		elif second_table_el==1:
			el_cohord['SRA_ORR']+=1
		elif second_table_el==2:
			el_cohord['SAA_ORR']+=1
		elif second_table_el=='NA':
			el_cohord['SNA_ORR']+=1
	elif first_table_el==1:
		if second_table_el==0:
			el_cohord['SRR_ORA']+=1
		elif second_table_el==1:
			el_cohord['SRA_ORA']+=1
		elif second_table_el==2:
			el_cohord['SAA_ORA']+=1
		elif second_table_el=='NA':
			el_cohord['SNA_ORA']+=1
	elif first_table_el==2:
		if second_table_el==0:
			el_cohord['SRR_OAA']+=1
		elif second_table_el==1:
			el_cohord['SRA_OAA']+=1
		elif second_table_el==2:
			el_cohord['SAA_OAA']+=1
		elif second_table_el=='NA':
			el_cohord['SNA_OAA']+=1
	elif first_table_el=='NA':
		if second_table_el==0:
			el_cohord['SRR_ONA']+=1
		elif second_table_el==1:
			el_cohord['SRA_ONA']+=1
		elif second_table_el==2:
			el_cohord['SAA_ONA']+=1
		elif second_table_el=='NA':
			el_cohord['SNA_ONA']+=1
	return el_cohord

#TODO:write down all command used to format all files involved in the script!
#OTHER_inputfile_name='ALL_VB_20151013_gwas_definitivo.new.s1_gwas_chr1.tped'
#SEQ_inputfile_name='ALL_VB_20151013_seq_definitivo.new.s2_seq.flipped_chr1.tped'
#ref_inputfile_name='ref_discordance_table.txt.single'
#indiv_inputfile_name='ALL_VB_20151013_seq_definitivo.new.s2_seq.flipped_chr1.tfam'
#chr=1
OTHER_inputfile_name=sys.argv[1]
SEQ_inputfile_name=sys.argv[2]
ref_inputfile_name=sys.argv[3]
chr=sys.argv[4]
indiv_inputfile_name=sys.argv[5]
#SAMPLE_geno_inputfile_name=sys.argv[6]
#SAMPLE_SEQ_inputfile_name=sys.argv[7]

#set the reference allele (the minor allele) using a list previously extracted from plink
ref_inputfile=open(ref_inputfile_name, 'r')
#read the file from 370K,previously converted and formatted: we want it in tped format and tab spaced!!!!so we use:
#plink --noweb --file <filename> --recode --transpose --tab --out <outname>
OTHER_inputfile=open(OTHER_inputfile_name, 'r')
#read the file from SEQ,previously converted and formatted: we want it in tped format and tab spaced!!!!so we use:
#plink --noweb --file <filename> --recode --transpose --tab --out <outname>
SEQ_inputfile=open(SEQ_inputfile_name, 'r')
#open the file with indiv list
indiv_inputfile=open(indiv_inputfile_name, 'r')
#read the file from 370k for calculation by sample,previously converted and formatted: we want it in ped format BUT tab spaced!!!!so we use:
#plink --noweb --file <filename> --recode --tab --out <outname>
#SAMPLE_geno_inputfile=open(SAMPLE_geno_inputfile_name, 'r')
##read the file from SEQ for calculation by sample,previously converted and formatted: we want it in ped format BUT tab spaced!!!!so we use:
##plink --noweb --file <filename> --recode --tab --out <outname>
#SAMPLE_SEQ_inputfile=open(SAMPLE_SEQ_inputfile_name, 'r')

#extract the indiv list
global indiv_list
indiv_list={}
pos=0
for indiv_line in indiv_inputfile:
	x=indiv_line.split()
	indiv_list[x[0]]=pos
	pos+=1
#set classes to calculate discordance
classes=['SRR_ORR','SRR_ORA','SRR_OAA','SRR_ONA','SRA_ORR','SRA_ORA','SRA_OAA','SRA_ONA','SAA_ORR','SAA_ORA','SAA_OAA','SAA_ONA','SNA_ORR','SNA_ORA','SNA_OAA','SNA_ONA']

#work with each chr
print chr
ref_all={}
#set the ref for each position for each chr, so we can count and translate them
for snp_line in ref_inputfile:
	x=snp_line.split()
	if x[0]==str(chr):
		# ref_all[x[1]]=x[3] , use the rs id instead of the position
		ref_all[x[2]]=x[3]

print len(ref_all)

###WORK BY SITE!!!######
#convert the SEQ files to RR/AA useful format
global seq_conversion
seq_conversion={}
rs_lines=SEQ_inputfile.readlines()

seq_RR_out=open('seq_RR_chr' + str(chr) + '.txt', 'w')
print >> seq_RR_out, "POS ALL1 ALL2 REF CONV"

seq_het_out=open('seq_het_chr' + str(chr) + '.txt', 'w')
print >> seq_het_out, "POS ALL1 ALL2 REF CONV"

seq_AA_out=open('seq_AA_chr' + str(chr) + '.txt', 'w')
print >> seq_AA_out, "POS ALL1 ALL2 REF CONV"

seq_NA_out=open('seq_NA_chr' + str(chr) + '.txt', 'w')
print >> seq_NA_out, "POS ALL1 ALL2 REF CONV"

#for rs_line in SEQ_inputfile:
for rs_line in rs_lines:
	#rs_line=rs_lines[0]
	#print (seq_conversion)
	snp=rs_line.strip().split('\t')
	# current_ref=ref_all[snp[3]] , switch to rsID instead of position checking 
	ref_snp = snp[1]
	current_ref=ref_all[ref_snp]
	global conversion_seq
	conversion=[]
	for i in xrange (4,len(snp)):
		if current_ref in snp[i] or current_ref in comp_geno(snp[i]):
			if snp[i]==" ".join([current_ref,current_ref]) or snp[i]==" ".join([comp_ref(current_ref),comp_ref(current_ref)]):
				#major allele monomorphic
				conv_value=0
				print >> seq_RR_out, "%s %s %s %s" %(ref_snp,snp[i],current_ref,conv_value)
			else:
				#heterozygote
				conv_value=1
				print >> seq_het_out, "%s %s %s %s" %(ref_snp,snp[i],current_ref,conv_value)
		elif snp[i]== '0 0':
			#missing genotype
			conv_value='NA'
			print >> seq_NA_out, "%s %s %s %s" %(ref_snp,snp[i],current_ref,conv_value)
		elif current_ref not in snp[i] and current_ref not in comp_geno(snp[i]):
			#minor allele monomorphic
			conv_value=2
			print >> seq_AA_out, "%s %s %s %s" %(ref_snp,snp[i],current_ref,conv_value)
		conversion.append(conv_value)
		#print 'TEST ON SEQ'
		#print conversion
	seq_conversion[ref_snp]=conversion
seq_RR_out.close()
seq_het_out.close()
seq_AA_out.close()
#sys.stdout=sys.stdout
print "Sono alla conversione della sequenza"

#now the same for the other set genotypes
gwas_RR_out=open('gwas_RR_chr' + str(chr) + '.txt', 'w')
print >> gwas_RR_out, "POS ALL1 ALL2 REF CONV"

gwas_het_out=open('gwas_het_chr' + str(chr) + '.txt', 'w')
print >> gwas_het_out, "POS ALL1 ALL2 REF CONV"

gwas_AA_out=open('gwas_AA_chr' + str(chr) + '.txt', 'w')
print >> gwas_AA_out, "POS ALL1 ALL2 REF CONV"

gwas_NA_out=open('gwas_NA_chr' + str(chr) + '.txt', 'w')
print >> gwas_NA_out, "POS ALL1 ALL2 REF CONV"

global gwas_conversion
gwas_conversion={}
gwas_lines=OTHER_inputfile.readlines()
for gwas_line in gwas_lines:
	snp=gwas_line.strip().split('\t')
	ref_snp = snp[1] #switch to rsId istead of position for cheking: change this back to ref_snp = snp[3] to use position again
	current_ref=ref_all[ref_snp]
	global conversion_gwas
	conversion_gwas=[]
	for i in xrange (4,len(snp)):
		if current_ref in snp[i] or current_ref in comp_geno(snp[i]):
			if snp[i]==" ".join([current_ref,current_ref]) or snp[i]==" ".join([comp_ref(current_ref),comp_ref(current_ref)]):
				#major allele monomorphic
				g_conv_value=0
				print >> gwas_RR_out, "%s %s %s %s" %(ref_snp,snp[i],current_ref,g_conv_value)
			else:
				#heterozygote
				g_conv_value=1
				print >> gwas_het_out, "%s %s %s %s" %(ref_snp,snp[i],current_ref,g_conv_value)
		elif snp[i]== '0 0':
			#missing genotype
			g_conv_value='NA'
			print >> gwas_NA_out, "%s %s %s %s" %(ref_snp,snp[i],current_ref,g_conv_value)
		elif current_ref not in snp[i] and current_ref not in comp_geno(snp[i]):
			#minor allele monomorphic
			g_conv_value=2
			print >> gwas_AA_out, "%s %s %s %s" %(ref_snp,snp[i],current_ref,g_conv_value)
		conversion_gwas.append(g_conv_value)
#		print 'TEST ON GWAS'
#		print conversion
	gwas_conversion[ref_snp]=conversion_gwas
#print gwas_conversion
gwas_RR_out.close()
gwas_het_out.close()
gwas_AA_out.close()
#NOW take the two converted dictionaries and scan foreach snp to count how many RR/RA/AA combination we have
global ref_count
ref_count={}

#count for each postion
for current_pos in ref_all.keys():
	global snp_cohord
	snp_cohord={}
	#initialize the cohord dictionary.We need to do it for each snp
	for cohord in classes: #63799
		snp_cohord[cohord]=0
#	print current_pos
	current_gwas=gwas_conversion[current_pos]
	current_seq=seq_conversion[current_pos]
	#first check if the two list are equal by size
	if len(current_gwas)==len(current_seq):
		# print "Same length check passed!!"
		#FIXME:use a more elegant way to do it.
		for code_index in xrange (0,len(current_seq)):
#			snp_cohord=table_counter(current_gwas[code_index],current_seq[code_index],classes)
			if current_gwas[code_index]==0:
				if current_seq[code_index]==0:
					snp_cohord['SRR_ORR']+=1
				elif current_seq[code_index]==1:
					snp_cohord['SRA_ORR']+=1
				elif current_seq[code_index]==2:
					snp_cohord['SAA_ORR']+=1
				elif current_seq[code_index]=='NA':
					snp_cohord['SNA_ORR']+=1
			elif current_gwas[code_index]==1:
				if current_seq[code_index]==0:
					snp_cohord['SRR_ORA']+=1
				elif current_seq[code_index]==1:
					snp_cohord['SRA_ORA']+=1
				elif current_seq[code_index]==2:
					snp_cohord['SAA_ORA']+=1
				elif current_seq[code_index]=='NA':
					snp_cohord['SNA_ORA']+=1
			elif current_gwas[code_index]==2:
				if current_seq[code_index]==0:
					snp_cohord['SRR_OAA']+=1
				elif current_seq[code_index]==1:
					snp_cohord['SRA_OAA']+=1
				elif current_seq[code_index]==2:
					snp_cohord['SAA_OAA']+=1
				elif current_seq[code_index]=='NA':
					snp_cohord['SNA_OAA']+=1
			elif current_gwas[code_index]=='NA':
				if current_seq[code_index]==0:
					snp_cohord['SRR_ONA']+=1
				elif current_seq[code_index]==1:
					snp_cohord['SRA_ONA']+=1
				elif current_seq[code_index]==2:
					snp_cohord['SAA_ONA']+=1
				elif current_seq[code_index]=='NA':
					snp_cohord['SNA_ONA']+=1
		ref_count[current_pos]=snp_cohord
	else:
		print "ERROR!!Different length for pos %s in chr %s" (current_pos, chr)


#Extract info for each site
snp_table_file=open('site_concordance_discordance_chr' + str(chr) + '.txt', 'w')
print >> snp_table_file, '\tChr %s\n' %(chr)
allsnp_table_file=open('site_concordance_discordance_table_chr' + str(chr) + '.txt', 'w')
print >> allsnp_table_file, 'CHR\tPOS\tOGC\tNRD'

for pos in ref_count.keys():
	snp_table=ref_count[pos]
	print >> snp_table_file, '\n\tSet1 vs Set2 (pos: %s)\n' %(pos)
	print >> snp_table_file, ' \t|RR\t|RA\t|AA\t|NA'
	print >> snp_table_file, '---\t---\t---\t---\t---\t'
	print >> snp_table_file, 'RR\t|%s\t|%s\t|%s\t|%s' %(snp_table['SRR_ORR'],snp_table['SRR_ORA'],snp_table['SRR_OAA'],snp_table['SRR_ONA'])
	print >> snp_table_file, '---\t---\t---\t---\t---\t'
	print >> snp_table_file, 'RA\t|%s\t|%s\t|%s\t|%s' %(snp_table['SRA_ORR'],snp_table['SRA_ORA'],snp_table['SRA_OAA'],snp_table['SRA_ONA'])
	print >> snp_table_file, '---\t---\t---\t---\t---\t'
	print >> snp_table_file, 'AA\t|%s\t|%s\t|%s\t|%s' %(snp_table['SAA_ORR'],snp_table['SAA_ORA'],snp_table['SAA_OAA'],snp_table['SAA_ONA'])
	print >> snp_table_file, '---\t---\t---\t---\t---\t'
	print >> snp_table_file, 'NA\t|%s\t|%s\t|%s\t|%s' %(snp_table['SNA_ORR'],snp_table['SNA_ORA'],snp_table['SNA_OAA'],snp_table['SNA_ONA'])
	print >> snp_table_file, '---\t---\t---\t---\t---\t'
	print >> snp_table_file, '\n'
	print >> allsnp_table_file, '%s\t%s\t%s\t%s' %(chr,pos,g_concordance(snp_table),nr_discordance(snp_table))
#####DONE BY SITE#####

#############################################
#now calculate the discordance per individual
global by_sample_count
by_sample_count={}
global actual_sample_gwas
actual_sample_gwas={}
global actual_sample_seq
actual_sample_seq={}

#count for each sample
print indiv_list

for current_sample in sorted(indiv_list.keys()):
#	pdb.set_trace()
	global ind_cohord
	ind_cohord={}
	#initialize the cohord dictionary.We need to do it for each snp
	for cohord in classes: #63799
		ind_cohord[cohord]=0
	print 'Current sample'
	print current_sample
	#create the genotype list for the current sample
	#the lists have to be in the same order
	current_sample_gwas=[]
	for current_pos in sorted(gwas_conversion.keys(),key=str.lower):
#		pdb.set_trace()
		sample_gwas=gwas_conversion[current_pos][indiv_list[current_sample]]
		current_sample_gwas.append(sample_gwas)

	actual_sample_gwas[current_sample]=current_sample_gwas
#	print 'sample GWAS'
#	print actual_sample_gwas[current_sample]
	current_sample_seq=[]
	for current_pos in sorted(seq_conversion.keys(),key=str.lower):
#		pdb.set_trace()
		sample_seq=seq_conversion[current_pos][indiv_list[current_sample]]
		current_sample_seq.append(sample_seq)
		
	actual_sample_seq[current_sample]=current_sample_seq
#	print 'sample SEQ'
#	print actual_sample_seq[current_sample]

 	#extract the actual sample genotype
	current_gwas=actual_sample_gwas[current_sample]
	current_seq=actual_sample_seq[current_sample]
#	pdb.set_trace()
	#first check if the two list are equal by size
	if len(current_gwas)==len(current_seq):
		print "Same length check passed!!"
		#FIXME:use a more elegant way to do it.
		for code_index in range (0,len(current_seq)):
			if current_gwas[code_index]==0:
				if current_seq[code_index]==0:
					ind_cohord['SRR_ORR']+=1
				elif current_seq[code_index]==1:
					ind_cohord['SRA_ORR']+=1
				elif current_seq[code_index]==2:
					ind_cohord['SAA_ORR']+=1
				elif current_seq[code_index]=='NA':
					ind_cohord['SNA_ORR']+=1
			elif current_gwas[code_index]==1:
				if current_seq[code_index]==0:
					ind_cohord['SRR_ORA']+=1
				elif current_seq[code_index]==1:
					ind_cohord['SRA_ORA']+=1
				elif current_seq[code_index]==2:
					ind_cohord['SAA_ORA']+=1
				elif current_seq[code_index]=='NA':
					ind_cohord['SNA_ORA']+=1
			elif current_gwas[code_index]==2:
				if current_seq[code_index]==0:
					ind_cohord['SRR_OAA']+=1
				elif current_seq[code_index]==1:
					ind_cohord['SRA_OAA']+=1
				elif current_seq[code_index]==2:
					ind_cohord['SAA_OAA']+=1
				elif current_seq[code_index]=='NA':
					ind_cohord['SNA_OAA']+=1
			elif current_gwas[code_index]=='NA':
				if current_seq[code_index]==0:
					ind_cohord['SRR_ONA']+=1
				elif current_seq[code_index]==1:
					ind_cohord['SRA_ONA']+=1
				elif current_seq[code_index]==2:
					ind_cohord['SAA_ONA']+=1
				elif current_seq[code_index]=='NA':
					ind_cohord['SNA_ONA']+=1
		by_sample_count[current_sample]=ind_cohord
	else:
		print "ERROR!!Different length for pos %s in chr %s" (current_pos, chr)
#	pdb.set_trace()
#Extract info for each site
indiv_table_file=open('sample_concordance_discordance_chr' + str(chr) + '.txt', 'w')
print >> indiv_table_file, '\tChr %s\n' %(chr)
allindiv_table_file=open('sample_concordance_discordance_table_chr' + str(chr) + '.txt', 'w')
print >> allindiv_table_file, 'CHR\tSAMPLE\tOGC\tNRD\tSRR_ORR\tSRR_ORA\tSRR_OAA\tSRR_ONA\tSRA_ORR\tSRA_ORA\tSRA_OAA\tSRA_ONA\tSAA_ORR\tSAA_ORA\tSAA_OAA\tSAA_ONA\tSNA_ORR\tSNA_ORA\tSNA_OAA\tSNA_ONA'

for sample in by_sample_count.keys():
	sample_table=by_sample_count[sample]
	print sample_table
	print >> indiv_table_file, '\n\tSet1 vs Set2 (sample: %s)\n' %(sample)
	print >> indiv_table_file, ' \t|RR\t|RA\t|AA\t|NA'
	print >> indiv_table_file, '---\t---\t---\t---\t---\t'
	print >> indiv_table_file, 'RR\t|%s\t|%s\t|%s\t|%s' %(sample_table['SRR_ORR'],sample_table['SRR_ORA'],sample_table['SRR_OAA'],sample_table['SRR_ONA'])
	print >> indiv_table_file, '---\t---\t---\t---\t---\t'
	print >> indiv_table_file, 'RA\t|%s\t|%s\t|%s\t|%s' %(sample_table['SRA_ORR'],sample_table['SRA_ORA'],sample_table['SRA_OAA'],sample_table['SRA_ONA'])
	print >> indiv_table_file, '---\t---\t---\t---\t---\t'
	print >> indiv_table_file, 'AA\t|%s\t|%s\t|%s\t|%s' %(sample_table['SAA_ORR'],sample_table['SAA_ORA'],sample_table['SAA_OAA'],sample_table['SAA_ONA'])
	print >> indiv_table_file, '---\t---\t---\t---\t---\t'
	print >> indiv_table_file, 'NA\t|%s\t|%s\t|%s\t|%s' %(sample_table['SNA_ORR'],sample_table['SNA_ORA'],sample_table['SNA_OAA'],sample_table['SNA_ONA'])
	print >> indiv_table_file, '---\t---\t---\t---\t---\t'
	print >> indiv_table_file, '\n'
	#1)Overall genotype concordance
	sample_overall_concordance= g_concordance(sample_table)
	#2)Non reference discordance rate
	sample_nonref_disc= nr_discordance(sample_table)
	print >> indiv_table_file,'Overall concordance: %s' %(sample_overall_concordance)
	print >> indiv_table_file,'Non ref-discordance: %s' %(sample_nonref_disc)
	print >> allindiv_table_file, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(chr,sample,g_concordance(sample_table),nr_discordance(sample_table),sample_table['SRR_ORR'],sample_table['SRR_ORA'],sample_table['SRR_OAA'],sample_table['SRR_ONA'],sample_table['SRA_ORR'],sample_table['SRA_ORA'],sample_table['SRA_OAA'],sample_table['SRA_ONA'],sample_table['SAA_ORR'],sample_table['SAA_ORA'],sample_table['SAA_OAA'],sample_table['SAA_ONA'],sample_table['SNA_ORR'],sample_table['SNA_ORA'],sample_table['SNA_OAA'],sample_table['SNA_ONA'])


#calculate the TOTAL discordance for the chr based on per site count
global total_count
total_count={}
for disc_class in classes:
	total_count[disc_class]=0
	for position in ref_count:
		#now update the counter for the current chr
		total_count[disc_class]+=ref_count[position][disc_class]

#now calculate the formula for:
#1)Overall genotype concordance
overall_concordance= g_concordance(total_count)
#2)Non reference discordance rate
nonref_disc= nr_discordance(total_count)
#3)Overall genotype concordance
#overall_NA= sum([total_count['SRR_ORR'],total_count['SRA_ORA'],total_count['SAA_OAA']])/(sum([total_count['SRR_ORR'],total_count['SRA_ORA'],total_count['SAA_OAA'],total_count['SRA_ORR'],total_count['SRA_ORA'],total_count['SRA_OAA'],total_count['SAA_ORR'],total_count['SAA_ORA'],total_count['SAA_OAA']]))
#4)Alleles missed in seq
#nonref_disc= sum([total_count['SRR_ORR'],total_count['SRA_ORA'],total_count['SAA_OAA']])/(sum([total_count['SRR_ORR'],total_count['SRA_ORA'],total_count['SAA_OAA'],total_count['SRA_ORR'],total_count['SRA_ORA'],total_count['SRA_OAA'],total_count['SAA_ORR'],total_count['SAA_ORA'],total_count['SAA_OAA']]))
#5)Alleles overcalled in seq

#now print everything in a fancy way!
out=open('overall_discordance_chr' + str(chr) + '.txt', 'w')
sys.stdout=out
print '\tVBI seq vs GWAS\n'
print_table(total_count)
print 'Overall concordance: %s' %(overall_concordance)
print 'Non ref-discordance: %s' %(nonref_disc)


