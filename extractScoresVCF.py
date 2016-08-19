#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher [at] uw.edu
:Date: *23.03.2014

:Description:

This script reads VCF from standard input and writes the retrieved 
scores to standard output, while generating a file with variants 
that could not be found in the downloaded file. Use it with a 
gzip compressed VCF and a downloaded score file (incl. 
downloaded index) as follows:

gunzip -c input.vcf.gz | python extractScoresVCF.py -p CADDscores.tsv.gz | gzip -c > scores.tsv.gz

Or if you gave execute permission to the script (chmod +x extractScoresVCF.py):

gunzip -c input.vcf.gz | ./extractScoresVCF.py -p CADDscores.tsv.gz | gzip -c > scores.tsv.gz

In case your input.vcf is not compressed by gzip replace 
"gunzip -c input.vcf.gz" by "cat input.vcf".

It will create an uncompressed indels_out.vcf in the same folder
which can then be submitted to the CADD webserver for scoring 
the remaining variants.

"""

import sys, os
import pysam
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p","--path", dest="path", help="Path to scored genome (default './whole_genome_SNVs_anno.tsv.gz')",default="./whole_genome_SNVs_anno.tsv.gz")
parser.add_option("--indels_out", dest="indels_out", help="Write indels and other not found variants to file (default 'indels_out.vcf')",default="indels_out.vcf")
(options, args) = parser.parse_args()

# VCF FIELDS
fchr = 0
fpos = 1
fdbsnp = 2
fref_allele = 3
falt_allele = 4
fgeno_qual = 5
fflag = 6
finfo = 7
fformat = 8
fvalues = 9

indel_out = open(options.indels_out,'w')
indel_out.write("#CHROM\tPOS\tID\tREF\tALT\n")

full_header = ["Chrom","Pos","Ref","Anc","Alt","Type","Length","isTv","isDerived","AnnoType","Consequence","ConsScore","ConsDetail","GC","CpG","mapAbility20bp","mapAbility35bp","scoreSegDup","priPhCons","mamPhCons","verPhCons","priPhyloP","mamPhyloP","verPhyloP","GerpN","GerpS","GerpRS","GerpRSpval","bStatistic","EncExp","EncH3K27Ac","EncH3K4Me1","EncH3K4Me3","EncNucleo","EncOCC","EncOCCombPVal","EncOCDNasePVal","EncOCFairePVal","EncOCpolIIPVal","EncOCctcfPVal","EncOCmycPVal","EncOCDNaseSig","EncOCFaireSig","EncOCpolIISig","EncOCctcfSig","EncOCmycSig","Segway","tOverlapMotifs","motifDist","motifECount","motifEName","motifEHIPos","motifEScoreChng","TFBS","TFBSPeaks","TFBSPeaksMax","isKnownVariant","ESP_AF","ESP_AFR","ESP_EUR","TG_AF","TG_ASN","TG_AMR","TG_AFR","TG_EUR","minDistTSS","minDistTSE","GeneID","FeatureID","CCDS","GeneName","cDNApos","relcDNApos","CDSpos","relCDSpos","protPos","relProtPos","Dst2Splice","Dst2SplType","Exon","Intron","oAA","nAA","Grantham","PolyPhenCat","PolyPhenVal","SIFTcat","SIFTval","RawScore","PHRED"]
selfields = -2

short_header = ["Chrom","Pos","Ref","Alt","RawScore","PHRED"]

header = False

if os.path.exists(options.path) and os.path.exists(options.path+".tbi"):
  filename = options.path
  sys.stderr.write("Opening %s...\n"%(filename))
  regionTabix = pysam.Tabixfile(filename,'r')
else:
  sys.stderr.write("Require valid file with scored variants.\n")
  sys.exit()


for line in sys.stdin:
  if line.startswith('#'): continue
  fields = line.rstrip().split('\t')
  chrom = fields[fchr]
  pos = int(fields[fpos])
  lref = fields[fref_allele]
  present_alleles = set(fields[falt_allele].split(','))

  for allele in present_alleles:
    found = False
    for VEPline in regionTabix.fetch(chrom,pos-1,pos):
      vfields = VEPline.rstrip().split('\t')
      if len(vfields) < len(short_header): continue
      if len(vfields) == len(short_header):
        if (vfields[2] == lref) and (vfields[3] == allele):
          if not header:
            sys.stdout.write("## (c) University of Washington and Hudson-Alpha Institute for Biotechnology 2013. All rights reserved.\n")
            sys.stdout.write("#"+"\t".join(short_header)+"\n")
            header = True
          sys.stdout.write(VEPline+"\n")
          found = True 
          break
      else:
        if (vfields[2] == lref) and (vfields[4] == allele):
          if not header:
            sys.stdout.write("## (c) University of Washington and Hudson-Alpha Institute for Biotechnology 2013. All rights reserved.\n")
            sys.stdout.write("#"+"\t".join(full_header)+"\n")
            header = True           
          sys.stdout.write(VEPline+"\n")
          found = True 
          break
    if not found: 
      fields[fdbsnp]='.'
      indel_out.write("\t".join(fields[:5])+"\n")

indel_out.close()