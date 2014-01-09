#!/software/bin/python-2.7
import os

outerrdir='~/Work/bash_scripts/tempouterr'
codedir='~/Work/bash_scripts/'
scriptname='split_vcf_in_chr.py'


#for chr in range ( 1, 23 ): 
chr="X"
cmdl='bsub -e %s/split.%s.err -o %s/split.%s.out \'/software/bin/python-2.7 %s/%s   %s \'' %( outerrdir,  chr, outerrdir,   chr,codedir, scriptname,  chr )
		
#print cmdl 
os.system(cmdl)  
