#!/software/bin/python-2.7
import os
import sys 

outerrdir=os.getcwd()
codedir='/nfs/users/nfs_m/mc14/Work/bash_scripts/'
scriptname='count_functional_annotations.py'
conseq_list_path=sys.argv[1]



for chr in range (1, 23):
# chr=22
  if chr==23:
  	chr="X"

cmdl='bsub -R"select[mem>2500] rusage[mem=2500]" -M2500  -e %s/format.%s.err -o %s/format.%s.out -q yesterday \'/software/bin/python-2.7 %s/%s %s %s\'' %( outerrdir, chr,  outerrdir,  chr,  codedir, scriptname, chr, conseq_list_path)

#print cmdl 
os.system(cmdl)

  
