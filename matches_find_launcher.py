#!/software/bin/python-2.7
import os




outerrdir='/nfs/users/nfs_m/mc14/Work/bash_scripts/tempouterr'
codedir='/nfs/users/nfs_m/mc14/Work/bash_scripts/'
scriptname='matches_find.py'



for chr in range (1, 24):
        cmdl='bsub -R"select[mem>2500] rusage[mem=2500]" -M2500000  -e %s/format.%s.err -o %s/format.%s.out \'/software/bin/python-2.7 %s/%s  %s  \'' %( outerrdir, chr,  outerrdir,  chr,  codedir, scriptname, chr)

        #print cmdl
        os.system(cmdl)

