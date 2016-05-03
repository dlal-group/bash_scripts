#!/usr/local/bin/bash

#Script to prepare data for plotting after imputation
    # for pan in uk10k1kg.ref
# for pan in INGI.shapeit CARL.shapeit FVG.shapeit VBI.shapeit 1000Gph1.shapeit 1000GP_Phase3.shapeit INGI_1000GPh3.shapeit uk10k1kg.ref
    # for pan in 1000GP_Phase3.shapeit

# pop=CARL
    # for pan in INGI.shapeit CARL.shapeit FVG.shapeit VBI.shapeit 1000Gph1.shapeit INGI_1000GPh3.shapeit uk10k1kg.ref
    # pan=1000Gph1.shapeit
    # pan=uk10k1kg.ref
    # do
basefolder2="/lustre/scratch113/projects/esgi-vbseq/31032016_IMPUTATION"
basefolder="/lustre/scratch113/teams/soranzo/users/jh21/imputed"

# for pop in CARL FVG INCIPE2 VBI
# for pop in fvg vbi
# for pop in carl fvg incipe2 vbi
# for pop in INCIPE2
# pop=CARL
# pop=FVG
# pop=VBI
# pop=INCIPE2
# for pan in TGP3_ALL.shapeit EUR.shapeit
# for pan in CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit
for pop in carl
do
for pan in uk10k1kg.ref
# for pan in CARL_FVG_VBI.shapeit CARL_FVG_VBI_TSI.shapeit CARL_FVG_VBI_TGP3_ALL.shapeit uk10k1kg.ref TGP3_ALL.shapeit EUR.shapeit CARL.shapeit
# pan=TGP3_ALL.shapeit
# pan=CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit
do

        