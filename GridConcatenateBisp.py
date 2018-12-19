#!/usr/bin/env python

from __future__ import print_function

import numpy as np


halfnum=0
run=0

gridname = 'LightConeHectorPatchyWideHDvFFT_IR06'

tableBisptest=np.load('/exports/pierre/EFTofBOSS/output/textfiles/TableBispSC%shalf%srun%sstep0.npy'%(gridname,halfnum,250*run))

#Ntot=70*48*72
Ntot = 100*48*48
tableBisp=np.zeros((Ntot,tableBisptest.shape[0],tableBisptest.shape[1]))

print("ready")


for nrank in range(240):
    for i in range(Ntot/240):
        
        tableBisp[i+nrank*Ntot/240]=np.load('/exports/pierre/EFTofBOSS/output/textfiles/TableBispSC%shalf%srun%sstep%s.npy'%(gridname,halfnum,nrank,i))
  
    print("Finished nrank of ", nrank)

#gridname_out = 'ChallengeHDvFFT_IR06'
#gridname_out = 'LightConeHectorPatchyHDvFFT_IR06'

np.save('/exports/pierre/EFTofBOSS/input/GridsEFT/TableBisp%s'%(gridname),tableBisp)
