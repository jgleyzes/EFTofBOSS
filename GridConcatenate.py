#!/usr/bin/env python

from __future__ import print_function

import numpy as np


halfnum=0
run=0
#gridname = 'PTchallengeCMASS2v1'
#gridname = 'conLightConeHectorv1.13'
#gridname = 'LightConeHectorHDvFFT'
#gridname = 'ChallengeHDvCUBA_IR06'
#gridname = 'ChallengeHDvFFT_IR06'
#gridname = 'LightConeHectorPatchyHDvFFT_IR2'
#gridname = 'LightConeDidaHDvFFT_IR06'

#gridname = 'ChallengeJapanHDvFFT_IR06'

gridname = 'ChallengeFullHDvFFT_IR06'


#tablepkstest=np.load('/scratch/users/kokron/textfiles/Tablepks0.npy')
tablePlintest=np.load('/exports/pierre/EFTofBOSS/output/textfiles/TablePlinSC%shalf%srun%sstep0.npy'%(gridname,halfnum,250*run))
tablePlooptest=np.load('/exports/pierre/EFTofBOSS/output/textfiles/TablePloopSC%shalf%srun%sstep0.npy'%(gridname,halfnum,250*run))
#tables8test=np.load('/exports/pierre/EFTofBOSS/output/textfiles/Tables8SC%shalf%srun%sstep0.npy'%(gridname,halfnum,250*run))
#tableftest=np.load('/scratch/users/kokron/textfiles/TablefSCPatchyhalf%srun%sstep0.npy'%(gridname,halfnum,250*run))
#tablePktest=np.load('/scratch/users/kokron/textfiles/TablepkSPBwidelassChalf%srun%sstep0.npy'%(gridname,halfnum,250*run+missing[0]))

Ntot = 150*48*72

#Ntot=70*48*72

#tablef=np.zeros(Ntot)
#tables8=np.zeros(Ntot)
tablePlin=np.zeros((Ntot,tablePlintest.shape[0],tablePlintest.shape[1]))
tablecoord=np.zeros((Ntot,3))
tablePloop=np.zeros((Ntot,tablePlooptest.shape[0],tablePlooptest.shape[1]))
#tablePk=np.zeros((Ntot,tablePktest.shape[0],tablePktest.shape[1]))

print("ready")


for nrank in range(240):
    for i in range(Ntot/240):
        #tables8[i+nrank*Ntot/240]=np.load('/exports/pierre/EFTofBOSS/output/textfiles/Tables8SC%shalf%srun%sstep%s.npy'%(gridname,halfnum,nrank,i))
        #tablesigsq[i]=np.load('/exports/pierre/EFTofBOSS/output/TablesigsqSC%shalf%srun%sstep%s.npy'%(gridname,halfnum,nrank,i))
        tablePlin[i+nrank*Ntot/240]=np.load('/exports/pierre/EFTofBOSS/output/textfiles/TablePlinSC%shalf%srun%sstep%s.npy'%(gridname,halfnum,nrank,i))
        tablePloop[i+nrank*Ntot/240]=np.load('/exports/pierre/EFTofBOSS/output/textfiles/TablePloopSC%shalf%srun%sstep%s.npy'%(gridname,halfnum,nrank,i))
        tablecoord[i+nrank*Ntot/240]=np.load('/exports/pierre/EFTofBOSS/output/textfiles/TablecoordSC%shalf%srun%sstep%s.npy'%(gridname,halfnum,nrank,i))
  
    print("Finished nrank of ", nrank)

#gridname_out = 'ChallengeHDvFFT_IR06'
#gridname_out = 'LightConeHectorPatchyHDvFFT_IR06'

np.save('/exports/pierre/EFTofBOSS/input/GridsEFT/TablePloop%s'%(gridname),tablePloop)
np.save('/exports/pierre/EFTofBOSS/input/GridsEFT/TablePlin%s'%(gridname),tablePlin)
#np.save('/exports/pierre/EFTofBOSS/output/Tables8%s'%(gridname),tables8)
#np.save('/exports/pierre/EFTofBOSS/output/Tablesigsq%s'%(gridname),tablesigsq)
np.save('/exports/pierre/EFTofBOSS/input/GridsEFT/Tablecoord%s'%(gridname),tablecoord)
