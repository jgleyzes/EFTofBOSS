#!/usr/bin/env python

from __future__ import print_function

import numpy as np

#Bisp=''
Bisp = 'withBispkmax0.1'
Planck = ''
#Planck = 'withPlanck'

chain=np.load('/exports/pierre/EFTofBOSS/output_fsigma8/textfiles/fsigma8chain%s%srank_0.npy'%(Bisp,Planck))

#chain = np.zeros((120,test.shape[0],test.shape[1]))

for nrank in range(240):
    chain=np.concatenate((chain,np.load('/exports/pierre/EFTofBOSS/output_fsigma8/textfiles/fsigma8chain%s%srank_%s.npy'%(Bisp,Planck,nrank))))
    print("Finished nrank of ", nrank)

print(np.shape(chain))

np.save('/exports/pierre/EFTofBOSS/output_fsigma8/fsigma8%s%s.npy'%(Bisp,Planck),chain)
