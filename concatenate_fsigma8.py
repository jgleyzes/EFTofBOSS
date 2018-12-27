#!/usr/bin/env python

from __future__ import print_function

import numpy as np

chain=np.load('/exports/pierre/EFTofBOSS/output_fsigma8/textfiles/fsigma8chainrank_0.npy')

#chain = np.zeros((120,test.shape[0],test.shape[1]))

for nrank in range(120):
    chain=np.concatenate((chain,np.load('/exports/pierre/EFTofBOSS/output_fsigma8/textfiles/fsigma8chainrank_%s.npy'%(nrank))))
    print("Finished nrank of ", nrank)

print(np.shape(chain))

np.save('/exports/pierre/EFTofBOSS/output_fsigma8/fsigma8.npy',chain)
