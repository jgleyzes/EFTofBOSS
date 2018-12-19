#!/usr/bin/env python

import numpy as np

#spectra = 'ps1D'
#sim = 'LightConeHectorNGC'
sim = 'LightConeHector_NGC'
spectra = 'Bispred'

Nb = 10
mean = np.loadtxt('/home/pchonje/Documents/EFTofLSS/EFTofBOSS/input/DataSims/'+spectra+'_'+sim+'_1.dat')
for i in range(Nb-1):
    box = np.loadtxt('/home/pchonje/Documents/EFTofLSS/EFTofBOSS/input/DataSims/'+spectra+'_'+sim+'_'+str(i+2)+'.dat')
    mean[:,3] += box[:,3]

mean[:,3] /= Nb

np.savetxt('/home/pchonje/Documents/EFTofLSS/EFTofBOSS/input/DataSims/'+spectra+'_'+sim+'_mean.dat', mean)