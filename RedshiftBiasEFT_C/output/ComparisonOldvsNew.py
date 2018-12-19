#!/usr/bin/python
# -*- coding: latin-1 -*-

import matplotlib.pyplot as plt 
import numpy as np

#############################################################################
### EFT
#############################################################################

N = 10

PATH_to_MODEL = '/home/pchonje/Documents/EFTofLSS/dev/RedshiftBiasEFT_C/output/Precise/PowerSpectra'
Nk = 100

ks = np.array([ 0.01, 0.0104566, 0.010934, 0.0114332, 0.0119553, 0.0125011, 0.0130719, 0.0136688, 0.0142929, 0.0149454, 0.0156278, 0.0163414, 0.0170875, 0.0178677, 0.0186835, 0.0195366, 0.0204286, 0.0213613, 0.0223366, 0.0233565, 0.0244229, 0.025538, 0.0267041, 0.0279233, 0.0291983, 0.0305314, 0.0319254, 0.0333831, 0.0349073, 0.0365011, 0.0381677, 0.0399104, 0.0417327, 0.0436381, 0.0456306, 0.047714, 0.0498925, 0.0521706, 0.0545526, 0.0570434, 0.0596479, 0.0623713, 0.0652191, 0.0681969, 0.0713107, 0.0745666, 0.0779712, 0.0815313, 0.0852539, 0.0891464, 0.0932167, 0.0974729, 0.101923, 0.106577, 0.111443, 0.116531, 0.121852, 0.127416, 0.133233, 0.139317, 0.145678, 0.152329, 0.159284, 0.166557, 0.174162, 0.182113, 0.190428, 0.199123, 0.208215, 0.217722, 0.227662, 0.238057, 0.248927, 0.260292, 0.272177, 0.284604, 0.297598, 0.311186, 0.325395, 0.340252, 0.355787, 0.372032, 0.389018, 0.40678, 0.425353, 0.444774, 0.465082, 0.486317, 0.508521, 0.53174, 0.556018, 0.581405, 0.607951, 0.635709, 0.664735, 0.695086, 0.726822, 0.760008, 0.794709, 0.830994 ])

EftP0 = []
EftP1 = []
EftP2 = []

Pole = [ '0', '2' , '4']

for n in Pole:

    DataLin = np.loadtxt(PATH_to_MODEL+ 'Linear_l'+ n +'.dat', comments='#', skiprows=0)
    DataLoop = np.loadtxt(PATH_to_MODEL+ '1loop_l'+ n +'.dat', comments='#', skiprows=0)

    # Power spectra

    P0 = 1.*(DataLin[:,1]+DataLoop[:,1])
    P1 = np.zeros((Nk,N))
    P2 = np.zeros((Nk,N,N))

    P1[:,0] = DataLin[:,2] + DataLoop[:,2]
    P1[:,1] = DataLoop[:,3]
    P1[:,2] = DataLoop[:,4]
    P1[:,3] = DataLoop[:,5]
    P1[:,4] = DataLoop[:,16]
    P1[:,5] = DataLoop[:,17]
    P1[:,6] = DataLoop[:,18]
    P1[:,7] = DataLoop[:,19]
    P1[:,8] = DataLoop[:,20]
    P1[:,9] = DataLoop[:,21]
    
    P2[:,0,0] = DataLin[:,3] + DataLoop[:,6]
    P2[:,0,1] = DataLoop[:,7]
    P2[:,0,2] = DataLoop[:,8]
    P2[:,0,3] = DataLoop[:,9]
    P2[:,0,4] = DataLoop[:,13]
    P2[:,0,5] = DataLoop[:,14]
    P2[:,0,6] = DataLoop[:,15]

    P2[:,1,1] = DataLoop[:,10]
    P2[:,1,3] = DataLoop[:,11]
    P2[:,3,3] = DataLoop[:,12]
    

    EftP0.append(P0)
    EftP1.append(P1)
    EftP2.append(P2)

def f(k, l, params):
    params_tensor = np.tensordot(params,params,axes=0)

    Phh = np.interp( k, ks, EftP0[l] )

    for i in range(N):
        Phh += params[i] * np.interp( k, ks, (EftP1[l])[:,i] )
        for j in range(N):
            Phh += params_tensor[i,j]  * np.interp( k, ks, (EftP2[l])[:,i,j] )
    return  Phh

##############
PATH_to_MODEL2 = '/home/pchonje/Documents/EFTofLSS/dev/RedshiftBiasEFT_C/output/PowerSpectra'

EftP0o = []
EftP1o = []
EftP2o = []

for n in Pole:

    DataLin = np.loadtxt(PATH_to_MODEL2+ 'Linear_l'+ n +'.dat', comments='#', skiprows=0)
    DataLoop = np.loadtxt(PATH_to_MODEL2+ '1loop_l'+ n +'.dat', comments='#', skiprows=0)

    # Power spectra
    P0 = 1.*(DataLin[:,1]+DataLoop[:,1])
    P1 = np.zeros((Nk,N))
    P2 = np.zeros((Nk,N,N))

    P1[:,0] = DataLin[:,2] + DataLoop[:,2]
    P1[:,1] = DataLoop[:,3]
    P1[:,2] = DataLoop[:,4]
    P1[:,3] = DataLoop[:,5]
    P1[:,4] = DataLoop[:,16]
    P1[:,5] = DataLoop[:,17]
    P1[:,6] = DataLoop[:,18]
    P1[:,7] = DataLoop[:,19]
    P1[:,8] = DataLoop[:,20]
    P1[:,9] = DataLoop[:,21]
    
    P2[:,0,0] = DataLin[:,3] + DataLoop[:,6]
    P2[:,0,1] = DataLoop[:,7]
    P2[:,0,2] = DataLoop[:,8]
    P2[:,0,3] = DataLoop[:,9]
    P2[:,0,4] = DataLoop[:,13]
    P2[:,0,5] = DataLoop[:,14]
    P2[:,0,6] = DataLoop[:,15]

    P2[:,1,1] = DataLoop[:,10]
    P2[:,1,3] = DataLoop[:,11]
    P2[:,3,3] = DataLoop[:,12]
    
    EftP0o.append(P0)
    EftP1o.append(P1)
    EftP2o.append(P2)

def g(k, l, params):
    params_tensor = np.tensordot(params,params,axes=0)

    Phh = np.interp( k, ks, EftP0o[l] )

    for i in range(N):
        Phh += params[i] * np.interp( k, ks, (EftP1o[l])[:,i] )
        for j in range(N):
            Phh += params_tensor[i,j]  * np.interp( k, ks, (EftP2o[l])[:,i,j] )
    return  Phh


params = np.array([1.92059074, 2.5996115 ,  -5.50871094,  -1.08271503,  -1.93388186, -14.56153162,  -2.88083332,   0.,  29.71909245, 0.])

plt.figure(1)
plt.semilogx(ks,f(ks,0,params),color='b')
plt.semilogx(ks,g(ks,0,params),color='g')
plt.grid()

error = np.loadtxt('/home/pchonje/Documents/EFTofLSS/dev/RedshiftBiasEFT_TestResum/output/ps1D_ChallengeQuarter_A.dat', comments='#', skiprows=0)

from scipy import interpolate
A = 0
Ne = len(error[:,0])
tck = interpolate.splrep(error[Ne/3*(A/2):Ne/3*(A/2+1),0], error[Ne/3*(A/2):Ne/3*(A/2+1),2], s=0)
Perr = interpolate.splev(ks, tck, der=0)

plt.figure(2)
plt.plot(ks,(f(ks,0,params)-g(ks,0,params))/Perr,color='b')
plt.grid()

plt.show()