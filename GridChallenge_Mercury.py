#!/usr/bin/env python

from __future__ import print_function

import time

t0 = time.time()

#import triangle
import numpy as np
import scipy as sp
import scipy.stats
from numpy.linalg import inv
import scipy.optimize as op
import os.path as opa
import scipy.interpolate as sp
import sys
import os
import pandas as pd
os.chdir('/'+os.path.dirname(__file__))


from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#print(rank)

###########################################
###  Globals ##############################
###########################################

# PATH GLOBALS & tests to check they are as expected

# Save this path (expect data and other needed executables present in relative paths!!!)
SCRATCH_PATH = '/exports/pierre/EFTofBOSS/'

# Data paths
OUTPATH = opa.abspath('/exports/pierre/EFTofBOSS/output/')
#print(OUTPATH)
if not opa.isdir(OUTPATH): raise Exception(OUTPATH + ' not there!')

# The CLASS Boltzmann code installation location
CLASSPATH = opa.abspath(opa.join(SCRATCH_PATH, 'class/')) # Expects the CLASS code compiled here!!!

# How to find the EFT model files and tools
EFT_PATH = opa.abspath(opa.join(SCRATCH_PATH,'RedshiftBiasEFTwithFFT/')) # Expects the EFT code compiled here!!!
#EFT_PATH = opa.abspath(opa.join(SCRATCH_PATH,'RedshiftBiasEFT_C/'))


if opa.isdir(EFT_PATH):
    import sys
    if EFT_PATH not in sys.path: 
        sys.path.append(EFT_PATH)
    import zbeft
else:
    raise Exception('Module not found at ' + EFT_PATH)


# COSMOLOGICAL GLOBALS: fiducial model (should match input sim data!)
dfcosmo = pd.read_csv(opa.join(SCRATCH_PATH,'input/DataFrameCosmosims.csv'),index_col=0)

simtype = 'LightConeDida'
#simtype = 'ChallengeA'
#simtype = 'PTchallengeCMASS2'
series_cosmo = dfcosmo.loc[simtype] 



#gridname = series_cosmo.loc['gridname']
gridname = 'LightConeDidaHDvFFT_IR06'
#gridname = 'ChallengeHDvFFT_IR06'
ocfid =  series_cosmo.loc['Omega_m']*series_cosmo.loc['h']**2-series_cosmo.loc['omega_b']
obfid = series_cosmo.loc['omega_b']

#fb = series_cosmo.loc['fb']
fb = obfid/ocfid


Omega_mfid = dfcosmo.loc[simtype,'Omega_m']
#h_TRUE = dfcosmo.loc[simtype,'h']
#lnAs_TRUE = dfcosmo.loc[simtype,'lnAs']
z_pk = dfcosmo.loc[simtype,'z_pk']
nsfid = dfcosmo.loc[simtype,'ns']



###########################################
###  Grid ##############################
###########################################

'''
### GRID FOR CHALLENGE
lnAsmin=2.5
lnAsmax=3.5
Ommin=0.23
Ommax=0.37
hmin=0.55
hmax=0.75
'''

'''
### GRID FOR DATA CONFIG
lnAsmin=2.0
lnAsmax=3.3
Ommin=0.25
Ommax=0.35
hmin=0.6
hmax=0.9
'''

'''
### GRID FOR PATCHY
lnAsmin=2.5
lnAsmax=3.8
Ommin=0.2
Ommax=0.4
hmin=0.5
hmax=0.9
'''

### GRID FOR Dida
lnAsmin=2.0
lnAsmax=3.8
Ommin=0.2
Ommax=0.4
hmin=0.5
hmax=0.9

thetatab=[]
nbinsAs=70
nbinsOm=48 #number of realisation of random function for each (aT0,aB0)
nbinsh=72 #number of realisation of random function for each (aT0,aB0)
for i in range(0,nbinsAs): #creates a list of couple (alphaT,alphaB) to use the parallelisation
    for j in range(0,nbinsOm):
        for k in range(0,nbinsh):
            lnAs=np.linspace(lnAsmin,lnAsmax,nbinsAs)[i] #Value of alphaT today
            Om=np.linspace(Ommin,Ommax,nbinsOm)[j] #Value of alphaT today
            h=np.linspace(hmin,hmax,nbinsh)[k] #Value of alphaT today
            thetatab.append([lnAs,Om,h])
thetatab=np.array(thetatab)

###########################################
###  Functions  ###########################
###########################################

#The pre-computed functions f(\Omega_m) and \Omega_m(f) (capital \Omega, not \omega which depends on h)
def cH(Om,a):
    return np.sqrt(Om/a+a**2*(1-Om))
def DgN(Om,a):
    
    return 5./2*Om*cH(Om,a)/a*scipy.integrate.quad(lambda x:cH(Om,x)**-3,0,a)[0]
def fN(Om,a):
    return (Om*(5*a - 3*DgN(Om,a)))/(2.*(a**3*(1 - Om) + Om)*DgN(Om,a))
    

    
def get_EFT_pars(lnAs, Om, h, iniposb,nrank):
    
    
    omega_m=Om*h**2


    # preparing the config dictionary from the globals on top (this should be in a config file!)
    Outpath=opa.abspath(opa.join(OUTPATH, 'intermediary/'))
    zbeftpath=EFT_PATH
    classpath=CLASSPATH
        
    keys = ['b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','DM','kren',
            'outpath','zbEFT_path','CLASS_path',
            'ln10^{10}A_s','omega_cdm','omega_b','h','z_pk','pid','n_s']
    valuesbs = [iniposb[0], iniposb[1], iniposb[2], iniposb[3], iniposb[4], iniposb[5], 
                iniposb[6], iniposb[7], iniposb[8], iniposb[9], False, 0.001,
                Outpath, zbeftpath, classpath,
                lnAs, omega_m/(1+fb), omega_m*fb/(1+fb),h,z_pk,str(nrank),nsfid]

    return dict(zip(keys, valuesbs))
iniposb = np.array([  2.03395135,  -6.15044179,   1.21410315,   7.19087139,
        11.61423533, -33.46605767,   1.58262629, -44.64033227,
        57.43130091,  26.44292187])

def wR(x):
    return (3*(-(x*np.cos(x)) + np.sin(x)))/x**3 

 
def CompPterms(theta,nrank,i,halfnum,run):

    lnAs, Om, h=theta
    # Get the EFT parameters
    pars = get_EFT_pars(lnAs, Om, h, iniposb,nrank)
    # Runs Pierre's code and save output to folder in path.
    path = zbeft.run_zbEFT(pars)
    
    # Get the k-values
    P0lin = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectraLinear_l0.dat'))))
    P2lin = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectraLinear_l2.dat'))))
    P4lin = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectraLinear_l4.dat'))))
    
    P0 = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loop_l0.dat'))))
    P2 = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loop_l2.dat'))))
    P4 = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loop_l4.dat'))))
    
    Ploop=np.array([P0,P2,P4]).reshape((P0.shape[0]*3,P0.shape[1]))
    Plin=np.array([P0lin,P2lin,P4lin]).reshape((P0lin.shape[0]*3,P0lin.shape[1]))
    #Pkclass=(np.loadtxt(opa.abspath(opa.join(path,'class_pk.dat')))).T
    #kclass,pkclass=Pkclass
    #wRtab=np.array(map(lambda x: wR(8*x),kclass))

    #s8=(1./(2*np.pi**2)*np.trapz(kclass**2*pkclass*wRtab**2,dx=kclass[1:]-kclass[:-1]))**0.5
    #print(path)
    np.save(opa.join(OUTPATH,'textfiles/TablePloopSC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),Ploop)
    np.save(opa.join(OUTPATH,'textfiles/TablePlinSC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),Plin)
    #np.save(opa.join(OUTPATH,'textfiles/Tables8SC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),s8)
    #np.save(opa.join(OUTPATH,'textfiles/TablepkClassSChallAhalf%srun%sstep%s'%(halfnum,250*run+nrank,i)),respkC)
    np.save(opa.join(OUTPATH,'textfiles/TablecoordSC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),theta)


###########################################
###  Data  ###########################
###########################################


halfnum = 0
ncore=240
run=0
Ntot=len(thetatab)
sizered=Ntot/ncore

def testcs(arraytheta,nrank,sizered):
    arrayred=arraytheta[nrank*sizered:(nrank+1)*sizered]
    counter = 0
    for i in range(sizered):
        CompPterms(arrayred[i],nrank,i,halfnum,run)


testcs(thetatab,rank,sizered)   
