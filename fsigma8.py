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

###########################################
###  Globals ##############################
###########################################

# PATH GLOBALS & tests to check they are as expected

# Save this path (expect data and other needed executables present in relative paths!!!)
SCRATCH_PATH = '/exports/pierre/EFTofBOSS/'

# Data paths
OUTPATH2 = opa.abspath('/exports/pierre/EFTofBOSS/output/')
OUTPATH = opa.abspath('/exports/pierre/EFTofBOSS/output_fsigma8/')
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

simtype = 'LightConeHector'
#simtype = 'ChallengeA'
#simtype = 'PTchallengeCMASS2'
series_cosmo = dfcosmo.loc[simtype] 


ocfid =  series_cosmo.loc['Omega_m']*series_cosmo.loc['h']**2-series_cosmo.loc['omega_b']
obfid = series_cosmo.loc['omega_b']

fb = obfid/ocfid

Omega_mfid = dfcosmo.loc[simtype,'Omega_m']
#h_TRUE = dfcosmo.loc[simtype,'h']
#lnAs_TRUE = dfcosmo.loc[simtype,'lnAs']
z_pk = dfcosmo.loc[simtype,'z_pk']
nsfid = dfcosmo.loc[simtype,'ns']


#The pre-computed functions f(\Omega_m) and \Omega_m(f) (capital \Omega, not \omega which depends on h)
def cH(Om,a):
    return np.sqrt(Om/a+a**2*(1-Om))
def DgN(Om,a):
    return 5./2*Om*cH(Om,a)/a*scipy.integrate.quad(lambda x:cH(Om,x)**-3,0,a)[0]
def fN(Om,a):
    return (Om*(5*a - 3*DgN(Om,a)))/(2.*(a**3*(1 - Om) + Om)*DgN(Om,a))
    

    
def get_EFT_pars(lnAs, Om, h):
    
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
                lnAs, omega_m/(1+fb), omega_m*fb/(1+fb),h,z_pk,0,nsfid]

    return dict(zip(keys, valuesbs))

iniposb = np.array([  2.03395135,  -6.15044179,   1.21410315,   7.19087139,
        11.61423533, -33.46605767,   1.58262629, -44.64033227,
        57.43130091,  26.44292187])

def wR(x):
    return (3*(-(x*np.cos(x)) + np.sin(x)))/x**3 

 
def sigma8(lnAs, Om, h):
    # Get the EFT parameters
    pars = get_EFT_pars(lnAs, Om, h)
    # Runs Pierre's code and save output to folder in path.
    path = zbeft.run_zbEFTClassonly(pars)
    
    Pkclass=(np.loadtxt(opa.abspath(opa.join(path,'class_pk.dat')))).T
    kclass,pkclass=Pkclass
    wRtab=np.array(map(lambda x: wR(8*x),kclass))

    s8=(1./(2*np.pi**2)*np.trapz(kclass**2*pkclass*wRtab**2,dx=kclass[1:]-kclass[:-1]))**0.5

    return s8


###########################################
###  Data  ###########################
###########################################

kmax=0.25
run = OUTPATH2 + '/samplerchainLightConeHectorNGCprior20.0gaussMarg'
boxnumber = 'data'
Bisp = ''

a=1./(1.+z_pk)

def getchain():
    samplerchaintab0 = np.array(np.load(run+"box_%skmax_%s%srun_0.npy"%(boxnumber,kmax,Bisp)))
    samplerchaintab1 = np.array(np.load(run+"box_%skmax_%s%srun_1.npy"%(boxnumber,kmax,Bisp)))
    samplerchaintab2 = np.array(np.load(run+"box_%skmax_%s%srun_2.npy"%(boxnumber,kmax,Bisp)))
    samplerchaintab3 = np.array(np.load(run+"box_%skmax_%s%srun_3.npy"%(boxnumber,kmax,Bisp)))

    nparam = samplerchaintab3.shape[-1]
    
    samplestot = np.array([samplerchaintab0[:,:,:].reshape((-1,nparam)),samplerchaintab1[:,:,:].reshape((-1,nparam)),samplerchaintab2[:,:,:].reshape((-1,nparam)),samplerchaintab3[:,:,:].reshape((-1,nparam))])
    samples = (samplestot[:,samplestot.shape[1]/2::10,:]).reshape((-1,nparam))

    print (np.shape(samples))

    return samples

def fsigma8chain(chain, nrank):
    fs8chain=np.zeros((chain.shape[0],5))
    i=0
    for theta in chain:
        print (theta)
        lnAs,Om,h,b1,b2,b4=theta
        fs8chain[i,:]=np.array([fN(Om,a), sigma8(lnAs,Om,h), b1, b2, b4])
        i=i+1
    np.save(OUTPATH+"/textfiles/fsigma8chain%srank_%s.npy"%(Bisp,nrank), fs8chain)


samples=getchain()
ncore=240
sizered=samples.shape[0]/ncore
arrayred=samples[rank*sizered:(rank+1)*sizered]

fsigma8samples=fsigma8chain(arrayred, rank)

