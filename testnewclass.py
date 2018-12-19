
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

# PATH GLOBALS & tests to check they are as expected

# Save this path (expect data and other needed executables present in relative paths!!!)
SCRATCH_PATH = '/home/pchonje/Documents/EFTofLSS/EFTofBOSS/'

# Data paths
OUTPATH = opa.abspath('/home/pchonje/Documents/EFTofLSS/EFTofBOSS/output/')
#print(OUTPATH)
if not opa.isdir(OUTPATH): raise Exception(OUTPATH + ' not there!')

# The CLASS Boltzmann code installation location
CLASSPATH = opa.abspath(opa.join(SCRATCH_PATH, 'class/')) # Expects the CLASS code compiled here!!!

# How to find the EFT model files and tools
EFT_PATH = opa.abspath(opa.join(SCRATCH_PATH,'RedshiftBiasEFTwithFFT/')) # Expects the EFT code compiled here!!!


if opa.isdir(EFT_PATH):
    import sys
    if EFT_PATH not in sys.path: 
        sys.path.append(EFT_PATH)
    import zbeft
else:
    raise Exception('Module not found at ' + EFT_PATH)

# COSMOLOGICAL GLOBALS: fiducial model (should match input sim data!)
dfcosmo = pd.read_csv(opa.join(SCRATCH_PATH,'input/DataFrameCosmosims.csv'),index_col=0)

#simtype = 'LightConeHector'

simtype = 'PTchallengeCMASS2'
series_cosmo = dfcosmo.loc[simtype] 

fb = series_cosmo.loc['fb']
#fb = obfid/ocfid


gridname = series_cosmo.loc['gridname']
ocfid =  series_cosmo.loc['Omega_m']*series_cosmo.loc['h']**2-series_cosmo.loc['omega_b']
obfid = series_cosmo.loc['omega_b']


Omega_mfid = dfcosmo.loc[simtype,'Omega_m']
#h_TRUE = dfcosmo.loc[simtype,'h']
#lnAs_TRUE = dfcosmo.loc[simtype,'lnAs']
z_pk = dfcosmo.loc[simtype,'z_pk']
nsfid = dfcosmo.loc[simtype,'ns'] 


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

import zbeft

pars = get_EFT_pars(3., 0.3, 0.7, iniposb, 0)
zbeft.run_zbEFT(pars)