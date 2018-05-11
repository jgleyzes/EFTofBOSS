#!/usr/bin/env python

from __future__ import print_function

import time

t0 = time.time()

#import emcee
#import triangle
import numpy as np
import scipy as sp
import scipy.stats
from numpy.linalg import inv
import scipy.optimize as op
import os.path as opa
import matplotlib.pyplot as pl
import scipy.interpolate
import sys
#import Cmodules_runPB
import os
import pandas as pd
THIS_PATH = opa.dirname(__file__)
INPATH = opa.abspath(opa.join(THIS_PATH,'input')) 
EFT_PATH = opa.join(THIS_PATH,'RedshiftBiasEFT_v2.3c/')# Expects the EFT code compiled here!!!
import zbeft




Outpath = opa.join(EFT_PATH,'output/')
zbeftpath = EFT_PATH
classpath = opa.join(THIS_PATH,'class/')
resumpath = opa.join(EFT_PATH,'resum_data/' )   
cosmorefpath = opa.join(EFT_PATH,'cosmo_ref/')  


dfcosmo = pd.read_csv(opa.join(INPATH,'DataFrameCosmosims.csv'),index_col=0)

Seriescosmo = dfcosmo.loc['LightConeDida']
    
paramfid = {'z_pk' : Seriescosmo['z_pk'],
    'format' : 'camb',
    'omega_cdm' :Seriescosmo['Omega_m']*Seriescosmo['h']**2-Seriescosmo['omega_b'], 
    'omega_b' :Seriescosmo['omega_b'] ,
    'n_s' : Seriescosmo['ns'],
    'P_k_max_h/Mpc' : 20,
    'ln10^{10}A_s' : Seriescosmo['lnAs'],
    'h' : Seriescosmo['h'],
    'output' : 'mPk'}    
            
EPSAbs_NoCosmoRef = 10
EPSRel_NoCosmoRef = 10
EPSAbs_YesCosmoRef = 10
EPSRel_YesCosmoRef = 10


              
keys = ['outpath','zbEFT_path','CLASS_path','PathToFolderRD','PathToFolderCosmoRef',
            'UseCosmoRef','ExportResummationMatrix','ImportResummationMatrix','ComputePowerSpectrum','ComputeBispectrum',
            'EpsAbs_NoCosmoRef', 'EpsRel_NoCosmoRef', 'EpsAbs_YesCosmoRef', 'EpsRel_YesCosmoRef',
            'ln10^{10}A_s','omega_cdm','omega_b','h','z_pk','n_s']
            
            
valuesbs = [Outpath, zbeftpath, classpath, resumpath, cosmorefpath,
                'yes','no','yes','yes','no',
                EPSAbs_NoCosmoRef, EPSRel_NoCosmoRef, EPSAbs_YesCosmoRef, EPSRel_YesCosmoRef,
                Seriescosmo['lnAs'], paramfid['omega_cdm'], Seriescosmo['omega_b'],Seriescosmo['h'],Seriescosmo['z_pk'],Seriescosmo['ns']]

pars = dict(zip(keys, valuesbs))


path = zbeft.run_zbEFT(pars)
 

P0lin = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectraLinear_l0.dat'))))
P2lin = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectraLinear_l2.dat'))))
P4lin = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectraLinear_l4.dat'))))
    
P0 = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loop_l0.dat'))))
P2 = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loop_l2.dat'))))
P4 = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loop_l4.dat'))))
   

P0nores = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loopNoResum_l0.dat'))))    
    
Ploop=np.array([P0,P2,P4]).reshape((P0.shape[0]*3,P0.shape[1]))
Plin=np.array([P0lin,P2lin,P4lin]).reshape((P0lin.shape[0]*3,P0lin.shape[1]))

    
