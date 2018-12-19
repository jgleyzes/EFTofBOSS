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
EFT_PATH = opa.abspath(opa.join(SCRATCH_PATH,'RedshiftBiasEFT_Bispectrum/')) # Expects the EFT code compiled here!!!


if opa.isdir(EFT_PATH):
    import sys
    if EFT_PATH not in sys.path: 
        sys.path.append(EFT_PATH)
    import zbeftWindow
else:
    raise Exception('Module not found at ' + EFT_PATH)


# COSMOLOGICAL GLOBALS: fiducial model (should match input sim data!)
dfcosmo = pd.read_csv(opa.join(SCRATCH_PATH,'input/DataFrameCosmosims.csv'),index_col=0)

simtype = 'LightConeHector'
#simtype = 'ChallengeA'
#simtype = 'PTchallengeCMASS2'
series_cosmo = dfcosmo.loc[simtype] 

gridname = 'LightConeHectorPatchyWideHDvFFT_IR06'
#gridname = 'LightConeHectorDataHDvFFT_IR06'


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
### GRID FOR DATA CONFIG
lnAsmin=2.0
lnAsmax=3.3
Ommin=0.25
Ommax=0.35
hmin=0.6
hmax=0.9

nbinsAs=100
nbinsOm=48 
nbinsh=48 


'''
### GRID FOR PATCHY
lnAsmin=2.5
lnAsmax=3.8
Ommin=0.2
Ommax=0.4
hmin=0.5
hmax=0.9

nbinsAs=70
nbinsOm=48 
nbinsh=72 


thetatab=[]
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

def get_EFT_pars(theta, pathtriangle,nrank):
    lnAs, Om, h=theta
    def Hfid(z):
            omegam=Omega_mfid
            return ((omegam)*(1+z)**3.+(1-omegam))**0.5
    def H(z):
            return ((Om)*(1+z)**3.+(1-Om))**0.5
    def DAfid(z):
            r=scipy.integrate.quad(lambda x:1./Hfid(x), 0, z)[0]
    	    return r/(1+z)
    def DA(z):
            r=scipy.integrate.quad(lambda x:1./H(x), 0, z)[0]
    	    return r/(1+z)        

        
    qperp = DA(z_pk)/DAfid(z_pk)
    qpar = Hfid(z_pk)/H(z_pk)
    omega_m=Om*h**2


    # preparing the config dictionary from the globals on top (this should be in a config file!)
    Outpath=opa.abspath(opa.join(OUTPATH, 'intermediary/'))
    zbeftpath=EFT_PATH
    classpath=CLASSPATH

        
    keys = [
            'outpath','zbEFT_path','CLASS_path','PathToTriangles' ,
            'ln10^{10}A_s','omega_cdm','omega_b','h','aperp','apar','pid','n_s']
    valuesbs = [
                Outpath, zbeftpath, classpath,pathtriangle,
                lnAs, omega_m/(1+fb), omega_m*fb/(1+fb),h,qperp,qpar,str(nrank),nsfid]

    return dict(zip(keys, valuesbs))
    
iniposb = np.array([  2.03395135,  -6.15044179,   1.21410315,   7.19087139,
        11.61423533, -33.46605767,   1.58262629, -44.64033227,
        57.43130091,  26.44292187])

def wR(x):
    return (3*(-(x*np.cos(x)) + np.sin(x)))/x**3 

KMAXBISP = 0.15
dataQ=np.loadtxt(opa.join(SCRATCH_PATH,'input/Window_functions/W_v4_NGC_mask_DR12cmass_10pc.txt')).T 

def CompPterms(theta,nrank,i,half):

    lnAs, Om, h=theta
    bispecconf=np.loadtxt(opa.join(SCRATCH_PATH,'input/DataSims/Bispred_LightConeHector_NGC_data.dat')).T
    Q1,Q2,Q3=bispecconf[:3,:]
    maskBisp = (Q1<KMAXBISP)&(Q2<KMAXBISP)&(Q3<KMAXBISP)
    np.savetxt(opa.join(OUTPATH,'intermediary/TrianglesConfigurationBispgridrank%s.dat')%nrank,np.array([Q1[maskBisp],Q2[maskBisp],Q3[maskBisp]]).T)
    # Get the EFT parameters
    pars = get_EFT_pars([lnAs, Om, h], opa.join(OUTPATH,'intermediary/TrianglesConfigurationBispgridrank%s.dat')%nrank,nrank)

    # Runs Pierre's code and save output to folder in path.
    #path = zbeftWindow.run_zbEFT(pars,dataQ)
    #path = zbeft.run_zbEFT(pars)
    path = zbeftWindow.run_zbEFT(pars,dataQ)
    
    # Get the k-values

    Coefftermscenter=(np.loadtxt(opa.abspath(opa.join(path,'BispectrumTreeMonopole.dat')))).T

    np.save(opa.join(OUTPATH,'textfiles/TableBispSC%shalf%srun%sstep%s'%(gridname,half,nrank,i)),Coefftermscenter)
    #shutil.rmtree(path)
    print('Doing pretty good at rank %s step %s'%(nrank,i))


###########################################
###  Data  ###########################
###########################################

halfnum = 0
ncore=240
Ntot=len(thetatab)
sizered=Ntot/ncore

def testcs(arraytheta,nrank,sizered):
    arrayred=arraytheta[nrank*sizered:(nrank+1)*sizered]
    counter = 0
    for i in range(sizered):
        CompPterms(arrayred[i],nrank,i,halfnum)

testcs(thetatab,rank,sizered)

