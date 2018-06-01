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
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import scipy.interpolate as sp
import sys
#import Cmodules_runPB
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
THIS_PATH = '/home/jgleyzes/EFTofLSS'#opa.dirname(opa.dirname(__file__))



# Data paths
OUTPATH = opa.abspath(opa.join(THIS_PATH,'output/')) # Need to test it's existence:
#print(OUTPATH)
if not opa.isdir(OUTPATH): raise Exception(OUTPATH + ' not there!')


# The CLASS Boltzmann code installation location
CLASSPATH = opa.abspath(opa.join(THIS_PATH, 'class/')) # Expects the CLASS code compiled here!!!

# How to find the EFT model files and tools
EFT_PATH = opa.abspath(opa.join(THIS_PATH,'RedshiftBiasEFT_v2.4/')) # Expects the EFT code compiled here!!!
FOLDERRDPATH = opa.abspath(opa.join(EFT_PATH,'resum_data/')) # Let's not use this
FOLDERCRPATH = opa.abspath(opa.join(EFT_PATH,'cosmo_ref/')) # Should get created if not there

if opa.isdir(EFT_PATH):
    import sys
    if EFT_PATH not in sys.path: 
        sys.path.append(EFT_PATH)
    import zbeft
else:
    raise Exception('Module not found at ' + EFT_PATH)

# COSMOLOGICAL GLOBALS: fiducial model (should match input sim data!)
dfcosmo = pd.read_csv(opa.join(OUTPATH,'DataFrameCosmosims.csv'),index_col=0)
simtype = 'ChallengeA'
series_cosmo = dfcosmo.loc[simtype]
    
gridname = series_cosmo.loc['gridname']

gridname = 'ChallengeA'




ocfid =  series_cosmo.loc['Omega_m']*series_cosmo.loc['h']**2-series_cosmo.loc['omega_b']
obfid = series_cosmo.loc['omega_b']
fb = obfid/ocfid
Omega_mfid = dfcosmo.loc[simtype,'Omega_m']
h_TRUE = dfcosmo.loc[simtype,'h']
lnAs_TRUE = dfcosmo.loc[simtype,'lnAs'] 
z_pk = dfcosmo.loc[simtype,'z_pk'] 
nsfid = dfcosmo.loc[simtype,'ns'] 

# EFT PRECISION GLOBALS
EPSAbs_NoCosmoRef = 0.1
EPSRel_NoCosmoRef = 1e-5
EPSAbs_YesCosmoRef = 0.3
EPSRel_YesCosmoRef = 5e-2


###########################################
###  Grid ##############################
###########################################

lnAsmin=2.5
lnAsmax=3.5
Ommin=0.2
Ommax=0.4
hmin=0.63
hmax=0.77

thetatab=[]

nbinsAs=100
nbins=50#number of realisation of random function for each (aT0,aB0)
for i in range(0,nbinsAs): #creates a list of couple (alphaT,alphaB) to use the parallelisation
    for j in range(0,nbins):
        for k in range(0,nbins):
            lnAs=np.linspace(lnAsmin,lnAsmax,nbinsAs)[i] #Value of alphaT today
            Om=np.linspace(Ommin,Ommax,nbins)[j] #Value of alphaT today
            h=np.linspace(hmin,hmax,nbins)[k] #Value of alphaT today
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
    Outpath=OUTPATH
    zbeftpath=EFT_PATH
    classpath=CLASSPATH
    FolderRDpath=FOLDERRDPATH
    FolderCRpath=FOLDERCRPATH
        
    keys = ['b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','DM','kren',
            'outpath','zbEFT_path','CLASS_path','PathToFolderRD','PathToFolderCosmoRef',
            'UseCosmoRef',
            'EpsAbs_NoCosmoRef', 'EpsRel_NoCosmoRef', 'EpsAbs_YesCosmoRef', 'EpsRel_YesCosmoRef',
            'ln10^{10}A_s','omega_cdm','omega_b','h','z_pk','pid','n_s']
    valuesbs = [iniposb[0], iniposb[1], iniposb[2], iniposb[3], iniposb[4], iniposb[5], 
                iniposb[6], iniposb[7], iniposb[8], iniposb[9], False, 0.001,
                Outpath, zbeftpath, classpath, FolderRDpath, FolderCRpath,
                'yes',
                EPSAbs_NoCosmoRef, EPSRel_NoCosmoRef, EPSAbs_YesCosmoRef, EPSRel_YesCosmoRef,
                lnAs, omega_m/(1+fb), omega_m*fb/(1+fb),h,z_pk,str(nrank),nsfid]

    return dict(zip(keys, valuesbs))
iniposb = np.array([  2.03395135,  -6.15044179,   1.21410315,   7.19087139,
        11.61423533, -33.46605767,   1.58262629, -44.64033227,
        57.43130091,  26.44292187])

def wR(x):
    return (3*(-(x*np.cos(x)) + np.sin(x)))/x**3 

 
def CompPterms(theta,nrank,i,halfnum,run):
    
    t0 = time.time()
    

    
    # initial position for the b_i, coming from fitting the sim with the right, fixed cosmology.


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
    Pkclass=(np.loadtxt(opa.abspath(opa.join(path,'class_planck2015_pk.dat')))).T
    kclass,pkclass=Pkclass
    wRtab=np.array(map(lambda x: wR(8*x),kclass))


    s8=(1./(2*np.pi**2)*np.trapz(kclass**2*pkclass*wRtab**2,dx=kclass[1:]-kclass[:-1]))**0.5

    np.save(opa.join(OUTPATH,'textfiles/TablePloopSC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),Ploop)
    np.save(opa.join(OUTPATH,'textfiles/TablePlinSC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),Plin)
    np.save(opa.join(OUTPATH,'textfiles/Tables8SC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),s8)
    #np.save(opa.join(OUTPATH,'textfiles/TablepkClassSChallAhalf%srun%sstep%s'%(halfnum,250*run+nrank,i)),respkC)
    np.save(opa.join(OUTPATH,'textfiles/TablecoordSC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),theta)
    print(time.time()-t0,Ploop.shape)


###########################################
###  Data  ###########################
###########################################



halfnum = int(sys.argv[1])
thetahalf=thetatab[halfnum*len(thetatab)/2:(halfnum+1)*len(thetatab)/2]
ncore=500
run=int(sys.argv[2])
Ntot=len(thetahalf)
sizered=Ntot/ncore
thetarun=thetahalf[run*len(thetahalf)/2:(run+1)*len(thetahalf)/2]


def testcs(arraytheta,nrank,sizered):

    arrayred=arraytheta[nrank*sizered:(nrank+1)*sizered]
    for i in range(sizered):
        CompPterms(arrayred[i],nrank,i,halfnum,run)


testcs(thetarun,rank,sizered)   

