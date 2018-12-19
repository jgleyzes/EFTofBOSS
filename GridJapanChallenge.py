#!/usr/bin/env python

from __future__ import print_function

import time
t0 = time.time()

import numpy as np
import os.path as opa
import sys
import os
import subprocess as sp
import configobj as cfg

import camb
from camb import model, initialpower
import pandas as pd
os.chdir('/'+os.path.dirname(__file__))

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#print(rank)

###########################################
###  Globals ##############################
###########################################
#SCRATCH_PATH = '/exports/pierre/EFTofBOSS/'
SCRATCH_PATH = '/home/pchonje/Documents/EFTofLSS/EFTofBOSS'

OUT_PATH = opa.abspath(opa.join(SCRATCH_PATH, 'output/'))
if not opa.isdir(OUT_PATH): raise Exception(OUT_PATH + ' not there!')

EFT_PATH = opa.abspath(opa.join(SCRATCH_PATH,'RedshiftBiasEFTwithFFT'))
if opa.isdir(EFT_PATH):
    import sys
    if EFT_PATH not in sys.path: 
        sys.path.append(EFT_PATH)
else:
    raise Exception('Module not found at ' + EFT_PATH)


###########################################
### Cosmologies ###########################
###########################################

simtype = 'PTchallengeCMASS2'

# COSMOLOGICAL GLOBALS: fiducial model (should match input sim data!)
dfcosmo = pd.read_csv(opa.join(SCRATCH_PATH,'input/DataFrameCosmosims.csv'),index_col=0)
series_cosmo = dfcosmo.loc[simtype] 
gridname = series_cosmo.loc['gridname']
fb = series_cosmo.loc['fb']
Omega_mfid = dfcosmo.loc[simtype,'Omega_m']
z_pk = dfcosmo.loc[simtype,'z_pk']
nsfid = dfcosmo.loc[simtype,'ns'] 


###########################################
###  Functions  ###########################
###########################################

DEFAULTCONFIG_zbEFT = cfg.ConfigObj({
                'z_pk':1.,
                'omega_b':0.022,
                'omega_cdm':0.118,
                'h':0.7,
                'knl':1.,
                'km':1.,
                'nbar':0.00952380952,
                'PathToLinearPowerSpectrum':'camb_pk.dat',
                'PathToFolderOutput':'output',
                'ComputePowerSpectrum':'yes',
                'ImportResummationMatrix':'no',
                'ExportResummationMatrix':'no',
                #'ComputeBispectrum':'no',
                #'EpsRel_IntegrBispectrumAP':1e-3,
                #'PathToTriangles' : ,
                #'aperp':1,
                #'apar':1
                })

def run_zbEFT(lnAs, Om, h, nrank, i):

    #
    omega_m = Om*h**2
    omega_cdm = omega_m/(1.+fb)
    omega_b = omega_m*fb/(1.+fb)

    # Run CAMB
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=h*100, ombh2= omega_b, omch2= omega_cdm, tau=0.09, num_massive_neutrinos=0, mnu=0.0, standard_neutrino_neff=3.046)
    pars.InitPower.set_params(As = np.exp(lnAs)*1e-10, ns=0.9649)
    pars.YHe = 0.2454
    pars.set_accuracy(AccuracyBoost=2)
    pars.set_matter_power(redshifts=[z_pk], kmax=2.1)

    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-3, maxkh=2.1, npoints = 150)
    s8 = np.array(results.get_sigma8())
    
    # Make a directory to put inputs outputs for zbEFT
    outpath = opa.abspath( opa.join( OUT_PATH, 'intermediary/rank%sstep%s' % (nrank, i) ) )
    if not opa.isdir(outpath):
        os.makedirs(outpath)


    # Save CAMB power spectrum in a file
    pathcambfile = opa.join( outpath, DEFAULTCONFIG_zbEFT['PathToLinearPowerSpectrum'] )
    np.savetxt( pathcambfile, zip(kh, pk[0,:]) )

    # Make configuration file
    configfile = opa.join(outpath,'zbeft.ini')
    keys = ['PathToLinearPowerSpectrum', 'PathToFolderOutput', 'omega_cdm','omega_b','h','z_pk']
    values = [pathcambfile, outpath, omega_cdm, omega_b, h, z_pk]
    config = dict(zip(keys, values))

    czbEFT = cfg.ConfigObj()
    czbEFT.filename = configfile

    for key in DEFAULTCONFIG_zbEFT.keys():
        try:
            czbEFT[key] = config[key]
        except KeyError:
            czbEFT[key] = DEFAULTCONFIG_zbEFT[key]
    
    czbEFT.write()

    # Run zbEFT
    zbeftpath = opa.join(EFT_PATH, 'RedshiftBiasEFT')
    logfile = opa.join(outpath, 'zbeft.log')
    with open(logfile,"wb") as out:
        process = sp.Popen( [zbeftpath, configfile], stdout=out, stderr=out )
        try:
            process.wait()
        except KeyboardInterrupt as e:
            process.kill()
            raise e

    # return output path
    return outpath


 
def CompPterms(theta,nrank,i,halfnum,run):

    lnAs, Om, h = theta

    # Compute power spectrum
    path = run_zbEFT(lnAs, Om, h, nrank, i)
    
    # Get the power spectrum
    P0lin = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectraLinear_l0.dat'))))
    P2lin = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectraLinear_l2.dat'))))
    P4lin = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectraLinear_l4.dat'))))
    
    P0 = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loop_l0.dat'))))
    P2 = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loop_l2.dat'))))
    P4 = (np.loadtxt(opa.abspath(opa.join(path,'PowerSpectra1loop_l4.dat'))))
    
    Ploop = np.array([P0,P2,P4]).reshape((P0.shape[0]*3,P0.shape[1]))
    Plin = np.array([P0lin,P2lin,P4lin]).reshape((P0lin.shape[0]*3,P0lin.shape[1]))

    np.save(opa.join(OUT_PATH,'textfiles/TablePloopSC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),Ploop)
    np.save(opa.join(OUT_PATH,'textfiles/TablePlinSC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),Plin)
    np.save(opa.join(OUT_PATH,'textfiles/Tables8SC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),s8)
    np.save(opa.join(OUT_PATH,'textfiles/TablecoordSC%shalf%srun%sstep%s'%(gridname,halfnum,250*run+nrank,i)),theta)


###########################################
###  Produce power spectrum grid  #########
###########################################


run_zbEFT(3., 0.2, 0.7, 1, 2)


'''
# Grid
lnAsmin=2.0
lnAsmax=3.5
Ommin=0.2
Ommax=0.4
hmin=0.5
hmax=1

thetatab=[]
nbinsAs=100
nbins=48
for i in range(0,nbinsAs):
    for j in range(0,nbins):
        for k in range(0,nbins):
            lnAs=np.linspace(lnAsmin,lnAsmax,nbinsAs)[i]
            Om=np.linspace(Ommin,Ommax,nbins)[j] 
            h=np.linspace(hmin,hmax,nbins)[k] 
            thetatab.append([lnAs,Om,h])
thetatab=np.array(thetatab)

halfnum = 0
run=0
ncore=240
Ntot=len(thetatab)
sizered=Ntot/ncore ### should be an integer !!!

def testcs(arraytheta, nrank, sizered):
    arrayred=arraytheta[nrank*sizered:(nrank+1)*sizered]
    counter = 0
    for i in range(sizered):
        CompPterms(arrayred[i],nrank,i,halfnum,run)


#testcs(thetatab,rank,sizered)

'''