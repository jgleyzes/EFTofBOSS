#!/usr/bin/env python

###########################################
###  Imports  #############################
###########################################

from __future__ import print_function

import os

import numpy as np
import scipy.stats

import os.path as opa
import time

###########################################
###  Globals ##############################
###########################################

# PATH GLOBALS & tests to check they are as expected

# Save this path (expect data and other needed executables present in relative paths!!!)
THIS_PATH  =  opa.dirname(__file__)



# Data paths
INPATH = opa.abspath(opa.join(THIS_PATH,'input')) 
OUTPATH = opa.abspath(opa.join(THIS_PATH,'output')) 


def check_if_multipoles_k_array(setk):
    return setk[len(setk)/3] == setk[0]  
    
    
def apply_window_PS(setPS,PS,setkout,withmask=True,windowk=0.1):
    
    """
    Apply the window function to the power spectrum by doing a convolution directly in fourier space, encoded in Qll.
    
    Qll is an array of shape l,l',k',k where so that P_l(k) = \int dk' \sum_l' Q_{l l'}(k',k) P_l'(k')
    
    The original k on which Qll is evaluated is given by setk_or and k' by setkp_or.
    
    Inputs:
    ------
    setPS: the array of k on which PS is evaluated (ideally, the full array from Pierre's code)
    PS: the multipoles of the PS (non concatenated), shape (3,len(setPS))
    setkout: the array of k on which the results should be evaluated.
    
    withmask: whether to only do the convolution over a small window around k
    windowk: the size of said window
    
    Output:
    ------
    PStransformed: the multipoles of power spectrum evaluted on setkout with the window function applied
    
    """
    
    if check_if_multipoles_k_array(setkout):
        setkout = setkout[:len(setkout)/3]


    # Load window matrices
    Qll = np.load(opa.join(INPATH,'Window_functions/Qll_LightConeHectorNGC.npy'))  
    setk_or = np.loadtxt(opa.join(INPATH,'Window_functions/kp_LightConeHectorNGC.txt'))  
    setkp_or = np.loadtxt(opa.join(INPATH,'Window_functions/k.dat'))

    

    # Apply masking centered around the value of k
    if withmask:
        kpgrid,kgrid = np.meshgrid(setkp_or,setk_or,indexing='ij')
        mask = (kpgrid<kgrid+windowk)&(kpgrid>kgrid-windowk)
        Qll = np.einsum('lpkn,kn->lpkn',Qll,mask)
    
    # the spacing (needed to do the convolution as a sum)
    deltak = setkp_or[1:] - setkp_or[:-1]
    deltak = np.concatenate([[0],deltak])
    Qll_weighted = np.einsum('lpkn,k->lpkn',Qll,deltak)
    
    # Only keep value of setkp_or in the relevant range
    maskred = ((setkp_or>setkout.min()-0.1*windowk)&(setkp_or<setkout.max()+windowk))
    kpred = setkp_or[maskred]
    
    Qll_weighted_red = Qll_weighted[:,:,maskred,:]
    
    # Interpolate Qll(k) on setkout
    Qll_data = scipy.interpolate.interp1d(setk_or,Qll_weighted_red,axis=-1)(setkout)
    
    
    PS_red = scipy.interpolate.interp1d(setPS,PS,axis=-1,bounds_error=False,fill_value='extrapolate')(kpred)
    
    PStransformed = np.einsum('lpkm,pk->lm',Qll_data,PS_red)
    
    
    
    return PStransformed


def apply_window_covariance(Cinv,setk,withmask=True,windowk=0.1):
    """
    Apply the window function to the inverse covariance by doing a 2 convolutions directly in fourier space, encoded in Qll.
    
    Qll is an array of shape l,l',k',k where so that P_l(k) = \int dk' \sum_l' Q_{l l'}(k',k) P_l'(k')
    
    The original k on which Qll is evaluated is given by setk_or and k' by setkp_or.
    
    Then we can write chi2 = (P^window-P^data) Cinv (P^window-P^data) as
    chi2 = P^{no window} Cinvww P^{no window} -2  P^data Cinvw P^{no window} + P^data Cinv P^data
    
    Inputs:
    ------
    Cinv: the array of k on which PS is evaluated (ideally, the full array from Pierre's code)
    PS: the multipoles of the PS (non concatenated), shape (3,len(setPS))
    setk: the array of k on which Cinv is evaluated
    withmask: whether to only do the convolution over a small window around k
    windowk: the size of said window
    
    Outputs:
    ------
    Cinvpw: Cinv convoluted once for the P_model Cinv Pdata term 
    Cinvpww: Cinv convoluted twice for the P_model Cinv Pmodel term   
    """
    if check_if_multipoles_k_array(setk):
        setk = setk[:len(setk)/3]
    nkin = len(setk)

    
    
    if Cinv.shape[0]/3 != nkin:
        raise Exception('The setk needs to match the array of k for Cinv')
    
    #Put the inverse covariance in shape (l,k,l',k')    
    Cinvllp = np.swapaxes(Cinv.reshape((3,nkin,3,nkin)),axis1=1,axis2=2)
    
    # Load window matrices
    Qll = np.load(opa.join(INPATH,'Window_functions/Qll_LightConeHectorNGC.npy'))  
    setk_or = np.loadtxt(opa.join(INPATH,'Window_functions/kp_LightConeHectorNGC.txt'))  
    setkp_or = np.loadtxt(opa.join(INPATH,'Window_functions/k.dat'))

    
    if withmask:
        kpgrid,kgrid = np.meshgrid(setkp_or,setk_or,indexing='ij')
        mask = (kpgrid<kgrid+windowk)&(kpgrid>kgrid-windowk)
        Qll = np.einsum('lpkn,kn->lpkn',Qll,mask)
    
    
    # the spacing (needed to do the convolution as a sum)
    deltak = setkp_or[1:] - setkp_or[:-1]
    deltak = np.concatenate([[0],deltak])
    Qll_weighted = np.einsum('lpkn,k->lpkn',Qll,deltak)
    
    # Only keep value of setkp_or in the relevant range
    maskred = ((setkp_or>setk.min()-0.05*windowk)&(setkp_or<setk.max()+windowk))
    kpred = setkp_or[maskred]
    
    Qll_weighted_red = Qll_weighted[:,:,maskred,:]
    
    # Put the Qll(k) on the same array as Cinv for the matrix multiplication
    Qll_out = scipy.interpolate.interp1d(setk_or,Qll_weighted_red,axis=-1)(setk)
    
    nkout = sum(maskred)
    
    
    # Cinv convoluted once for the P_model Cinv Pdata term
    Cinvllpw = np.einsum('likp,imnp->lmkn', Cinvllp,Qll_out)
    
    # Cinv convoluted twice for the P_model Cinv Pmodel term
    Cinvllpww = np.einsum('imnk,ilkp->mlnp', Qll_out,Cinvllpw)
    
    # Standard form for the matrices (concatenated multipoles)
    Cinvpw = np.swapaxes(Cinvllpw,axis1=1,axis2=2).reshape((3*nkin,3*nkout)) 
    Cinvpww = np.swapaxes(Cinvllpww,axis1=1,axis2=2).reshape((3*nkout,3*nkout))
    
    return kpred,Cinvpw,Cinvpww