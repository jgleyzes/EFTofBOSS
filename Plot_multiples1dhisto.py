#!/usr/bin/env python

###########################################
###  Imports  #############################
###########################################

from __future__ import print_function

import numpy as np
import scipy
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.interpolate as sp
import os.path as opa
import scipy.ndimage as sn
import scipy.optimize as so
import pandas as pd
import brewer2mpl
themap0 = brewer2mpl.get_map('Greens','Sequential', 3)
lcol0 = themap0.mpl_colors

themap1 = brewer2mpl.get_map('Reds','Sequential', 3)
lcol1 = themap1.mpl_colors

themap2 = brewer2mpl.get_map('Purples','Sequential', 3)
lcol2 = themap2.mpl_colors

themap3 = brewer2mpl.get_map('Blues','Sequential', 3)
lcol3 = themap3.mpl_colors

themap4 = brewer2mpl.get_map('PuRd','Sequential', 3)
lcol4 = themap4.mpl_colors

lcoltab = [lcol0,lcol1,lcol2,lcol3,lcol4]
namecol = ['Green','Red','Purple','Blue',"PuRd"]

###########################################
###  Globals ##############################
###########################################

# PATH GLOBALS & tests to check they are as expected

# Save this path (expect data and other needed executables present in relative paths!!!)
THIS_PATH = opa.dirname(__file__)

# Data paths
INPATH = opa.abspath(opa.join(THIS_PATH,'input')) 
OUTPATH = opa.abspath(opa.join(THIS_PATH,'output')) 


# The dataframe with the different cosmologies

dfcosmo = pd.read_csv(opa.join(INPATH,'DataFrameCosmosims.csv'),index_col=0)


###########################################
###  Functions ########
###########################################


def interpolate_pdf_common(paramtab,pdf_with_support):
    
    """
    Return pdf evaluated on the support paramtab
    """
    
    param_mid,pdf = pdf_with_support
    return np.array([paramtab,sp.interpolate.interp1d(param_mid,pdf,bounds_error=False,fill_value = 0)(paramtab)])
    
def get_mean_pdf(array_pdfs,axis = 0):
    
    """
    Compute the mean of a list of pdfs (for different boxes) as their product. Pdfs must have same support (use function interpolate_pdf_common beforehand)
    """
    
    return np.prod(array_pdfs,axis = axis)

def get_bias(pdf_with_support,true_params):
    """
    Compute the bias (ie deviation of the mean of pdf from true parameter)
    
    Inputs
    ----------
    pdf_with_support: the pdf with its support, ie array of values for the parameter (eg [x,pdf(x)])
    true_params: the list of the true values of the parameters
    
    Outputs
    ----------
    The bias
    """
    tab, pdf = pdf_with_support
    
    param_mean = sum(tab*pdf)/sum(pdf)

    bias = param_mean-true_params
    
    return bias

def norm_max_pdf(pdf):
    """
    Normalizes a pdf to have a maximum of one
    """
    return pdf/pdf.max()

def getminmax(kmax,runtype,boxnumber):
    
    """ Gets the min and max values for each parameters in the chains
    Inputs
    ----------
    kmax: the kmax of the chains to be studied
    runtype: the name of the run studied
    boxnumber: the number of the box studied
    Outputs
    ----------
    tuple of lists with min values and max values
    """
    samplerchaintab0 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_0.npy")%(boxnumber,kmax,Bisp)))
    samplerchaintab1 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_1.npy")%(boxnumber,kmax,Bisp)))
    samplerchaintab2 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_2.npy")%(boxnumber,kmax,Bisp)))
    samplerchaintab3 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_3.npy")%(boxnumber,kmax,Bisp)))

    nparam = samplerchaintab3.shape[-1]
    
    samplestot = np.array([samplerchaintab0[:,:,:].reshape((-1,nparam)),samplerchaintab1[:,:,:].reshape((-1,nparam)),samplerchaintab2[:,:,:].reshape((-1,nparam)),samplerchaintab3[:,:,:].reshape((-1,nparam))])
    samples = (samplestot[:,samplestot.shape[1]/2:,:]).reshape((-1,nparam))
    return np.min(samples,axis=0),np.max(samples,axis=0)

def get_true_params(dfcosmo,simtype,nparammax):
    """ Gets the true value of the parameters
    Inputs
    ----------
    dfcosmo: the dataframe with the info about the simulations
    simtype: the name of the simulation studied
    nparammax: the maximum number of parameters varied (should be 13 for power spectrum only or 14 for Bispectrum)
    Outputs
    ----------
    List with the values of the true parameters (for b_i with i>1, returns 0)
    """
    
    Omega_mfid = dfcosmo.loc[simtype,'Omega_m']
    h_TRUE = dfcosmo.loc[simtype,'h']
    lnAs_TRUE = dfcosmo.loc[simtype,'lnAs'] 
    b_fid = 2.
    
    return [lnAs_TRUE,Omega_mfid,h_TRUE,b_fid]+[0]*(nparammax-4)
        

def getpdf(name,kmax,runtype,boxnumber,dictparam,nbins = 75, smoothness=0.01):
    """ Gets the pdf for the desired paramaters
    Inputs
    ----------
    params : string
        The name of the parameter
    kmax : float 
        The maximum k, which determines the chains
    runtype : string
         Describes the type of run for the chains
    boxnumber : int or string
        The box for the data on which the MCMC ran
    nmeas : int
        The number of points in the data (for the p-value)
    dictparam : dict
        Related the name of each parameter to its index in the chains
    
    nbins : int
        Number of bins in the histograms
    smoothness : float
        The smoothing factor for Gaussian smoothing of the pdf.
    Outputs
    ----------
    pdfparam: 2d array
        The array of parameter values and smoothed pdf for the desired parameters
    meanparam: float
        The mean value
    maxlikeparam: float
       The maximum of the 1D pdf
    onesigmaparam: float
        The one sigma value for the parameter
    """
    samplerchaintab0 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_0.npy")%(boxnumber,kmax,Bisp)))
    samplerchaintab1 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_1.npy")%(boxnumber,kmax,Bisp)))
    samplerchaintab2 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_2.npy")%(boxnumber,kmax,Bisp)))
    samplerchaintab3 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_3.npy")%(boxnumber,kmax,Bisp)))

    nparam = samplerchaintab3.shape[-1]
    
    samplestot = np.array([samplerchaintab0[:,:,:].reshape((-1,nparam)),samplerchaintab1[:,:,:].reshape((-1,nparam)),samplerchaintab2[:,:,:].reshape((-1,nparam)),samplerchaintab3[:,:,:].reshape((-1,nparam))])
    samples = (samplestot[:,samplestot.shape[1]/2:,:]).reshape((-1,nparam))

    mcmcarray = np.array(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0], v[4]-v[1], v[1]-v[3]), zip(*np.percentile(samples, [15.86555, 50, 84.13445, 2.2775, 97.7225], axis = 0))))    
    
    index_name = dictparam[name]
    
    #Compute histogram
    Hist, edges = np.histogram(samples[:,index_name],  bins = (nbins), normed = True)
    param_bin = edges[1]-edges[0]
    param_mid = 0.5*(edges[1:]+edges[:-1])
    
    #Compute pdf from histogram (by multiplying by binning) and then apply Gaussian smoothing
    pdf_raw = Hist*param_bin
    pdf_smooth = sn.filters.gaussian_filter(pdf_raw,nbins*smoothness)/np.sum(pdf_raw)
        
    meanparam = sum(pdf_smooth*param_mid)/sum(pdf_smooth)
    maxlikeparam = param_mid[Hist == Hist.max()][0]        
    
    pdfparam = np.array([param_mid,pdf_smooth])

    onesigmaparam = 0.5*(mcmcarray[index_name,1]+mcmcarray[index_name,2])
        
    return pdfparam,meanparam,maxlikeparam,onesigmaparam



                                                                ########################
                                                                ###  Main Program ########
                                                                #########################


# The parameters 


boxnumber = 6#'G' 
KMAX = 0.1
runtype = 'PatchyDidaNGCsigsq'
simtype = 'LightConeDida'

namefile = ''
Bisp = ''
kmaxtab = [0.2]
boxtab = np.array([2,3,5,7,8,10])#np.arange(1,11)#
nbins = 75
smoothness = 0.01

#The min and max values over all the samples
minparamaters = np.min(np.array([getminmax(kmax,runtype,box)[0] for kmax in kmaxtab for box in boxtab]),axis=0)
maxparamaters = np.max(np.array([getminmax(kmax,runtype,box)[1] for kmax in kmaxtab for box in boxtab]),axis=0)

nparam = len(maxparamaters)

labelnames = [r'$\ln(10^{10}A_s)$',r'$\Omega_{\rm m}$','h']+[r'$b_{%s}$' %i for i in range(1,nparam-2)]

true_valueparam = get_true_params(dfcosmo,simtype,nparam)


#The dictionaries between parameter names and the chains
dictparam = {}
dict_true = {}
dicttab = {}
for i,name in enumerate(labelnames):
    dictparam[name] = i
    dict_true[name] = true_valueparam[i]
    dicttab[name] = np.linspace(minparamaters[i],maxparamaters[i],nbins)


###########################################
###  Get the pdf, bias and one sigma ########
###########################################


dict_pdf_box_kmax = {}
dict_mean_box_kmax = {}
dict_bias_box_kmax = {}
dict_onesigma_box_kmax = {}



params = [r'$\ln(10^{10}A_s)$',r'$\Omega_{\rm m}$','h',r'$b_{1}$']

for name in  params: 
    pdf_box_kmax = []
    mean_box_kmax = []
    bias_box_kmax = []
    onesigma_box_kmax = []
    for boxnumber in boxtab:
        pdf_box = [] 
        mean_box = []
        bias_box = []
        onesigma_box = []
        for kmax in kmaxtab:
            pdfparam,meanparam,maxlikeparam,onesigmaparam = getpdf(name,kmax,runtype,boxnumber,dictparam,smoothness=smoothness,nbins=nbins)
            
            pdf_common = interpolate_pdf_common(dicttab[name],pdfparam)
            
            pdf_box.append(pdf_common)
            mean_box.append(meanparam)
            bias_box.append(get_bias(pdfparam,dict_true[name]))
            onesigma_box.append(onesigmaparam)    
        
        pdf_box_kmax.append(pdf_box)
        mean_box_kmax.append(mean_box)
        bias_box_kmax.append(bias_box)
        onesigma_box_kmax.append(onesigma_box)
    
    dict_pdf_box_kmax[name] = np.array(pdf_box_kmax)
    dict_mean_box_kmax[name] = np.array(mean_box_kmax)
    dict_bias_box_kmax[name] = np.array(bias_box_kmax)
    dict_onesigma_box_kmax[name] = np.array(onesigma_box_kmax)  
        

  
###########################################
###  Plotting ########
###########################################    
      
          
plt.figure(figsize = (5*len(params),4+4*len(kmaxtab)))
gs2 = gridspec.GridSpec(len(kmaxtab),len(params) )

plt.suptitle('Comparison of the pdf for run ' + runtype,fontsize=20)
for i in range(len(kmaxtab))[::-1]:
    for iname,name in enumerate(params):
        
        ax = plt.subplot(gs2[i, iname])
            
        paramtab = dict_pdf_box_kmax[name][0,0,0]
        pdfbox = np.array([norm_max_pdf(pdf) for pdf in dict_pdf_box_kmax[name][:,i,1]])
        paramfid = dict_true[name]
        paramonesigma = dict_onesigma_box_kmax[name][:,i]
        onesigmamean = np.mean(paramonesigma**2,axis = 0)**0.5
        biasbox = dict_bias_box_kmax[name][:,i]
        biasboxavg = np.mean(biasbox,axis=0)
            
        pdfmean = get_mean_pdf(pdfbox)
        plt.plot(paramtab,pdfbox.T,color = lcoltab[i][-2])
        plt.plot(paramtab,pdfmean/pdfmean.max(),color = lcoltab[i][-1],lw = 3)
        
        plt.axvline(paramfid,color = 'grey',alpha = 0.7)
        plt.axvline(paramfid-onesigmamean,color = 'grey',ls = '--')
        plt.axvline(paramfid+onesigmamean,color = 'grey',ls = '--')
        
            
        plt.text(paramfid-4*onesigmamean,1.75,'Avg Bias  = {:.3f}'.format(biasboxavg))
        plt.text(paramfid-4*onesigmamean,1.40,r'$\sqrt{\sigma_{\rm syst}^2+\sigma_{\rm stat}^2}$'+'  = {:.3f}'.format(np.sqrt(biasboxavg**2+onesigmamean**2)))
        plt.text(paramfid-4*onesigmamean,1.05,r'$<\sigma_{\rm stat}>$'+'  = {:.3f}'.format(onesigmamean))
            
        plt.tick_params(
            axis = 'y',
                which = 'both',      # both major and minor ticks are affected
                left = 'off',      # ticks along the bottom edge are off
    	        right = 'off',         # ticks along the top edge are off
                labelleft = 'off')
        
        if i == len(kmaxtab)-1:
            plt.xlabel(name,fontsize = 15)
        else:
            plt.tick_params(axis = 'x',          # changes apply to the x-axis
                which = 'both',      # both major and minor ticks are affected
                bottom = 'on',      # ticks along the bottom edge are off
                top = 'on',         # ticks along the top edge are off
                labelbottom = 'off')
        if iname == len(params) - 1:
            plt.text(paramfid+2*onesigmamean,1.4,r'$k_{\rm max}$'+' ={:.2f}'.format(kmaxtab[i]),fontsize=16)
        plt.xlim(paramfid-4*onesigmamean,paramfid+4*onesigmamean)
        plt.ylim(0,2)
gs2.update(wspace = 0.1, hspace = 0.05)  
plt.show()

plt.savefig(opa.join(OUTPATH,'Figs/ComparisonPDF_%s%suptokmax_%s.pdf')%(runtype,namefile,max(kmaxtab)))