#!/usr/bin/env python

###########################################
###  Imports  #############################
###########################################

from __future__ import print_function

import numpy as np
import scipy as sp
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

# Change directory for this run... (is this really neccessary?)
#os.chdir(opa.abspath(opa.join(THIS_PATH,'../MCMC')))

# Data paths
INPATH = opa.abspath(opa.join(THIS_PATH,'input')) 
OUTPATH = opa.abspath(opa.join(THIS_PATH,'output')) 



# The dataframe with the different cosmologies

dfcosmo = pd.read_csv(opa.join(INPATH,'DataFrameCosmosims.csv'),index_col=0)



###########################################
###  Functions ##############################
##########################################
    
def cH(Om,a):
    return np.sqrt(Om/a+a**2*(1-Om))
def DgN(Om,a):
    
    return 5./2*Om*cH(Om,a)/a*scipy.integrate.quad(lambda x:cH(Om,x)**-3,0,a)[0]
def fN(Om,a):
    return (Om*(5*a - 3*DgN(Om,a)))/(2.*(a**3*(1 - Om) + Om)*DgN(Om,a))    
       
def s8(theta,s8interp):
    lnAs,Om,h = theta
    return s8interp((lnAs,Om,h))

def H(Om,h,z):
        return h*((Om)*(1+z)**3.+(1-Om))**0.5
def DA(Om,h,z):
        r = scipy.integrate.quad(lambda x:1./H(Om,h,x), 0, z)[0]
    	return r/(1+z)           

def AParam(theta,theta_true,z_pk):
    lnAs,Om,h = theta
    lnAsfid,Omfid,hfid = theta_true
    
    qperp  =  DA(Om,h,z_pk)/DA(Omfid,hfid,z_pk)
    qpar  =  H(Omfid,hfid,z_pk)/H(Om,h,z_pk)
    return [qperp,qpar]
   
def find_confidence_intervalor(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def getsigma(name,histoparams):
    return np.array([np.percentile(histoparams[name],50),np.percentile(histoparams[name],84.13445)-np.percentile(histoparams[name],50),np.percentile(histoparams[name],50)-np.percentile(histoparams[name],15.86555)])

def get_true_values(dfcosmo,simtype,nparammax):
    Omega_mfid = dfcosmo.loc[simtype,'Omega_m']
    h_TRUE = dfcosmo.loc[simtype,'h']
    lnAs_TRUE = dfcosmo.loc[simtype,'lnAs'] 
    b_fid = 2.
    
    return [lnAs_TRUE,Omega_mfid,h_TRUE,b_fid]+[0]*(nparammax-4)    



def get_grid(gridname,nbinsAs=100,nbins=50):
    
    """ Computes the power spectra given the b_i and the EFT power spectra

        Inputs
        ------
        gridname : The name of grid associated to the sim
        nbinsAs : number of bins for As (default is 100)
        nbinsAs : number of bins for Om and h (default is 50)

        Outputs
        ------
        The min,max values for the three parameters as well as the interpolation for the linear and loop power spectra
    """
    
    thetatab = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/Tablecoord%s.npy'%gridname)))

    theta3D = thetatab.reshape((nbinsAs,nbins,nbins,3))

    lnAstab = theta3D[:,0,0,0]
    Omtab = theta3D[0,:,0,1]
    htab = theta3D[0,0,:,2]


    lnAsmin = lnAstab.min()
    lnAsmax = lnAstab.max()
    Ommin = Omtab.min()
    Ommax = Omtab.max()
    hmin = htab.min()
    hmax = htab.max()

    Tables8 = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/Tables8%s.npy'%gridname)))

    s8interp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),Tables8.reshape((nbinsAs,nbins,nbins)),bounds_error=False,fill_value=-100.)
    return lnAsmin,lnAsmax,Ommin,Ommax,hmin,hmax,s8interp    
            
def density_contoursingle(xdata,ydata,nbinsx,nbinsy,smoothness,lcol,ax = None, **contour_kwargs):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """

    # Compute the confidence levels
    H, xedges, yedges = np.histogram2d(xdata, ydata, bins = (nbinsx,nbinsy), normed = True)
    x_bin_sizes = (xedges[1:] - xedges[:-1])[0]
    y_bin_sizes = (yedges[1:] - yedges[:-1])[0]

    pdfor = (H*(x_bin_sizes*y_bin_sizes))

    pdf = sn.filters.gaussian_filter(pdfor,nbinsx*smoothness)/np.sum(pdfor)

    zero_sigma = so.brentq(find_confidence_intervalor, 0., 1., args = (pdf, 0.))
    one_sigma = so.brentq(find_confidence_intervalor, 0., 1., args = (pdf, 0.68))
    two_sigma = so.brentq(find_confidence_intervalor, 0., 1., args = (pdf, 0.95))
    three_sigma = so.brentq(find_confidence_intervalor, 0., 1., args = (pdf, 0.99))
    levels = [three_sigma,two_sigma,one_sigma,zero_sigma ]

    # The smoothed histogram

    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])


    #Zint = scipy.interpolate.interp2d(X, Y, pdf, kind = 'cubic')

    if ax ==  None:
        contour = plt.contourf(X, Y, pdf.T, levels = levels,colors = (lcol[0], lcol[1], lcol[2]),origin = "lower",linewidths = 3,alpha = 0.5, **contour_kwargs)

    else: 
        contour = ax.contourf(X, Y, pdf.T, levels = levels,colors = (lcol[0], lcol[1], lcol[2]),origin = "lower",linewidths = 3,alpha = 0.5, **contour_kwargs)

    
    return contour



def getHisto(kmax,runtype,boxnumber,nmeas,dictparam,params = [r'$\ln(10^{10}A_s)$',r'$\Omega_{\rm m}$',r'$h$'],Bisp='',s8interp=None,z_pk=0.5,theta_true=None):
    """ Gets the chains for the desired paramaters
    Inputs
    ----------
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
    params : list of string
        The desired paramters for the output. Default is cosmology only
    Outputs
    ----------
    histoparams: dict
        The chains for the desired parameters
    pvalue: float
        the p-value for the fit
    mcmcarray: numpy array
       The central and percentiles for the whole chains 
    """

    samplerchaintab0 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_0.npy")%(boxnumber,kmax,Bisp)))
    samplerchaintab1 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_1.npy")%(boxnumber,kmax,Bisp)))
    samplerchaintab2 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_2.npy")%(boxnumber,kmax,Bisp)))
    samplerchaintab3 = np.array(np.load(opa.join(OUTPATH,"Chains/samplerchain"+runtype+"box_%skmax_%s%srun_3.npy")%(boxnumber,kmax,Bisp)))

    nparam = samplerchaintab3.shape[-1]
    
    samplestot = np.array([samplerchaintab0[:,:,:].reshape((-1,nparam)),samplerchaintab1[:,:,:].reshape((-1,nparam)),samplerchaintab2[:,:,:].reshape((-1,nparam)),samplerchaintab3[:,:,:].reshape((-1,nparam))])
    samples = (samplestot[:,samplestot.shape[1]/2:,:]).reshape((-1,nparam))
    minchi2 = np.array(np.loadtxt(opa.join(OUTPATH,"Chains/minchi2"+runtype+"box_%skmax_%s%s.txt")%(boxnumber,kmax,Bisp)))
    
    dof = nmeas-nparam
    pvalue = stats.distributions.chi2.sf(minchi2[-1],dof)

    
    histoparams = {}
    for name in params:
        if name in dictparam.keys():
            histoparams[name] = samples[:,dictparam[name]]
        if name == r'$f$':
            histof = np.array([fN(Om,1./(1+z_pk)) for Om in samples[:,dictparam[r'$\Omega_{\rm m}$']]])
            name = r'$f$'
            histoparams[r'$f$'] = histof

        if name == r'$\sigma_8$':
            if type(s8interp) == type(None):
                raise Exception('You want to compute sigma_8 but did not provide a interpolated function s8interp') 
            else :
                histos8 =  np.array([s8(theta,s8interp) for theta in samples[:,[dictparam[r'$\ln(10^{10}A_s)$'],dictparam[r'$\Omega_{\rm m}$'],dictparam['$h$']]]])
                histoparams[r'$\sigma_8$'] = histos8
                name = r'$f$'
        
        if name == r'$q_{\rm perp}$' or name == r'$q_{\rm par}$':
            histoAP =  np.array([AParam(theta,theta_true,z_pk)for theta in samples[:,[dictparam[r'$\ln(10^{10}A_s)$'],dictparam[r'$\Omega_{\rm m}$'],dictparam['$h$']]]])
            histoparams[r'$q_{\rm perp}$'] = histoAP[:,0]
            histoparams[r'$q_{\rm par}$'] = histoAP[:,0]
        
        if name == r'$f\sigma_8$':
            if type(s8interp) == type(None):
                raise Exception('You want to compute sigma_8 but did not provide a interpolated function s8interp') 
            else :
                histof = np.array([fN(Om,1./(1+z_pk)) for Om in samples[:,dictparam[r'$\Omega_{\rm m}$']]])
                histos8 =  np.array([s8(theta,s8interp) for theta in samples[:,[dictparam[r'$\ln(10^{10}A_s)$'],dictparam[r'$\Omega_{\rm m}$'],dictparam['$h$']]]])
                histoparams[name] = histof*histos8
        
        if name == r'$b_1\sigma_8$':
            if type(s8interp) == type(None):
                raise Exception('You want to compute sigma_8 but did not provide a interpolated function s8interp') 
            else :
                histob1 =  samples[:,dictparam[r'$b_{1}$']]
                histos8 =  np.array([s8(theta,s8interp) for theta in samples[:,[dictparam[r'$\ln(10^{10}A_s)$'],dictparam[r'$\Omega_{\rm m}$'],dictparam['$h$']]]])
                histoparams[name] = histob1*histos8

        
    
    return histoparams,pvalue
    

    
###########################################
###  Data ##############################
###########################################


#### Choice for the data #####

boxnumber = 5#'G' 
KMAX = 0.2
runtype = 'PatchyDidaNGCsigsq'
simtype = 'LightConeDida'

if 'NGC' in runtype:
    ZONE = 'NGC'
elif 'SGC' in runtype:
    ZONE = 'SGC'
else :
    ZONE = ''
    
plottype = 'Original'

# Three differents type of plots: In terms of the CLASS parameters ('Original'), with f, \sigma_8, h ('fs8h') or usual LSS-type (b1\sigma_8,f_\sigma_8,q_perp,qpar

paramsplottype = {'Original':[r'$\ln(10^{10}A_s)$',r'$\Omega_{\rm m}$',r'$h$'],
                  'fs8h' : [r'$f$',r'$\sigma_8$',r'$h$'],
                   'LSS':[r'$b_1\sigma_8$',r'$f\sigma_8$',r'$q_{\rm perp}$',r'$q_{\rm par}$']}

gridname = dfcosmo.loc[simtype,'gridname']
nparammax = 20
z_pk = dfcosmo.loc[simtype,'z_pk']
kmin = 0.01
kPS,PSdata,_ = np.loadtxt(opa.join(INPATH,'DataSims/ps1D_%s%s_%s.dat'%(simtype,ZONE,boxnumber))).T

#We need the min max values, as well as an interpolation of s8(lnAs,Om,h) on the grid.

lnAsmin,lnAsmax,Ommin,Ommax,hmin,hmax,s8interp = get_grid(gridname)


 
                                                  #######################################################
                                                  ###  True value of parameters (put in a dictionary) ###
                                                  #######################################################

    
labelnames = [r'$\ln(10^{10}A_s)$',r'$\Omega_{\rm m}$',r'$h$']+[r'$b_{%s}$' %i for i in range(1,nparammax-2)]

truevalueparam = get_true_values(dfcosmo,simtype,nparammax)
dictparam = {}
dict_true_values = {}
for i,name in enumerate(labelnames):
    dictparam[name] = i
    dict_true_values[name] = truevalueparam[i]

theta_true = np.array([dict_true_values[r'$\ln(10^{10}A_s)$'],  dict_true_values[r'$\Omega_{\rm m}$'],dict_true_values[r'$h$']])      


dict_true_values[r'$f$'] = fN(dict_true_values[r'$\Omega_{\rm m}$'],1./(1+z_pk))
dict_true_values[r'$\sigma_8$'] = s8(theta_true,s8interp)
dict_true_values[r'$q_{\rm perp}$'] = 1.
dict_true_values[r'$q_{\rm par}$'] = 1.
dict_true_values[r'$b_1\sigma_8$'] = 2*s8(theta_true,s8interp)
dict_true_values[r'$f\sigma_8$'] = fN(dict_true_values[r'$\Omega_{\rm m}$'],1./(1+z_pk))*s8(theta_true,s8interp)


                                                  ###############################
                                                  ###  Getting the histograms ###
                                                  ###############################

histoparams = []
pvalue = []
mcmcarray = []
  
params = paramsplottype[plottype]

#Can choose to plot multiples runs at the same time by adding the corresponding KMAX and runtype in configtab
configtab = [[KMAX ,runtype]]
kmaxtab = [float(x[0]) for x in configtab]

#Get number of points to compute p-value of chi^2
nmeas = {}
for i in kmaxtab:
    nmeas[str(i)] =len(kPS[(kPS <= i)&(kPS >= kmin)])

for i,config in enumerate(configtab):
    kmax,runtype = config
    res = getHisto(kmax,runtype,boxnumber,nmeas[str(kmax)],dictparam,params=params,z_pk=z_pk,theta_true=theta_true)
    histoparams.append(res[0])
    pvalue.append(res[1])

                                                                                            
                                                    ################################
                                                    ###  Plotting the 1-d histos ###
                                                    ################################

#parameters for the histograms
nbhist1D = 30
nbhist2D = 50
smooth1D = 0.05
smooth2D = 0.03



gs = gridspec.GridSpec(len(configtab),len(params) )
fig1 = plt.figure(figsize = (18,10))
nametab = [ 'EFT at 1-loop']
for i in range(len(configtab))[::-1]:
  if i == len(configtab)-1:
    
    for iname,name in enumerate(params):
        ax00 = plt.subplot(gs[ i,iname])
        medval,sigmaplus,sigmaminus = getsigma(name,histoparams[i])

        plt.hist(histoparams[i][name],bins = nbhist1D,color = lcoltab[i][-1],label = r'$\sigma$({}) = ${:.3f}$'.format(name,0.5*(sigmaplus+sigmaminus)),alpha = 0.5,normed = True)
        
        plt.axvline(medval+sigmaplus,ls = '--',color = lcoltab[i][-1],lw=2)
        plt.axvline(medval-sigmaminus,ls = '--',color = lcoltab[i][-1],lw=2)
        plt.axvline(medval,color = lcoltab[i][-1],lw=2)
        plt.axvline(dict_true_values[name],color = 'black',lw=2)
        plt.xlabel(name,fontsize=15)
        plt.legend(loc = 1)
    
        plt.tick_params(
        axis = 'y',
        which = 'both',      
        left = 'off',      
        right = 'off',      
        labelleft = 'off')
        if iname == 0:
            plt.ylabel(r'%s with $k_{\rm max} = $%s'%(nametab[i],kmaxtab[i]),rotation = 'horizontal',fontsize=15)

    
  else:
    for iname,name in enumerate(params):
        axOr = plt.subplot(gs[ -1,iname])
        ax0 = plt.subplot(gs[i,iname],sharex = axOr,sharey = axOr)
        medval,sigmaplus,sigmaminus = getsigma(name,histoparams[i])

        plt.hist(histoparams[i][name],bins = nbhist1D,color = lcoltab[i][-1],label = r'$\sigma$({}) = ${:.3f}$'.format(name,0.5*(sigmaplus+sigmaminus)),alpha = 0.5,normed = True)
        plt.axvline(medval+sigmaplus,ls = '--',color = lcoltab[i][-1],lw=2)
        plt.axvline(medval-sigmaminus,ls = '--',color = lcoltab[i][-1],lw=2)
        plt.axvline(medval,color = lcoltab[i][-1],lw=2)
        plt.axvline(dict_true_values[name],color = 'black',lw=2)
        plt.xlabel(name,fontsize=15)
        plt.legend(loc = 1,fontsize=15)
    
        plt.tick_params(
        axis = 'y',
        which = 'both',     
        left = 'off',      
        right = 'off',         
        labelleft = 'off')
   
        plt.tick_params(axis = 'x',         
        which = 'both',      
        bottom = 'off',      
        top = 'off',         
        labelbottom = 'off')
        if iname == 0:
            plt.ylabel(r'%s with $k_{\rm max} = $%s'%(nametab[i],kmaxtab[i]),rotation = 'horizontal',fontsize=15)


gs.update(wspace = 0.25, hspace = 0.05) 
plt.show()
plt.savefig(opa.join(OUTPATH,'Figs/1D_histo_%s_%s_box%s.pdf')%(plottype,runtype,boxnumber))



                                                    ################################
                                                    ### Plotting the corner plot ###
                                                    ################################
                                                                                            
gs2 = gridspec.GridSpec(len(params), len(params))
plt.figure(figsize = (18,10))

stringtitle = ''
for i in range(len(kmaxtab)):    
    
    title = 'In %s \n'%namecol[i]+r'%s with $k_{\rm max} = %s$'%(nametab[i],kmaxtab[i]) +r', $p$'+'-value best fit : {:.2f} \n \n '.format(pvalue[i])
    stringtitle+= title
    plt.suptitle(stringtitle,x=0.6,fontsize = 15)

    for inamex,namex in enumerate(params):
        axdiag = plt.subplot(gs2[inamex, inamex])
        medval,sigmaplus,sigmaminus = getsigma(namex,histoparams[i])

        plt.hist(histoparams[i][namex],bins = nbhist1D,color = lcoltab[i][-1],label = r'$\sigma$('+namex+') = ${:.3f}$'.format(0.5*(sigmaplus+sigmaminus)),alpha = 0.5,normed = True)
        
        plt.axvline(medval+sigmaplus,ls = '--',color = lcoltab[i][-1],lw=2)
        plt.axvline(medval-sigmaminus,ls = '--',color = lcoltab[i][-1],lw=2)
        plt.axvline(medval,color = lcoltab[i][-1],lw=2)
        plt.axvline(dict_true_values[namex],color = 'black')

        plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.3),
          fancybox = True, shadow = True, ncol = 2,fontsize=15)
        plt.tick_params(
        axis = 'y',
        which = 'both',      
        left = 'on',      
        right = 'on',         
        labelleft = 'off')
        if inamex == len(params)-1:
            plt.xlabel(namex,fontsize=15)
        else :
            plt.tick_params(
                        axis = 'y',
                        which = 'both',      
                        left = 'on',      
                        right = 'on',        
                        labelleft = 'off')
            plt.tick_params(axis = 'x',         
                    which = 'both',      
                    bottom = 'on',      
                    top = 'on',         
                    labelbottom = 'off')

        for jnamey in range(len(params))[::-1]:
            namey = params[jnamey]
            if jnamey > inamex:
                ax2d = plt.subplot(gs2[jnamey,inamex])
                density_contoursingle(histoparams[i][namex],histoparams[i][namey],nbhist2D,nbhist2D,smooth2D,lcoltab[i],ax = ax2d)
                plt.axvline(dict_true_values[namex],color = 'black')
                plt.axhline(dict_true_values[namey],color = 'black')
                
                if inamex == 0:
                    plt.ylabel(namey,fontsize=15)
                else :
                     plt.tick_params(
                        axis = 'y',
                        which = 'both',      
                        left = 'on',      
                        right = 'on',         
                        labelleft = 'off')
                if jnamey == len(params) -1:
                    plt.xlabel(namex,fontsize=15)
                else :
                    plt.tick_params(axis = 'x',          
                    which = 'both',      
                    bottom = 'on',      
                    top = 'on',         
                    labelbottom = 'off')
                      
gs2.update(wspace = 0.025, hspace = 0.05)  
plt.show()
plt.savefig(opa.join(OUTPATH,'Figs/CornerPlot%s_%s_box%s.pdf')%(plottype,runtype,boxnumber))


