import numpy as np
import scipy as sp
import scipy.stats
import scipy.optimize as op
import os.path as opa
import scipy.interpolate
import scipy.interpolate as sp
import sys
import matplotlib.pyplot as plt
import mcfit
from classy import Class
import pandas as pd
import fftlog
# Save this path (expect data and other needed executables present in relative paths!!!)
THIS_PATH = opa.dirname(opa.dirname(__file__))

# Data paths
OUTPATH = opa.abspath(opa.join(THIS_PATH,'MCMC/output/'))

import APpowerspectraNkmu

def check_if_multipoles_k_array(setk):
    return setk[len(setk)/3] == setk[0]

    
def cH(Om,a):
    return np.sqrt(Om/a+a**2*(1-Om))
def DgN(Om,a):
    
    return 5./2*Om*cH(Om,a)/a*scipy.integrate.quad(lambda x:cH(Om,x)**-3,0,a)[0]
def fN(Om,a):
    return (Om*(5*a - 3*DgN(Om,a)))/(2.*(a**3*(1 - Om) + Om)*DgN(Om,a))
    
            
def Classtomultipoles(setk,Pkclass,Om,zpk,b1=1):
    kmulti = np.concatenate([setk,setk,setk])
    f = fN(Om,1./(1+zpk))
    P0k = (b1**2 + (2*f)*b1/3. + f**2/5.)*Pkclass
    P2k = ((4*f*(7*b1 + 3*f))/21.)*Pkclass
    P4k = ((8*f**2)/35.)*Pkclass
    return np.array([kmulti,np.concatenate([P0k,P2k,P4k])])
        

def get_powerlaw_junc(k_junc,Pfunc,dx,damp=False,adamp=1):
    
    """ Get the coefficient for the power law using (b * (k/k_junc)**c) = P(k) and d/dk (b * (k/k_junc)**c) = P(k)) at k_junc.

        Inputs
        ------
        k_junc : (float) The value of k at the junction
        Pfunc : (function) The interpolated function of P(k)
        dx : (float) The step for the numerical derivative of P(k)
        damp : (boolean) whether to damp the power law at high k with with exp(-a k/k_junc)
        adamp : (float) the value of in for the damping exp(-a k/k_junc)

        Outputs
        ------
        The coefficient b,c (and a if damping)
    """
    
    Pk_junc = Pfunc(k_junc)
    Pk1_junc = scipy.misc.derivative(Pfunc, k_junc, dx=dx, n=1, args=(), order=3)
    
    a = adamp
    if damp :
            b = np.exp(a)*Pk_junc
            c = (Pk1_junc*k_junc+a*Pk_junc)/Pk_junc
            return b,c,a
    else :
            b = Pk_junc
            c = Pk1_junc*k_junc/b
            return b,c
    

def ExtrapolationPk(Pk,setk,setkextrap,k_junc_low = 0.01,k_junc_high=0.4,damp = False,adamp=1):
    
    """ Performs an extrapolation with power laws at low and high k, using matching value and slope at the junction.

        Inputs
        ------
        Pk : (array) The values for the power spectrum (one mutlipole at a time) to be extrapolated
        setk : (array) The k-range over which power spectra is evaluated
        setkextrap : (array) the range for the extrapolation
        k_junc_low : (float) the point at which the junction is made for low k (should be within setk)
        k_junc_high : (float) the point at which the junction is made for high k (should be within setk)
        damp : (boolean) whether to damp the power law at high k with with exp(-a k/k_junc)
        adamp : (float) the value of in for the damping exp(-a k/k_junc)

        Outputs
        ------
        The extrapolated power spectrum
    """
    
    kmin = setk.min()
    kmax = setk.max()
    
    if k_junc_low < kmin:
        raise Exception('The junction point at low k must be within initial range ')
    elif k_junc_high > kmax:
        raise Exception('The junction point at high k must be within initial range. kmax = {} and k_junc = {} '.format(kmax,k_junc_high))
   
    else :
        Pkfunc = sp.interpolate.interp1d(setk,Pk,kind='cubic',bounds_error = False, fill_value = 'extrapolate')
        dx = np.mean(setk[1:]-setk[:-1])
        
        b_low,c_low = get_powerlaw_junc(k_junc_low,Pkfunc,dx)
        lowk = setkextrap[setkextrap < k_junc_low]
        Plowk = b_low*(lowk/k_junc_low)**c_low
    

        
        
        highk = setkextrap[setkextrap > k_junc_high]
        if damp :
            b_high,c_high,a = get_powerlaw_junc(k_junc_high,Pkfunc,dx,damp=damp,adamp=adamp)
            Phighk = b_high*(highk/k_junc_high)**c_high*np.exp(-a*(highk/k_junc_high))
        else :
            b_high,c_high = get_powerlaw_junc(k_junc_high,Pkfunc,dx)
            Phighk = b_high*(highk/k_junc_high)**c_high
    
        kforextrap = np.concatenate([lowk,setk[(setk>=k_junc_low)&(setk<=k_junc_high)], highk])
        
        Pforextrap = np.concatenate([Plowk,Pk[(setk>=k_junc_low)&(setk<=k_junc_high)],Phighk])
    
        ExtraPk = sp.interpolate.interp1d(kforextrap,Pforextrap,kind='cubic',bounds_error = False, fill_value = 'extrapolate')(setkextrap)
    
        return ExtraPk 
          
def computePS(bvals,datalin,dataloop,setkin,setkout):
    """ Computes the power spectra given the b_i and the EFT power spectra

        Inputs
        ------
        bvals : The values for the b_i
        datalin : the linear power spectra from the EFT, with shape (multipoles, b_i, k)
        dataloop : the loop power spectra from the EFT, with shape (multipoles, b_i, k)
        setkin : the values of k for the intput power spectra
        setkout : the values of k for the output power spectra

        Outputs
        ------
        The power spectra multipoles, non-concatenated
    """
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop
    b1,b2,b3,b4,b5,b6,b7,b8,b9,b10 = bvals
    cvals = np.array([1,b1,b2,b3,b4,b1*b1,b1*b2,b1*b3,b1*b4,b2*b2,b2*b4,b4*b4,b1*b5,b1*b6,b1*b7,b5,b6 ,b7,b8 ,b9,b10])
    
    nkin = len(setkin)/3
    nkout = len(setkout)/3
    
    P0 = (sp.interp1d(setkin[:nkin],np.dot(cvals,data0)+datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2],kind=3,bounds_error = False,fill_value = 'extrapolate'))(setkout[:nkout])
    P2 = (sp.interp1d(setkin[:nkin],np.dot(cvals,data2)+datalin2[0]+b1*datalin2[1]+b1*b1*datalin2[2],kind=3,bounds_error = False,fill_value = 'extrapolate'))(setkout[:nkout])
    P4 = (sp.interp1d(setkin[:nkin],np.dot(cvals,data4)+datalin4[0]+b1*datalin4[1]+b1*b1*datalin4[2],kind=3,bounds_error = False,fill_value = 'extrapolate'))(setkout[:nkout])
    return np.array([P0,P2,P4])

def transformQ(Pkin,setkin,setkout,dataQ,n=2**12,kr=1,extrap=True,setkextrap = 10**(np.linspace(-5,1.2,400)),dampl4=True,**kwargs):
    stab,dummy,Q0,Q2,Q4,Q6,Q8,mueff=dataQ
    
    if not check_if_multipoles_k_array(setkin):
        setkin = np.concatenate([setkin,setkin,setkin])
    nkin = len(Pkin)/3
    
    if extrap:
        P0k = ExtrapolationPk(Pkin[:nkin],setkin[:nkin],setkextrap,**kwargs)
        P2k = ExtrapolationPk(Pkin[nkin:2*nkin],setkin[:nkin],setkextrap,**kwargs)
        P4k = ExtrapolationPk(Pkin[2*nkin:3*nkin],setkin[:nkin],setkextrap,damp=True,**kwargs)
        setkin = np.concatenate([setkextrap,setkextrap,setkextrap])
        nkin =  len(setkin)/3
    else:
        P0k = Pkin[:nkin]
        P2k = Pkin[nkin:2*nkin]
        P4k = Pkin[2*nkin:3*nkin]
    
    
    P0int = scipy.interpolate.interp1d(setkin[:nkin],P0k,kind=3,bounds_error = False, fill_value = 'extrapolate')
    P2int = scipy.interpolate.interp1d(setkin[:nkin],P2k,kind=3,bounds_error = False, fill_value = 'extrapolate')
    P4int = scipy.interpolate.interp1d(setkin[:nkin],P4k,kind=3,bounds_error = False, fill_value = 'extrapolate')
    
    logkmin = np.log(setkin.min())
    logkmax = np.log(setkin.max())
    

    # Bias exponent: q = 0 is unbiased
    q = 0


    # Tell fhti to change kr to low-ringing value
    kropt = 1
    
    # Central point log10(r_c) of periodic interval
    logkc = (logkmin + logkmax)/2

    # Central index (1/2 integral if n is even)
    nc = (n + 1)/2.0

    dlnk=float(logkmax - logkmin)/(n-1)

    k = np.exp(np.linspace(logkmin,logkmax,n))

    P0k=P0int(k)
    P2k=P2int(k)
    P4k=P4int(k)

    kr0, wsave0, ok0 = fftlog.fhti(n, 0.5, dlnk, q, kr, kropt)
    kr2, wsave2, ok2 = fftlog.fhti(n, 2+0.5, dlnk, q, kr, kropt)
    kr4, wsave4, ok4 = fftlog.fhti(n, 4+0.5, dlnk, q, kr, kropt)

    logrc = np.log(kr) - logkc
    r = np.exp(logrc + (np.arange(1, n+1) - nc)*dlnk)


    ak0 = k**1.5*P0k/(2*np.pi)**1.5
    xi0r = fftlog.fht(ak0.copy(), wsave0, -1)/r**1.5
    
    ak2 = k**1.5*P2k/(2*np.pi)**1.5
    xi2r = fftlog.fht(ak2.copy(), wsave2, -1)/r**1.5
    
    ak4 = k**1.5*P4k/(2*np.pi)**1.5
    xi4r = fftlog.fht(ak4.copy(), wsave4, -1)/r**1.5
    
    Q0r=(sp.interp1d(stab,Q0/Q0[0],kind='cubic',bounds_error=False,fill_value='extrapolate'))(r)
    Q2r=(sp.interp1d(stab,Q2/Q0[0],kind='cubic',bounds_error=False,fill_value=0))(r)
    Q4r=(sp.interp1d(stab,Q4/Q0[0],kind='cubic',bounds_error=False,fill_value=0))(r)
    Q6r=(sp.interp1d(stab,Q6/Q0[0],kind='cubic',bounds_error=False,fill_value=0))(r)
    Q8r=(sp.interp1d(stab,Q8/Q0[0],kind='cubic',bounds_error=False,fill_value=0))(r)
    
    
    xihat0=xi0r*Q0r+1./5*xi2r*Q2r+1./9*xi4r*Q4r
    xihat2=xi0r*Q2r+xi2r*(Q0r+2./7*Q2r+2./7*Q4r)+xi4r*(2./7*Q2r+100./693*Q4r+25./143*Q6r)
    xihat4=xi0r*Q4r+xi2r*(18./35*Q2r+20./77*Q4r+45./143*Q6r)+xi4r*(Q0r+20./77*Q2r+162./1001*Q4r+20./143*Q6r+490./2431*Q8r)
    
    ar0=xihat0*r**1.5*(2*np.pi)**1.5
    ar2=xihat2*r**1.5*(2*np.pi)**1.5
    ar4=xihat4*r**1.5*(2*np.pi)**1.5
    
    P0hat=fftlog.fht(ar0.copy(), wsave0, 1)/k**1.5
    P2hat=fftlog.fht(ar2.copy(), wsave2, 1)/k**1.5
    P4hat=fftlog.fht(ar4.copy(), wsave4, 1)/k**1.5
    
    nkout = len(setkout)/3
    P0hatout = (sp.interp1d(k,P0hat,bounds_error = False, fill_value = 'extrapolate'))(setkout[:nkout])
    P2hatout = (sp.interp1d(k,P2hat,bounds_error = False, fill_value = 'extrapolate'))(setkout[:nkout])
    P4hatout = (sp.interp1d(k,P4hat,bounds_error = False, fill_value = 'extrapolate'))(setkout[:nkout])
    
    
    
    return np.array([P0hatout,P2hatout,P4hatout])
     
if __name__ ==  "__main__":
    



    # Model
    kmodel = np.load(opa.join(OUTPATH,'Ploop_fid_PatchyDida.npy'))[0,0]#kclass#[maskk]#
    kmodel3  = np.concatenate([kmodel ,kmodel ,kmodel ])
    Ploopfid = np.load(opa.join(OUTPATH,'Ploop_fid_PatchyDida.npy'))[:,1:]
    Plinfid = np.load(opa.join(OUTPATH,'Plin_fid_PatchyDida.npy'))[:,1:]
    inipos = np.array([2.03395135,  -6.15044179,   1.21410315,   7.19087139,
        11.61423533, -33.46605767,   1.58262629, -44.64033227,
        57.43130091,  26.44292187])##np.array([2]+[0]*9)#

    kmodelfine = np.exp(np.linspace(np.log(kmodel.min()),np.log(kmodel.max()),100))
    kmodelfine3 = np.concatenate([kmodelfine,kmodelfine,kmodelfine])
    kmin = kmodel.min()
    kmax = kmodel.max()
    P0model,P2model,P4model = APpowerspectraNkmu.changetoAPnoNK(computePS(inipos,Plinfid,Ploopfid,kmodel3,kmodelfine3),kmodelfine3,kmodelfine3,1,1)

    
    
    dataQ = np.loadtxt(opa.join(OUTPATH,'DataBispec/W_v4_NGC_mask_DR12cmass_10pc.txt')).T

    Pkin = np.concatenate([P0model,P2model,P4model])
    setkin = kmodelfine3
    setkout = np.concatenate([10**(np.linspace(-5,1.2,400))]*3)
    
    PStransformed = transformQ(Pkin,setkin,setkout,dataQ,n=64*4)
    
    nkout = len(setkout)/3
    nkin = len(setkin)/3
    
    dictcolor = {'0' : 'blue', '1' : 'red', '2' : 'green'}
    plt.figure()
    for l in range(3):
         
        plt.loglog(setkin[l*nkin:(l+1)*nkin],Pkin[l*nkin:(l+1)*nkin], color=dictcolor[str(l)],label='l = ' +str(l) + ' Original')
        plt.loglog(setkout[l*nkout:(l+1)*nkout],PStransformed[l*nkout:(l+1)*nkout], color=dictcolor[str(l)], ls='--',label='l = ' + str(l)+ r' with $Q_\ell$')
    #plt.axvline((2*np.pi)/Lbox,color='k')
    plt.axvline(kmin,color='grey',ls='--')
    plt.axvline(kmax,color='grey',ls='--')
    #plt.text(1,1e5,'kr= ' + str(krtest))
    plt.show()

     