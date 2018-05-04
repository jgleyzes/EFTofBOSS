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
# Save this path (expect data and other needed executables present in relative paths!!!)
THIS_PATH = opa.dirname(__file__)

# Change directory for this run... (is this really neccessary?)
#os.chdir(opa.abspath(opa.join(THIS_PATH,'../MCMC')))

# Data paths
INPATH = opa.abspath(opa.join(THIS_PATH,'input')) 
OUTPATH = opa.abspath(opa.join(THIS_PATH,'output')) 
import APpowerspectraNkmu
FFTLOG_PATH = '/Users/jgleyzes/Documents/Projets/EFTofLSS/Python/fftlog-master/' # Expects the EFT code compiled here!!!
if opa.isdir(FFTLOG_PATH):
    import sys
    if FFTLOG_PATH not in sys.path: 
        sys.path.append(FFTLOG_PATH)
    #import fftlog
else:
    raise Exception('Module not found at ' + FFTLOG_PATH)
    
def cH(Om,a):
    return np.sqrt(Om/a+a**2*(1-Om))
def DgN(Om,a):
    
    return 5./2*Om*cH(Om,a)/a*scipy.integrate.quad(lambda x:cH(Om,x)**-3,0,a)[0]
def fN(Om,a):
    return (Om*(5*a - 3*DgN(Om,a)))/(2.*(a**3*(1 - Om) + Om)*DgN(Om,a))
    
def check_if_multipoles_k_array(setk):
    return setk[len(setk)/3] == setk[0]
                
def Classtomultipoles(setk,Pkclass,Om,zpk,b1=1):
    kmulti = np.concatenate([setk,setk,setk])
    f = fN(Om,1./(1+zpk))
    #f = 0.51
    P0k = (b1**2 + (2*f)*b1/3. + f**2/5.)*Pkclass
    P2k = ((4*f*(7*b1 + 3*f))/21.)*Pkclass
    P4k = ((8*f**2)/35.)*Pkclass
    return np.array([kmulti,np.concatenate([P0k,P2k,P4k])])

def get_powerlaw_junc(kjunc,Pfunc,dlnx=0.1,damp=False):
    
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
    
    dx = dlnx*kjunc
    P = Pfunc(kjunc)
    P1 = scipy.misc.derivative(Pfunc, kjunc, dx=dx, n=1, args=(), order=3)
    P2 = scipy.misc.derivative(Pfunc, kjunc, dx=dx, n=2, args=(), order=5)
    
    
    if damp :
            a = (kjunc*(kjunc*P1**2 - P*(P1 + kjunc*P2)))/P**2
            b = float(P)
            c = (kjunc**2*(P1**2 - P*P2))/P**2
            return a,b,c
    else :
            b = P
            c = P1*kjunc/b
            return b,c
            

 
def damptanh(k,ktr,sig):
    return ((1 + np.tanh((-k + ktr)/sig))/2.)
            
def damptanhlog(k,ktr,sig):
    return ((1 + np.tanh((-np.log(k) + np.log(ktr))/sig))/2.)
            
                
def TBBKS(k,a):
    return (0.08986547601495114*np.log(1 + 11.12774387166911*a*k))/(a*k*(1 + 18.481165260549865*a*k + 5896.387755102042*a**2*k**2 + 17504.580955408142*a**3*k**3 + 1.0337805625000012e6*a**4*k**4)**0.25)
    

def fitBBKSk(k_junc,Pfunc,dx):
    P = Pfunc(k_junc)
    P1 = scipy.misc.derivative(Pfunc, k_junc, dx=dx, n=1, args=(), order=5)
    P2 = scipy.misc.derivative(Pfunc, k_junc, dx=dx, n=2, args=(), order=7)
    def bc(a):
        T = TBBKS(k_junc,a)
        T1 = scipy.misc.derivative(lambda x:TBBKS(x,a), k_junc, dx=dx, n=1, args=(), order=5)
        b = P/T**2
        c = k_junc*(P1/P - (2*T1)/T)
        return b,c
        
    def functoroot(a):
        T = TBBKS(k_junc,a)
        T1 = scipy.misc.derivative(lambda x:TBBKS(x,a), k_junc, dx=dx, n=1, args=(), order=5)
        T2 = scipy.misc.derivative(lambda x:TBBKS(x,a), k_junc, dx=dx, n=2, args=(), order=7)
        
        return -1 + (k_junc*(P1**2/P + (2*P*(-(k_junc*T1**2) + T*(T1 + k_junc*T2)))/(k_junc*T**2)))/(P1 + k_junc*P2)
    result = scipy.optimize.root(functoroot,1.1)
    a = result.x[0]
    b,c = bc(a)
    return a,b,c
    

def ExtrapolationPk(Pk,setk,setkextrap,k_junc_low = 0.01,k_junc_high=0.4,ktr=2,sig=0.5,withlog=False,damp = False):
    
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
        raise Exception('The junction point at low k must be within initial range. kmin = {} and k_junc = {} '.format(kmin,k_junc_low))
    elif k_junc_high > kmax:
        raise Exception('The junction point at high k must be within initial range. kmax = {} and k_junc = {} '.format(kmax,k_junc_high))
   
    else :
        Pkfunc = sp.interpolate.interp1d(setk,Pk,kind='cubic',bounds_error = False, fill_value = 'extrapolate')
        
        
        
        lowk = setkextrap[setkextrap < k_junc_low]
        dx = 0.01*k_junc_low
        a_low,b_low,c_low = fitBBKSk(k_junc_low,Pkfunc,dx)
        Plowk = np.array([TBBKS(k,a_low)**2 for k in lowk])*b_low*(lowk/k_junc_low)**c_low#
        
        highk = setkextrap[setkextrap > k_junc_high]
        
        if damp :
            a_high,b_high,c_high = get_powerlaw_junc(k_junc_high,Pkfunc,damp=damp)
            Phighk = b_high*(highk/k_junc_high)**(c_high)*np.exp(-a_high*(highk-k_junc_high)/k_junc_high)
            print(a_high,b_high,c_high)
        else :
            a_high,b_high,c_high = fitBBKSk(k_junc_high,Pkfunc,dx)#b_high,c_high = get_powerlaw_junc(k_junc_high,Pkfunc,dx)
            Phighk = np.array([TBBKS(k,a_high)**2 for k in highk])*b_high*(highk/k_junc_high)**c_high#
            
        kforextrap = np.concatenate([lowk,setk[(setk>=k_junc_low)&(setk<=k_junc_high)], highk])
        
        Pforextrap = np.concatenate([Plowk,Pk[(setk>=k_junc_low)&(setk<=k_junc_high)],Phighk])
        
        if withlog:
            dampfunction = np.array([damptanhlog(k,ktr,sig) for k in kforextrap])
        else:
            dampfunction = np.array([damptanh(k,ktr,sig) for k in kforextrap])
        ExtraPk = sp.interpolate.interp1d(kforextrap,Pforextrap*dampfunction,kind='cubic',bounds_error = False, fill_value = 'extrapolate')(setkextrap)
    
        return ExtraPk 
          
def computePS(bvals,datalin,dataloop,setkin,setkout,withsq=1):
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
    
    P11int =  sp.interpolate.interp1d(setkin[:nkin], datalin[0,-1])
    sigsq = withsq*2*scipy.integrate.quad(lambda x: x**2/(2*np.pi)**2*(P11int(x))**2,setkin.min(),setkin.max())[0]
    
    P0 = (sp.interp1d(setkin[:nkin],np.dot(cvals,data0)+datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2]-2*(-b1 + b2 + b4)**2*sigsq,kind=2,bounds_error = False,fill_value = 'extrapolate'))(setkout[:nkout])
    P2 = (sp.interp1d(setkin[:nkin],np.dot(cvals,data2)+datalin2[0]+b1*datalin2[1]+b1*b1*datalin2[2],kind=2,bounds_error = False,fill_value = 'extrapolate'))(setkout[:nkout])
    P4 = (sp.interp1d(setkin[:nkin],np.dot(cvals,data4)+datalin4[0]+b1*datalin4[1]+b1*b1*datalin4[2],kind=2,bounds_error = False,fill_value = 'extrapolate'))(setkout[:nkout])
    return np.array([P0,P2,P4])
   
    
    
def Ptoxi(Pkin,setkin,n=64*64,extrap=True,dampl4=False,setkextrap = 10**(np.linspace(-5,1.2,400)),**kwargs):
    
    """ Computes the correlation function corresponding to the power spectra

        Inputs
        ------
        Pkin : (array) The input power spectra (concatenated multipoles)
        setkin : (array) the range of k on which the input power spectra in computed
        dataQ : (array) contains the real space multipoles of the window functions (see eq 14 of https://arxiv.org/pdf/1511.07799.pdf)
        n : (integer) the number of step of the internal log spaced k-range for the FFTlog (default is max, 4096)
        extrap : (boolean) whether to extrapolate with power laws at high and small k. Done by matching value and slope at k_junc_low and k_junc_high.
        setkextrap : (log-spaced array) the k range for the extrapolation
        dampL4 : whether to damp the extrapolation of L4 (otherwise it can create artifacts in the FFTlog)
        kwargs : keyword arguments to be passed for the extrapolation (such as k_junc_low and k_junc_high)
        
        Outputs
        ------
        The power spectra multipoles with window functon applied, non-concatenated
    """
    
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    
    nkin = len(setkin)
    if extrap:
        P0k = ExtrapolationPk(Pkin[:nkin],setkin[:nkin],setkextrap,damp=dampl4,**kwargs)
        P2k = ExtrapolationPk(Pkin[nkin:2*nkin],setkin[:nkin],setkextrap,damp=dampl4,**kwargs)
        P4k = ExtrapolationPk(Pkin[2*nkin:3*nkin],setkin[:nkin],setkextrap,damp=dampl4,**kwargs)
        setkin = np.concatenate([setkextrap,setkextrap,setkextrap])
        nkin =  len(setkin)/3
    else:
        P0k = Pkin[:nkin]
        P2k = Pkin[nkin:2*nkin]
        P4k = Pkin[2*nkin:3*nkin]
    
    P0int = scipy.interpolate.interp1d(setkin[:nkin],P0k,kind=2,bounds_error = False, fill_value = 'extrapolate')
    P2int = scipy.interpolate.interp1d(setkin[:nkin],P2k,kind=2,bounds_error = False, fill_value = 'extrapolate')
    P4int = scipy.interpolate.interp1d(setkin[:nkin],P4k,kind=2,bounds_error = False, fill_value = 'extrapolate')
    
    logkmin = np.log(setkin.min())
    logkmax = np.log(setkin.max())
    logrmin = np.log(1)-logkmax
    logrmax = np.log(1)-logkmin
    k = np.exp(np.linspace(logkmin,logkmax,n))
    r = np.exp(np.linspace(logrmin,logrmax,n))
    
    P0_logspaced = P0int(k)
    P2_logspaced = P2int(k)
    P4_logspaced = P4int(k)
    
    r, xi0r = mcfit.cosmology.P2xi(k, l=0, q=1.5, N=None, lowring=True)(P0_logspaced)
    r, xi2r = mcfit.cosmology.P2xi(k, l=2, q=1.5, N=None, lowring=True)(P2_logspaced)
    r, xi4r = mcfit.cosmology.P2xi(k, l=4, q=1.5, N=None, lowring=True)(P4_logspaced)
    

    
    return np.array([r, xi0r,xi2r,xi4r])
    
if __name__ ==  "__main__":
    
    dfcosmo = pd.read_csv(opa.join(INPATH,'DataFrameCosmosims.csv'),index_col=0)
    simtype = 'LightConeHectorTrue'
    ZONE = 'NGC'
    Seriescosmo = dfcosmo.loc[simtype]
    gridname = Seriescosmo.loc['gridname']
    cosmo = Class()
    
    paramfid = {'z_pk' : Seriescosmo['z_pk'],
    'format' : 'camb',
    'omega_cdm' :Seriescosmo['Omega_m']*Seriescosmo['h']**2-Seriescosmo['omega_b'], 
    'omega_b' :Seriescosmo['omega_b'] ,
    'n_s' : Seriescosmo['ns'],
    'P_k_max_h/Mpc' : 30,
    'ln10^{10}A_s' : Seriescosmo['lnAs'],
    'h' : Seriescosmo['h'],
    'output' : 'mPk'}
    cosmo.set(paramfid)
    cosmo.compute()
    setkclass = 10**(np.linspace(-5,np.log10(30),1000))
    pkclass= np.array(map(lambda x: cosmo.pk(x,paramfid['z_pk']),setkclass*paramfid['h']))*paramfid['h']**3
    cosmo.struct_cleanup()
    cosmo.empty()
    
    simtype = 'LightConeHector'
    
    kclass,PSclass = Classtomultipoles(setkclass,pkclass,(paramfid['omega_cdm']+paramfid['omega_b'])/paramfid['h']**2,paramfid['z_pk'],b1=2)
    PSextra0,PSextra2,PSextra4 = PSclass.reshape((3,len(setkclass)))
    setkextrap = 10**(np.linspace(-5,1.2,800))
    # Model
    kmodelor = np.load(opa.join(OUTPATH,'TablePloopTrue_highkZ_%s.npy'%(gridname)))[0,0]#kclass#[maskk]#'Ploop_fidmoreZhigh_PatchyDida.npy'
    
    Ploopfidor = np.load(opa.join(OUTPATH,'TablePloopTrue_highkZ_%s.npy'%(gridname)))[:,1:]
    Plinfidor = np.load(opa.join(OUTPATH,'TablePlinTrue_highkZ_%s.npy'%(gridname)))[:,1:]

    kmodel = kmodelor
    kmodel3  = np.concatenate([kmodel ,kmodel ,kmodel ])
    
    Ploopfid = Ploopfidor
    Plinfid = Plinfidor
    

    
    
    
    inipos = np.array([   1.92188101,   -3.92210966,   -7.55391499,    4.89261009,
        -13.79322164,  -18.77971407,   -2.39509908,   0,
        0,0.      ])
    #inipos = [1.89466591,   -9.01241666,    2.71013798,    9.04969125,
    #     18.2062663 , -177.19343897,  199.99999695] + [0]*3
    kmodelfine = np.exp(np.linspace(np.log(kmodel.min()),np.log(kmodel.max()),400))
    kmodelfine3 = np.concatenate([kmodelfine,kmodelfine,kmodelfine])
    kmin = kmodel.min()
    kmax = kmodel.max()
    P0model,P2model,P4model = computePS(inipos,Plinfid,Ploopfid,kmodel3,kmodel3)
    
    
    kjunhigh = 0.8
    ktr = 2
    sig = 0.5
    withlog = False
    dampl4 = False
    
    
    Pkin = np.concatenate([P0model,P2model,P4model])#PSclass[(kclass<=kmax)]#np.concatenate([P0model,P2model,P4model])#np.concatenate([P0model,P2model,P4model])#PSclass#np.concatenate([P0model,P2model,P4model])##
    setkin = kmodel3#kclass[(kclass<=kmax)]#kmodel3#kmodel3#kPS#[(kclass>=kmin)]#&(kclass<=kmax)
    
    r,xi0,xi2,xi4 = Ptoxi(Pkin,setkin,n=2**12,extrap=True,setkextrap = setkextrap,k_junc_low=setkin.min(),k_junc_high=kjunhigh,ktr=ktr,sig=sig,withlog=withlog,dampl4=dampl4)#,k_junc_low=setkin.min(),k_junc_high=setkin.max())

    xil = np.array([xi0,xi2,xi4])
    
    rclass,xi0class,xi2class,xi4class = np.loadtxt(opa.join(OUTPATH,'xiclass_%s.dat'%simtype))
    xilclass = np.array([xi0class,xi2class,xi4class])
    
    
    dictcolor = {'0' : 'blue', '1' : 'red', '2' : 'green'}
    plt.figure()
    plt.suptitle(simtype)
    for l in range(3):
         
        plt.plot(r,r**2*xil[l], color=dictcolor[str(l)],label=r'$\ell$ = ' +str(2*l) ,lw=2)
        plt.plot(rclass,rclass**2*xilclass[l], color=dictcolor[str(l)],ls='--',label = r'$\ell$ = %s (CLASS)' %(2*l) )
    
    rD,xi0D = np.loadtxt(opa.join(OUTPATH,'CFdata/dr12v6c_rmu_cmass_NS_z43z75000%s.mono'%1)).T[:,1:]
    xi0av = np.zeros(xi0D.shape)
    xi2av = np.zeros(xi0D.shape)
    xi4av = np.zeros(xi0D.shape)
    for i in range(1,10):
     rD,xi0D = np.loadtxt(opa.join(OUTPATH,'CFdata/dr12v6c_rmu_cmass_NS_z43z75000%s.mono'%i)).T[:,1:]
     rD,xi2D = np.loadtxt(opa.join(OUTPATH,'CFdata/dr12v6c_rmu_cmass_NS_z43z75000%s.quad'%i)).T[:,1:]
     rD,xi4D = np.loadtxt(opa.join(OUTPATH,'CFdata/dr12v6c_rmu_cmass_NS_z43z75000%s.four'%i)).T[:,1:]
     xi0av += xi0D
     xi2av += xi2D
     xi4av += xi4D
     xilD = np.array([xi0D,xi2D,xi4D]) 
     xilav = np.array([xi0av,xi2av,xi4av])/9   
    for l in range(3):
        #pass 
        plt.plot(rD,rD**2*xilav[l], label=r'$\ell$ = ' +str(2*l) + ' Data avg',lw=2,color='cyan' )
        
        
        #plt.plot(rb1,rb1**2*xilb1[l], color=dictcolor[str(l)],ls='-.',label = r'$\ell$ = %s (b1)' %(2*l) )
    plt.xlim(0.01,200)
    plt.ylim(-100,200)
    plt.xlabel(r'r',fontsize=15)
    plt.ylabel(r'$r^2\xi_\ell$',fontsize=15)
    plt.legend()
    
    plt.axvline(2*np.pi/kmin,color='grey',ls='--')
    plt.axvline(2*np.pi/0.2,color='grey',ls='--')
    #plt.text(1,1e5,'kr= ' + str(krtest))
    plt.show()

     