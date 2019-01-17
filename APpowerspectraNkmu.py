import numpy as np
import scipy
import scipy.interpolate as sp

def check_if_multipoles_k_array(setk):
    return setk[len(setk)/3] == setk[0]
    

def changetoAPbinning(Pk,setkin,setkout,qperp,qpar,TableNkmu,l68 = None): #setk is the k array corresponding to Pk, while x is the output array, the one corresponding to the data.
    """ Applies the Alcock Paczinski effect as well as binning effect

        Inputs
        ------
        Pk : The power spectra multipoles concatenated [P0,P2,P4]
        setkin : the values of k for which each multipoles is evaluated
        setkout : the values of k for which to evaluate the final power spectra (x should be contained in setkin)
        qperp, qpar : the value of the AP parameters
        Tablemu : table computed from the data for each (k,mu) bin with kmean, mucentral, Number of modes
        
        Outputs
        ------
        The transformed power spectra, non-concatenated
    """
    _,kmean,mucent,nkmu = TableNkmu #import the data from the sims. mucent are central values, kmean the mean.
    

    
    
    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkout[:len(setkout)/3]


    # Add l=6,8 contribution
    if type(l68) != type(None):
        Pkloc = np.concatenate([Pk,l68])
    else:
        Pkloc = Pk
    #print (Pkloc.shape)


    # Interpolate the multipoles
    Pkint = scipy.interpolate.interp1d(setkin,Pkloc,axis=-1,kind='cubic',bounds_error = False,fill_value = 'extrapolate')

    
    # Define the grid with the right kmax and kmin and reshape into (k,mu)
    
    kmin = setkout.min()
    kmax = setkout.max()
    kmeanx = kmean[(kmean>= kmin)&(kmean <= kmax)]
    mucentx = mucent[(kmean>= kmin)&(kmean <= kmax)]
   
    Nbink = len(kmeanx)/100
    Nbinmu = 100
    
    kgrid = kmeanx.reshape((Nbink,Nbinmu))   
    mugrid = mucentx.reshape((Nbink,Nbinmu))
    
    # Reshape N(k,mu) on the grid with right kmin and kmax
    
    nkmux = nkmu[(kmean>= kmin)&(kmean <= kmax)]
    nkgrid = nkmux.reshape((Nbink,Nbinmu)) 
    
    # Interpolate the mu part of N(k,mu)
    
    nkgridint = sp.interp1d(mugrid[0,:],nkgrid,axis = 1,kind = 'nearest',bounds_error = False,fill_value = 'extrapolate')
    
    # New array of mu with more points (better precision for the integration)
    
    muacc = np.linspace(0.,1.,1000)
    
    mugrid,kgrid = np.meshgrid(muacc,np.unique(kmeanx))  
    
    # AP factors
    F = float(qpar/qperp)
    k = kgrid/qperp*(1+mugrid**2*(F**-2-1))**0.5
    mup = mugrid/F*(1+mugrid**2*(F**-2-1))**-0.5
    


    # Goes from the multipoles back to P(k,mu) and apply AP
    if type(l68) == type(None):
        arrayLegendremup = nkgridint(muacc)*np.array([scipy.special.legendre(0)(mup),
                                scipy.special.legendre(2)(mup),
                                scipy.special.legendre(4)(mup)])
    else:
        #print ('l68 Legendre')
        arrayLegendremup = nkgridint(muacc)*np.array([scipy.special.legendre(0)(mup),
                                    scipy.special.legendre(2)(mup),
                                    scipy.special.legendre(4)(mup),
                                    scipy.special.legendre(6)(mup),
                                    scipy.special.legendre(8)(mup)])

    
    arrayLegendremugrid = np.array([2*(2*0+1.)/(2*qperp**2*qpar)*scipy.special.legendre(0)(mugrid),
                              2*(2*2.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(2)(mugrid),
                              2*(2*4.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(4)(mugrid)])

    #print (Pkint.shape, arrayLegendremup.shape)
    Pkmu = np.einsum('lkm,lkm->km',Pkint(k),arrayLegendremup)


    # Normalization for N(k,mu)dmu
    
    nk = np.trapz(nkgridint(muacc),x = muacc,axis = 1)
    
    # Back to multipoles (factor of 2 because we integrate an even function from 0 to 1 instead of -1 to 1)
    
    Integrandmu = np.einsum('km,lkm->lkm',Pkmu,arrayLegendremugrid)
    Pk_AP = np.trapz(Integrandmu,x = mugrid,axis = -1)/nk

    
    # interpolate on the wanted k-array for output
    Pk_AP_out =(scipy.interpolate.interp1d(np.unique(kmeanx),Pk_AP,axis=-1,bounds_error = False,fill_value = 'extrapolate'))(setkout) 

    
    return Pk_AP_out
    

def changetoAPnobinning(Pk,setkin,setkout,qperp,qpar,nbinsmu = 500,l68 = None): 
    """ Applies the Alcock Paczinski effect but no binning effect to an list of Pi (for the marginalization).

        Inputs
        ------
        Pik : The list of power spectra multipoles, shape (n_i,3,n_k)
        setkin : the values of k for which each multipoles is evaluated
        setkout : the values of k for which to evaluate the final power spectra (x should be contained in setk)
        qperp, qpar : the value of the AP parameters
        nbinsmu : number of mu-bins for the integration (default 500)
        
        Outputs
        ------
        The transformed power spectra, non-concatenated 
        """
    
    muacc = np.linspace(0.,1.,nbinsmu)

    
    
    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkout[:len(setkout)/3]
    #print(Pk.shape, " is pk shape") 


    # Add l=6,8 contribution
    if type(l68) != type(None):
        Pkloc = np.concatenate([Pk,l68])
    else:
        Pkloc = Pk


    # Interpolate the multipoles
    Pkint = scipy.interpolate.interp1d(setkin,Pkloc,axis=-1,kind='cubic',bounds_error = False,fill_value = 'extrapolate')

    
    # Define the grid with the right kmax and kmin and reshape into (k,mu)

    kgrid,mugrid = np.meshgrid(setkout,muacc,indexing='ij')
    
    
    # AP factors
    F = float(qpar/qperp)
    k = kgrid/qperp*(1+mugrid**2*(F**-2-1))**0.5
    mup = mugrid/F*(1+mugrid**2*(F**-2-1))**-0.5

    
   # Goes from the multipoles back to P(k,mu) and apply AP
    if type(l68) == type(None):
        arrayLegendremup = np.array([scipy.special.legendre(0)(mup),
                                scipy.special.legendre(2)(mup),
                                scipy.special.legendre(4)(mup)]) 
    else:
        arrayLegendremup = np.array([scipy.special.legendre(0)(mup),
                                    scipy.special.legendre(2)(mup),
                                    scipy.special.legendre(4)(mup),
                                    scipy.special.legendre(6)(mup),
                                    scipy.special.legendre(8)(mup)])
        
    arrayLegendremugrid = np.array([2*(2*0+1.)/(2*qperp**2*qpar)*scipy.special.legendre(0)(mugrid),
                                  2*(2*2.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(2)(mugrid),
                                  2*(2*4.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(4)(mugrid)])


    Pkmu = np.einsum('lkm,lkm->km',Pkint(k),arrayLegendremup)
    #(P0int(k)*scipy.special.legendre(0)(mup)+P2int(k)*scipy.special.legendre(2)(mup)+P4int(k)*scipy.special.legendre(4)(mup))

    
    # Back to multipoles (factor of 2 because we integrate an even function from 0 to 1 instead of -1 to 1)
    #print(Pkmu.shape, arrayLegendremugrid.shape) 
    Integrandmu = np.einsum('km,lkm->lkm',Pkmu,arrayLegendremugrid)
    
    Pk_AP = np.trapz(Integrandmu,x = mugrid,axis = -1)

    return Pk_AP




def changetoAPnobinningPi(Pik,setkin,setkout,qperp,qpar,nbinsmu = 500): 
    """ Applies the Alcock Paczinski effect but no binning effect to an list of Pi (for the marginalization).

        Inputs
        ------
        Pik : The list of power spectra multipoles, shape (n_i,3,n_k)
        setkin : the values of k for which each multipoles is evaluated
        setkout : the values of k for which to evaluate the final power spectra (x should be contained in setk)
        qperp, qpar : the value of the AP parameters
        nbinsmu : number of mu-bins for the integration (default 500)
        
        Outputs
        ------
        The transformed power spectra, non-concatenated 
        """
    
    muacc = np.linspace(0.,1.,nbinsmu)

    
    
    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkout[:len(setkout)/3]
    
    # Interpolate the multipoles
    
    Piint = scipy.interpolate.interp1d(setkin,Pik,axis=-1,kind='cubic',bounds_error = False,fill_value = 'extrapolate')

    
    # Define the grid with the right kmax and kmin and reshape into (k,mu)

    
    kgrid,mugrid = np.meshgrid(setkout,muacc,indexing='ij')
    
    
    # AP factors
    F = float(qpar/qperp)
    k = kgrid/qperp*(1+mugrid**2*(F**-2-1))**0.5
    mup = mugrid/F*(1+mugrid**2*(F**-2-1))**-0.5

    
    # Goes from the multipoles back to P(k,mu) and apply AP
    
    arrayLegendremup = np.array([scipy.special.legendre(0)(mup),
                                scipy.special.legendre(2)(mup),
                                scipy.special.legendre(4)(mup)])
    
    arrayLegendremugrid = np.array([2*(2*0+1.)/(2*qperp**2*qpar)*scipy.special.legendre(0)(mugrid),
                              2*(2*2.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(2)(mugrid),
                              2*(2*4.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(4)(mugrid)])

    Pikmu = np.einsum('ilkm,lkm->ikm',Piint(k),arrayLegendremup)
    #(P0int(k)*scipy.special.legendre(0)(mup)+P2int(k)*scipy.special.legendre(2)(mup)+P4int(k)*scipy.special.legendre(4)(mup))

    
    # Back to multipoles (factor of 2 because we integrate an even function from 0 to 1 instead of -1 to 1)
    
    Integrandmu = np.einsum('ikm,lkm->ilkm',Pikmu,arrayLegendremugrid)
    
    Pi_AP = np.trapz(Integrandmu,x = mugrid,axis = -1)

    

    return Pi_AP


def changetoAPbinningPi(Pik,setkin,setkout,qperp,qpar,TableNkmu): #setk is the k array corresponding to Pk, while x is the output array, the one corresponding to the data.
    """ Applies the Alcock Paczinski effect as well as binning effect

        Inputs
        ------
        Pk : The list of power spectra multipoles, shape (n_i,3,n_k)
        setkin : the values of k for which each multipoles is evaluated
        setkout : the values of k for which to evaluate the final power spectra (x should be contained in setkin)
        qperp, qpar : the value of the AP parameters
        Tablemu : table computed from the data for each (k,mu) bin with kmean, mucentral, Number of modes
        
        Outputs
        ------
        The transformed power spectra, non-concatenated
    """
    _,kmean,mucent,nkmu = TableNkmu #import the data from the sims. mucent are central values, kmean the mean.
    

    
    
    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkout[:len(setkout)/3]
    
    # Interpolate the multipoles
    
    Pkint = scipy.interpolate.interp1d(setkin,Pik,axis=-1,kind='cubic',bounds_error = False,fill_value = 'extrapolate')
    
    # Define the grid with the right kmax and kmin and reshape into (k,mu)
    
    kmin = setkout.min()
    kmax = setkout.max()
    kmeanx = kmean[(kmean>= kmin)&(kmean <= kmax)]
    mucentx = mucent[(kmean>= kmin)&(kmean <= kmax)]
   
    Nbink = len(kmeanx)/100
    Nbinmu = 100
    
    kgrid = kmeanx.reshape((Nbink,Nbinmu))   
    mugrid = mucentx.reshape((Nbink,Nbinmu))
    
    # Reshape N(k,mu) on the grid with right kmin and kmax
    
    nkmux = nkmu[(kmean>= kmin)&(kmean <= kmax)]
    nkgrid = nkmux.reshape((Nbink,Nbinmu)) 
    
    # Interpolate the mu part of N(k,mu)
    
    nkgridint = sp.interp1d(mugrid[0,:],nkgrid,axis = 1,kind = 'nearest',bounds_error = False,fill_value = 'extrapolate')
    
    # New array of mu with more points (better precision for the integration)
    
    muacc = np.linspace(0.,1.,1000)
    
    mugrid,kgrid = np.meshgrid(muacc,np.unique(kmeanx))  
    
    # AP factors
    F = float(qpar/qperp)
    k = kgrid/qperp*(1+mugrid**2*(F**-2-1))**0.5
    mup = mugrid/F*(1+mugrid**2*(F**-2-1))**-0.5
    
    # Goes from the multipoles back to P(k,mu) and apply AP
    
    arrayLegendremup = nkgridint(muacc)*np.array([scipy.special.legendre(0)(mup),
                                scipy.special.legendre(2)(mup),
                                scipy.special.legendre(4)(mup)])
    
    arrayLegendremugrid = np.array([2*(2*0+1.)/(2*qperp**2*qpar)*scipy.special.legendre(0)(mugrid),
                              2*(2*2.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(2)(mugrid),
                              2*(2*4.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(4)(mugrid)])
    
    Pkmu = np.einsum('ilkm,lkm->ikm',Pkint(k),arrayLegendremup)


    # Normalization for N(k,mu)dmu
    
    nk = np.trapz(nkgridint(muacc),x = muacc,axis = 1)
    
    # Back to multipoles (factor of 2 because we integrate an even function from 0 to 1 instead of -1 to 1)
    
    Integrandmu = np.einsum('ikm,lkm->ilkm',Pkmu,arrayLegendremugrid)
    Pk_AP = np.trapz(Integrandmu,x = mugrid,axis = -1)/nk

    
    # interpolate on the wanted k-array for output
    Pk_AP_out =(scipy.interpolate.interp1d(np.unique(kmeanx),Pk_AP,axis=-1,bounds_error = False,fill_value = 'extrapolate'))(setkout) 

    
    return Pk_AP_out