import numpy as np
import os.path as opa
import subprocess as sp
import configobj as cfg

import camb
from camb import model, initialpower


DEFAULTCONFIG_zbEFT = cfg.ConfigObj({
                'z_pk':1.,
                'omega_b':0.022,
                'omega_cdm':0.118,
                'h':0.7,
                'knl':1.,
                'km':1.,
                'nbar':0.00952380952,
                'PathToLinearPowerSpectrum':'camb_pk.dat',
                'outpath':'output',
                'ComputePowerSpectrum':'yes',
                'ImportResummationMatrix':'no',
                'ExportResummationMatrix':'no',
                'ComputeBispectrum':'no',
                'EpsRel_IntegrBispectrumAP':1e-3,
                'PathToTriangles' : '',
                'aperp':1,
                'apar':1
                })

def run_zbEFT(OUT_PATH, EFT_PATH, lnAs, Om, h, fb, z_pk):

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
    pars.set_matter_power(redshifts=[z_pk], kmax=2)

    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-3, maxkh=2, npoints = 128)
    s8 = np.array(results.get_sigma8())
    
    # Make a directory to put inputs outputs for zbEFT
    outpath = opa.abspath( opa.join( OUT_PATH, 'intermediary/lnAs%sOm%sh%s' % (lnAs, Om, h) ) )
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
    zbeftpath = EFT_PATH
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