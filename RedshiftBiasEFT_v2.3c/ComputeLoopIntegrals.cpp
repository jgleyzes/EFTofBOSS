#include "RedshiftBiasEFT.h"
#include <cuba.h>

#define NDIM 2
#define NCOMP 1
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 10000000
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define KEY 0

#define NVEC 1
#define VERBOSE 0
#define LAST 4
#define SPIN NULL

int ComputeLoopIntegrals (const PrecisionIntegr & Eps, const string & PathToFolder,const string & PathToFolderRD,const string & PathToFolderCosmoRef, const ParametersCosmology & Target, const redshift & z0, const ParamsP11 & params, const bool & UseCosmoRef, PowerSpectraNoResum * Ps) {
	
	double epsabs = Eps[0][0] ;
	double epsrel = Eps[0][1] ;

	if (UseCosmoRef == 1) { 
		epsabs = Eps[1][0] ;
		epsrel = Eps[1][1] ;
	}

	double Integrals[NI][Nk] ; 

	ParamsIntegr integr ;
	integr.p11 = params ;
	integr.UseCosmoRef = UseCosmoRef ;
	integr.RescaleFactor = 1. ;



	if (UseCosmoRef == 1) {

		// We get the cosmological parameters and the redshift of the cosmo_ref
		redshift z_ref ; ParametersCosmology cosmo_ref ;
		
		//plenty of dummy thingy
		double d1,d2,d3 ; string d4,d5,d5bis,d5ter, d5quater ; YesNo d6,d6bis,d7,d8,d9 ; PrecisionIntegr d10 ; double d11,d12,d13 ;

		string cosmo_ref_config_str = PathToFolderCosmoRef + "/cosmo_ref.ini" ; 
		char *cosmo_ref_config = new char[cosmo_ref_config_str.length() + 1] ; 
		strcpy(cosmo_ref_config, cosmo_ref_config_str.c_str()) ;
		LoadConfigFile (cosmo_ref_config, d1, d2, d3, z_ref, cosmo_ref, d4, d5,d5bis,d5ter, d5quater, d6, d6bis,d7, d8, d9, d10, d11, d12, d13) ;
		
		// If it happens that the cosmology and redshift under evaluation are the same as the ones of cosmo_ref, we are done.
		if ( cosmo_ref[0] == Target[0] && cosmo_ref[1] == Target[1] && cosmo_ref[2] == Target[2] && cosmo_ref[3] == Target[3] && cosmo_ref[4] == Target[4] && z_ref == z0) {
			LoadCosmoRef(PathToFolderCosmoRef, Ps) ; 
			return 0 ;
		}

		string cosmo_ref_pk_data = PathToFolderCosmoRef + "/cosmo_ref_pk.dat" ; 
		LoadP11 ( cosmo_ref_pk_data, cosmo_ref, z_ref, integr.p11ref) ;
		integr.RescaleFactor = pow( P11(0.2,integr.p11)/P11(0.2,integr.p11ref) , 2 ) ;
	}


	

	int neval, fail, nregions ;

	///// 

	// Loop on the ks defined in "klist.h"
	unsigned int j = 0 ;
	while ( j < Nk) {

		integr.k = klist[j] ;

		for (unsigned int id = 0 ; id < NI ; id++) {
			integr.id = id ;
			double res, err, prob ; 

			double epsrel_scaled = epsrel ;
			if (UseCosmoRef == 1) epsrel_scaled *= ps1D_errrel[j]*ps1D_inverrelmax ;			

			// CUBA v3
			//Cuhre(NDIM, NCOMP, Integrand_LoopIntegrals_CUBA, &integr, epsrel, EPSABS, 0, MINEVAL, MAXEVAL, KEY, &nregions, &neval, &fail, &res, &err, &prob) ;
			//Vegas(NDIM, NCOMP, Integrand_LoopIntegrals_CUBA, &integr, epsrel, EPSABS, 0, SEED, MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE, &neval, &fail, &res, &err, &prob) ;
			
			// CUBA v4.2
			Cuhre (NDIM,NCOMP,Integrand_LoopIntegrals_CUBA, &integr, NVEC, epsrel_scaled, epsabs, VERBOSE, MINEVAL, MAXEVAL, KEY, STATEFILE, SPIN, &nregions, &neval, &fail, &res, &err, &prob) ;

			Integrals[id][j] = res ;
		}

		j++ ;
	}

	///// CosmoRef ////////////////////////////
	PowerSpectraNoResum Pref ;
	if (UseCosmoRef == true) {
		LoadCosmoRef(PathToFolderCosmoRef, &Pref) ; 
	}
	///////////////////////////////////////////


	for (unsigned int i = 0 ; i < Nl ; i++) {

    	for (unsigned int m = 0 ; m < Nk ; m++) {
    		(*Ps)[i][0][m] = (*MultipoleExpansion[19])[i]*Integrals[19][m] + (*MultipoleExpansion[20])[i]*Integrals[20][m] + (*MultipoleExpansion[21])[i]*Integrals[21][m] ;	// *1
    		(*Ps)[i][1][m] = (*MultipoleExpansion[11])[i]*Integrals[11][m] + (*MultipoleExpansion[12])[i]*Integrals[12][m] + (*MultipoleExpansion[13])[i]*Integrals[13][m] ;	// *b1
    		(*Ps)[i][2][m] = (*MultipoleExpansion[14])[i]*Integrals[14][m] + (*MultipoleExpansion[15])[i]*Integrals[15][m] ;													// *b2
    		(*Ps)[i][3][m] = (*MultipoleExpansion[16])[i]*Integrals[16][m] ;																									// *b3
    		(*Ps)[i][4][m] = (*MultipoleExpansion[17])[i]*Integrals[17][m] + (*MultipoleExpansion[18])[i]*Integrals[18][m] ;													// *b4
    		(*Ps)[i][5][m] = (*MultipoleExpansion[0])[i]*Integrals[0][m] + (*MultipoleExpansion[3])[i]*Integrals[3][m] + (*MultipoleExpansion[4])[i]*Integrals[4][m] ;			// *b1*b1
    		(*Ps)[i][6][m] = (*MultipoleExpansion[1])[i]*Integrals[1][m] + (*MultipoleExpansion[5])[i]*Integrals[5][m] ;														// *b1*b2
    		(*Ps)[i][7][m] = (*MultipoleExpansion[6])[i]*Integrals[6][m] ;																										// *b1*b3
    		(*Ps)[i][8][m] = (*MultipoleExpansion[2])[i]*Integrals[2][m] + (*MultipoleExpansion[7])[i]*Integrals[7][m] ;														// *b1*b4
    		(*Ps)[i][9][m] = (*MultipoleExpansion[8])[i]*Integrals[8][m] ;																										// *b2*b2
    		(*Ps)[i][10][m] = (*MultipoleExpansion[9])[i]*Integrals[9][m] ;																										// b2*b4
    		(*Ps)[i][11][m] = (*MultipoleExpansion[10])[i]*Integrals[10][m] ;																									// b4*b4
    	
    		///// CosmoRef ////////////////////////////
    		if (UseCosmoRef == true) for (unsigned int n = 0 ; n < 12 ; n++) (*Ps)[i][n][m] += integr.RescaleFactor * Pref[i][n][m] ;
    		///////////////////////////////////////////
    	}
	}

	return 0 ;
}


static int Integrand_LoopIntegrals_CUBA (const int *ndim, const double a[], const int *ncomp, double ff[], void *params) {

	ParamsIntegr p = *(ParamsIntegr *) params ;

	double q = CutIR + (CutUV-CutIR)*a[0] ;
	double x = -1.+2.*a[1] ;

	if (p.UseCosmoRef == 0) ff[0] = (CutUV-CutIR)*2. * pow(q,2)/pow(2.*M_PI,3) * (*Integrands[p.id]) (p.k, q, x, p.p11) ;
	else if (p.UseCosmoRef == 1) ff[0] = (CutUV-CutIR)*2. * pow(q,2)/pow(2.*M_PI,3) * ( (*Integrands[p.id]) (p.k, q, x, p.p11) - p.RescaleFactor * (*Integrands[p.id]) (p.k, q, x, p.p11ref) ) ;

	return 0 ;
}


