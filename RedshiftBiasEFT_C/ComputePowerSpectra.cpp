#include "RedshiftBiasEFT.h"

void ComputePowerSpectra1LoopNoResum (const PrecisionIntegr & Eps, const string & PathToFolder, const string & PathToFolderCosmoRef, 
	const ParametersCosmology & Target, const redshift & z0, const ParamsP11 & params, 
	const YesNo & UseCosmoRef, const YesNo & UVsubP22, const YesNo & MoreKs,
	const double & nbar, const double & km, const double & knl, PowerSpectraNoResum * Ps) {
	ComputeLoopIntegrals (Eps, PathToFolder, PathToFolderCosmoRef, Target, z0, params, UseCosmoRef, UVsubP22, MoreKs, Ps) ;
	ComputeCounterTerms (params, nbar, km, knl, Ps) ;
}

///////

void ComputePowerSpectraLinearNoResum (const ParamsP11 & params, PowerSpectraNoResum * Ps, const YesNo & MoreKs = 0) {
	
	double f1 = params.f ;

	size_t Nk_local = Nk ;
	if (MoreKs == true) Nk_local = Nh ;

	for (unsigned int i = 0 ; i < 3 ; i++) {
		for (unsigned int m = 0 ; m < Nk_local; m++) {
			double k_local = klist[m] ;
			if (MoreKs == true) k_local = khigh[m] ;
			double P11k = P11(k_local,params) ;

			(*Ps)[i][0][m] = m4[i] *f1*f1 *P11k ; 	// 1
			(*Ps)[i][1][m] = m2[i] *2.*f1 *P11k ; 	// b1
			(*Ps)[i][2][m] = m0[i] *P11k ; 			// b1*b1
		}
	}
}