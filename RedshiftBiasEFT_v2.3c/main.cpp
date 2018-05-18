#include "RedshiftBiasEFT.h"
#include "ResumEFT.h"

#include <ctime>

int main(int argc, char *argv[]) {

	if (argc <= 1) {
		cerr << "Error: no configuration file specified." << endl ;
		exit(EXIT_FAILURE) ;
	}
	
	else {

		// Default values
		double nbar = 1./105. , km = 1. , knl = 1. ;
		redshift z0 = 0.67 ;
		ParametersCosmology cosmo ; for (unsigned int i = 0 ; i < Nc ; i++) cosmo[i] = Reference[i] ;
		YesNo UseRef = 0, ImportM = 1, ExportM = 0, ComputePowerSpectrum = 1 , ComputeBispectrum = 1 ;
		string PathToFolder = "./" ;
		string PathToFolderRD ;
		string PathToFolderCosmoRef ;
		string PathToLinearPowerSpectrum ;
		string PathToTriangles ; 
		PrecisionIntegr Eps = { { 0.1,1e-3 } , { 0.1,1e-2 } } ;
		double aperp = 1., apar = 1. ;

		double EpsRel_IntegrBispectrumAP = 1e-3 ; 

		LoadConfigFile (argv[1], nbar, km, knl, z0, cosmo, \
			PathToFolder,PathToFolderRD,PathToFolderCosmoRef, PathToLinearPowerSpectrum, PathToTriangles, 
			ComputePowerSpectrum, ComputeBispectrum, UseRef, ImportM, ExportM, 
			Eps, EpsRel_IntegrBispectrumAP, aperp, apar) ;

		///////////////////////////
		//
		///////////////////////////

		ParamsP11 paramsP11 ;
		LoadP11 (PathToLinearPowerSpectrum, cosmo, z0, paramsP11) ;


		PowerSpectraNoResum Ps1Loop ;
		PowerSpectraNoResum PsLinear ;

		static StoreM TableM ;


		int start_s=clock() ;

		///////////// ONE-LOOP POWER SPECTRUM WITH RESUMMATION IN REDSHIFT SPACE (3 FIRST MULTIPOLES) ///////
		if (ComputePowerSpectrum == true) {

			if (ImportM == false) ResumM (PathToFolderRD, paramsP11, ExportM, &TableM) ;

			ComputePowerSpectraLinearNoResum  (paramsP11, &PsLinear) ;
			ComputePowerSpectra1LoopNoResum (Eps, PathToFolder,PathToFolderRD,PathToFolderCosmoRef, cosmo, z0, paramsP11, UseRef, nbar, km, knl, &Ps1Loop) ;

			//ExportPowerSpectraNoResum (PathToFolder, 0, &PsLinear) ;
			ExportPowerSpectraNoResum (PathToFolder, 1, &Ps1Loop) ;

			ResumPowerSpectra (PathToFolder,PathToFolderRD, paramsP11, &PsLinear, &Ps1Loop, ImportM, ExportM, &TableM) ;
		}


		///////////// TREE-LEVEL BISPECTRUM WITH AP EFFECT /////////////////////////////////////////////////
		if (ComputeBispectrum == true) {
			ComputeBispectrumMonopole (EpsRel_IntegrBispectrumAP, PathToFolder, PathToTriangles, paramsP11, nbar, aperp, apar) ;
		}

		int stop_s=clock() ;

		cout << "RedshiftBiasEFT ran in " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " seconds." << endl ;

		
		///////////////////////////
		///
		///////////////////////////
	}

	return 0 ;

}
