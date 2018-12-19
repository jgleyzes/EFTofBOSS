#include "RedshiftBiasEFT.h"
#include "ResumEFT.h"

/*
#include <ctime>
#include <ratio>
#include <chrono>


using namespace std::chrono;
*/

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
		YesNo UseRef = 0, ImportM = 0, ExportM = 0, ComputePowerSpectrum = 1, UVsubP22 = 1, MoreKs = 0 ;
		string PathToFolder = "./" ;
		string PathToFolderCosmoRef ;
		string PathToLinearPowerSpectrum ;
		PrecisionIntegr Eps = { { 0.1,1e-3 } , { 0.1,1e-2 } } ;

		LoadConfigFile (argv[1], nbar, km, knl, z0, cosmo, PathToFolder, PathToFolderCosmoRef, PathToLinearPowerSpectrum, ComputePowerSpectrum, UseRef, ImportM, ExportM, UVsubP22, MoreKs, Eps) ;

		///////////////////////////
		//
		///////////////////////////

		ParamsP11 paramsP11 ;
		LoadP11 (PathToLinearPowerSpectrum, cosmo, z0, paramsP11) ;

		PowerSpectraNoResum Ps1Loop ;
		PowerSpectraNoResum PsLinear ;

		static StoreM TableM ;

		//steady_clock::time_point start = steady_clock::now();

		if (ImportM == false) ResumM (PathToFolder, paramsP11, ExportM, &TableM) ;

		if (ComputePowerSpectrum == true) {
			ComputePowerSpectraLinearNoResum  (paramsP11, &PsLinear, false) ;
			ComputePowerSpectra1LoopNoResum (Eps, PathToFolder, PathToFolderCosmoRef, cosmo, z0, paramsP11, UseRef, UVsubP22, false, nbar, km, knl, &Ps1Loop) ;

			//ExportPowerSpectraNoResum (PathToFolder, 0, &PsLinear, false) ;
			//ExportPowerSpectraNoResum (PathToFolder, 1, &Ps1Loop, false) ;

			ResumPowerSpectra (PathToFolder, paramsP11, &PsLinear, &Ps1Loop, ImportM, ExportM, &TableM) ;

			if (MoreKs == true) {

				cerr << "Evaluating at more k points..." << endl ;

				PowerSpectraNoResum Ph1Loop ;
				PowerSpectraNoResum PhLinear ;

				ComputePowerSpectraLinearNoResum  (paramsP11, &PhLinear, true) ;
				ComputePowerSpectra1LoopNoResum (Eps, PathToFolder, PathToFolderCosmoRef, cosmo, z0, paramsP11, UseRef, UVsubP22, true, nbar, km, knl, &Ph1Loop) ;

				ExportPowerSpectraNoResum (PathToFolder, 0, &PhLinear, true) ;
				ExportPowerSpectraNoResum (PathToFolder, 1, &Ph1Loop, true) ;
			}

		}
		
		/*
		steady_clock::time_point stop = steady_clock::now();
		duration<double> runtime = duration_cast<duration<double>>(stop-start);
		cout << "RedshiftBiasEFT ran in " << runtime.count() << " seconds." << endl ;
		*/
		
		///////////////////////////
		///
		///////////////////////////
	}

	return 0 ;

}