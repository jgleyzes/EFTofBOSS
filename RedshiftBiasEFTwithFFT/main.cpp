#include "RedshiftBiasEFT.h"
#include "ResumEFT.h"


#include <ctime>
#include <ratio>
#include <chrono>

using namespace std::chrono;


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
		YesNo ImportM = 0, ExportM = 0, ComputePowerSpectrum = 1 ;
		string PathToFolder = "./" ;
		string PathToLinearPowerSpectrum ;

		LoadConfigFile (argv[1], nbar, km, knl, z0, cosmo, PathToFolder, PathToLinearPowerSpectrum, ComputePowerSpectrum, ImportM, ExportM) ;

		///////////////////////////
		//
		///////////////////////////

		ParamsP11 paramsP11 ;
		LoadP11 (PathToLinearPowerSpectrum, cosmo, z0, paramsP11) ;

		PowerSpectraNoResum Ps1Loop ;
		PowerSpectraNoResum PsLinear ;

		static StoreM TableM ;

		steady_clock::time_point start = steady_clock::now();

		if (ImportM == false) ResumM (PathToFolder, paramsP11, ExportM, &TableM) ;

		steady_clock::time_point stop2 = steady_clock::now();
		duration<double> runtime2 = duration_cast<duration<double>>(stop2-start);
		cout << "Resummation matrices computed in " << runtime2.count() << " seconds." << endl ;

		if (ComputePowerSpectrum == true) {
			ComputePowerSpectraLinearNoResum  (paramsP11, &PsLinear) ;
			ComputePowerSpectra1LoopNoResum (paramsP11, nbar, km, knl, &Ps1Loop) ;

			ExportPowerSpectraNoResum (PathToFolder, 0, &PsLinear) ;
			ExportPowerSpectraNoResum (PathToFolder, 1, &Ps1Loop) ;

			ResumPowerSpectra (PathToFolder, paramsP11, &PsLinear, &Ps1Loop, ImportM, ExportM, &TableM) ;

		}

		
		steady_clock::time_point stop = steady_clock::now();
		duration<double> runtime = duration_cast<duration<double>>(stop-start);
		cout << "RedshiftBiasEFTwithFFT ran in " << runtime.count() << " seconds." << endl ;
		
		
		///////////////////////////
		///
		///////////////////////////
	}

	return 0 ;

}