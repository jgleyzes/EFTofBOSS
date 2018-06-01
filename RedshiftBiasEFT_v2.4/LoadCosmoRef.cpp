#include "RedshiftBiasEFT.h"

void LoadCosmoRef (const string & PathToFolderCosmoRef, PowerSpectraNoResum * Ps) { ;

	for (unsigned int i = 0 ; i < Nl ; i++) {

		ostringstream filename ;
		filename << PathToFolderCosmoRef << "/PowerSpectra1loopNoResum_l" << 2*i << ".dat" ; 
		ifstream ref(filename.str(), ios::in) ;

		if (!ref.is_open()) {
			cerr << "There was a problem opening the reference cosmology loop integrals files." <<PathToFolderCosmoRef<< endl ;
			exit(EXIT_FAILURE) ;
		}
		
		else {
			ref.clear() ;
			ref.seekg(0, ios::beg) ;

			for (unsigned int m = 0 ; m < Nk ; m++) {
				for (unsigned int n = 0 ; n < N1-9 ; n++) {
					ref >> (*Ps)[i][n][m] ;
				}
			}

			ref.close() ;
		}

	}
    
}

void ExportPowerSpectraNoResum (const string & PathToFolder, const unsigned & order, PowerSpectraNoResum * Ps) {

	size_t Np = 12 ; 
	if (order == 0) Np = N0 ;
	if (order == 1) Np = N1-9 ;  // -9 counterterms

	for (unsigned int i = 0 ; i < 5 ; i++) {

		ostringstream filename ;

		if (order == 1) filename << PathToFolder << "/PowerSpectra1loopNoResum_l" << 2*i << ".dat" ;
		if (order == 0) filename << PathToFolder << "/PowerSpectraLinearNoResum_l" << 2*i << ".dat" ;

		ofstream write (filename.str(), ios::out | ios::trunc) ;

		for (unsigned int m = 0 ; m < Nk ; m++) {
			for (unsigned int n = 0 ; n < Np ; n++) {
				write << setw(12) << (*Ps)[i][n][m] << " " ;
			}

			write << endl ;
		}

		write.close() ;
	}
}