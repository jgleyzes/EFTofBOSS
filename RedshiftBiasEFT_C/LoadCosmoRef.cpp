#include "RedshiftBiasEFT.h"

void LoadCosmoRef (const string & PathToFolderCosmoRef, PowerSpectraNoResum * Ps) { ;

	for (unsigned int i = 0 ; i < Nl ; i++) {

		ostringstream filename ;
		filename << PathToFolderCosmoRef << "/PowerSpectra1loopNoResum_l" << 2*i << ".dat" ; 
		ifstream ref(filename.str(), ios::in) ;

		if (!ref.is_open()) {
			cerr << "There was a problem opening the reference cosmology loop integrals files." << endl ;
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