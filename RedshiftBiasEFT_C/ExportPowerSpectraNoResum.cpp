#include "RedshiftBiasEFT.h"

void ExportPowerSpectraNoResum (const string & PathToFolder, const unsigned & order, PowerSpectraNoResum * Ps, const YesNo & MoreKs) {

	size_t Np ;
	if (order == 0) Np = N0 ;
	if (order == 1) Np = N1 ; 

	int pre = 15 ; // white spacing

	for (unsigned int i = 0 ; i < 5 ; i++) {

		ostringstream filename ;
		if (order == 1) filename << PathToFolder << "/PowerSpectra1loopNoResum_l" << 2*i << ".dat" ;
		if (order == 0) filename << PathToFolder << "/PowerSpectraLinearNoResum_l" << 2*i << ".dat" ;

		ofstream write (filename.str(), ios::out | ios::trunc) ;
		/*
		if (order == 0)
    		write << "#" << setw(pre-1) << "k" << setw(pre) << "1" << setw(pre) << "b1" << setw(pre) << "b1*b1" << endl ;

		if (order == 1)
			 write << "#" << setw(pre-1) << "k" <<  setw(pre) << "1" << setw(pre) << "b1"  << setw(pre) << "b2"  << setw(pre) << "b3"  << setw(pre) << "b4"  
			<< setw(pre) << "b1*b1"  << setw(pre) << "b1*b2"  << setw(pre) << "b1*b3"  << setw(pre) << "b1*b4"  << setw(pre) << "b2*b2"  << setw(pre) << "b2*b4"  << setw(pre) << "b4*b4" 
			<< setw(pre) << "b1*b5"  << setw(pre) << "b1*b6"  << setw(pre) << "b1*b7"  
			<< setw(pre) << "b5"  << setw(pre) << "b6"  << setw(pre) << "b7"  << setw(pre) << "b8"  << setw(pre) << "b9"  << setw(pre) << "b10" << endl ;
		*/
		for (unsigned int m = 0 ; m < Nk ; m++) {
			//write << setw(pre) << klist[m] ;
			for (unsigned int n = 0 ; n < Np ; n++) write << setw(pre) << (*Ps)[i][n][m] ;
			write << endl ;
		}
		write.close() ;
	}
}