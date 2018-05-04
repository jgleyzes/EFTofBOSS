#include "ResumEFT.h"

void LoadM (const string & PathToFolderRD, const unsigned int & order, const double & k, const unsigned int & l, const unsigned int & lp, InterpFunc & InterpM) {

	ostringstream filename ;

    filename << PathToFolderRD << "M" << order << "_" << k << "_" << l << "_" << lp << ".dat" ;
    
	ifstream data(filename.str(), ios::in) ;

	if (!data.is_open()) { 
		cerr << "Problem loading M file:" << filename.str() << endl ; 
		exit(EXIT_FAILURE) ; 
	}

	else {
		//cout << "Loading M file:" << filename.str() << endl ; 
		double qdata[Nkp] ;
		double Mdata[Nkp] ;

		for (unsigned int i = 0 ; i < Nkp ; i++) data >> qdata[i] >> Mdata[i] ;

		// Interpolation
		InterpM.accel = gsl_interp_accel_alloc () ;
		InterpM.interp = gsl_spline_alloc (gsl_interp_cspline, Nkp) ;
		gsl_spline_init (InterpM.interp, qdata, Mdata, Nkp) ;

		data.close() ;
	}
}

void UnloadInterp (InterpFunc & params) {
	gsl_spline_free (params.interp) ;
    gsl_interp_accel_free (params.accel) ;
}