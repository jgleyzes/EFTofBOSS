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

#define EPSREL 0.0001
#define EPSABS 0.1

void P22UV (const string & PathToFolder, const ParamsP11 & params) {

	ParamsP11 local = params ;

	int neval, fail, nregions ;
	double res, err, prob ; 
	Cuhre (NDIM,NCOMP,Integrand_P22UV_CUBA, &local, NVEC, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, KEY, STATEFILE, SPIN, &nregions, &neval, &fail, &res, &err, &prob) ;

    ostringstream filename ;
	filename << PathToFolder << "/P22UV.dat" ;
    ofstream write (filename.str(), ios::out | ios::trunc) ;
    write << "#   *2*(-b1+b2+b4)^2" << endl ;
    write << res ;
    write.close() ;
}

static int Integrand_P22UV_CUBA (const int *ndim, const double a[], const int *ncomp, double ff[], void *params) {

	ParamsP11 p = *(ParamsP11 *) params ;

	double q = CutIR + (CutUV-CutIR)*a[0] ;
	double x = -1.+2.*a[1] ;

	ff[0] = (CutUV-CutIR)*2. * pow(q,2)/pow(2.*M_PI,2) * pow(P11(q,p),2) ;

	return 0 ;
}