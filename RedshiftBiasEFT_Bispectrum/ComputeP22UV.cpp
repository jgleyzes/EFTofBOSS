#include "RedshiftBiasEFT.h"
#include <cuba.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#define EPSABS 0.001
#define EPSREL 1e-5

double P22UV (const ParamsP11 & params) {
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;
	gsl_function F ;
  	F.function = &Integrand_P22UV ;
  	ParamsP11 local = params ;
  	F.params = &local ;
	double res, err ;
	gsl_integration_qag (&F, CutIR, CutUV, EPSABS, EPSREL, 1000, GSL_INTEG_GAUSS41, w, &res, &err) ; 
	gsl_integration_workspace_free (w) ;
	return res ; 
}

double Integrand_P22UV (double q, void * params) {
	return pow(q,2)/pow(M_PI,2) * pow(P11(q,*(ParamsP11 *) params),2) ;
}