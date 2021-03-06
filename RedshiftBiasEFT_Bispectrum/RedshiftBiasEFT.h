#ifndef REDSHIFTBIASEFT_H
#define REDSHIFTBIASEFT_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring> 
#include <sstream>
#include <fstream>
#include <vector>
#include <complex>

#include <gsl/gsl_spline.h>

#include <cuba.h>

#include "klist.h"

using namespace std ;

// Math stuff
double Heaviside (const double & a) ;
const double Pi = M_PI ;
const double E = exp(1.) ;



// Struct declarations
struct InterpFunc {
	gsl_interp_accel * accel ;
	gsl_spline * interp ;
} ;

struct ParamsP11 {
	double f ;
	gsl_interp_accel * accel ;
	gsl_spline * interp ;
} ;

struct ParamsIntegr {
	double k ;
	ParamsP11 p11 ;
	int id ;
	bool UseCosmoRef ;
	ParamsP11 p11ref ;
	double RescaleFactor ;
} ;

typedef double redshift ;

// Precision of the evalution of the loop integrals: 2 for UseCosmoRef = Yes/No, 2 for EpsAbs,EpsRel 
typedef double PrecisionIntegr[2][2] ; 


/** Cosmology **/
const size_t Nc = 5 ; // Number of cosmological parameters 
typedef double ParametersCosmology[Nc] ; // cosmological parameters
const string ParametersCosmologyNames[Nc] = { "A_s", "n_s", "h", "omega_b", "omega_cdm" } ; // A_s refers to ln(10^10 A_s)

/* Reference cosmology: Planck2015 */
const ParametersCosmology Reference = { 3.094, 0.9645, 0.6727, 0.02225, 0.1198 } ;

// Linear Growth rate f: GrowthFunction.cpp
double LinearGrowthRate (const ParametersCosmology & p, const redshift & z);

// LinearPowerSpectrum.cpp
void LoadP11 (const string & LinearPowerSpectrumData, const ParametersCosmology & cosmo, const redshift & z, ParamsP11 & p) ;
double P11 (const double & q, const ParamsP11 & params) ;
void UnloadP11 (const ParamsP11 & params) ;


// Integrands.cpp
double Integrands (const int j, const double & k, const double & q, const double & x, const ParamsP11 & p) ;
const size_t NI = 22 ; // Number of Integrals
const double qUV = 1. ; // the UV term propto (k/q)^2 are subtracted up to this lower bound. 

// Multipole expansion
const size_t Nl = 5 ;
typedef double MultipoleMoments[Nl] ; // l = 0, 2, 4, 6, 8
// We call Mi with i: power of mu
const MultipoleMoments m0 = { 1., 0., 0., 0., 0. } ;
const MultipoleMoments m2 = { 1./3., 2./3., 0., 0., 0. } ;
const MultipoleMoments m4 = { 1./5., 4./7., 8./35., 0., 0. } ;
const MultipoleMoments m6 = { 1./7., 10./21., 24./77., 16./231., 0. } ;
const MultipoleMoments m8 = { 1./9., 40./99., 48./148., 64./495., 128./6435. } ;

// For each integral we associate a pointer to multipole moments.
const MultipoleMoments * const MultipoleExpansion[NI] = { &m2, &m2, &m2, &m4, &m0, &m0, &m0, &m0, &m0, &m0, &m0, &m2, &m4, &m6, &m2, &m4, &m2, &m2, &m4, &m4, &m6, &m8 } ;

// ComputePowerSpectra.cpp
const size_t N0 = 3 ; 	// 3 linear terms
const size_t N1 = 21 ; // 12 1-Loop PowerSpectra + 9 CounterTerms
typedef double PowerSpectraNoResum[Nl][N1][Nk] ;

void ComputePowerSpectraLinearNoResum (const ParamsP11 & params, PowerSpectraNoResum * Ps) ;
void ComputePowerSpectra1LoopNoResum (const PrecisionIntegr & Eps, const string & PathToFolder,const string & PathToFolderRD,const string & PathToFolderCosmoRef, const ParametersCosmology & Target, const redshift & z0, const ParamsP11 & params, const bool & UseCosmoRef, const double & Nbar, const double & kM, const double & kNL, PowerSpectraNoResum *) ;


// ComputeLoopIntegrals.cpp
int ComputeLoopIntegrals (const PrecisionIntegr & Eps, const string & PathToFolder,const string & PathToFolderRD,const string & PathToFolderCosmoRef, const ParametersCosmology & Target, const redshift & z0, const ParamsP11 & params, const bool & UseCosmoRef, PowerSpectraNoResum *) ;
static int Integrand_LoopIntegrals_CUBA (const int *ndim, const double a[], const int *ncomp, double ff[], void *params) ;

// ComputeCounterTerms.cpp
void ComputeCounterTerms (const ParamsP11 & params, const double & Nbar, const double & kM, const double & kNL, PowerSpectraNoResum * Ps) ;

// ComputeP22UV.cpp
double P22UV (const ParamsP11 & params) ;
double Integrand_P22UV (double q, void * params) ;

// LoadCosmoRef.cpp
void ExportPowerSpectraNoResum (const string & PathToFolder, const unsigned & order, PowerSpectraNoResum * Ps) ;
void LoadCosmoRef (const string & PathToFolderCosmoRef, PowerSpectraNoResum * Ps) ;

// LoadConfigFile.cpp
typedef bool YesNo ;
void LoadConfigFile (char * ConfigFile, double & nbar, double & km, double & knl, redshift & z0, ParametersCosmology & cosmo, 
	string & PathToFolder, string & PathToFolderRD,string & PathToFolderCosmoRef, string & PathToLinearPowerSpectrum, string & PathToTriangles, 
	YesNo & ComputePowerSpectrum, YesNo & ComputeBispectrum, YesNo & UseRef, YesNo & ImportM, YesNo & ExportM, 
	PrecisionIntegr & Eps, double & EpsRel_IntegrBispectrumAP, double & aperp, double & apar) ;


///////////// TREE-LEVEL BISPECTRUM WITH AP EFFECT //////////////////////


struct ParamsIntegrBispectrumAP {
	double k1,k2,k3 ;
	ParamsP11 p11 ;
	double nbar ;
	double aperp, apar ;
	int id ;
} ;

// IntegrandBispectrumAP.cpp
double IntegrandBispectrumAP (const int & id, const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, const ParamsP11 & p, const double & nbar, const double & aperp, const double & apar) ;


void ScoccimarroTransform (const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, 
	double & q1, double & q2, double & q3, double & nu1, double & nu2, double & nu3) ;

void ScoccimarroTransformWithAP (const double & aperp, const double & apar, 
	const double & k1, const double & k2, const double & k3, const double & mu1, const double & phi, 
	double & q1, double & q2, double & q3, double & nu1, double & nu2, double & nu3) ;

// ComputeBispectrum.cpp
void ComputeBispectrumMonopole (const double & EpsRel, const string & PathToFolder, const string & PathToTriangles, const ParamsP11 & params, const double & nbar, const double & aperp, const double & apar) ;
static int IntegrandBispectrumAP_CUBA (const int *ndim, const double a[], const int *ncomp, double ff[], void *params) ;



#endif