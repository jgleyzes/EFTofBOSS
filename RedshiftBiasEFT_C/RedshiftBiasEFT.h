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
const ParametersCosmology Reference = { 3.094, 0.9645, 0.6727, 0.04917, 0.2647 } ;

// Linear Growth rate f: GrowthFunction.cpp
double LinearGrowthRate (const ParametersCosmology & p, const redshift & z) ;

// LoadConfigFile.cpp
typedef bool YesNo ;
void LoadConfigFile (char * ConfigFile, double & nbar, double & km, double & knl, redshift & z0, ParametersCosmology & cosmo, 
	string & PathToFolder, string & PathToFolderCosmoRef, string & PathToLinearPowerSpectrum, 
	YesNo & ComputePowerSpectrum, YesNo & UseRef, YesNo & ImportM, YesNo & ExportM, YesNo & UVsubP22, YesNo & MoreKs,
	PrecisionIntegr & Eps) ;

// LinearPowerSpectrum.cpp
void LoadP11 (const string & LinearPowerSpectrumData, const ParametersCosmology & cosmo, const redshift & z, ParamsP11 & p) ;
double P11 (const double & q, const ParamsP11 & params) ;
void UnloadP11 (const ParamsP11 & params) ;


// Integrands.cpp
typedef double Integrand (const double & k, const double & q, const double & x, const ParamsP11 & p) ;

Integrand I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12, I13, I14, I15, I16, I17, I18, I19, I20, I21, I22, I23, I24, I25, I26, I27, I28, I29, I30, I31, I32, I33, I34, I35, I36, I37, I38, I39, I40, I41, I42, I43 ;
// We found out that many are zero: only 22 among the original 43 are used.

const size_t NI = 22 ; 
Integrand * const Integrands[NI] = { I1, I2, I3, I4, I5, I6, I7, I8, I10, I12, I17, I20, I21, I22, I24, I25, I28, I32, I33, I41, I42, I43 } ;


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

void ComputePowerSpectraLinearNoResum (const ParamsP11 & params, PowerSpectraNoResum * Ps, const YesNo & MoreKs) ;
void ComputePowerSpectra1LoopNoResum (const PrecisionIntegr & Eps, const string & PathToFolder, const string & PathToFolderCosmoRef, 
	const ParametersCosmology & Target, const redshift & z0, const ParamsP11 & params, 
	const YesNo & UseCosmoRef, const YesNo & UVsubP22, const YesNo & MoreKs,
	const double & Nbar, const double & kM, const double & kNL, PowerSpectraNoResum *) ;


// ComputeLoopIntegrals.cpp
int ComputeLoopIntegrals (const PrecisionIntegr & Eps, const string & PathToFolder, const string & PathToFolderCosmoRef, 
	const ParametersCosmology & Target, const redshift & z0, const ParamsP11 & params, 
	const YesNo & UseCosmoRef, const YesNo & UVsubP22, const YesNo & MoreKs,
	PowerSpectraNoResum *) ;
static int Integrand_LoopIntegrals_CUBA (const int *ndim, const double a[], const int *ncomp, double ff[], void *params) ;


// ComputeCounterTerms.cpp
void ComputeCounterTerms (const ParamsP11 & params, const double & Nbar, const double & kM, const double & kNL, PowerSpectraNoResum * Ps) ;

// ComputeP22UV.cpp
double P22UV (const ParamsP11 & params) ;
double Integrand_P22UV (double q, void * params) ;

// LoadCosmoRef.cpp
void ExportPowerSpectraNoResum (const string & PathToFolder, const unsigned & order, PowerSpectraNoResum * Ps, const YesNo & MoreKs) ;
void LoadCosmoRef (const string & PathToFolder, PowerSpectraNoResum * Ps) ;

#endif
