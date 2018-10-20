#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <gsl/gsl_spline.h>
#include "fftlog.h"
#include <gsl/gsl_sf_bessel.h>

using namespace std;

int main(int argc, char** argv){

	string box = argv[1] ;
	//string box = "ChallengeQuarter_A" ;


	///////////////////////////////////////////////////////////////
	// Read Window Functions in configuration spaces 
	////////////////////////////////////////////////////////////////
	ifstream read("Qs_"+box+".txt", ios::in) ;

	size_t N = 0 ;
	double si, dummy, Q0i, Q2i, Q4i, Q6i, Q8i ;
	
	while (read >> si >> dummy >> Q0i >> Q2i >> Q4i >> Q6i >> Q8i >> dummy) N++ ;

	// reset ifstream buffer at the beginning of the document
	read.clear() ;
	read.seekg(0, ios::beg) ;

	double s[N], W[5][N] ;

	for (unsigned int i = 0 ; i < N ; i++) 
		read >> s[i] >> dummy >> W[0][i] >> W[1][i] >> W[2][i] >> W[3][i] >> W[4][i] >> dummy ;

	read.close() ;
	
	///////////////////////////////////////////////////////////////
	// Read the array of k's of the data
	////////////////////////////////////////////////////////////////

	ifstream readkp("kp_"+box+".txt") ;

	size_t Nk = 0 ;
	
	while (readkp >> dummy) Nk++ ;

	readkp.clear() ;
	readkp.seekg(0, ios::beg) ;

	double kp[Nk] ;

	for (unsigned int i = 0 ; i < Nk ; i++) {
		readkp >> kp[i] ;
	}

	readkp.close() ;

	/*
	size_t Nk = 100 ;
	double kmin = 0.01 ;
	double kmax = 0.4 ;

	double dk = (kmax-kmin)/((double)Nk-1.) ;

	double kp[Nk] ;

	for (unsigned int i = 0 ; i < Nk ; i++) {
		kp[i] = kmin + i*dk ;
	}
	*/

	//////

	size_t Nfft = 1024*4 ; 
	double ds = log(s[N-2]/s[0]) / ((double)Nfft-1.) ;
	double slog[Nfft] ;
	for (unsigned int i = 0 ; i < Nfft ; i++) slog[i] = s[0] * exp(i*ds) ;

	//////
	double k[Nfft] ;

	///////////////////////////////////////////////////////////////
	// 3-j symbol squared times 2A+1 for A = 0,2,4, l = 0,2,4, L = 0,2,4,6,8 (n = 0,2)
	////////////////////////////////////////////////////////////////
	double C0[3][3][5] = { 
		{ {1., 0., 0., 0., 0.}, {0., 1./5., 0., 0., 0.}, {0., 0., 1./9., 0., 0.} },
		{ {0., 1., 0., 0., 0.}, {1., 2./7., 2./7., 0., 0.}, {0., 2./7., 100./693., 25./143., 0.} },
		{ {0., 0., 1., 0., 0.}, {0., 18./35., 20./77., 45./143., 0.}, {1., 20./77., 162./1001., 20./143., 490./2431.} } } ;

	////////

	int pre = 15 ;
	//cout << N << ' ' << Nk << endl ;

	double il, ilp ;

	///////////////////////////////////////////////////////////////
	// Compute Q_l,l' in Fourier space
	////////////////////////////////////////////////////////////////

	for (unsigned int l = 0 ; l < 3 ; l++) { // A = 0,2,4 : Pi_A

		if (l == 1) il = -1. ;
		else il = 1. ;

		for (unsigned int p = 0 ; p < 3 ; p++) { // l = 0,2,4 : Xi_l

			if (p == 1) ilp = -1. ;
			else ilp = 1. ;

			cout << 2*l << ", " << 2*p << endl ;

			double Q[Nk][Nfft] ;

			for (unsigned int j = 0 ; j < Nk ; j++) { // k'

				dcomplex IntegrandPi[Nfft] ;
				dcomplex Pi[Nfft] ;
				for (unsigned int i = 0 ; i < Nfft ; i++) IntegrandPi[i] = 0. ;
				
				for (unsigned int q = 0 ; q < 5 ; q++) { // L = 0,2,4,6,8 : Q_L

					gsl_interp_accel * acc = gsl_interp_accel_alloc () ;
					gsl_spline * spline = gsl_spline_alloc (gsl_interp_cspline, N) ;
					gsl_spline_init (spline, s, W[q], N) ;

					for (unsigned int i = 0 ; i < Nfft ; i++) 
						IntegrandPi[i] += pow(slog[i],1.5) * C0[l][p][q] * gsl_spline_eval (spline, slog[i], acc) * gsl_sf_bessel_jl(2*p, slog[i]*kp[j]) ; // Integrand to Hankel transform

					gsl_spline_free (spline) ;
    				gsl_interp_accel_free (acc) ;
    			}

				fht(Nfft, slog, IntegrandPi, k, Pi, 2*l+0.5, 0.) ; // Hankel transform by fftlog

				for (unsigned int i = 0 ; i < Nfft ; i++) Q[j][i] = il * ilp * pow(k[i],2) * pow(2./M_PI,0.5) * pow(k[i],-1.5) * real(Pi[i]) ;

			}
			/////// Write to file Q_l,l' /////////////////////////

			ostringstream filename ; 
			filename << "Q" << to_string(2*l) << to_string(2*p) << "_" << box << ".dat" ;
			
			ofstream write(filename.str(), ios::out | ios::trunc) ;

			//write << "#" << setw(pre-1) << "(k,kp)" ;
			//for (unsigned int j = 0 ; j < Nk ; j++) write << setw(pre) << kp[j] ;
			//write << endl ;

			for (unsigned int i = 0 ; i < Nfft ; i++) {
				//write << setw(pre) << k[i] ;
				for (unsigned int j = 0 ; j < Nk ; j++) write << setw(pre) << Q[j][i] ;
				write << endl ;
			}

			write.close() ;
		}
	}

	ofstream writek("k.dat", ios::out | ios::trunc) ;
	for (unsigned int i = 0 ; i < Nfft ; i++) writek << k[i] << endl ;
	writek.close() ;



	////////////////////////////////////////////////////////////////

	
	return 0 ;
}