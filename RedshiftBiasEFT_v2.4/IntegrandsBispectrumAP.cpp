#include "RedshiftBiasEFT.h"

double IntegrandBispectrumAP (const int & id, const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, const ParamsP11 & p, const double & nbar, const double & aperp, const double & apar) {
	
	double f1 = p.f ; 

	double q1, q2, q3, nu1, nu2, nu3 ;
	
	//ScoccimarroTransform (k1, k2, k3, mu, phi, q1, q2, q3, nu1, nu2, nu3) ;
	//double APfactor = 1. ;

	ScoccimarroTransformWithAP (aperp, apar, k1, k2, k3, mu, phi, q1, q2, q3, nu1, nu2, nu3) ;
	double APfactor = 1./(pow(apar,2)*pow(aperp,4)) ;

	switch(id) {
		case 0:
			return APfactor * ( (pow(f1,3)*pow(q1,-2)*pow(q2,-2)*pow(q3,-1)*(P11(q1,p)*pow(nu1,2)*(q3*P11(q2,p)*pow(nu2,2)*pow(nu1*q1 + nu2*q2,2)*((q2 + q3)*pow(q1,2) + q2*(q2*q3 - pow(q2,2) + pow(q3,2)))*(28*f1*nu1*nu2*q2*pow(q1,3) + 11*pow(q1,4) + 11*pow(q2,4) + pow(q1,2)*(34*pow(q2,2) - 15*pow(q3,2)) + 14*f1*nu1*nu2*q1*q2*(2*pow(q2,2) - pow(q3,2)) - 15*pow(q2,2)*pow(q3,2) + 4*pow(q3,4)) + q2*P11(q3,p)*pow(nu3,2)*(2*pow(q1,2) + 2*pow(q2,2) - pow(q3,2))*(14*f1*nu1*nu3*q2*(q2 + q3)*pow(q1,3) + (4*q2 + 7*q3)*pow(q1,4) - 14*f1*nu1*nu3*q1*pow(q2,2)*(-(q2*q3) + pow(q2,2) - pow(q3,2)) + pow(q1,2)*(-8*pow(q2,3) + 20*q2*pow(q3,2) + 7*pow(q3,3)) + q2*(-7*q3*pow(q2,3) + 4*pow(q2,4) - 8*pow(q2,2)*pow(q3,2) + 7*q2*pow(q3,3) + 4*pow(q3,4)))*pow(nu1*q1 + nu3*q3,2))*pow(2*pow(q1,2) + 2*pow(q2,2) - pow(q3,2),-1)*pow((q2 + q3)*pow(q1,2) + q2*(q2*q3 - pow(q2,2) + pow(q3,2)),-1) - q1*P11(q2,p)*P11(q3,p)*pow(nu2,2)*pow(nu3,2)*(-7*(2*f1*nu2*nu3*q2 + q3)*pow(q1,4) + 4*pow(q1,5) + 7*q3*pow(q2,2)*(pow(q2,2) + pow(q3,2)) - 2*pow(q1,3)*(-7*f1*nu2*nu3*q2*q3 + 4*pow(q2,2) + 4*pow(q3,2)) + 7*pow(q1,2)*(2*f1*nu2*nu3*pow(q2,3) + 2*f1*nu2*nu3*q2*pow(q3,2) + pow(q3,3)) + 2*q1*(7*f1*nu2*nu3*q3*pow(q2,3) + 2*pow(q2,4) + 10*pow(q2,2)*pow(q3,2) + 2*pow(q3,4)))*pow(nu2*q2 + nu3*q3,2)*pow(-(q3*pow(q1,2)) + pow(q1,3) - q3*pow(q2,2) - q1*(pow(q2,2) + pow(q3,2)),-1)))/14. ) ;
			break ;

		case 1:
			return APfactor * ( (pow(f1,2)*pow(q1,-2)*pow(q2,-2)*pow(q3,-1)*(P11(q1,p)*(q3*P11(q2,p)*((q2 + q3)*pow(q1,2) + q2*(q2*q3 - pow(q2,2) + pow(q3,2)))*(2*nu1*nu2*q2*(14*f1*pow(nu1,4) + 11*pow(nu2,2) + pow(nu1,2)*(11 + 28*f1*pow(nu2,2)))*pow(q1,5) + (11*pow(nu1,4) + 25*pow(nu1,2)*pow(nu2,2))*pow(q1,6) + pow(nu2,2)*pow(q2,2)*(pow(q2,2) - pow(q3,2))*((25*pow(nu1,2) + 11*pow(nu2,2))*pow(q2,2) - (11*pow(nu1,2) + 4*pow(nu2,2))*pow(q3,2)) + pow(q1,4)*((pow(nu1,4)*(34 + 84*f1*pow(nu2,2)) + 11*pow(nu2,4) + pow(nu1,2)*(87*pow(nu2,2) + 84*f1*pow(nu2,4)))*pow(q2,2) - 3*pow(nu1,2)*(5*pow(nu1,2) + 12*pow(nu2,2))*pow(q3,2)) + 2*nu1*nu2*q2*pow(q1,3)*(2*(7*f1*pow(nu1,4) + pow(nu2,2)*(17 + 7*f1*pow(nu2,2)) + pow(nu1,2)*(17 + 28*f1*pow(nu2,2)))*pow(q2,2) - (7*f1*pow(nu1,4) + 15*pow(nu2,2) + pow(nu1,2)*(15 + 14*f1*pow(nu2,2)))*pow(q3,2)) + 2*nu1*nu2*q1*q2*((pow(nu2,2)*(11 + 14*f1*pow(nu2,2)) + pow(nu1,2)*(11 + 28*f1*pow(nu2,2)))*pow(q2,4) - (pow(nu2,2)*(15 + 7*f1*pow(nu2,2)) + pow(nu1,2)*(15 + 14*f1*pow(nu2,2)))*pow(q2,2)*pow(q3,2) + 4*(pow(nu1,2) + pow(nu2,2))*pow(q3,4)) + pow(q1,2)*((pow(nu1,4)*(11 + 84*f1*pow(nu2,2)) + 34*pow(nu2,4) + pow(nu1,2)*(87*pow(nu2,2) + 84*f1*pow(nu2,4)))*pow(q2,4) - 3*(2*pow(nu1,2)*pow(nu2,2)*(12 + 7*f1*pow(nu2,2)) + pow(nu1,4)*(5 + 14*f1*pow(nu2,2)) + 5*pow(nu2,4))*pow(q2,2)*pow(q3,2) + pow(nu1,2)*(4*pow(nu1,2) + 11*pow(nu2,2))*pow(q3,4))) + q2*P11(q3,p)*(2*pow(q1,2) + 2*pow(q2,2) - pow(q3,2))*(pow(nu1,2)*(7*q3*(pow(nu1,2) + 2*pow(nu3,2)) + q2*(4*pow(nu1,2) + 11*pow(nu3,2)))*pow(q1,6) + 2*nu1*nu3*pow(q1,5)*(q2*q3*(7*f1*pow(nu1,4) + 4*pow(nu3,2) + 2*pow(nu1,2)*(2 + 7*f1*pow(nu3,2))) + 7*f1*pow(nu1,2)*(pow(nu1,2) + 2*pow(nu3,2))*pow(q2,2) + 7*(pow(nu1,2) + pow(nu3,2))*pow(q3,2)) + pow(q1,4)*(7*q3*pow(nu1,2)*pow(nu3,2)*(1 + 6*f1*(pow(nu1,2) + pow(nu3,2)))*pow(q2,2) - (8*pow(nu1,4) + 15*pow(nu1,2)*pow(nu3,2))*pow(q2,3) + 2*q2*(pow(nu1,4)*(10 + 21*f1*pow(nu3,2)) + pow(nu1,2)*pow(nu3,2)*(19 + 21*f1*pow(nu3,2)) + 2*pow(nu3,4))*pow(q3,2) + 7*(pow(nu1,4) + 3*pow(nu1,2)*pow(nu3,2) + pow(nu3,4))*pow(q3,3)) - 2*nu1*nu3*pow(q1,3)*(q3*(-7*f1*pow(nu1,4) + 8*pow(nu3,2) + pow(nu1,2)*(8 - 14*f1*pow(nu3,2)))*pow(q2,3) + 7*f1*pow(nu1,2)*(pow(nu1,2) + 2*pow(nu3,2))*pow(q2,4) - 7*f1*(pow(nu1,4) + 4*pow(nu1,2)*pow(nu3,2) + pow(nu3,4))*pow(q2,2)*pow(q3,2) - q2*(2*pow(nu1,2)*(10 + 7*f1*pow(nu3,2)) + pow(nu3,2)*(20 + 7*f1*pow(nu3,2)))*pow(q3,3) - 7*(pow(nu1,2) + pow(nu3,2))*pow(q3,4)) + q2*pow(nu3,2)*(pow(q2,2) - pow(q3,2))*(-7*q3*pow(nu1,2)*pow(q2,3) + 7*pow(nu1,2)*pow(q2,4) + (-3*pow(nu1,2) + 4*pow(nu3,2))*pow(q2,2)*pow(q3,2) - 7*q2*(pow(nu1,2) + pow(nu3,2))*pow(q3,3) - 4*(pow(nu1,2) + pow(nu3,2))*pow(q3,4)) + 2*nu1*nu3*q1*q2*q3*(-7*q3*(pow(nu3,2) + pow(nu1,2)*(1 + 2*f1*pow(nu3,2)) + f1*pow(nu3,4))*pow(q2,3) + 4*(pow(nu1,2) + pow(nu3,2))*pow(q2,4) + (pow(nu3,2)*(-8 + 7*f1*pow(nu3,2)) + 2*pow(nu1,2)*(-4 + 7*f1*pow(nu3,2)))*pow(q2,2)*pow(q3,2) + 7*q2*(pow(nu3,2) + pow(nu1,2)*(1 + 2*f1*pow(nu3,2)) + f1*pow(nu3,4))*pow(q3,3) + 4*(pow(nu1,2) + pow(nu3,2))*pow(q3,4)) + pow(q1,2)*(-7*q3*pow(nu1,2)*(pow(nu1,2) + 2*pow(nu3,2) + 6*f1*pow(nu1,2)*pow(nu3,2) + 6*f1*pow(nu3,4))*pow(q2,4) + (4*pow(nu1,4) - 3*pow(nu1,2)*pow(nu3,2))*pow(q2,5) + 2*(pow(nu1,2) + pow(nu3,2))*(-4*pow(nu3,2) + pow(nu1,2)*(-4 + 21*f1*pow(nu3,2)))*pow(q2,3)*pow(q3,2) + 7*pow(nu1,2)*(pow(nu1,2) + 3*pow(nu3,2) + 6*f1*pow(nu1,2)*pow(nu3,2) + 6*f1*pow(nu3,4))*pow(q2,2)*pow(q3,3) + q2*(4*pow(nu1,4) + 31*pow(nu1,2)*pow(nu3,2) + 20*pow(nu3,4))*pow(q3,4) + 7*pow(nu3,2)*(pow(nu1,2) + pow(nu3,2))*pow(q3,5))))*pow(2*pow(q1,2) + 2*pow(q2,2) - pow(q3,2),-1)*pow((q2 + q3)*pow(q1,2) + q2*(q2*q3 - pow(q2,2) + pow(q3,2)),-1) - q1*P11(q2,p)*P11(q3,p)*(-7*q3*pow(nu2,2)*pow(nu3,2)*pow(q1,6) + 7*pow(nu2,2)*pow(nu3,2)*pow(q1,7) + pow(q1,5)*(8*nu2*nu3*q2*q3*(pow(nu2,2) + pow(nu3,2)) + (4*pow(nu2,4) - 3*pow(nu2,2)*pow(nu3,2))*pow(q2,2) + 2*pow(nu3,2)*(-5*pow(nu2,2) + 2*pow(nu3,2))*pow(q3,2)) - 7*pow(q1,4)*(q3*pow(nu2,2)*(pow(nu2,2) + 2*pow(nu3,2) + 6*f1*pow(nu2,2)*pow(nu3,2) + 6*f1*pow(nu3,4))*pow(q2,2) + 2*f1*nu3*pow(nu2,3)*(pow(nu2,2) + 2*pow(nu3,2))*pow(q2,3) + 2*nu2*nu3*q2*(pow(nu2,2) + pow(nu3,2) + 2*f1*pow(nu2,2)*pow(nu3,2) + f1*pow(nu3,4))*pow(q3,2) + pow(nu3,4)*pow(q3,3)) + 7*q3*pow(q2,2)*(2*nu2*nu3*q3*(pow(nu2,2) + pow(nu3,2))*pow(q2,3) + (pow(nu2,4) + 2*pow(nu2,2)*pow(nu3,2))*pow(q2,4) + (pow(nu2,4) + 3*pow(nu2,2)*pow(nu3,2) + pow(nu3,4))*pow(q2,2)*pow(q3,2) + 2*nu2*nu3*q2*(pow(nu2,2) + pow(nu3,2))*pow(q3,3) + pow(nu3,2)*(pow(nu2,2) + pow(nu3,2))*pow(q3,4)) - pow(q1,3)*(-2*nu2*nu3*q3*(7*f1*pow(nu2,4) - 8*pow(nu3,2) + 2*pow(nu2,2)*(-4 + 7*f1*pow(nu3,2)))*pow(q2,3) + (8*pow(nu2,4) + 15*pow(nu2,2)*pow(nu3,2))*pow(q2,4) - 2*(pow(nu2,2) + pow(nu3,2))*(-4*pow(nu3,2) + pow(nu2,2)*(-4 + 21*f1*pow(nu3,2)))*pow(q2,2)*pow(q3,2) - 2*nu2*nu3*q2*(pow(nu3,2)*(-8 + 7*f1*pow(nu3,2)) + 2*pow(nu2,2)*(-4 + 7*f1*pow(nu3,2)))*pow(q3,3) + pow(nu3,2)*(pow(nu2,2) + 8*pow(nu3,2))*pow(q3,4)) + 7*pow(q1,2)*(q3*pow(nu2,2)*pow(nu3,2)*(1 + 6*f1*(pow(nu2,2) + pow(nu3,2)))*pow(q2,4) + 2*f1*nu3*pow(nu2,3)*(pow(nu2,2) + 2*pow(nu3,2))*pow(q2,5) + 2*f1*nu2*nu3*(pow(nu2,4) + 4*pow(nu2,2)*pow(nu3,2) + pow(nu3,4))*pow(q2,3)*pow(q3,2) + pow(nu2,2)*(pow(nu2,2) + 3*pow(nu3,2) + 6*f1*pow(nu2,2)*pow(nu3,2) + 6*f1*pow(nu3,4))*pow(q2,2)*pow(q3,3) + 2*nu2*nu3*q2*(pow(nu3,2) + pow(nu2,2)*(1 + 2*f1*pow(nu3,2)) + f1*pow(nu3,4))*pow(q3,4) + pow(nu3,2)*(pow(nu2,2) + pow(nu3,2))*pow(q3,5)) + q1*(2*nu2*nu3*q3*(7*f1*pow(nu2,4) + 4*pow(nu3,2) + 2*pow(nu2,2)*(2 + 7*f1*pow(nu3,2)))*pow(q2,5) + (4*pow(nu2,4) + 11*pow(nu2,2)*pow(nu3,2))*pow(q2,6) + 2*(pow(nu2,4)*(10 + 21*f1*pow(nu3,2)) + pow(nu2,2)*pow(nu3,2)*(19 + 21*f1*pow(nu3,2)) + 2*pow(nu3,4))*pow(q2,4)*pow(q3,2) + 2*nu2*nu3*(2*pow(nu2,2)*(10 + 7*f1*pow(nu3,2)) + pow(nu3,2)*(20 + 7*f1*pow(nu3,2)))*pow(q2,3)*pow(q3,3) + (4*pow(nu2,4) + 31*pow(nu2,2)*pow(nu3,2) + 20*pow(nu3,4))*pow(q2,2)*pow(q3,4) + 8*nu2*nu3*q2*(pow(nu2,2) + pow(nu3,2))*pow(q3,5) + 4*pow(nu3,2)*(pow(nu2,2) + pow(nu3,2))*pow(q3,6)))*pow(-(q3*pow(q1,2)) + pow(q1,3) - q3*pow(q2,2) - q1*(pow(q2,2) + pow(q3,2)),-1)))/14. ) ;
			break ;

		case 2:
			return APfactor * ( (pow(f1,2)*pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(P11(q2,p)*P11(q3,p)*pow(nu2,2)*pow(nu3,2)*pow(q1,2)*(pow(q1,4) + pow(q2,4) + 12*pow(q2,2)*pow(q3,2) - 2*pow(q1,2)*(pow(q2,2) + pow(q3,2)) + pow(q3,4)) + P11(q1,p)*pow(nu1,2)*(P11(q3,p)*pow(nu3,2)*pow(q2,2)*(pow(q1,4) - 2*pow(q1,2)*(pow(q2,2) - 6*pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)) + P11(q2,p)*pow(nu2,2)*pow(q3,2)*(pow(q1,4) + 2*pow(q1,2)*(6*pow(q2,2) - pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)))))/7. ) ;
			break ;

		case 3:
			return APfactor * ( 2*pow(f1,2)*(P11(q2,p)*P11(q3,p)*pow(nu2,2)*pow(nu3,2) + P11(q1,p)*pow(nu1,2)*(P11(q2,p)*pow(nu2,2) + P11(q3,p)*pow(nu3,2))) ) ;
			break ;

		case 4:
			return APfactor * ( (P11(q1,p) + P11(q2,p) + P11(q3,p))*pow(nbar,-1) ) ;
			break ;

		case 5:
			return APfactor * ( (f1*pow(q1,-2)*pow(q2,-2)*pow(q3,-1)*(P11(q1,p)*(q3*P11(q2,p)*((q2 + q3)*pow(q1,2) + q2*(q2*q3 - pow(q2,2) + pow(q3,2)))*(2*nu1*nu2*q2*(11 + 14*f1*(2*pow(nu1,2) + pow(nu2,2)))*pow(q1,5) + (25*pow(nu1,2) + 14*pow(nu2,2))*pow(q1,6) + pow(q1,4)*((28*f1*pow(nu1,4) + 4*pow(nu1,2)*(19 + 28*f1*pow(nu2,2)) + pow(nu2,2)*(53 + 28*f1*pow(nu2,2)))*pow(q2,2) - 3*(12*pow(nu1,2) + 7*pow(nu2,2))*pow(q3,2)) + pow(q2,2)*(pow(q2,2) - pow(q3,2))*((14*pow(nu1,2) + 25*pow(nu2,2))*pow(q2,2) - (7*pow(nu1,2) + 11*pow(nu2,2))*pow(q3,2)) + 2*nu1*nu2*q2*pow(q1,3)*((34 + 42*f1*(pow(nu1,2) + pow(nu2,2)))*pow(q2,2) - (15 + 7*f1*(2*pow(nu1,2) + pow(nu2,2)))*pow(q3,2)) + 2*nu1*nu2*q1*q2*((11 + 14*f1*(pow(nu1,2) + 2*pow(nu2,2)))*pow(q2,4) - (15 + 7*f1*(pow(nu1,2) + 2*pow(nu2,2)))*pow(q2,2)*pow(q3,2) + 4*pow(q3,4)) + pow(q1,2)*((28*f1*pow(nu1,4) + 4*pow(nu2,2)*(19 + 7*f1*pow(nu2,2)) + pow(nu1,2)*(53 + 112*f1*pow(nu2,2)))*pow(q2,4) - (14*f1*pow(nu1,4) + pow(nu2,2)*(57 + 14*f1*pow(nu2,2)) + pow(nu1,2)*(57 + 56*f1*pow(nu2,2)))*pow(q2,2)*pow(q3,2) + (11*pow(nu1,2) + 7*pow(nu2,2))*pow(q3,4))) + q2*P11(q3,p)*(2*pow(q1,2) + 2*pow(q2,2) - pow(q3,2))*((7*q3*(2*pow(nu1,2) + pow(nu3,2)) + q2*(11*pow(nu1,2) + 7*pow(nu3,2)))*pow(q1,6) + 2*nu1*nu3*pow(q1,5)*(q2*q3*(4 + 7*f1*(2*pow(nu1,2) + pow(nu3,2))) + 7*f1*(2*pow(nu1,2) + pow(nu3,2))*pow(q2,2) + 7*pow(q3,2)) + pow(q1,4)*(7*q3*(2*f1*pow(nu1,4) + pow(nu3,2) + pow(nu1,2)*(1 + 8*f1*pow(nu3,2)) + 2*f1*pow(nu3,4))*pow(q2,2) - (15*pow(nu1,2) + 7*pow(nu3,2))*pow(q2,3) + 2*q2*(7*f1*pow(nu1,4) + pow(nu3,2)*(9 + 7*f1*pow(nu3,2)) + pow(nu1,2)*(17 + 28*f1*pow(nu3,2)))*pow(q3,2) + 14*(pow(nu1,2) + pow(nu3,2))*pow(q3,3)) + 2*nu1*nu3*q1*q2*q3*(-7*q3*(1 + f1*(pow(nu1,2) + 2*pow(nu3,2)))*pow(q2,3) + 4*pow(q2,4) + (-8 + 7*f1*(pow(nu1,2) + 2*pow(nu3,2)))*pow(q2,2)*pow(q3,2) + 7*q2*(1 + f1*(pow(nu1,2) + 2*pow(nu3,2)))*pow(q3,3) + 4*pow(q3,4)) + 2*nu1*nu3*pow(q1,3)*(q3*(-8 + 7*f1*(2*pow(nu1,2) + pow(nu3,2)))*pow(q2,3) - 7*f1*(2*pow(nu1,2) + pow(nu3,2))*pow(q2,4) + 21*f1*(pow(nu1,2) + pow(nu3,2))*pow(q2,2)*pow(q3,2) + q2*(20 + 7*f1*(pow(nu1,2) + 2*pow(nu3,2)))*pow(q3,3) + 7*pow(q3,4)) + q2*(pow(q2,2) - pow(q3,2))*(-7*q3*(pow(nu1,2) + pow(nu3,2))*pow(q2,3) + 7*(pow(nu1,2) + pow(nu3,2))*pow(q2,4) - (7*pow(nu1,2) + 3*pow(nu3,2))*pow(q2,2)*pow(q3,2) - 7*q2*pow(nu3,2)*pow(q3,3) - 4*pow(nu3,2)*pow(q3,4)) + pow(q1,2)*(-7*q3*(2*f1*pow(nu1,4) + pow(nu3,2) + pow(nu1,2)*(2 + 8*f1*pow(nu3,2)) + 2*f1*pow(nu3,4))*pow(q2,4) - (3*pow(nu1,2) + 7*pow(nu3,2))*pow(q2,5) + 2*(7*f1*pow(nu1,4) + pow(nu3,2)*(-4 + 7*f1*pow(nu3,2)) + 4*pow(nu1,2)*(-1 + 7*f1*pow(nu3,2)))*pow(q2,3)*pow(q3,2) + 7*(2*f1*pow(nu1,4) + pow(nu1,2)*(3 + 8*f1*pow(nu3,2)) + 2*(pow(nu3,2) + f1*pow(nu3,4)))*pow(q2,2)*pow(q3,3) + q2*(11*pow(nu1,2) + 27*pow(nu3,2))*pow(q3,4) + 7*pow(nu3,2)*pow(q3,5))))*pow(2*pow(q1,2) + 2*pow(q2,2) - pow(q3,2),-1)*pow((q2 + q3)*pow(q1,2) + q2*(q2*q3 - pow(q2,2) + pow(q3,2)),-1) - q1*P11(q2,p)*P11(q3,p)*(-7*q3*(pow(nu2,2) + pow(nu3,2))*pow(q1,6) + 7*(pow(nu2,2) + pow(nu3,2))*pow(q1,7) + 7*q3*pow(q2,2)*(pow(q2,2) + pow(q3,2))*(2*nu2*nu3*q2*q3 + (2*pow(nu2,2) + pow(nu3,2))*pow(q2,2) + pow(nu3,2)*pow(q3,2)) - pow(q1,5)*(-8*nu2*nu3*q2*q3 + (3*pow(nu2,2) + 7*pow(nu3,2))*pow(q2,2) + 2*(7*pow(nu2,2) + 5*pow(nu3,2))*pow(q3,2)) - 7*pow(q1,4)*(q3*(2*f1*pow(nu2,4) + pow(nu3,2) + pow(nu2,2)*(2 + 8*f1*pow(nu3,2)) + 2*f1*pow(nu3,4))*pow(q2,2) + 2*f1*nu2*nu3*(2*pow(nu2,2) + pow(nu3,2))*pow(q2,3) + 2*nu2*nu3*q2*(1 + f1*(pow(nu2,2) + 2*pow(nu3,2)))*pow(q3,2) - pow(nu2,2)*pow(q3,3)) + pow(q1,3)*(2*nu2*nu3*q3*(-8 + 7*f1*(2*pow(nu2,2) + pow(nu3,2)))*pow(q2,3) - (15*pow(nu2,2) + 7*pow(nu3,2))*pow(q2,4) + 2*(7*f1*pow(nu2,4) + pow(nu3,2)*(-4 + 7*f1*pow(nu3,2)) + 4*pow(nu2,2)*(-1 + 7*f1*pow(nu3,2)))*pow(q2,2)*pow(q3,2) + 2*nu2*nu3*q2*(-8 + 7*f1*(pow(nu2,2) + 2*pow(nu3,2)))*pow(q3,3) + (7*pow(nu2,2) - pow(nu3,2))*pow(q3,4)) + 7*pow(q1,2)*(q3*(2*f1*pow(nu2,4) + pow(nu3,2) + pow(nu2,2)*(1 + 8*f1*pow(nu3,2)) + 2*f1*pow(nu3,4))*pow(q2,4) + 2*f1*nu2*nu3*(2*pow(nu2,2) + pow(nu3,2))*pow(q2,5) + 6*f1*nu2*nu3*(pow(nu2,2) + pow(nu3,2))*pow(q2,3)*pow(q3,2) + (2*f1*pow(nu2,4) + pow(nu2,2)*(3 + 8*f1*pow(nu3,2)) + 2*(pow(nu3,2) + f1*pow(nu3,4)))*pow(q2,2)*pow(q3,3) + 2*nu2*nu3*q2*(1 + f1*(pow(nu2,2) + 2*pow(nu3,2)))*pow(q3,4) + pow(nu3,2)*pow(q3,5)) + q1*(2*nu2*nu3*q3*(4 + 7*f1*(2*pow(nu2,2) + pow(nu3,2)))*pow(q2,5) + (11*pow(nu2,2) + 7*pow(nu3,2))*pow(q2,6) + 2*(7*f1*pow(nu2,4) + pow(nu3,2)*(9 + 7*f1*pow(nu3,2)) + pow(nu2,2)*(17 + 28*f1*pow(nu3,2)))*pow(q2,4)*pow(q3,2) + 2*nu2*nu3*(20 + 7*f1*(pow(nu2,2) + 2*pow(nu3,2)))*pow(q2,3)*pow(q3,3) + (11*pow(nu2,2) + 27*pow(nu3,2))*pow(q2,2)*pow(q3,4) + 8*nu2*nu3*q2*pow(q3,5) + 4*pow(nu3,2)*pow(q3,6)))*pow(-(q3*pow(q1,2)) + pow(q1,3) - q3*pow(q2,2) - q1*(pow(q2,2) + pow(q3,2)),-1)))/14. ) ;
			break ;

		case 6:
			return APfactor * ( (f1*pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(P11(q2,p)*P11(q3,p)*(pow(nu2,2) + pow(nu3,2))*pow(q1,2)*(pow(q1,4) + pow(q2,4) + 12*pow(q2,2)*pow(q3,2) - 2*pow(q1,2)*(pow(q2,2) + pow(q3,2)) + pow(q3,4)) + P11(q1,p)*(P11(q3,p)*(pow(nu1,2) + pow(nu3,2))*pow(q2,2)*(pow(q1,4) - 2*pow(q1,2)*(pow(q2,2) - 6*pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)) + P11(q2,p)*(pow(nu1,2) + pow(nu2,2))*pow(q3,2)*(pow(q1,4) + 2*pow(q1,2)*(6*pow(q2,2) - pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)))))/7. ) ;
			break ;

		case 7:
			return APfactor * ( 2*f1*(P11(q2,p)*P11(q3,p)*(pow(nu2,2) + pow(nu3,2)) + P11(q1,p)*(P11(q2,p)*(pow(nu1,2) + pow(nu2,2)) + P11(q3,p)*(pow(nu1,2) + pow(nu3,2)))) ) ;
			break ;

		case 8:
			return APfactor * ( (pow(q1,-2)*pow(q2,-2)*pow(q3,-1)*(q1*P11(q2,p)*P11(q3,p)*(-pow(q1,4) + pow(q1,2)*pow(q3,2) + pow(q2,2)*(pow(q2,2) + pow(q3,2)) + 2*f1*q1*q2*(q2*q3*(pow(nu2,2) + pow(nu3,2)) + nu2*nu3*pow(q2,2) + nu2*nu3*pow(q3,2))) + P11(q1,p)*(q3*P11(q2,p)*(2*f1*nu1*nu2*q2*pow(q1,3) + pow(q1,4) + 2*f1*nu1*nu2*q1*pow(q2,3) + pow(q2,4) + pow(q1,2)*(2*(1 + f1*(pow(nu1,2) + pow(nu2,2)))*pow(q2,2) - pow(q3,2)) - pow(q2,2)*pow(q3,2)) + q2*P11(q3,p)*(q3*(q3 + 2*f1*q2*(pow(nu1,2) + pow(nu3,2)))*pow(q1,2) + 2*f1*nu1*nu3*q2*pow(q1,3) + pow(q1,4) - pow(q2,4) + 2*f1*nu1*nu3*q1*q2*pow(q3,2) + pow(q2,2)*pow(q3,2)))))/2. ) ;
			break ;

		case 9:
			return APfactor * ( (pow(q1,-2)*pow(q2,-2)*pow(q3,-2)*(P11(q2,p)*P11(q3,p)*pow(q1,2)*(pow(q1,4) + pow(q2,4) + 12*pow(q2,2)*pow(q3,2) - 2*pow(q1,2)*(pow(q2,2) + pow(q3,2)) + pow(q3,4)) + P11(q1,p)*(P11(q3,p)*pow(q2,2)*(pow(q1,4) - 2*pow(q1,2)*(pow(q2,2) - 6*pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)) + P11(q2,p)*pow(q3,2)*(pow(q1,4) + 2*pow(q1,2)*(6*pow(q2,2) - pow(q3,2)) + pow(pow(q2,2) - pow(q3,2),2)))))/7. ) ;
			break ;

		case 10:
			return APfactor * ( 2*(P11(q2,p)*P11(q3,p) + P11(q1,p)*(P11(q2,p) + P11(q3,p))) ) ;
			break ;
	}
}

void ScoccimarroTransform (const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, 
	double & q1, double & q2, double & q3, double & nu1, double & nu2, double & nu3) {

	double cos12 = (k1*k1 + k2*k2 - k3*k3) / (2.*k1*k2) ;

	nu1 = mu ;
	nu2 = mu * cos12 - sqrt(1-mu*mu)*sqrt(1.-cos12*cos12)*cos(phi) ;
	nu3 = - k1/k3 * nu1 - k2/k3 * nu2 ;

	q1 = k1 ;
	q2 = k2 ;
	q3 = k3 ;
}


// If the condition of closed triangle is broken by the AP effect, return 0 and the integrands are set to 0.
void ScoccimarroTransformWithAP (const double & aperp, const double & apar, 
	const double & k1, const double & k2, const double & k3, const double & mu, const double & phi, 
	double & q1, double & q2, double & q3, double & nu1, double & nu2, double & nu3) {

	double F = apar/aperp ;

	double cos12 = (k1*k1 + k2*k2 - k3*k3) / (2.*k1*k2) ;

	nu1 = mu ;
	nu2 = mu * cos12 - sqrt(1-mu*mu)*sqrt(1.-cos12*cos12)*cos(phi) ;
	nu3 = - k1/k3 * nu1 - k2/k3 * nu2 ;

	q1 = k1/aperp * sqrt(1.+ nu1*nu1 * (pow(F,-2)-1.)) ;
	q2 = k2/aperp * sqrt(1.+ nu2*nu2 * (pow(F,-2)-1.)) ;
	q3 = k3/aperp * sqrt(1.+ nu3*nu3 * (pow(F,-2)-1.)) ;

	nu1 = nu1/F * 1./sqrt(1.+ nu1*nu1 * (pow(F,-2)-1.)) ;
	nu2 = nu2/F * 1./sqrt(1.+ nu2*nu2 * (pow(F,-2)-1.)) ;
	nu3 = nu3/F * 1./sqrt(1.+ nu3*nu3 * (pow(F,-2)-1.)) ;
	
}