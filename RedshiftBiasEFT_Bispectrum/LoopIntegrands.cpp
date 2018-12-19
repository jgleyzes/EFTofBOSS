#include "RedshiftBiasEFT.h"

double Integrands (const int j, const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ; 

	switch(j) {
		case 0:
			return (f1*Pi*P11(q,params)*pow(k,2)*(2*P11(k,params)*(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2))*(-28*pow(k,4)*pow(x,2) + 4*pow(q,4)*(6 - 25*pow(x,2) + 12*pow(x,4)) + 8*pow(k,2)*pow(q,2)*(3 - 10*pow(x,2) + 14*pow(x,4)) + 7*f1*(-1 + pow(x,2))*(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2))) - 7*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(-4*(q + k*x)*(2*k*q*x + pow(k,2) + 2*pow(q,2))*(k*x + q*(-1 + 2*pow(x,2))) + f1*pow(k,2)*(-1 + pow(x,2))*pow(k + 2*q*x,2))*pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(-4*(-q + k*x)*(-2*k*q*x + pow(k,2) + 2*pow(q,2))*(q + k*x - 2*q*pow(x,2)) + f1*pow(k,2)*(-1 + pow(x,2))*pow(k - 2*q*x,2))*pow(2*k*q*x + pow(k,2) + pow(q,2),2)))*pow(-4*pow(k,2)*pow(q,3)*pow(x,2) + q*pow(pow(k,2) + pow(q,2),2),-2))/14. ;
			break ; 
		case 1:
			return (4*f1*Pi*P11(q,params)*pow(k,2)*pow(q,-1)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(q + k*x - 2*q*pow(x,2))*(-14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2) - Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(k*x + q*(-1 + 2*pow(x,2)))*(14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*pow(2*k*q*x + pow(k,2) + pow(q,2),-2)))/7. ;
			break ; 
		case 2:
			return 4*f1*Pi*P11(q,params)*pow(k,2)*pow(q,-1)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(q + k*x - 2*q*pow(x,2))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-1) + Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(q - k*x - 2*q*pow(x,2))*pow(2*k*q*x + pow(k,2) + pow(q,2),-1)) ;
			break ; 
		case 3:
			return (Pi*P11(q,params)*pow(f1,2)*pow(k,2)*(pow(k,2)*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(2*pow(q,2) + (k + 2*q*x)*(6*q*x*(-1 + pow(x,2)) + k*(-1 + 3*pow(x,2))))*pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(2*pow(q,2) + (k - 2*q*x)*(-k + 6*q*x + 3*k*pow(x,2) - 6*q*pow(x,3)))*pow(2*k*q*x + pow(k,2) + pow(q,2),2)) - 2*P11(k,params)*(-1 + 3*pow(x,2))*pow(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2),2))*pow(-4*pow(k,2)*pow(q,3)*pow(x,2) + q*pow(pow(k,2) + pow(q,2),2),-2))/2. ;
			break ; 
		case 4:
			return (Pi*P11(q,params)*(-2*P11(k,params)*(-2*k*q*x + pow(k,2) + pow(q,2))*(2*k*q*x + pow(k,2) + pow(q,2))*(21*pow(k,6)*pow(x,2) + 8*pow(q,6)*(-1 + 3*pow(x,2)) + pow(k,2)*pow(q,4)*(-34 + 155*pow(x,2) - 132*pow(x,4)) - 2*pow(k,4)*pow(q,2)*(13 - 42*pow(x,2) + 42*pow(x,4))) + 21*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*pow(q - k*x,2)*pow(2*k*q*x + pow(k,2) + pow(q,2),2)*pow(-2*k*q*x + pow(k,2) + 2*pow(q,2),2) + Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*pow(q + k*x,2)*pow(-2*k*q*x + pow(k,2) + pow(q,2),2)*pow(2*k*q*x + pow(k,2) + 2*pow(q,2),2)))*pow(-4*pow(k,2)*pow(q,3)*pow(x,2) + q*pow(pow(k,2) + pow(q,2),2),-2))/21. ;
			break ; 
		case 5:
			return (4*Pi*P11(q,params)*pow(q,-1)*(-((q - k*x)*Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(-2*k*q*x + pow(k,2) + 2*pow(q,2))*(-14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)) - (q + k*x)*Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(2*k*q*x + pow(k,2) + 2*pow(q,2))*(14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*pow(2*k*q*x + pow(k,2) + pow(q,2),-2)))/7. ;
			break ; 
		case 6:
			return (-16*Pi*P11(k,params)*P11(q,params)*(2*pow(q,4)*(1 - 3*pow(x,2)) + pow(k,4)*(3 - 8*pow(x,2) + pow(x,4)) + pow(k,2)*pow(q,2)*(5 - 22*pow(x,2) + 25*pow(x,4)))*pow(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2),-1))/21. ;
			break ; 
		case 7:
			return 4*Pi*P11(q,params)*pow(q,-1)*(-((q - k*x)*Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(-2*k*q*x + pow(k,2) + 2*pow(q,2))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-1)) - (q + k*x)*Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(2*k*q*x + pow(k,2) + 2*pow(q,2))*pow(2*k*q*x + pow(k,2) + pow(q,2),-1)) ;
			break ; 
		case 8:
			return (4*Pi*P11(q,params)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*pow(-14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)),2) + Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*pow(14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)),2)))/49. ;
			break ; 
		case 9:
			return (8*Pi*P11(q,params)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(-14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-1) + Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*pow(2*k*q*x + pow(k,2) + pow(q,2),-1)))/7. ;
			break ; 
		case 10:
			return 4*Pi*P11(q,params)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params) + Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)) ;
			break ; 
		case 11:
			return (f1*Pi*P11(q,params)*(3*pow(k,2)*((q + k*x)*Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(2*k*q*x + pow(k,2) + 2*pow(q,2))*((2 + 7*f1)*q + 14*k*x + (12 - 7*f1)*q*pow(x,2))*pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + (-q + k*x)*Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(-2*k*q*x + pow(k,2) + 2*pow(q,2))*(14*k*x + q*(-2 - 7*f1 + (-12 + 7*f1)*pow(x,2)))*pow(2*k*q*x + pow(k,2) + pow(q,2),2)) - 4*P11(k,params)*(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2))*(21*pow(k,6)*pow(x,2) + 4*pow(q,6)*(-1 + 3*pow(x,2)) + pow(k,4)*pow(q,2)*(-10 + 9*f1 + 6*(8 - 3*f1)*pow(x,2) + 9*(-8 + f1)*pow(x,4)) + pow(k,2)*pow(q,4)*(-14 + 91*pow(x,2) - 72*pow(x,4) + 9*f1*pow(-1 + pow(x,2),2))))*pow(-4*pow(k,2)*pow(q,3)*pow(x,2) + q*pow(pow(k,2) + pow(q,2),2),-2))/21. ;
			break ; 
		case 12:
			return (Pi*P11(q,params)*pow(f1,2)*pow(k,2)*(2*P11(k,params)*(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2))*(-28*pow(k,4)*pow(x,2) + 2*pow(q,4)*(9 - 44*pow(x,2) + 21*pow(x,4)) + 2*pow(k,2)*pow(q,2)*(9 - 34*pow(x,2) + 53*pow(x,4)) + 7*f1*(-1 + pow(x,2))*(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2))) - Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(q*x*pow(k,3)*(5 - 89*pow(x,2)) + 14*pow(q,4)*(1 - 3*pow(x,2)) + 7*f1*pow(k,2)*((k + 2*q*x)*(k + 3*q*x) - pow(q,2))*(-1 + pow(x,2)) - 28*pow(k,4)*pow(x,2) - 84*k*pow(q,3)*pow(x,3) + pow(k,2)*pow(q,2)*(9 - 55*pow(x,2) - 66*pow(x,4)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),2) - Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(14*pow(q,4)*(1 - 3*pow(x,2)) - 28*pow(k,4)*pow(x,2) + q*x*pow(k,3)*(-5 + 89*pow(x,2)) + 7*f1*pow(k,2)*(-1 + pow(x,2))*(-5*k*q*x + pow(k,2) + pow(q,2)*(-1 + 6*pow(x,2))) + 84*k*pow(q,3)*pow(x,3) + pow(k,2)*pow(q,2)*(9 - 55*pow(x,2) - 66*pow(x,4)))*pow(2*k*q*x + pow(k,2) + pow(q,2),2))*pow(-4*pow(k,2)*pow(q,3)*pow(x,2) + q*pow(pow(k,2) + pow(q,2),2),-2))/7. ;
			break ; 
		case 13:
			return Pi*P11(q,params)*pow(f1,3)*pow(k,2)*(pow(k,2)*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(pow(k,2)*(-1 + 3*pow(x,2)) + k*q*x*(-7 + 11*pow(x,2)) + pow(q,2)*(1 - 9*pow(x,2) + 10*pow(x,4)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(k*q*x*(7 - 11*pow(x,2)) + pow(k,2)*(-1 + 3*pow(x,2)) + pow(q,2)*(1 - 9*pow(x,2) + 10*pow(x,4)))*pow(2*k*q*x + pow(k,2) + pow(q,2),2)) - 2*P11(k,params)*(-1 + 3*pow(x,2))*pow(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2),2))*pow(-4*pow(k,2)*pow(q,3)*pow(x,2) + q*pow(pow(k,2) + pow(q,2),2),-2) ;
			break ; 
		case 14:
			return (2*f1*Pi*P11(q,params)*pow(k,2)*pow(q,-1)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(-14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*(14*k*x + q*(-2 - 7*f1 + (-12 + 7*f1)*pow(x,2)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2) - Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*((2 + 7*f1)*q + 14*k*x + (12 - 7*f1)*q*pow(x,2))*(14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*pow(2*k*q*x + pow(k,2) + pow(q,2),-2)))/49. ;
			break ; 
		case 15:
			return (2*Pi*P11(q,params)*pow(f1,2)*pow(k,2)*pow(q,-1)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(q + 2*k*x - 3*q*pow(x,2))*(-14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2) - Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*(2*k*x + q*(-1 + 3*pow(x,2)))*pow(2*k*q*x + pow(k,2) + pow(q,2),-2)))/7. ;
			break ; 
		case 16:
			return (-16*f1*Pi*P11(k,params)*P11(q,params)*(2*pow(q,4)*(1 - 3*pow(x,2)) + pow(k,4)*(3 - 8*pow(x,2) + pow(x,4)) + pow(k,2)*pow(q,2)*(5 - 22*pow(x,2) + 25*pow(x,4)))*pow(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2),-1))/21. ;
			break ; 
		case 17:
			return (2*f1*Pi*P11(q,params)*pow(k,2)*pow(q,-1)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(14*k*x + q*(-2 - 7*f1 + (-12 + 7*f1)*pow(x,2)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-1) - Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*((2 + 7*f1)*q + 14*k*x + (12 - 7*f1)*q*pow(x,2))*pow(2*k*q*x + pow(k,2) + pow(q,2),-1)))/7. ;
			break ; 
		case 18:
			return 2*Pi*P11(q,params)*pow(f1,2)*pow(k,2)*pow(q,-1)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(q + 2*k*x - 3*q*pow(x,2))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-1) + Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(q - 2*k*x - 3*q*pow(x,2))*pow(2*k*q*x + pow(k,2) + pow(q,2),-1)) ;
			break ; 
		case 19:
			return (Pi*P11(q,params)*pow(f1,2)*pow(k,2)*pow(q,-2)*(pow(k,2)*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*(392*pow(k,2)*pow(x,2) - 56*k*q*x*(-2 - 7*f1 + (-12 + 7*f1)*pow(x,2)) + pow(q,2)*(56*f1*(1 + 5*pow(x,2) - 6*pow(x,4)) + 147*pow(f1,2)*pow(-1 + pow(x,2),2) + 8*pow(1 + 6*pow(x,2),2))) + Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*(392*pow(k,2)*pow(x,2) + 56*k*q*x*(-2 - 7*f1 + (-12 + 7*f1)*pow(x,2)) + pow(q,2)*(56*f1*(1 + 5*pow(x,2) - 6*pow(x,4)) + 147*pow(f1,2)*pow(-1 + pow(x,2),2) + 8*pow(1 + 6*pow(x,2),2)))) - 112*P11(k,params)*(7*pow(k,4)*pow(x,2) + 2*pow(k,2)*pow(q,2)*(1 + 2*pow(x,2) - 10*pow(x,4) + 3*f1*pow(-1 + pow(x,2),2)) + pow(q,4)*(2 + 9*pow(x,2) - 4*pow(x,4) + 6*f1*pow(-1 + pow(x,2),2)))*pow(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2),-1)))/392. ;
			break ; 
		case 20:
			return (Pi*P11(q,params)*pow(f1,3)*pow(k,2)*(4*P11(k,params)*(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2))*(4*(-7*pow(k,4)*pow(x,2) + pow(q,4)*(3 - 19*pow(x,2) + 9*pow(x,4)) + pow(k,2)*pow(q,2)*(3 - 14*pow(x,2) + 25*pow(x,4))) + 7*f1*(-1 + pow(x,2))*(-4*pow(k,2)*pow(q,2)*pow(x,2) + pow(pow(k,2) + pow(q,2),2))) - pow(k,2)*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(-4*(q + 7*k*x + 6*q*pow(x,2))*(2*k*x + q*(-1 + 3*pow(x,2))) + 7*f1*(-1 + pow(x,2))*(12*k*q*x + 2*pow(k,2) + 3*pow(q,2)*(-1 + 5*pow(x,2))))*pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(-4*(q - 7*k*x + 6*q*pow(x,2))*(-2*k*x + q*(-1 + 3*pow(x,2))) + 7*f1*(-1 + pow(x,2))*(-12*k*q*x + 2*pow(k,2) + 3*pow(q,2)*(-1 + 5*pow(x,2))))*pow(2*k*q*x + pow(k,2) + pow(q,2),2)))*pow(-4*pow(k,2)*pow(q,3)*pow(x,2) + q*pow(pow(k,2) + pow(q,2),2),-2))/28. ;
			break ; 
		case 21:
			return (Pi*P11(q,params)*pow(f1,4)*pow(k,2)*pow(q,-2)*(P11(k,params)*(8 - 24*pow(x,2)) + pow(k,2)*(Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(8*k*q*x*(3 - 5*pow(x,2)) + 4*pow(k,2)*(-1 + 3*pow(x,2)) + pow(q,2)*(3 - 30*pow(x,2) + 35*pow(x,4)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2) + Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*(4*pow(k,2)*(-1 + 3*pow(x,2)) + 8*k*q*x*(-3 + 5*pow(x,2)) + pow(q,2)*(3 - 30*pow(x,2) + 35*pow(x,4)))*pow(2*k*q*x + pow(k,2) + pow(q,2),-2))))/8. ;
			break ; 
	}
}
