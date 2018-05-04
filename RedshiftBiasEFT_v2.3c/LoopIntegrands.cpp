#include "RedshiftBiasEFT.h"

double I1 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (f1*Pi*P11(q,params)*pow(k,2)*pow(q,-2)*
     pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
     (2*P11(k,params)*(pow(k,4) + pow(q,4) + 
          2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)))*
        (7*f1*(pow(k,4) + pow(q,4) + 
             2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)))*(-1 + pow(x,2)) - 
          4*(7*pow(k,4)*pow(x,2) + 
             pow(q,4)*(-6 + 25*pow(x,2) - 12*pow(x,4)) - 
             2*pow(k,2)*pow(q,2)*(3 - 10*pow(x,2) + 14*pow(x,4)))) - 
       7*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (-4*(pow(k,4)*pow(x,2) + 2*pow(q,4)*(-1 + 2*pow(x,2)) + 
                2*k*x*pow(q,3)*(-1 + 4*pow(x,2)) + 
                4*q*pow(k,3)*pow(x,3) + 
                pow(k,2)*pow(q,2)*(-1 + 4*pow(x,2) + 4*pow(x,4))) + 
             f1*pow(k,2)*(-1 + pow(x,2))*pow(k + 2*q*x,2))*
           pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + 
          Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (-4*(2*k*x*pow(q,3)*(1 - 4*pow(x,2)) + pow(k,4)*pow(x,2) + 
                2*pow(q,4)*(-1 + 2*pow(x,2)) - 4*q*pow(k,3)*pow(x,3) + 
                pow(k,2)*pow(q,2)*(-1 + 4*pow(x,2) + 4*pow(x,4))) + 
             f1*pow(k,2)*(-1 + pow(x,2))*pow(k - 2*q*x,2))*
           pow(2*k*q*x + pow(k,2) + pow(q,2),2))))/14. ;
}

double I2 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (4*f1*Pi*P11(q,params)*pow(k,2)*pow(q,-1)*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
     (-(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
          (7*pow(q,3)*(-1 + 2*pow(x,2)) + x*pow(k,3)*(5 + 2*pow(x,2)) + 7*k*x*pow(q,2)*(-1 + 4*pow(x,2)) + 
            q*pow(k,2)*(-5 + 22*pow(x,2) + 4*pow(x,4)))*pow(-2*k*q*x + pow(k,2) + pow(q,2),2)) + 
       Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (7*pow(q,3)*(1 - 2*pow(x,2)) + x*pow(k,3)*(5 + 2*pow(x,2)) + 7*k*x*pow(q,2)*(-1 + 4*pow(x,2)) + 
          q*pow(k,2)*(5 - 22*pow(x,2) - 4*pow(x,4)))*pow(2*k*q*x + pow(k,2) + pow(q,2),2)))/7. ;
}

double I3 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return  4*f1*Pi*P11(q,params)*pow(k,2)*pow(q,-1)*
   (Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
      P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
      (q*pow(k,2) + x*pow(k,3) + k*x*pow(q,2)*(3 - 4*pow(x,2)) + 
        pow(q,3)*(1 - 2*pow(x,2))) + 
     Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
      P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
      (q*pow(k,2) - x*pow(k,3) + pow(q,3)*(1 - 2*pow(x,2)) + 
        k*x*pow(q,2)*(-3 + 4*pow(x,2))))*
   pow(-2*k*q*x + pow(k,2) + pow(q,2),-1)*
   pow(2*k*q*x + pow(k,2) + pow(q,2),-1) ;
}

double I4 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (Pi*P11(q,params)*pow(f1,2)*pow(k,2)*pow(q,-2)*
     pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
     (pow(k,2)*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (4*k*q*x*(-2 + 3*pow(x,2)) + pow(k,2)*(-1 + 3*pow(x,2)) + 
             2*pow(q,2)*(1 - 6*pow(x,2) + 6*pow(x,4)))*
           pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + 
          Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (4*k*q*x*(2 - 3*pow(x,2)) + pow(k,2)*(-1 + 3*pow(x,2)) + 
             2*pow(q,2)*(1 - 6*pow(x,2) + 6*pow(x,4)))*
           pow(2*k*q*x + pow(k,2) + pow(q,2),2)) - 
       2*P11(k,params)*(-1 + 3*pow(x,2))*
        pow(pow(k,4) + pow(q,4) + 2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)),
         2)))/2. ;
}

double I5 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (Pi*P11(q,params)*pow(q,-2)*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
     (-2*P11(k,params)*(21*pow(k,10)*pow(x,2) + 
          8*pow(q,10)*(-1 + 3*pow(x,2)) + 
          pow(k,2)*pow(q,8)*(-50 + 235*pow(x,2) - 228*pow(x,4)) - 
          2*pow(k,8)*pow(q,2)*(13 - 63*pow(x,2) + 84*pow(x,4)) + 
          2*pow(k,6)*pow(q,4)*
           (-43 + 224*pow(x,2) - 318*pow(x,4) + 168*pow(x,6)) + 
          2*pow(k,4)*pow(q,6)*
           (-51 + 277*pow(x,2) - 484*pow(x,4) + 264*pow(x,6))) + 
       21*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           pow(q*pow(k,4) + x*pow(k,5) + 2*pow(q,5) + 
             x*pow(k,3)*pow(q,2)*(3 - 4*pow(x,2)) + 
             3*pow(k,2)*pow(q,3)*(1 - 2*pow(x,2)),2) + 
          Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           pow(q*pow(k,4) - x*pow(k,5) + 2*pow(q,5) + 
             3*pow(k,2)*pow(q,3)*(1 - 2*pow(x,2)) + 
             x*pow(k,3)*pow(q,2)*(-3 + 4*pow(x,2)),2))))/21. ;
}

double I6 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (-4*Pi*P11(q,params)*pow(q,-1)*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
     (Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (56*k*x*pow(q,4) + 14*pow(q,5) + x*pow(k,5)*(5 + 2*pow(x,2)) + 
          x*pow(k,3)*pow(q,2)*(41 + 36*pow(x,2)) + 
          pow(k,2)*pow(q,3)*(17 + 74*pow(x,2)) + 
          q*pow(k,4)*(5 + 26*pow(x,2) + 4*pow(x,4)))*
        pow(-2*k*q*x + pow(k,2) + pow(q,2),2) - 
       Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (56*k*x*pow(q,4) - 14*pow(q,5) + x*pow(k,5)*(5 + 2*pow(x,2)) + 
          x*pow(k,3)*pow(q,2)*(41 + 36*pow(x,2)) - 
          pow(k,2)*pow(q,3)*(17 + 74*pow(x,2)) - 
          q*pow(k,4)*(5 + 26*pow(x,2) + 4*pow(x,4)))*
        pow(2*k*q*x + pow(k,2) + pow(q,2),2)))/7. ;
}

double I7 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (-16*Pi*P11(k,params)*P11(q,params)*
     (2*pow(q,4)*(1 - 3*pow(x,2)) + 
       pow(k,4)*(3 - 8*pow(x,2) + pow(x,4)) + 
       pow(k,2)*pow(q,2)*(5 - 22*pow(x,2) + 25*pow(x,4)))*
     pow(-2*k*q*x + pow(k,2) + pow(q,2),-1)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-1))/21. ;
}

double I8 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return -4*Pi*P11(q,params)*pow(q,-1)*
   (Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
      P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
      (q*pow(k,4) + x*pow(k,5) + 2*pow(q,5) + 
        x*pow(k,3)*pow(q,2)*(3 - 4*pow(x,2)) + 
        3*pow(k,2)*pow(q,3)*(1 - 2*pow(x,2))) + 
     Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
      P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
      (q*pow(k,4) - x*pow(k,5) + 2*pow(q,5) + 
        3*pow(k,2)*pow(q,3)*(1 - 2*pow(x,2)) + 
        x*pow(k,3)*pow(q,2)*(-3 + 4*pow(x,2))))*
   pow(-2*k*q*x + pow(k,2) + pow(q,2),-1)*
   pow(2*k*q*x + pow(k,2) + pow(q,2),-1) ;
}

double I9 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I10 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (4*Pi*P11(q,params)*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
     (Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        pow(7*pow(q,4) + 2*pow(k,2)*pow(q,2)*(6 - 13*pow(x,2)) - 
          4*q*x*pow(k,3)*(-1 + pow(x,2)) + pow(k,4)*(5 + 2*pow(x,2)),2)\
        + Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        pow(7*pow(q,4) + 2*pow(k,2)*pow(q,2)*(6 - 13*pow(x,2)) + 
          4*q*x*pow(k,3)*(-1 + pow(x,2)) + pow(k,4)*(5 + 2*pow(x,2)),2))
     )/49. ;
}

double I11 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I12 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (8*Pi*P11(q,params)*(Heaviside(-q + 
          pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (7*pow(q,4) + 2*pow(k,2)*pow(q,2)*(6 - 13*pow(x,2)) - 
          4*q*x*pow(k,3)*(-1 + pow(x,2)) + pow(k,4)*(5 + 2*pow(x,2))) + 
       Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (7*pow(q,4) + 2*pow(k,2)*pow(q,2)*(6 - 13*pow(x,2)) + 
          4*q*x*pow(k,3)*(-1 + pow(x,2)) + pow(k,4)*(5 + 2*pow(x,2))))*
     pow(-2*k*q*x + pow(k,2) + pow(q,2),-1)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-1))/7. ;
}

double I13 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I14 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I15 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I16 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I17 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 4*Pi*P11(q,params)*(Heaviside(-q + 
        pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
      P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params) + 
     Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
      P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)) ;
}

double I18 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I19 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I20 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (f1*Pi*P11(q,params)*pow(q,-2)*pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
     (3*pow(k,2)*(Heaviside(-q + 
             pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (14*k*x + q*(2 - 7*f1*(-1 + pow(x,2)) + 12*pow(x,2)))*
           (x*pow(k,3) + 4*k*x*pow(q,2) + 2*pow(q,3) + 
             pow(k,2)*(q + 2*q*pow(x,2)))*
           pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + 
          Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (x*pow(k,3) + 4*k*x*pow(q,2) - 2*pow(q,3) - 
             pow(k,2)*(q + 2*q*pow(x,2)))*
           (14*k*x + q*(7*f1*(-1 + pow(x,2)) - 2*(1 + 6*pow(x,2))))*
           pow(2*k*q*x + pow(k,2) + pow(q,2),2)) - 
       4*P11(k,params)*(pow(k,4) + pow(q,4) + 
          2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)))*
        (21*pow(k,6)*pow(x,2) + 4*pow(q,6)*(-1 + 3*pow(x,2)) + 
          pow(k,4)*pow(q,2)*(-10 + 48*pow(x,2) - 72*pow(x,4) + 
             9*f1*pow(-1 + pow(x,2),2)) + 
          pow(k,2)*pow(q,4)*(-14 + 91*pow(x,2) - 72*pow(x,4) + 
             9*f1*pow(-1 + pow(x,2),2)))))/21. ;
}

double I21 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (Pi*P11(q,params)*pow(f1,2)*pow(k,2)*pow(q,-2)*
     pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
     (2*P11(k,params)*(pow(k,4) + pow(q,4) + 
          2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)))*
        (7*f1*(pow(k,4) + pow(q,4) + 
             2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)))*(-1 + pow(x,2)) + 
          2*(-14*pow(k,4)*pow(x,2) + 
             pow(q,4)*(9 - 44*pow(x,2) + 21*pow(x,4)) + 
             pow(k,2)*pow(q,2)*(9 - 34*pow(x,2) + 53*pow(x,4)))) - 
       Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (q*x*pow(k,3)*(5 - 89*pow(x,2)) + 
          14*pow(q,4)*(1 - 3*pow(x,2)) - 28*pow(k,4)*pow(x,2) + 
          7*f1*pow(k,2)*(-1 + pow(x,2))*
           (5*k*q*x + pow(k,2) + pow(q,2)*(-1 + 6*pow(x,2))) - 
          84*k*pow(q,3)*pow(x,3) + 
          pow(k,2)*pow(q,2)*(9 - 55*pow(x,2) - 66*pow(x,4)))*
        pow(-2*k*q*x + pow(k,2) + pow(q,2),2) - 
       Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (14*pow(q,4)*(1 - 3*pow(x,2)) - 28*pow(k,4)*pow(x,2) + 
          q*x*pow(k,3)*(-5 + 89*pow(x,2)) + 
          7*f1*pow(k,2)*(-1 + pow(x,2))*
           (-5*k*q*x + pow(k,2) + pow(q,2)*(-1 + 6*pow(x,2))) + 
          84*k*pow(q,3)*pow(x,3) + 
          pow(k,2)*pow(q,2)*(9 - 55*pow(x,2) - 66*pow(x,4)))*
        pow(2*k*q*x + pow(k,2) + pow(q,2),2)))/7. ;
}

double I22 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return Pi*P11(q,params)*pow(f1,3)*pow(k,2)*pow(q,-2)*
   pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
   pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
   (pow(k,2)*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
         P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
         (pow(k,2)*(-1 + 3*pow(x,2)) + k*q*x*(-7 + 11*pow(x,2)) + 
           pow(q,2)*(1 - 9*pow(x,2) + 10*pow(x,4)))*
         pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + 
        Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
         P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
         (k*q*x*(7 - 11*pow(x,2)) + pow(k,2)*(-1 + 3*pow(x,2)) + 
           pow(q,2)*(1 - 9*pow(x,2) + 10*pow(x,4)))*
         pow(2*k*q*x + pow(k,2) + pow(q,2),2)) - 
     2*P11(k,params)*(-1 + 3*pow(x,2))*
      pow(pow(k,4) + pow(q,4) + 2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)),2)) ;
}

double I23 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I24 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (2*f1*Pi*P11(q,params)*pow(k,2)*pow(q,-1)*
     (Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (-14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*
        (14*k*x + q*(7*f1*(-1 + pow(x,2)) - 2*(1 + 6*pow(x,2))))*
        pow(-2*k*q*x + pow(k,2) + pow(q,2),-2) - 
       Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (14*k*q*x + 7*pow(q,2) + pow(k,2)*(5 + 2*pow(x,2)))*
        (14*k*x + q*(2 - 7*f1*(-1 + pow(x,2)) + 12*pow(x,2)))*
        pow(2*k*q*x + pow(k,2) + pow(q,2),-2)))/49. ;
}

double I25 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (2*Pi*P11(q,params)*pow(f1,2)*pow(k,2)*pow(q,-1)*
     (Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (7*pow(q,3)*(1 - 3*pow(x,2)) + 2*x*pow(k,3)*(5 + 2*pow(x,2)) + 
          42*k*pow(q,2)*pow(x,3) + 
          q*pow(k,2)*(5 - 41*pow(x,2) - 6*pow(x,4)))*
        pow(-2*k*q*x + pow(k,2) + pow(q,2),-2) - 
       Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (2*x*pow(k,3)*(5 + 2*pow(x,2)) + 7*pow(q,3)*(-1 + 3*pow(x,2)) + 
          42*k*pow(q,2)*pow(x,3) + 
          q*pow(k,2)*(-5 + 41*pow(x,2) + 6*pow(x,4)))*
        pow(2*k*q*x + pow(k,2) + pow(q,2),-2)))/7. ;
}

double I26 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I27 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I28 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (-16*f1*Pi*P11(k,params)*P11(q,params)*
     (2*pow(q,4)*(1 - 3*pow(x,2)) + 
       pow(k,4)*(3 - 8*pow(x,2) + pow(x,4)) + 
       pow(k,2)*pow(q,2)*(5 - 22*pow(x,2) + 25*pow(x,4)))*
     pow(-2*k*q*x + pow(k,2) + pow(q,2),-1)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-1))/21. ;
}

double I29 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I30 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I31 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I32 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (2*f1*Pi*P11(q,params)*pow(k,2)*pow(q,-1)*
     (-(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
          P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
          (-2*k*q*x + pow(k,2) + pow(q,2))*
          (14*k*x + q*(2 - 7*f1*(-1 + pow(x,2)) + 12*pow(x,2)))) + 
       Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
        P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
        (2*k*q*x + pow(k,2) + pow(q,2))*
        (14*k*x + q*(7*f1*(-1 + pow(x,2)) - 2*(1 + 6*pow(x,2)))))*
     pow(-2*k*q*x + pow(k,2) + pow(q,2),-1)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-1))/7. ;
}

double I33 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 2*Pi*P11(q,params)*pow(f1,2)*pow(k,2)*pow(q,-1)*
   (Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
      P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
      (2*x*pow(k,3) + pow(q,3)*(1 - 3*pow(x,2)) + 
        2*k*x*pow(q,2)*(2 - 3*pow(x,2)) + q*pow(k,2)*(1 + pow(x,2))) + 
     Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
      P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
      (-2*x*pow(k,3) + pow(q,3)*(1 - 3*pow(x,2)) + 
        q*pow(k,2)*(1 + pow(x,2)) + 2*k*x*pow(q,2)*(-2 + 3*pow(x,2))))*
   pow(-2*k*q*x + pow(k,2) + pow(q,2),-1)*
   pow(2*k*q*x + pow(k,2) + pow(q,2),-1) ;
}

double I34 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I35 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I36 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I37 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I38 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I39 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I40 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return 0. ;
}

double I41 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (Pi*P11(q,params)*pow(f1,2)*pow(k,2)*pow(q,-2)*
     (-112*P11(k,params)*pow(pow(k,4) + pow(q,4) + 
          2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)),-1)*
        (7*pow(k,4)*pow(x,2) + 
          2*pow(k,2)*pow(q,2)*
           (1 + 2*pow(x,2) - 10*pow(x,4) + 3*f1*pow(-1 + pow(x,2),2)) + 
          pow(q,4)*(2 + 9*pow(x,2) - 4*pow(x,4) + 
             6*f1*pow(-1 + pow(x,2),2))) + 
       pow(k,2)*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
           (392*pow(k,2)*pow(x,2) - 
             56*k*q*x*(7*f1*(-1 + pow(x,2)) - 2*(1 + 6*pow(x,2))) + 
             pow(q,2)*(-56*f1*(-1 - 5*pow(x,2) + 6*pow(x,4)) + 
                147*pow(f1,2)*pow(-1 + pow(x,2),2) + 
                8*pow(1 + 6*pow(x,2),2))) + 
          Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
           (392*pow(k,2)*pow(x,2) + 
             56*k*q*x*(7*f1*(-1 + pow(x,2)) - 2*(1 + 6*pow(x,2))) + 
             pow(q,2)*(-56*f1*(-1 - 5*pow(x,2) + 6*pow(x,4)) + 
                147*pow(f1,2)*pow(-1 + pow(x,2),2) + 
                8*pow(1 + 6*pow(x,2),2))))))/392. ;
}

double I42 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (Pi*P11(q,params)*pow(f1,3)*pow(k,2)*pow(q,-2)*
     pow(-2*k*q*x + pow(k,2) + pow(q,2),-2)*
     pow(2*k*q*x + pow(k,2) + pow(q,2),-2)*
     (4*P11(k,params)*(pow(k,4) + pow(q,4) + 
          2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)))*
        (7*f1*(pow(k,4) + pow(q,4) + 
             2*pow(k,2)*pow(q,2)*(1 - 2*pow(x,2)))*(-1 + pow(x,2)) - 
          4*(7*pow(k,4)*pow(x,2) + 
             pow(k,2)*pow(q,2)*(-3 + 14*pow(x,2) - 25*pow(x,4)) + 
             pow(q,4)*(-3 + 19*pow(x,2) - 9*pow(x,4)))) - 
       pow(k,2)*(Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (7*f1*(-1 + pow(x,2))*
              (12*k*q*x + 2*pow(k,2) + 3*pow(q,2)*(-1 + 5*pow(x,2))) + 
             4*(k*q*x*(5 - 33*pow(x,2)) - 14*pow(k,2)*pow(x,2) + 
                pow(q,2)*(1 + 3*pow(x,2) - 18*pow(x,4))))*
           pow(-2*k*q*x + pow(k,2) + pow(q,2),2) + 
          Heaviside(-q + pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (7*f1*(-1 + pow(x,2))*
              (-12*k*q*x + 2*pow(k,2) + 3*pow(q,2)*(-1 + 5*pow(x,2))) + 
             4*(-14*pow(k,2)*pow(x,2) + k*q*x*(-5 + 33*pow(x,2)) + 
                pow(q,2)*(1 + 3*pow(x,2) - 18*pow(x,4))))*
           pow(2*k*q*x + pow(k,2) + pow(q,2),2))))/28. ;
}

double I43 (const double & k, const double & q, const double & x, const ParamsP11 & params) {

	double f1 = params.f ;

	return (Pi*P11(q,params)*pow(f1,4)*pow(k,2)*pow(q,-2)*
     (-8*P11(k,params)*(-1 + 3*pow(x,2)) + 
       pow(k,2)*(Heaviside(-q + 
             pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(-2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (4*pow(k,2)*(-1 + 3*pow(x,2)) - 8*k*q*x*(-3 + 5*pow(x,2)) + 
             pow(q,2)*(3 - 30*pow(x,2) + 35*pow(x,4)))*
           pow(-2*k*q*x + pow(k,2) + pow(q,2),-2) + 
          Heaviside(-q + pow(2*k*q*x + pow(k,2) + pow(q,2),0.5))*
           P11(pow(2*k*q*x + pow(k,2) + pow(q,2),0.5),params)*
           (4*pow(k,2)*(-1 + 3*pow(x,2)) + 8*k*q*x*(-3 + 5*pow(x,2)) + 
             pow(q,2)*(3 - 30*pow(x,2) + 35*pow(x,4)))*
           pow(2*k*q*x + pow(k,2) + pow(q,2),-2))))/8. ;
}
