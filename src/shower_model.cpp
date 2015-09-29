#include <math.h>
#include <cmath>
//#include <vector>
#include <iostream>

#include "shower_model.h"
#include "utils.h"

using namespace std;

const double rM = 79.0; // m (EAS page 199).
//const double rM = 100.0; // m (EAS page 199).


// e is in eV
// units of ( km^2 * year * sr * eV) ^-1
double flux(double e)
{  
  double gamma1 = 3.23;
  double gamma2 = 2.63;
  double logehalf=19.63;
  double logwc = 0.15;  
  double norm_low = 1.1577e44*8.01e-3; // to match normaliztion at 18.72
  double norm_hi  = 5.45e30; // normaliztion from fig5 of 1307.5059v1
  if (log10(e)<18.72)
    return norm_low*pow(e,-gamma1);    
  return norm_hi*( pow(e,-gamma2)/(1.0+exp( (log10(e)-logehalf)/(logwc))));
}




double logn_gamma (double logE){
   // from Daniel's fit to Corsika...
   double a = -2.36;
   double b = 1.15;
   double logn = a + b*(logE - log(1E9));  // Daniel's params in terms of GeV!!!
   return logn;
}

double logn_mu(double logE){
  // old mu- only:
  //double a = -4.36;
  //double b = 0.93;
  // new mu- and mu+
  double a = -3.57;
  double b = 0.93;
  double logn = a + b*(logE - log(1E9));  // Daniel's params in terms of GeV!!!
  return logn;
}

double logE_gamma (double logn_gamma){
   // invert above:
   double a = -2.36;
   double b = 1.15;
   double logE = (logn_gamma - a) / b;
   return logE+ log(1E9);   // Daniel's params in terms of GeV!!!
}

double logE_mu (double logn_gamma){
   // invert above:
   double a = -4.36;
   double b = 0.93;
   double logE = (logn_gamma - a) / b;
   return logE + log(1E9);  // Daniel's params in terms of GeV!!!
}

double shower_pdf(double x, double y, double s_sin2theta, double s_phi){
  // Debug with Gaussian:
  //const double gnorm = 1.0 / sqrt(2 * M_PI);
  //double dx = x/100*s_theta;
  //double dy = y/100*s_phi;
  //return gnorm * ( exp(-0.5*dx*dx)/(100*s_theta) * exp(-0.5*dy*dy)/(100*s_phi));

  //cout << s_theta << "\n";
  //cout << s_phi << "\n";

  // don't understand this yet, just copying from Daniel:
  double R = sqrt( x*x + y*y );
  double det_phi = std::atan2(y,x);      
  //double r = 0.5*R * sqrt(3.0-cos(2*(det_phi-s_phi)) + 2*cos(det_phi-s_phi)*cos(det_phi-s_phi)*cos(2*s_sin2theta));

  double cos2_dphi  = sq(cos(det_phi - s_phi));
  // simplified version of above:
  double r = R * sqrt(1.0 - cos2_dphi * s_sin2theta);

  // from EAS, page 199, eq 4.162
  // with overall factor of N_e removed
  //const double s = 1.9;
  const double s = 1.9;

  double norm = (1.0/(2*3.14159*rM*rM));
  double r1 = pow(r/rM,s-2);   
  double r2 = (double) pow( 1.0 + (r/rM),s-4.5);
  double gamma = tgammaf(4.5-s)/(tgammaf(s)*tgammaf(4.5-(2*s)));
  double pdf = norm*r1*r2*gamma;
  //cout << "r:  " << r << "\n";
  //cout << "pdf:  " << pdf << "\n";
  return pdf;
} 
