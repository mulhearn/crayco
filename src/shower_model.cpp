#include <math.h>
#include <cmath>
//#include <vector>

#include "shower_model.h"
#include "utils.h"

const double rM = 79.0; // m (EAS page 199).

const int np=12;
float fluxd[np];

// eV
float energy[np] = {
   1e15,
   5e15,
   1e16,
   5e16,
   1e17,
   5e17,
   1e18,
   5e18,
   1e19,
   5e19,
   1e20,
    5e20
};

void fill_flux()
{
   //  E^{2.6} * flux, GeV^2.6/(GeV m^2 s sr)
   float Eflux_plot[np] = {
      1e4, // 1e15
      8e3,
      6e3, // 1e16
      3e3,
      2e3,// 1e17
      1e3,
      6e2,// 1e18
      2e2,
      1.5e2, // 1e19
      1e2,
      5e1, // 1e20
      5e0, // 5e20
   };

   for (int i=0;i<np;i++)
   {
      fluxd[i] = Eflux_plot[i] 
         * 1e6  //  m^2 /km^2
         * 3.155e7  //s/year
         / pow((energy[i]*1e-9),2.6); // remove E^2.6
   }

}

bool got_flux=false;

// e in eV
// flux in number / (km2 yr sr GeV)
double flux(double e)
{
   if (!got_flux)
   {
      fill_flux();
      got_flux=true;
   }
  
   if (e > 5E20) return 0.0;

   // look for our index
   int ilow=0;
   while ( energy[ilow+1] < e)
      ilow++;
   
   // linear interpolation of the logs
   double loge_low = log(energy[ilow]);
   double loge_hi = log(energy[ilow+1]);
   double logf_low = log(fluxd[ilow]);
   double logf_hi = log(fluxd[ilow+1]);
  
   double logf_int = logf_low + (log(e)-loge_low)*(logf_hi - logf_low)/(loge_hi - loge_low);
   return exp(logf_int);   
}

double logn_gamma (double logE){
   // from Daniel's fit to Corsika...
   double a = -2.36;
   double b = 1.15;
   double logn = a + b*(logE);
   return logn;
}

double logn_mu(double logE){
  // old mu- only:
  //double a = -4.36;
  //double b = 0.93;
  // new mu- and mu+
  double a = -3.57;
  double b = 0.93;
  double logn = a + b*(logE);
  return logn;
}

double logE_gamma (double logn_gamma){
   // invert above:
   double a = -2.36;
   double b = 1.15;
   double logE = (logn_gamma - a) / b;
   return logE;
}

double logE_mu (double logn_gamma){
   // invert above:
   double a = -4.36;
   double b = 0.93;
   double logE = (logn_gamma - a) / b;
   return logE;
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
