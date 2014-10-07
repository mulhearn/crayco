#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <math.h>
#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include <iomanip>

#include "TMinuit.h"


using namespace std;

const float max_theta = 1.22; // 70 degrees
const float energy = 1.0e17;
const float rM = 79.0;// m (EAS page 199).
const float H = 10e3; // 100 km up

//  log (E) = a + b log (N)
//  log(e) - a = b log N
//  (log(E)-a)/b = log N
// exp ( log (E) - a / b ) = N
// a = 1.28, b = 0.92 (proton, http://particle.astro.ru.nl/pub/icrc09-1044.pdf)

float nch (float E)
{
  float a = 1.28;
  float b = 0.92;
  float N = exp( (log(E) - a) / b );
  return N;
}


// from EAS, page 199, eq 4.162
// with overall factor of N_e removed
float ldf(float r, float s){
  float norm = (1.0/(2*3.14159*rM*rM));
  float r1 = pow(r/rM,s-2);
  float r2 = (float) pow( 1.0 + (r/rM),s-4.5);
  float gamma = tgammaf(4.5-s)/(tgammaf(s)*tgammaf(4.5-(2*s)));
  return norm
    *r1
    *r2
    *gamma;
}

// don't understand this yet, just copying from Daniel:
float calc_r(double x, double y, double shower_theta, double shower_phi){
   float R = sqrt( x*x + y*y );
   float det_phi = std::atan2(y,x);      
   float r = 0.5*R * sqrt(3.0-cos(2*(det_phi-shower_phi)) + 2*cos(det_phi-shower_phi)*cos(det_phi-shower_phi)*cos(2*shower_theta));
   return r;
}



class shower_fcn {
public:
   // needed for TMinuit fcn interface
   static shower_fcn & instance(){ 
      static shower_fcn x;
      return x;
   }

   // the actual fcn using data from the class shower_fcn:
   double calc_fcn(double * par){
      double logl = 0.0;
      double logn  = par[0];
      double theta = par[1];
      double phi   = par[2];
      //cout << "logn:   " << logn << "\n";
      //cout << "theta:  " << theta << "\n";
      //cout << "phi:    " << phi << "\n";
      
      // simple sanity check:
      //double x = (logn  - gen_s_logn) / 2.0;
      //double y = (theta - gen_s_theta) / 2.0;
      //double z = (phi   - gen_s_phi) / 2.0;
      //logl = -x*x -y*y -z*z;
      //nlogl = x*x + y*y + z*z;
      //cout << "nlogl:  " << logl << "\n";
      //return -logl;

      for (int i=0; i<d_x.size(); i++){         
         double r   = calc_r(d_x[i], d_y[i], theta, phi);
         double l   = ldf(r, gen_s_age);
         //double p   = exp(logn) * d_eff * l;         
         double p   = l;         
         if (p >= 1.0) continue;
         //cout << p << "\n";

         if (d_hit[i]){
            logl += log(p);
         } else {
            //logl += log(1-p);
         }
      }
      cout << "logl =   " << logl << "\n";
      return -2.0*logl;

   }

   // detector data:
   vector<double> d_x;
   vector<double> d_y;
   vector<int>    d_hit;
   // fixed parameters (for now)
   double d_eff;
   
   // parameters used for generating shower data:
   double gen_s_age;   
   double gen_s_logn;     
   double gen_s_theta; 
   double gen_s_phi;   


private:
   shower_fcn(){}
};

// TMinuit fcn interface:
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   f = shower_fcn::instance().calc_fcn(par);
}


void generate(TRandom & rng, 
              double s_age, double s_logn, double s_theta, double s_phi, int d_n, double d_eff){

   cout << "s_logn = " << s_logn << "\n";
   cout << "d_eff = " << d_eff << "\n";

   shower_fcn & sfcn = shower_fcn::instance();
   sfcn.d_x.clear();
   sfcn.d_y.clear();
   sfcn.d_hit.clear();
   sfcn.d_eff = d_eff;
   sfcn.gen_s_age   = s_age;
   sfcn.gen_s_logn  = s_logn;
   sfcn.gen_s_theta = s_theta;
   sfcn.gen_s_phi   = s_phi;

   for (int i=0; i<d_n; i++){
      double x = -500 + rng.Uniform() * 1000;
      double y = -500 + rng.Uniform() * 1000;
      double r   = calc_r(x,y,s_theta,s_phi);
      double l   = ldf(r,s_age);
      double p   = exp(s_logn) * d_eff * l;

      double hit = 0;
      double xp = rng.Uniform();
      if (xp < p) hit = 1;

      cout << "x = " << setw(8) << x << " ";
      cout << "y = " << setw(8) << y << " ";
      cout << "hit = " << setw(8) << hit << "\n";      
      cout << "r:     " << r << "\n";
      cout << "ldf:   " << l << "\n";
      cout << "prob:  " << p << "\n";
      //cout << "-log(p):    " << -log(p) << "\n";
      //cout << "-log(1-p):  " << -log(1-p) << "\n";      

      sfcn.d_x.push_back(x);
      sfcn.d_y.push_back(y);
      sfcn.d_hit.push_back(hit);      
   }
}

void do_fit(){
   shower_fcn & sfcn = shower_fcn::instance();
   
   cout << "s_logn = " << sfcn.gen_s_logn << "\n";
   cout << "d_eff = " << sfcn.d_eff << "\n";

   for (int i=0; i<sfcn.d_x.size(); i++){
      cout << "x = " << sfcn.d_x[i] << " ";
      cout << "y = " << sfcn.d_y[i] << " ";
      cout << "hit = " << sfcn.d_hit[i] << "\n";
   }


   TMinuit *gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
   gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
   
   // Set starting values and step sizes for parameters
   static Double_t vstart[3] = {sfcn.gen_s_logn, sfcn.gen_s_theta, sfcn.gen_s_phi};
   static Double_t step[3]   = {0.01,    0.01,  0.01};
   gMinuit->mnparm(0, "s_logn",   vstart[0], step[0], 0,0,ierflg);
   gMinuit->mnparm(1, "s_theta",  vstart[1], step[1], 0,0,ierflg);
   gMinuit->mnparm(2, "s_phi",    vstart[2], step[2], 0,0,ierflg);

   gMinuit->FixParameter(0);

   
   // Now ready for minimization step
   arglist[0] = 1000;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   
   // Print results
   //Double_t amin,edm,errdef;
   //Int_t nvpar,nparx,icstat;
   //gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   //gMinuit->mnprin(3,amin);   

}


int main(int argc, char *argv[])
{
   TRandom rng;

   double N   = nch(energy);
   double s_logn  = log(N);
   double s_age   = 1.0;
   double s_theta = 0.15; 
   double s_phi = 1.1; 
   int d_n = 10000; 
   double d_eff = 1E-11;
   vector<double> d_x; 
   vector<double> d_y; 
   vector<int> d_hit;

   cout << "d_eff = " << d_eff << "\n";

   
   generate(rng, s_age, s_logn, s_theta, s_phi, d_n, d_eff);   
   do_fit();

}

