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
#include <Math/ProbFuncMathCore.h>

using namespace std;

class shower_fcn {
public:
   // needed for TMinuit fcn interface
   static shower_fcn & instance(){ 
      static shower_fcn x;
      return x;
   }

   // the actual fcn using data from the class shower_fcn:
   double calc_fcn(double * par){
      const double norm = 1.0 / sqrt(2 * M_PI);
      double logl = 0.0;
      double logn = fabs(par[0]);
      double sigx = fabs(par[1]);
      double sigy = fabs(par[2]);
      //cout << "sigx:  " << sigx << "\n";
      //cout << "sigy:  " << sigy << "\n";

      //double x = sigx - gen_s_sigx;
      //double y = sigy - gen_s_sigy;
      //logl = -x*x - y*y;

      for (int i=0; i<d_x.size();i++){
         double dx = d_x[i]/sigx;
         double dy = d_y[i]/sigy;
         double p  = norm * ( exp(-0.5*dx*dx)/sigx * exp(-0.5*dy*dy)/sigy);
         double mu = exp(logn) * p * d_xs;
         double p0 = ROOT::Math::poisson_cdf(0,mu);
         double p1 = 1.0 - p0;
         //cout << "mu:  " << mu << "\n";
         //cout << "p0:  " << p0 << "\n";
         //cout << "p1:  " << p1 << "\n";

         if (d_h[i]) logl += log(p1);
         else        logl += log(p0);
      }
      //cout << "logl =   " << logl << "\n";
      return -2.0*logl;
   }

   // detector data:
   vector<double> d_x;
   vector<double> d_y;
   vector<int>    d_h;
   //fixed parameters (for now)
   double d_xs;
   
   // parameters used for generating shower data:
   double gen_s_logn;
   double gen_s_sigx;   
   double gen_s_sigy;

   // fitted parameters:
   double fit_s_logn;
   double fit_s_sigx;
   double fit_s_sigy;
   double unc_s_logn;
   double unc_s_sigx;
   double unc_s_sigy;


private:
   shower_fcn(){}
};

// TMinuit fcn interface:
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   f = shower_fcn::instance().calc_fcn(par);
}


void generate(TRandom & rng, int s_logn, double s_sigx, double s_sigy, int d_n, double d_xs){ 
   const double norm = 1.0 / sqrt(2 * M_PI);
   shower_fcn & sfcn = shower_fcn::instance();
   sfcn.d_x.clear();
   sfcn.d_y.clear();
   sfcn.d_h.clear();
   sfcn.d_xs = d_xs;
   sfcn.gen_s_logn      = s_logn;
   sfcn.gen_s_sigx   = fabs(s_sigx);
   sfcn.gen_s_sigy   = fabs(s_sigy);

   for (int i=0; i<d_n; i++){
      double x = -500 + rng.Uniform() * 1000;
      double y = -500 + rng.Uniform() * 1000;
      double dx = x/s_sigx;
      double dy = y/s_sigy;
      double p = norm * ( exp(-0.5*dx*dx)/s_sigx * exp(-0.5*dy*dy)/s_sigy);
      double mu = p * exp(s_logn) * d_xs;
      double p0 = ROOT::Math::poisson_cdf(0,mu);

      double hit = 0;
      double xp = rng.Uniform();
      if (xp > p0) hit = 1;
      sfcn.d_x.push_back(x);
      sfcn.d_y.push_back(y);
      sfcn.d_h.push_back(hit);
   }

   //shower_fcn & sfcn = shower_fcn::instance();
   //sfcn.d_x.clear();
   //sfcn.d_y.clear();
   //sfcn.gen_s_logn      = N;
   //sfcn.gen_s_sigx   = fabs(s_sigx);
   //sfcn.gen_s_sigy   = fabs(s_sigy);

   //for (int i=0; i<N; i++){
   //double x = rng.Gaus() * s_sigx;
   //  double y = rng.Gaus() * s_sigy;
   // cout << "x = " << setw(8) << x << " ";
   // cout << "y = " << setw(8) << y << "\n";
   // sfcn.d_x.push_back(x);
   // sfcn.d_y.push_back(y);
   //}
}

void do_fit(){
   shower_fcn & sfcn = shower_fcn::instance();
   
   //for (int i=0; i<sfcn.d_x.size(); i++){
      //cout << "x = " << sfcn.d_x[i] << " ";
      //cout << "y = " << sfcn.d_y[i] << " ";
      //cout << "h = " << sfcn.d_h[i] << "\n";
   //}

   TMinuit *gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
   gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
   
   // Set starting values and step sizes for parameters
   static Double_t vstart[3] = {7.0, 100.0, 100.0};
   static Double_t step[3]   = {0.01, 0.01,  0.01};
   gMinuit->mnparm(0, "s_logn",  vstart[0], step[0], 0,0,ierflg);
   gMinuit->mnparm(1, "s_sigx",  vstart[1], step[1], 0,0,ierflg);
   gMinuit->mnparm(2, "s_sigy",  vstart[2], step[2], 0,0,ierflg);

   //gMinuit->FixParameter(0);
   
   // Now ready for minimization step
   arglist[0] = 1000;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   
   gMinuit->GetParameter(0, sfcn.fit_s_logn, sfcn.unc_s_logn);
   gMinuit->GetParameter(1, sfcn.fit_s_sigx, sfcn.unc_s_sigx);
   gMinuit->GetParameter(2, sfcn.fit_s_sigy, sfcn.unc_s_sigy);

   cout << "generated logn:  " << setw(12) << sfcn.gen_s_logn << " ";
   cout << "fitted value:  " << setw(12) << sfcn.fit_s_logn << "+/-" << sfcn.unc_s_logn << "\n";
   cout << "generated sigx:  " << setw(12) << sfcn.gen_s_sigx << " ";
   cout << "fitted value:  " << setw(12) << sfcn.fit_s_sigx << "+/-" << sfcn.unc_s_sigx << "\n";
   cout << "generated sigy:  " << setw(12) << sfcn.gen_s_sigy << " ";
   cout << "fitted value:  " << setw(12) << sfcn.fit_s_sigy << "+/-" << sfcn.unc_s_sigy << "\n";  
}


int main(int argc, char *argv[])
{
   TRandom rng;
   double s_logn = 8;
   double s_sigx = 210.0;
   double s_sigy = 300.0;
   int    d_n    = 10000;
   double d_xs   = 1.0;

   cout << "s_n:  " << exp(s_logn) << "\n";

   generate(rng, s_logn, s_sigx, s_sigy, d_n, d_xs); 
   do_fit();

}

