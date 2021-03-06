#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <math.h>
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TCanvas.h"
#include <iomanip>

#include "TMinuit.h"
#include <Math/ProbFuncMathCore.h>

#include "shower_model.h"
#include "utils.h"
#include "shower_fit.h"
#include "shower_mc.h"

using namespace std;

const double shower_mc::sin2theta_max = sq(sin(SHOWER_MC_THETA_MAX));

void shower_mc::print(){
   cout << "INFO: phones / km^2         " << d_n << "\n";
   cout << "INFO: size of grid (km):    " << d_size << "\n";
   cout << "INFO: total phones:         " << tot_n() << "\n";
   cout << "INFO: xsec times gamma eff: " << d_xs_gamma << "\n";
   cout << "INFO: xsec times muon eff:  " << d_xs_mu << "\n";
   cout << "INFO: flat background:      " << d_flat << "\n";
}

void shower_mc::init_shower_fcn(){
  shower_fcn & sfcn = shower_fcn::instance();
  //wgt = 1.0; // this function produces  unweighted events
  if (! mode_manual_locations){
    sfcn.d_x.clear();
    sfcn.d_y.clear();
  }
  sfcn.d_h.clear();
  sfcn.d_xs_gamma       = d_xs_gamma;
  sfcn.d_xs_mu          = d_xs_mu;
  sfcn.d_flat           = d_flat; 
  sfcn.d_size           = d_size; 
  sfcn.gen_s_loge       = 0;
  sfcn.gen_s_sin2theta  = 0;
  sfcn.gen_s_phi        = 0;
}

void shower_mc::flat_in_circle(double r, double & x, double & y){
  double r2 = r*r;  
  do {
    x = (-1.0 + 2.0*rng.Uniform()) * r;
    y = (-1.0 + 2.0*rng.Uniform()) * r;    
  } while ((x*x + y*y) > r2);
}


int shower_mc::prescale(double f){
  double x = rng.Uniform()*f;
  if (x > 1.0){
    wgt = 0;
    return 0;
  } else {
    wgt = wgt*f;
    return 1;
  }
}

void shower_mc::generate_locations(){
  if (mode_manual_locations) return;
  shower_fcn & sfcn = shower_fcn::instance();
  int tot_n = this->tot_n();
  for (int i=0; i<tot_n; i++){
    double x = (-500.0 + rng.Uniform() * 1000) * d_size;
    double y = (-500.0 + rng.Uniform() * 1000) * d_size;
    sfcn.d_x.push_back(x);
    sfcn.d_y.push_back(y);
  }
}

void shower_mc::generate(double s_loge, double s_sin2theta, double s_phi, int verbose, double s_x, double s_y){ 
  const double norm = 1.0 / sqrt(2 * M_PI);
   shower_fcn & sfcn = shower_fcn::instance();
   init_shower_fcn();
   generate_locations();

   //wgt = 1.0; // this function produces  unweighted events

   sfcn.gen_s_loge       = s_loge;
   sfcn.gen_s_sin2theta  = s_sin2theta;
   sfcn.gen_s_phi        = fabs(s_phi);
   sfcn.gen_s_x          = s_x;
   sfcn.gen_s_y          = s_y;

   double s_logn_mu    = logn_mu(s_loge);
   double s_logn_gamma = logn_gamma(s_loge);

   if (verbose){
     cout << "INFO: generate\n";
     cout << "INFO:  s_loge:        " << s_loge << "\n";
     cout << "INFO:  s_logn_mu:     " << s_logn_mu << "\n";
     cout << "INFO:  s_logn_gamma:  " << s_logn_gamma << "\n";
     cout << "INFO:  s_x:           " << s_x << "\n";
     cout << "INFO:  s_y:           " << s_y << "\n";
     cout << "INFO:  d_n:     " << d_n << "\n";
     cout << "INFO:  d_xs_gamma:    " << d_xs_gamma << "\n";
     cout << "INFO:  d_xs_mu:    " << d_xs_mu << "\n";
     cout << "INFO:  d_size:  " << d_size << "\n";
     cout << "INFO:  d_flat:  " << d_flat << "\n";
   }


   int nhit = 0;
   sfcn.d_h.clear();
   for (int i=0; i<sfcn.d_x.size(); i++){
     double p = shower_pdf(sfcn.d_x[i - s_x], sfcn.d_y[i] - s_y, s_sin2theta, s_phi);
     double mu = p * (exp(s_logn_gamma) * d_xs_gamma + exp(s_logn_mu) * d_xs_mu)+ d_flat;
     double p0 = ROOT::Math::poisson_cdf(0,mu);
      
     double hit = 0;
     hit = rng.Poisson(mu);
     if (hit > 0.0) nhit+=1;
     sfcn.d_h.push_back(hit);
   }
   if (verbose){
     cout << "INFO:  generated hits:     " << nhit << "\n";
   }
}

void shower_mc::generate(double E, int verbose, double fix_sin2theta, double s_x, double s_y){
   double s_loge       = log(E);
   double s_sin2theta  = fix_sin2theta;
   if (s_sin2theta < 0.0)
      s_sin2theta = rng.Uniform() * sin2theta_max;
   double s_phi        = rng.Uniform() * 2 * M_PI;

   generate(s_loge, s_sin2theta, s_phi, verbose, s_x, s_y);

   shower_fcn & sfcn = shower_fcn::instance();   
   if (verbose) cout << "INFO:  number of hits:  " << sfcn.count_hits() << "\n";
}

double shower_mc::log_prob_noise(int k){
  if (d_flat >= 1.0) return 0;
  if (d_flat <= 0.0) return 0;

  int n = tot_n();  
  double sum = 0.0;
  for (int i=0; i<k; i++){
    double r = ((double) (n-i)) / ((double) (k-i));
    sum = sum + log(r);
  }    
  sum += k * log(d_flat);
  sum += (n-k) * log(1.0 - d_flat);
  return sum;
}

double shower_mc::prob_noise(int k, double min){
  int min_log = log(min);

  static const int MIN = -100;
  double lp = log_prob_noise(k);
  if (lp < min_log) return 0.0;
  return exp(lp);
}

void shower_mc::noise(int verbose){
   shower_fcn & sfcn = shower_fcn::instance();
   init_shower_fcn();
   generate_locations();

   wgt = 1.0; // this function produces  unweighted events
   
   if (verbose){
     cout << "INFO:  noise\n";
     cout << "INFO:  d_n:     " << d_n << "\n";
     cout << "INFO:  d_size:  " << d_size << "\n";
     cout << "INFO:  d_flat:  " << d_flat << "\n";
   }

   int nhit = 0;     
   for (int i=0; i<sfcn.d_x.size(); i++){
      double xp = rng.Uniform();
      double hit = 0;
      if (xp < d_flat) {
	hit = 1;
	nhit++;
      }
      sfcn.d_h.push_back(hit);
   }
   if (verbose){
     cout << "INFO:  generated hits:     " << nhit << "\n";
   }

}

void shower_mc::optimize_weighted_noise(double pmin, int verbose){
  kmin = kmax = -1;
  pmax = 0;

  for (int k=0; k<tot_n(); k++){
    double p = prob_noise(k, pmin);
    //cout << "k:  " << k << "\n";
    //cout << "p:  " << p << "\n";
    if (p > 0){
      if (kmin < 0) kmin = k;
      if (p > pmax) pmax = p;
      kmax = k;
    } else {
      if (kmin >= 0) {
	break;
      }
    }
  }
  if (verbose) {
    cout << "kmin = " << kmin << "\n";
    cout << "kmax = " << kmax << "\n";
    cout << "pmax = " << pmax << "\n";
  }
}

void shower_mc::weighted_noise(int verbose){
  shower_fcn & sfcn = shower_fcn::instance();
  init_shower_fcn();

  if (verbose==VERBOSE){
    cout << "INFO:  weighted_noise\n";
    cout << "INFO:  d_n:           " << d_n << "\n";
    cout << "INFO:  d_size:        " << d_size << "\n";
    cout << "INFO:  d_flat:        " << d_flat << "\n";
    cout << "INFO:  kmin:          " << kmin << "\n";
    cout << "INFO:  kmax:          " << kmax << "\n";
    cout << "INFO:  pmax:          " << pmax << "\n";
  }
  
  generate_locations();
  
  int nhit = kmin + rng.Uniform()*(kmax-kmin);
  wgt = prob_noise(nhit) / pmax;  
  
  for (int i=0; i<sfcn.d_x.size(); i++){
    sfcn.d_h.push_back(i<nhit);
  }
  if (verbose){
    cout << "INFO:  generated hits:     " << nhit << "\n";
  }
}


double shower_mc::weighted_shower_energy(double Emin, double Emax, double ngen, int verbose){
  double NORM = 2*3.14 / 31536000;  // 2 pi / seconds per year -- convert flux to Hz
  
  double lmin = log(Emin);
  double lrat = log(Emax/Emin);
  
  double l = lmin + lrat*rng.Uniform();
  double E = exp(l);
  
  wgt = flux(E) * E * lrat / ngen * NORM;

  if (verbose){
    cout << "E shower:  " << E << "  wgt:  " << wgt << "\n";
    cout << "flux:  " << flux(E) << "\n";
    cout << "ngen:  " << ngen << "\n";
    cout << "wgt:  " << wgt << "\n";
  }
  return E;

}
