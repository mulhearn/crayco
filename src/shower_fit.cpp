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

#include "utils.h"
#include "shower_model.h"
#include "shower_fit.h"

using namespace std;

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  f = shower_fcn::instance().calc_fcn(par);
}


shower_fcn & shower_fcn::instance(){ 
  static shower_fcn x;
  return x;
}

double shower_fcn::count_hits(){
  double count_hits = 0;
  for (int i=0; i<d_x.size(); i++){
    if (d_h[i] > 0.0){ 
      count_hits+=1.0;
    }  
  }
  return count_hits;
}

double shower_fcn::calc_fcn(double * par){
  const double norm = 1.0 / sqrt(2 * M_PI);
  double logl     = 0.0;
  double loge      = par[0];
  double s_logn_mu    = logn_mu(loge);
  double s_logn_gamma = logn_gamma(loge);
  double sin2theta = par[1];
  double phi       = par[2];
  double flat      = par[3];
  //cout << "theta:  " << theta << "\n";
  //cout << "phi:  " << phi << "\n";
  
  //double x = theta - gen_s_theta;
  //double y = phi - gen_s_phi;
  //logl = -x*x - y*y;
  
  for (int i=0; i<d_x.size();i++){
    double p  = shower_pdf(d_x[i], d_y[i], sin2theta, phi);
    double mu = (exp(s_logn_gamma) * d_xs_gamma + exp(s_logn_mu) * d_xs_mu)* p + fabs(flat);
    //double p0 = ROOT::Math::poisson_cdf(0,mu);
    double p0 = exp(-mu);  // equivalent to above...
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

int shower_fit(int verbose){
  shower_fcn & sfcn = shower_fcn::instance();
  double count_hits = sfcn.count_hits();
  double fhit = count_hits / sfcn.d_x.size();
  if (verbose == VERBOSE){
     cout << "INFO:  number of phone hits " << count_hits << "\n";
     cout << "INFO:  number of phones " << sfcn.d_x.size() << "\n";
     cout << "INFO:  fraction of phones hit " << fhit << "\n";
  }

  // experimental, to remove bias:
  //count_hits = count_hits - sfcn.d_x.size() * sfcn.d_flat;

  if (count_hits < 5){
     return 0;
  }
  if (fhit > 0.90){
     return 0;
  }

  TMinuit *gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;

  if (verbose == QUIET){
     arglist[0] = -1;
     gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);
     gMinuit->mnexcm("SET NOWarnings",arglist,0,ierflg);
  }

  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  // Set starting values and step sizes for parameters
  //static Double_t vstart[3] = {sfcn.gen_s_loge,  sfcn.gen_s_theta,   sfcn.gen_s_phi};

  // this is a bit of cheat, but we should be able to get rough estimate before fit:
  // for now these are set near, but not identical to, the correct values...
  Double_t vstart[] = {sfcn.gen_s_loge, sfcn.gen_s_sin2theta,  sfcn.gen_s_phi, sfcn.d_flat};  
  
  //static Double_t vstart[3] = {sfcn.gen_s_loge, sfcn.gen_s_sin2theta,  sfcn.gen_s_phi};  
  Double_t step[]   = {100.0, 0.1, 0.1, 0.0};
  gMinuit->mnparm(0, "s_loge",       vstart[0], step[0], 0, 0, ierflg);
  gMinuit->mnparm(1, "s_sin2theta",  vstart[1], step[1], 0, 0, ierflg);
  gMinuit->mnparm(2, "s_phi",        vstart[2], step[2], 0, 0, ierflg);
  gMinuit->mnparm(3, "flat",         vstart[3], step[3], 0, 0, ierflg);

  //gMinuit->FixParameter(0);
  //gMinuit->FixParameter(1);
  //gMinuit->FixParameter(2);
  gMinuit->FixParameter(3);
  
  // Now ready for minimization step


  arglist[0] = 10000;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  //cout << "Minuit Status after MIGRAD:  " << gMinuit->GetStatus() << "\n";

  //gMinuit->Release(1);
  //gMinuit->Release(2);
  //gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
  gMinuit->GetParameter(0, sfcn.fit_s_loge, sfcn.unc_s_loge);
  gMinuit->GetParameter(1, sfcn.fit_s_sin2theta, sfcn.unc_s_sin2theta);
  gMinuit->GetParameter(2, sfcn.fit_s_phi, sfcn.unc_s_phi);

  if (verbose == VERBOSE){
     cout << "generated logn:       " << setw(12) << sfcn.gen_s_loge << " ";
     cout << "fitted value:         " << setw(12) << sfcn.fit_s_loge << " +/- " << sfcn.unc_s_loge << "\n";
     cout << "generated sin2 theta:  " << setw(12) << sfcn.gen_s_sin2theta << " ";
     cout << "fitted value:         " << setw(12) << sfcn.fit_s_sin2theta << " +/- " << sfcn.unc_s_sin2theta << "\n";
     cout << "generated phi:        " << setw(12) << sfcn.gen_s_phi << " ";
     cout << "fitted value:         " << setw(12) << sfcn.fit_s_phi << " +/- " << sfcn.unc_s_phi << "\n";  
  }

  delete gMinuit;
  return 1;
}

