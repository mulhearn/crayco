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

shower_fcn::shower_fcn(){
  int mode_cheat = 0;
  int mode_fix_theta_phi = 0;
  int mode_fix_x_y = 0;
}

shower_fcn & shower_fcn::instance(){ 
  static shower_fcn x;
  return x;
}


void shower_fcn::cheat(){ 
  start_s_loge       =   gen_s_loge;	  
  start_s_sin2theta  =   gen_s_sin2theta; 
  start_s_phi        =   gen_s_phi;	  
  start_s_x          =   gen_s_x; 	  
  start_s_y          =   gen_s_y;         
}

void shower_fcn::estimate(){ 
  start_s_sin2theta  =   0;
  start_s_phi        =   0;

  start_s_loge       =   log(1E18);
  int best_count = 0;
  start_s_x = 0;
  start_s_y = 0;
  double a = -d_size*500;
  double b = d_size*500;
  double step = 20;
  double MIN_R2 = sq(100);

  for (double x = a; x<b; x+=step){ 
    for (double y = a; y<b; y+=step){       
      int count = 0;
      for (int i=0; i<d_x.size();i++){
	if (d_h[i]==0) continue;
	double dx = (d_x[i]-x);
	double dy = (d_y[i]-y);
	double r2 = dx*dx + dy*dy;
	if (r2 < MIN_R2) count++;
      }
      //cout << x << " " << y << " " << count << "\n";
      if (count > best_count){
	best_count = count;
	start_s_x = x;
	start_s_y = y;
      }
    }
  }  
  //cout << "INFO: best count is  " << best_count << " out of a total of " << count_hits() << "\n";
  //cout << "INFO: starting x position is:  " << start_s_x << "\n";
  //cout << "INFO: starting y position is:  " << start_s_y << "\n";
}

double shower_fcn::count_hits(int nmin){
  double count_hits = 0;
  for (int i=0; i<d_x.size(); i++){
    if (d_h[i] >= nmin){ 
      count_hits+=1.0;
    }  
  }
  return count_hits;
}


double shower_fcn::count_hits_core(double r){
  double r2 = r*r;
  double count_hits = 0;
  for (int i=0; i<d_x.size(); i++){
    if (d_h[i] == 0) continue;
    double d_r2 = d_x[i]*d_x[i] + d_y[i]*d_y[i];
    if (d_r2 > r2) continue;
    count_hits+=1.0;
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
  double s_x       = par[4];
  double s_y       = par[5];

  //cout << "theta:  " << theta << "\n";
  //cout << "phi:  " << phi << "\n";
  
  //double x = theta - gen_s_theta;
  //double y = phi - gen_s_phi;
  //logl = -x*x - y*y;
  
  for (int i=0; i<d_x.size();i++){
    double p  = shower_pdf(d_x[i]-s_x, d_y[i]-s_y, sin2theta, phi);
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


double shower_fcn::noise_only_fcn(){
  double logl = 0;
  
  for (int i=0; i<d_x.size();i++){
    double mu = fabs(d_flat);
    double p0 = exp(-mu);  
    double p1 = 1.0 - p0;
    //cout << "mu:  " << mu << "\n";
    //cout << "p0:  " << p0 << "\n";
    //cout << "p1:  " << p1 << "\n";    
    if (d_h[i]) logl += log(p1);
    else        logl += log(p0);
  }
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

  //if (count_hits < 2){
  //return 0;
  //}
  if (fhit > 0.90){
     return 0;
  }

  if (sfcn.mode_cheat){
    static int warned = 0;
    if (!warned){
      cout << "WARNING:  cheating for the starting values, results may be optimistic...\n";
      warned = 1;
    }
    sfcn.cheat();  
  } else {
    sfcn.estimate();
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
  Double_t vstart[] = {sfcn.start_s_loge, sfcn.start_s_sin2theta,  sfcn.start_s_phi, sfcn.d_flat, sfcn.start_s_x, sfcn.start_s_y};  
  Double_t step[]   = {10.0, 0.1, 0.1, 0.0, 10.0, 10.0};

  gMinuit->mnparm(0, "s_loge",       vstart[0], step[0], 0, 0, ierflg);
  gMinuit->mnparm(1, "s_sin2theta",  vstart[1], step[1], 0, 1.0, ierflg);
  gMinuit->mnparm(2, "s_phi",        vstart[2], step[2], 0, 0, ierflg);
  gMinuit->mnparm(3, "flat",         vstart[3], step[3], 0, 0, ierflg);
  gMinuit->mnparm(4, "s_x",          vstart[4], step[4], 0, 0, ierflg);
  gMinuit->mnparm(5, "s_y",          vstart[5], step[5], 0, 0, ierflg);

  //gMinuit->FixParameter(0);
  gMinuit->FixParameter(1);
  gMinuit->FixParameter(2);
  gMinuit->FixParameter(3);
  gMinuit->FixParameter(4);
  gMinuit->FixParameter(5);
  
  // Now ready for minimization step
  arglist[0] = 10000;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  if (gMinuit->GetStatus() != 0){
    cout << "WARNING:  Minuit Status after MIGRAD (step 1):  " << gMinuit->GetStatus() << "\n";
  }

  if (!sfcn.mode_fix_x_y){
    gMinuit->Release(4);
    gMinuit->Release(5);
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    if (gMinuit->GetStatus() != 0){
      cout << "WARNING:  Minuit Status after MIGRAD (step 2):  " << gMinuit->GetStatus() << "\n";
    } 
  }

  if (!sfcn.mode_fix_theta_phi){
    gMinuit->FixParameter(0);
    gMinuit->FixParameter(3);
    gMinuit->FixParameter(4);
    gMinuit->FixParameter(5);

    // first find phi for a faily extreme value of theta:
    if (!sfcn.mode_cheat){
      gMinuit->mnparm(1, "s_sin2theta",  0.75, step[1], 0, 1.0, ierflg);
    }
    gMinuit->Release(2);
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    if (gMinuit->GetStatus() != 0){
      cout << "WARNING:  Minuit Status after MIGRAD (step 4):  " << gMinuit->GetStatus() << "\n";
    }

    // then find optimal theta holding phi constant:
    gMinuit->FixParameter(1);
    gMinuit->Release(2);
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    if (gMinuit->GetStatus() != 0){
      cout << "WARNING:  Minuit Status after MIGRAD (step 5):  " << gMinuit->GetStatus() << "\n";
    }
    // release everything but phi:
    gMinuit->Release(0);
    gMinuit->Release(2);
    if (!sfcn.mode_fix_x_y){
      gMinuit->Release(4);
      gMinuit->Release(5);
    }
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    if (gMinuit->GetStatus() != 0){
      cout << "WARNING:  Minuit Status after MIGRAD (step 6):  " << gMinuit->GetStatus() << "\n";
    } 

    // release everything:
    gMinuit->Release(1);
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    if (gMinuit->GetStatus() != 0){
      cout << "WARNING:  Minuit Status after MIGRAD (step 7):  " << gMinuit->GetStatus() << "\n";
    } 
  }

  gMinuit->GetParameter(0, sfcn.fit_s_loge, sfcn.unc_s_loge);
  gMinuit->GetParameter(1, sfcn.fit_s_sin2theta, sfcn.unc_s_sin2theta);
  gMinuit->GetParameter(2, sfcn.fit_s_phi, sfcn.unc_s_phi);
  gMinuit->mnstat(sfcn.fmin, sfcn.fedm, sfcn.errdef, sfcn.npari, sfcn.nparx, sfcn.istat);

  if (verbose == VERBOSE){
     cout << "generated logn:       " << setw(12) << sfcn.gen_s_loge << " ";
     cout << "fitted value:         " << setw(12) << sfcn.fit_s_loge << " +/- " << sfcn.unc_s_loge << "\n";
     cout << "generated sin2 theta:  " << setw(12) << sfcn.gen_s_sin2theta << " ";
     cout << "fitted value:         " << setw(12) << sfcn.fit_s_sin2theta << " +/- " << sfcn.unc_s_sin2theta << "\n";
     cout << "generated phi:        " << setw(12) << sfcn.gen_s_phi << " ";
     cout << "fitted value:         " << setw(12) << sfcn.fit_s_phi << " +/- " << sfcn.unc_s_phi << "\n";
     cout << "fmin:                 " << setw(12) << sfcn.fmin << "\n";
     cout << "fedm:                 " << setw(12) << sfcn.fedm << "\n";
     cout << "errdef:               " << setw(12) << sfcn.errdef << "\n";
     cout << "istat:                " << setw(12) << sfcn.istat << "\n";
  }
  if (sfcn.istat != 3){
    static int warn_count = 0;
    warn_count++;
    if (warn_count < 5){
      cout << "WARNING:  Minuit fit returns istat: " << setw(12) << sfcn.istat << "\n";
    } else if (warn_count == 5){
      cout << "WARNING:  Supressing further warnings regarding istat...\n";
    }
  }

  delete gMinuit;
  return 1;
}

