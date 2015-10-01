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


int main(int argc, char * argv[]){
  shower_fcn & sfcn = shower_fcn::instance();
  shower_mc mc;

  if (argc < 10) { 
    cout << "usage:  resolution <tag> <nruns> <nperkm2> <xs_mu> <xs_gamma> <flat> <emin> <emax> <escl>\n";
    cout << "eg: resolution demo 100 1000 5E-5 1E19 1E21 10.0 \n";     
    return 0; 
  }

  const char * tag = argv[1];
  int nruns        = atoi(argv[2]);
  mc.d_n           = atoi(argv[3]);
  mc.d_xs_mu       = atof(argv[4]);
  mc.d_xs_gamma    = atof(argv[5]);
  mc.d_flat        = atof(argv[6]);
  double emin      = atof(argv[7]);
  double emax      = atof(argv[8]);
  double escl      = atof(argv[9]);

  mc.d_size = 1.0;
  mc.print();

  char fname[200];
  sprintf(fname, "resolution_%s.root", tag);
  TFile f(fname, "RECREATE"); 
   
  // data for the standard graphs.
  vector<double> pen, mhits;
  vector<double> pen_fit, erms, trms, prms, eunc, tunc, punc; 

  int ifit = 0;
  int fit_fails = 0;

  for (double E=emin; E<=emax; E *= escl){
    //cout << "INFO: evaluating primary energy " << E << "\n";
    calc_stats hits;
    int n = 1000;
    for (int i=0; i<n; i++){
      mc.generate(E, QUIET);
      double nhits = sfcn.count_hits();
      //cout << "INFO:  number of hits:  " << nhits << "\n";
      hits.add(nhits);
    }
    double mean_hits = hits.mean();
    cout << "INFO: E:  " << E << " phone hits mean:  " << mean_hits << " +/- " << hits.unc();
    cout << " rms: " << hits.rms() << "\n";
    pen.push_back(E);
    mhits.push_back(mean_hits);

    if (mean_hits >= 5.0){
      double cnt=0;
      calc_stats eres;
      calc_stats pres;
      calc_stats tres;
      for (int i=0; i<nruns; i++){
	mc.generate(E, QUIET, 0.75);
	if (shower_fit(QUIET)){      
	  ifit++;
	  double e_fit = exp(sfcn.fit_s_loge);
	  double e_gen = exp(sfcn.gen_s_loge);
	  double e_fres = (e_fit - e_gen) / e_gen;
	  //cout << "e_fres:  " << e_fres << "\n";
	  if (fabs(e_fres) < 5.0){
	    double dp = delta_phi(sfcn.fit_s_phi,sfcn.gen_s_phi);
	    double dt = delta_sin2theta(sfcn.fit_s_sin2theta, sfcn.gen_s_sin2theta);
	    eres.add(e_fres);
	    pres.add(dp);
	    tres.add(dt);
	    cnt++;
	  } else {
	    fit_fails++;
	  }
	} else {
	  fit_fails++;
	}      
      }	
      cout << " -->    frac energy resolution:  " << eres.rms() << " bias:  " << eres.mean() << "\n";
      cout << " -->    theta resolution:        " << tres.rms() << " bias:  " << tres.mean() << "\n";
      cout << " -->    phi resolution:          " << pres.rms() << " bias:  " << pres.mean() << "\n";
	
      pen_fit.push_back(E);
      erms.push_back(eres.rms());
      eunc.push_back(eres.rms() / sqrt(cnt));
      trms.push_back(tres.rms());
      tunc.push_back(tres.rms() / sqrt(cnt));
      prms.push_back(pres.rms());
      punc.push_back(pres.rms() / sqrt(cnt));
    }
  }

  cout << "INFO:  " << fit_fails << " failures out of " << ifit << " fits.\n";
 
  TGraphErrors *geres_err = new TGraphErrors(pen_fit.size(), &pen_fit[0], &erms[0], NULL, &eunc[0]);
  TGraphErrors *gtres_err = new TGraphErrors(pen_fit.size(), &pen_fit[0], &trms[0], NULL, &tunc[0]);
  TGraphErrors *gpres_err = new TGraphErrors(pen_fit.size(), &pen_fit[0], &prms[0], NULL, &punc[0]);
  TGraphErrors *gnhits_pen = new TGraphErrors(pen.size(), &pen[0], &mhits[0], NULL, NULL);

  f.cd();
  geres_err->SetName("eres");
  gtres_err->SetName("tres");
  gpres_err->SetName("pres");
  gnhits_pen->SetName("mhit");
  geres_err->Write();
  gtres_err->Write();
  gpres_err->Write();
  gnhits_pen->Write();
  f.Close();
   
  return 0;
}



