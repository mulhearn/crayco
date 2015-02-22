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

TH1F * hnhits;
TGraph * gnhits;
TH1F * hereco;


void nhits_study(shower_mc & mc, double nruns, const char * tag){
  shower_fcn & sfcn = shower_fcn::instance();
  static const int MAX = 1000;
  static double n[MAX];
  static double p[MAX];
  int cnt = 0;

  char name[100];
  static const int NMIN = 0;
  static const int NMAX = 500;
  static const int NBIN = 250;
  double binsize = ((float) (NMAX - NMIN))/ ((float) NBIN);

  for (int k=0; k<mc.tot_n(); k++){
    if (cnt >= MAX) break;

    n[cnt] = k;
    p[cnt] = nruns*binsize*mc.prob_noise(k,1E-15);
    if (p[cnt] == 0.0){
      if (cnt > 0) break; // we are done!
      continue;  // we haven't started yet...
    }
    //cout << "noise prob for " << k << " hits " << p[cnt] << "\n";
    cnt++;
   }
   gnhits = new TGraph(cnt, &n[0], &p[0]);
   sprintf(name, "g_%s_nhits", tag);
   gnhits->SetName(name);
   gnhits->Write();


   double tot_wgt;



  sprintf(name, "h_%s_nhits_unwgt", tag);
  hnhits      = new TH1F(name,"",NBIN, NMIN, NMAX);
  tot_wgt = 0;
  while (tot_wgt < nruns){
    mc.noise(QUIET);
    tot_wgt += mc.wgt;     
    hnhits->Fill(sfcn.count_hits(), mc.wgt);
  }
  hnhits->Write();
  
  sprintf(name, "h_%s_nhits_wgt", tag);
  hnhits      = new TH1F(name,"",NBIN, NMIN, NMAX);
  hnhits->Sumw2();
  mc.optimize_weighted_noise(1E-15, VERBOSE);
  tot_wgt = 0;  
  cnt = 0;
  while (tot_wgt < nruns){
    cnt++;
    mc.weighted_noise(QUIET);
    tot_wgt += mc.wgt;     
    hnhits->Fill(sfcn.count_hits(), mc.wgt);
  }
  cout << "INFO:  it took " << cnt << " throws to reach target of " << nruns << "\n";
  hnhits->Write();
}


void noise_study(shower_mc & mc, int nruns, double E, const char * htag){
  shower_fcn & sfcn = shower_fcn::instance();

  int noise_only = 0;
  if (E <= 0) {
    cout << "INFO:  Generating Noise Only" << "\n";
    noise_only = 1;
    mc.optimize_weighted_noise(1E-15, VERBOSE);
  } else {
    cout << "INFO:  Energy of Generated Shower:  " << E << "\n";
  }


  

  
  // setup the histograms
  char name[100];
  sprintf(name, "h_%s_ereco_pass", htag);
  hereco      = new TH1F(name,"",80, 9.0, 13.0);
  sprintf(name, "h_%s_nhits_pass", htag);
  hnhits      = new TH1F(name,"",40, 0, 400);  
  if (noise_only){
    hereco->Sumw2();
    hnhits->Sumw2();
  } 

  double tot_wgt = 0;
  int i=0;
  while (tot_wgt < nruns){
    if (i%100 == 0) cout << i << "\n";
    i++;
    
    if (noise_only) mc.weighted_noise(QUIET); 
    else            mc.generate(E, QUIET, 0.0);
    tot_wgt += mc.wgt;
    
    //cout << "actual generated s_loge:  " << sfcn.gen_s_loge << "\n";
    sfcn.gen_s_loge = 10.0;
     //sfcn.d_flat     = 0.0; 
    if (shower_fit(QUIET)){
      //cout << "fitted s_loge:  " << sfcn.fit_s_loge << "\n";
      //cout << "uncertainty s_loge:  " << sfcn.unc_s_loge << "\n";
      double e_fit = exp(sfcn.fit_s_loge);
      double sig = sfcn.fit_s_loge / sfcn.unc_s_loge;

      if (sig > 100.0) {
	//cout << "logn significance:  " << sig << "\n";
	//cout << "fitted energy:  " << e_fit << "\n";
	hereco->Fill(log(e_fit) / log(10.0), mc.wgt);
	hnhits->Fill(sfcn.count_hits(), mc.wgt);
      }
    }
  }
  if (noise_only){
    cout << "INFO:  threw " << i << " events with total weight of " << tot_wgt << "\n";
  }


  hereco->Write();
  hnhits->Write();
}


int main(int argc, char * argv[]){      
   shower_fcn & sfcn = shower_fcn::instance();
   shower_mc mc;

   if (argc < 2) { 
     cout << "usage:  noise_study <nruns> \n";
     cout << "eg: noise_study 100\n";
     return 0; 
   }
   int nruns     = atoi(argv[1]);

   mc.rng.SetSeed(2014);

   mc.d_flat     = 0.2;
   mc.d_size     = 1.0;
   mc.d_n = 1000;
   mc.d_xs_gamma       = 0;
   mc.d_xs_mu          = 1E-5; 
   mc.print();   


   TFile fout("noise.root", "RECREATE");
   fout.cd();

   mc.d_flat     = 0.02;
   nhits_study(mc,10000, "02_noise");
   
   mc.d_flat     = 0.2;   
   nhits_study(mc,10000, "2_noise");   
    
   mc.d_flat     = 0.02;
   noise_study(mc, nruns, 0, "02_noise");
   noise_study(mc, nruns, 1E11, "02_signal");


   //mc.d_flat     = 0.2;
   //noise_study(mc, nruns, 0, "2_noise");
   //noise_study(mc, nruns, 1E11, "2_signal");
   


   fout.Close();

   return 0;
}
