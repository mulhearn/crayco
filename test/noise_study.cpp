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
#include "TTree.h"
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
TH1F * hllr;
TH1F * hsigf;
TTree * tree;
//float CUT1, CUT2;
int nhit;
int ncore;
double ereco;
double eshower;

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
  mc.optimize_weighted_noise(1E-20, VERBOSE);
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


void noise_study(shower_mc & mc, int nruns, double Emin, double Emax, const char * htag, double rate_scale){
  shower_fcn & sfcn = shower_fcn::instance();

  int noise_only = 0;
  if (Emax <= 0) {
    cout << "INFO:  Generating Noise Only" << "\n";
    noise_only = 1;
    mc.optimize_weighted_noise(1E-15, VERBOSE);
  } else {
    cout << "INFO:  Generating showers in range:  " << Emin << " to " << Emax << "\n";
  }
  
  // setup the histograms
  char name[100];
  sprintf(name, "h_%s_ereco", htag);
  hereco      = new TH1F(name,"",80, 18.0, 22.0);
  sprintf(name, "h_%s_nhits", htag);
  hnhits      = new TH1F(name,"",500, 0, 1000);  
  sprintf(name, "h_%s_llr", htag);
  hllr      = new TH1F(name,"",100, 0.0, 1.0);  
  sprintf(name, "h_%s_sigf", htag);
  hsigf      = new TH1F(name,"",100, 0, 10);  

  if (noise_only){
    hereco->Sumw2();
    hnhits->Sumw2();
    hllr->Sumw2();
    hsigf->Sumw2();
  } 

  double count = 0;
  double tot_wgt = 0;
  double pass_wgt = 0;
  int i=0;
  while (count < nruns){

    if (i%100 == 0) cout << i << "\n";
    i++;
    
    if (noise_only) {
      mc.weighted_noise(QUIET); 
      tot_wgt += mc.wgt;
      count   += mc.wgt;
    } else {
      mc.generate(mc.weighted_shower_energy(Emin,Emax,nruns,QUIET), QUIET, 0.0);
      tot_wgt += mc.wgt;
      count   += 1.0;
    } 
    //cout << "actual generated s_loge:  " << sfcn.gen_s_loge << "\n";

    sfcn.gen_s_loge = log(1E10);
    //sfcn.d_flat     = 0.0; 
    if (shower_fit(QUIET)){
      //cout << "fitted s_loge:  " << sfcn.fit_s_loge << "\n";
      //cout << "uncertainty s_loge:  " << sfcn.unc_s_loge << "\n";      
      double sigf     = sfcn.fit_s_loge / sfcn.unc_s_loge / sqrt(mc.tot_n());
      double llr       = (sfcn.noise_only_fcn() - sfcn.fmin) / mc.tot_n();      
      double log10_e_fit = sfcn.fit_s_loge / log(10.0);      
      //if (noise_only && (log10_e_fit>11)){
      //cout << "FOUND A PROBLEMATIC GUY...\n";
      //cout << e_fit << "\n";
      //cout << sig << "\n";
      //cout << sfcn.fmin << "\n";
      //shower_fit(VERBOSE);
      //} 

      // fill remaining ntuple fields
      nhit    = sfcn.count_hits();
      //ncore   = sfcn.count_hits_core(200.0);
      ncore   = sfcn.count_hits_core(150.0);
      ereco   = sfcn.fit_s_loge / log(10.0);
      //cout << "ereco:  " << setw(12) << ereco << " wgt:  " << mc.wgt << "\n";
      eshower = 0.0;
      //if (Emax >0.0) eshower = log(Emax) / log(10.0);
      tree->Fill();

      hnhits->Fill(sfcn.count_hits(), mc.wgt*rate_scale);



      //if (! noise_only){
      //cout << sigf << " " << chisqndf << "\n";	
      //}
      //cout << "llr:  " << llr << "\n";
      
      //if ((! noise_only) || (log10_e_fit > 10)){ // plot only the bad guys!
      hllr->Fill(llr, mc.wgt*rate_scale);
      hsigf->Fill(sigf, mc.wgt*rate_scale);
      //}
      
      //if ((CUT1>0.0) && (llr < CUT1)) continue;
      //if ((CUT2>0.0) && (sigf < CUT2)) continue;	  
      hereco->Fill(sfcn.fit_s_loge / log(10.0), mc.wgt*rate_scale);
      pass_wgt += mc.wgt;
    }    
  }

  cout << "INFO:  threw " << i << " events with total weight of " << tot_wgt << "\n";
  cout << "INFO:  selection efficiency:  " << 100 * pass_wgt / tot_wgt << "\n";
  hereco->Write();
  hnhits->Write();
  hllr->Write();
  hsigf->Write();
}


int main(int argc, char * argv[]){      
   shower_fcn & sfcn = shower_fcn::instance();
   shower_mc mc;


   //cout << flux(1E20)*1E20 << "\n";
   //cout << flux(1E19)*1E19 << "\n";
   //return 0;

   if (argc < 6) { 
     cout << "usage:  noise_study <nruns> <nperkm2> <xs_mu> <xs_gamma> <tres>\n";
     cout << "eg: noise_study 100 1000 5E-5 0.0 0.2\n";     
     return 0; 
   }





   int nruns     = atoi(argv[1]);
   mc.d_n        = atoi(argv[2]);
   mc.d_xs_mu    = atof(argv[3]);
   mc.d_xs_gamma = atof(argv[4]);
   double tres   = atof(argv[5]);
   double zbrate = 1.0 / tres;
   double murate = 0.01; // Hz
   mc.d_flat = murate * tres;

   cout << "INFO: timing resolution:  " << tres << " seconds.\n";
   cout << "INFO: zero bias rate:  :  " << zbrate << " Hz.\n";
   cout << "INFO: assuming background muon rate is : " << murate << " Hz.\n";
   cout << "INFO: flat rate is :  " << mc.d_flat << "\n";

        
   //CUT1 = atof(argv[5]);
   //CUT2 = atof(argv[6]);
   //cout << "applying CUT1:  " << CUT1 << "\n";
   //cout << "applying CUT2:  " << CUT2 << "\n";

   mc.rng.SetSeed(2014);

   mc.d_size     = 0.5;

   TFile fout("noise.root", "RECREATE");
   fout.cd();

   tree = new TTree("noise","");

   std::vector<int> *p_d_h = &sfcn.d_h;
   std::vector<double> *p_d_x = &sfcn.d_x;
   std::vector<double> *p_d_y = &sfcn.d_y;
   tree->Branch("d_h","vector<int>",&p_d_h);
   tree->Branch("d_x","vector<double>",&p_d_x);
   tree->Branch("d_y","vector<double>",&p_d_y);
   tree->Branch("xs_mu",&mc.d_xs_mu);
   tree->Branch("xs_gamma",&mc.d_xs_gamma);
   tree->Branch("d_flat",&mc.d_flat);
   tree->Branch("wgt",&mc.wgt);
   tree->Branch("nhit",&nhit);
   tree->Branch("ncore",&ncore);
   tree->Branch("ereco",&ereco);
   tree->Branch("eshower",&eshower);


   mc.print();
   noise_study(mc, nruns, 0, 0, "noise", zbrate / nruns);
   noise_study(mc, 10*nruns, 1E19, 2E19, "signala", 1.0);
   noise_study(mc, nruns, 1E20, 2E20, "signalb", 1.0);

   tree->Write();
   fout.Close();

   return 0;
}
