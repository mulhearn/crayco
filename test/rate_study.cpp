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
#include <vector>

#include "shower_model.h"
#include "utils.h"
#include "shower_fit.h"
#include "shower_mc.h"

using namespace std;

TTree * tree;
int sample;
int istat;
int nhit;
int ncore;
int mcore;
double ereco;
double ureco;
double eshower;
double sigf, llr;

int mode_dofit;
int mode_colocate_phones;
double zbrate;  // zero-bias rate, for scaling...
double rshower=0; 
double noise_tol = 1E-18;


double min_hits = 5;

void nhits_study(shower_mc & mc, double nruns, const char * tag){
  TGraph * gnhits;
  TH1F * hnhits;
  shower_fcn & sfcn = shower_fcn::instance();

  static const int MAX = 1000;
  static double n[MAX];
  static double p[MAX];
  int cnt = 0;
  
  char name[100];
  static const double NMIN = 0;
  static const double NMAX = 100;
  static const int NBIN = 100;
  static const double nmin = double(NMIN) - 0.5;
  static const double nmax = double(NMAX) + 0.5;
  static const int nbin = 101;
  double binsize = ((float) (NMAX - NMIN))/ ((float) NBIN);
  
  for (int k=0; k<mc.tot_n(); k++){
    if (cnt >= MAX) break;
    
    n[cnt] = k;
    p[cnt] = nruns*binsize*mc.prob_noise(k,1E-18);
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
  hnhits      = new TH1F(name,"",nbin, nmin, nmax);
  tot_wgt = 0;
  while (tot_wgt < nruns){
    mc.noise(QUIET);
    tot_wgt += mc.wgt;     
    hnhits->Fill(sfcn.count_hits(), mc.wgt);
    cout << mc.wgt << "\n";
  }
  hnhits->Write();
  
  sprintf(name, "h_%s_nhits_wgt", tag);
  hnhits      = new TH1F(name,"", nbin, nmin, nmax);
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


void shower_rate_study(shower_mc & mc, int nruns, double Emin, double Emax, const char * tag){  
  char name[100];
  sprintf(name, "h_%s_eshower_nowgt", tag);
  TH1F * hn      = new TH1F(name,"",100, 8, 21);  
  sprintf(name, "h_%s_eshower_wgt", tag);
  TH1F * hw      = new TH1F(name,"",100, 8, 21);  
  hw->Sumw2();
  
  for (int i=0; i<nruns; i++){
    double egen = log(mc.weighted_shower_energy(Emin,Emax,nruns,QUIET)) / log(10.0);    
    double wgt  = mc.wgt;
    hn->Fill(egen);
    hw->Fill(egen, wgt);    
  }
  hn->Write();
  hw->Write();
}

int max_core(double d_size, double r_size){
  shower_fcn & sfcn = shower_fcn::instance();
  int best_count = 0;
  double a = -d_size*500;
  double b = d_size*500;
  double step = 20;
  double MIN_R2 = sq(r_size);

  for (double x = a; x<b; x+=step){ 
    for (double y = a; y<b; y+=step){       
      int count = 0;
      for (int i=0; i<sfcn.d_x.size();i++){
	if (sfcn.d_h[i]==0) continue;
	double dx = (sfcn.d_x[i]-x);
	double dy = (sfcn.d_y[i]-y);
	double r2 = dx*dx + dy*dy;
	if (r2 < MIN_R2) count++;
      }
      //cout << x << " " << y << " " << count << "\n";
      if (count > best_count){
	best_count = count;
      }
    }
  }  
  return best_count;
}

void noise_study(shower_mc & mc, int nruns, double Emin, double Emax){
  shower_fcn & sfcn = shower_fcn::instance();
  double rate_scale = 1.0;
  
  
  


  int noise_only = 0;
  if (Emax <= 0) {
    cout << "INFO:  Generating Noise Only" << "\n";
    noise_only = 1;
    mc.optimize_weighted_noise(noise_tol, VERBOSE);
    rate_scale = zbrate / ((double) nruns);
    cout << "INFO:  will generate weighted equivalent of " << nruns << "\n";
    cout << "INFO:  zero-bias rate (externally set) is " << zbrate << "\n";
    cout << "INFO:  will scale event weights by additional factor of " << rate_scale << "\n";
  } else {
    double Ad = mc.d_size * mc.d_size;
    rate_scale = Ad;
    cout << "INFO:  Generating showers in range:  " << Emin << " to " << Emax << "\n";    
    cout << "INFO:  Will generate " << nruns << " of events weighted to integrate to flux in Hz/km^2\n";
    cout << "INFO:  simulated grid area:  (km^2)  " << Ad << "\n";
    cout << "INFO:  will scale event weights by additional factor of " << rate_scale << "\n";
  }
  
  double count = 0;
  double tot_wgt = 0;
  double pass_wgt = 0;
  double s_x=0, s_y=0;
  int i=0;
  while (count < nruns){
    if (i%100 == 0) cout << i << "\n";
    i++;
    
    if (mode_colocate_phones){
      double x = (-500.0 + mc.rng.Uniform() * 1000) * mc.d_size;
      double y = (-500.0 + mc.rng.Uniform() * 1000) * mc.d_size;
      sfcn.d_x.clear();
      sfcn.d_y.clear();
      for (int i=0; i<mc.d_n;i++){
	sfcn.d_x.push_back(x);
	sfcn.d_y.push_back(y);
      }
    }

    if (noise_only) {
      mc.weighted_noise(QUIET); 
      int n = sfcn.count_hits();

      if (n < 10) {if (!mc.prescale(10.0)) continue;}

      tot_wgt += mc.wgt;
      count   += mc.wgt; // not a bug!  (we produce the weighted equivalent of nruns of noise data...)  
      eshower = 0.0;      
    } else {
      if (rshower > 0.0){
	mc.flat_in_circle(1000.0*rshower, s_x, s_y);
	cout << "s_x:  " << s_x << "s_y:  " << s_y << "\n";
      }
      double egen = mc.weighted_shower_energy(Emin,Emax,nruns,QUIET);
      eshower = log(egen)/log(10.0);
      mc.generate(egen, QUIET, 0.0, s_x, s_y);  // NOTE theta fixed at theta=0... valid?
      tot_wgt += mc.wgt;
      count   += 1.0;
    } 
    sfcn.gen_s_loge = log(1E10);
    nhit    = sfcn.count_hits();      
    ncore   = sfcn.count_hits_core(150.0);
    mcore   = max_core(1.0, 150.0);
    mc.wgt  = mc.wgt * rate_scale;

    sigf  = 0;
    llr   = 0;
    ereco = 0.0;
    ureco = 0.0;
    istat = -1;
    //double noise_hits = mc.d_size * mc.d_n * mc.d_flat;
    if (mode_dofit && (nhit > min_hits)){
      if (shower_fit(QUIET)){
	sigf    = sfcn.fit_s_loge / sfcn.unc_s_loge / sqrt(mc.tot_n());
	llr     = (sfcn.noise_only_fcn() - sfcn.fmin) / mc.tot_n();            
	ereco   = sfcn.fit_s_loge / log(10.0);
	ureco   = sfcn.unc_s_loge / log(10.0);
	istat   = sfcn.istat;
      }   
    }
    tree->Fill();

  }
  cout << "INFO:  threw " << count << " events with total weight of " << tot_wgt << "\n";
}


int main(int argc, char * argv[]){      
   shower_fcn & sfcn = shower_fcn::instance();
   shower_mc mc;
   char name[100];

   sfcn.mode_fix_theta_phi = 1;
   mode_colocate_phones = 0;


   //cout << flux(1E20)*1E20 << "\n";
   //cout << flux(1E19)*1E19 << "\n";
   //return 0;

   if (argc < 9) { 
     cout << "usage:  noise_study <tag> <mode> <fit> <nruns> <nperkm2> <xs_mu> <xs_gamma> <tres> [ENERGIES]\n";
     cout << "eg: rate_study demo 0 0 100 1000 5E-5 0.0 0.2 1E16 1E19 1E20\n";     
     cout << "mode = 0:  default mode\n";
     cout << "mode = 1:  weighted/unweighted nhit cross check\n";
     cout << "mode = 2:  shower energy rate cross check\n";
     cout << "mode = 3:  fixed location rate study\n";
     return 0; 
   }

   const char * tag = argv[1];
   int mode      = atoi(argv[2]);
   mode_dofit    = atoi(argv[3]);
   int nruns     = atoi(argv[4]);
   mc.d_n        = atoi(argv[5]);
   mc.d_xs_mu    = atof(argv[6]);
   mc.d_xs_gamma = atof(argv[7]);
   double tres   = atof(argv[8]);
   vector<double> ebenchmark;
   for (int i=9; i<argc; i++){
     double E = atof(argv[i]);
     cout << "INFO: adding benchmark at E=" << E << "\n";
     ebenchmark.push_back(E);
   }

   zbrate = 1.0 / tres;
   double murate = 0.01; // Hz
   mc.d_flat = murate * tres;

   mc.d_size     = 1.0;

   cout << "INFO: timing resolution:  " << tres << " seconds.\n";
   cout << "INFO: zero bias rate:  :  " << zbrate << " Hz.\n";
   cout << "INFO: assuming background muon rate is : " << murate << " Hz.\n";
   cout << "INFO: flat rate is :  " << mc.d_flat << "\n";

   if (mode_dofit){
     cout << "INFO: will perform fit!\n";
   }

   mc.rng.SetSeed(2014);

   if (mode == 1){
     TFile fout("root/noise_control.root", "RECREATE");
     fout.cd();
     nhits_study(mc, nruns, "noise_control");     
     return 0;
   }

   if (mode == 2){
     TFile fout("root/eshower_control.root", "RECREATE");
     fout.cd();     
     shower_rate_study(mc, nruns, 1E10, 1E20, "eshower_control");     
     return 0;
   }
   
   sprintf(name,"root/rate_%s.root", tag);
   cout << "INFO: writing to file " << name << "\n";
   TFile fout(name, "RECREATE");
   fout.cd();
   tree = new TTree("rate","");

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
   tree->Branch("mcore",&mcore);
   tree->Branch("ereco",&ereco);
   tree->Branch("ureco",&ureco);
   tree->Branch("eshower",&eshower);
   tree->Branch("sigf",&sigf);
   tree->Branch("llr",&llr);
   tree->Branch("sample",&sample);
   tree->Branch("istat",&istat);
   mc.print();

   sample = 0;

   if (mode == 4){  
     mc.mode_manual_locations=1;
     mode_colocate_phones=1;
   }

   noise_study(mc, 2*nruns,   0,    0);
   for (int i=0; i<ebenchmark.size(); i++){
     double E = ebenchmark[i];
     double f = log(1E21/E) / log(10.0);
     if (f< 1.0) f=1.0;
     sample++;
     cout << "INFO: sample " << sample << " E:  " << E << " stat factor:  " << f << "\n";
     noise_study(mc, f*nruns, E, 1E21);
   }
   tree->Write();
   fout.Close();
   return 0;
}
