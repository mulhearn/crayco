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

class mcarray {
public:
   TRandom rng;
   double d_n;     // number of detectors per square kilometer
   double d_size;  // array size, in kilometers
   double d_xs_gamma; // cross-section times gamma efficiency of detector
   double d_xs_mu;    // cross-section times muon efficiency of detector
   double d_flat;  // flat background rate (fakes)

   double tot_n() { return d_n * d_size * d_size; }

   // generate MC data for the phone array:
   void generate(double s_loge, double s_theta, double s_phi, int verbose=VERBOSE);    
   // throws a random phi, random theta if fix theta < 0.theta, phi
   void generate(double E, int verbose=VERBOSE, double fix_theta=-1.0);
  
   // print
   void print();

   // specific calculations:
   double mean_hits(double E);
   void resolutions(int nruns, double E, vector<double> &erms, vector<double> &eunc, vector<double> &trms, 
                    vector<double> &tunc, vector<double> &prms, vector<double> &punc);

   // histogram filling (leaks memory as currently written...)
   void fill_histograms(double E);
   TH1F * hereco;
   TH1F * hfereco;

   void shower_histograms();
   TH1F * hshower_true;
   TH1F * hshower_reco;
   TH1F * hshower_near;
};

void mcarray::print(){
   cout << "INFO: phones / km^2         " << d_n << "\n";
   cout << "INFO: size of grid (km):    " << d_size << "\n";
   cout << "INFO: total phones:         " << tot_n() << "\n";
   cout << "INFO: xsec times gamma eff: " << d_xs_gamma << "\n";
   cout << "INFO: xsec times muon eff:  " << d_xs_mu << "\n";
   cout << "INFO: flat background:      " << d_flat << "\n";
}

void mcarray::generate(double s_loge, double s_sin2theta, double s_phi, int verbose){ 
   const double norm = 1.0 / sqrt(2 * M_PI);
   shower_fcn & sfcn = shower_fcn::instance();

   sfcn.d_x.clear();
   sfcn.d_y.clear();
   sfcn.d_h.clear();
   sfcn.d_xs_gamma       = d_xs_gamma;
   sfcn.d_xs_mu          = d_xs_mu;
   sfcn.d_flat           = d_flat; 

   sfcn.gen_s_loge       = s_loge;
   sfcn.gen_s_sin2theta  = s_sin2theta;
   sfcn.gen_s_phi        = fabs(s_phi);

   double s_logn_mu    = logn_mu(s_loge);
   double s_logn_gamma = logn_gamma(s_loge);

   if (verbose){
      cout << "INFO:  s_loge:        " << s_loge << "\n";
      cout << "INFO:  s_logn_mu:     " << s_logn_mu << "\n";
      cout << "INFO:  s_logn_gamma:  " << s_logn_gamma << "\n";
      cout << "INFO:  d_n:     " << d_n << "\n";
      cout << "INFO:  d_xs_gamma:    " << d_xs_gamma << "\n";
      cout << "INFO:  d_xs_mu:    " << d_xs_mu << "\n";
      cout << "INFO:  d_size:  " << d_size << "\n";
      cout << "INFO:  d_flat:  " << d_flat << "\n";
   }



   int tot_n = d_n * d_size * d_size;
   int nhit = 0;
   for (int i=0; i<tot_n; i++){
      double x = (-500.0 + rng.Uniform() * 1000) * d_size;
      double y = (-500.0 + rng.Uniform() * 1000) * d_size;
      double p = shower_pdf(x, y, s_sin2theta, s_phi);
      double mu = p * (exp(s_logn_gamma) * d_xs_gamma + exp(s_logn_mu) * d_xs_mu)+ d_flat;
      double p0 = ROOT::Math::poisson_cdf(0,mu);
      
      double hit = 0;
      double xp = rng.Uniform();
      if (xp > p0) hit = 1;

      if (hit > 0.0) nhit+=1;
      sfcn.d_x.push_back(x);
      sfcn.d_y.push_back(y);
      sfcn.d_h.push_back(hit);
   }
   if (verbose){
      cout << "INFO:  generated hits:     " << nhit << "\n";
   }

}

void mcarray::generate(double E, int verbose, double fix_sin2theta){
  const double sin2theta_max = sq(sin(1.22));
  shower_fcn & sfcn = shower_fcn::instance();
  double s_loge       = log(E);
  
  //double s_sin2theta  = 0.75;   
  double s_sin2theta  = fix_sin2theta;
  if (s_sin2theta < 0.0)
    s_sin2theta = rng.Uniform() * sin2theta_max;
  double s_phi        = rng.Uniform() * 2 * M_PI;
  
  generate(s_loge, s_sin2theta, s_phi, verbose);
  
  if (verbose) cout << "INFO:  number of hits:  " << sfcn.count_hits() << "\n";
}

double mcarray::mean_hits(double E){
   calc_stats hits;
   shower_fcn & sfcn = shower_fcn::instance();
   int n = 10;
   for (int i=0; i<n; i++){
      generate(E, QUIET);
      double nhits = sfcn.count_hits();
      //cout << "INFO:  number of hits:  " << nhits << "\n";
      hits.add(nhits);
   }
   cout << "INFO: E:  " << E << " phone hits mean:  " << hits.mean() << " +/- " << hits.unc() << " rms: " << hits.rms() << "\n";
   return hits.mean();
}


void mcarray::resolutions(int nruns, double E, vector<double> &erms, vector<double> &eunc, vector<double> &trms, 
vector<double> &tunc, vector<double> &prms, vector<double> &punc){
   calc_stats eres;
   calc_stats pres;
   calc_stats tres;

   calc_stats epull;
   calc_stats ppull;
   calc_stats tpull;

   shower_fcn & sfcn = shower_fcn::instance();
   double cnt = 0.0;
   for (int i=0; i<nruns; i++){
      generate(E, QUIET, 0.75);
      if (shower_fit(QUIET)){
         double e_fit = exp(sfcn.fit_s_loge);
         double e_gen = exp(sfcn.gen_s_loge);
         double e_fres = (e_fit - e_gen) / e_gen;
         double pull_loge  = pull(sfcn.gen_s_loge, sfcn.fit_s_loge, sfcn.unc_s_loge);
         double pull_sin2theta = pull(sfcn.gen_s_sin2theta, sfcn.fit_s_sin2theta, sfcn.unc_s_sin2theta);
         double pull_phi   = pull(sfcn.gen_s_phi, sfcn.fit_s_phi, sfcn.unc_s_phi);
         double de = (sfcn.fit_s_loge - sfcn.gen_s_loge);
         double dp = delta_phi(sfcn.fit_s_phi,sfcn.gen_s_phi);
         double dt = delta_sin2theta(sfcn.fit_s_sin2theta, sfcn.gen_s_sin2theta);
         //cout << "pull logn:        " << pull_logn << "\n";
         //cout << "pull sin2theta:   " << pull_sin2theta << "\n";
         //cout << "pull phi:         " << pull_phi << "\n";
         //cout << "delta logn:       " << dn << "\n";
         //cout << "delta theta:      " << dt << "\n";
         //cout << "delta phi:        " << dp << "\n";

         //cout << "e_fit:  " << e_fit << " e_gen:  " << e_gen << "\n";
         //cout << "fractional energy resolution:  " << e_fres << "\n";
 
         if (fabs(e_fres) < 5.0){
               eres.add(e_fres);
               cnt += 1.0;
               pres.add(dp);
               tres.add(dt);
               epull.add(pull_loge);
               ppull.add(pull_phi);
               tpull.add(pull_sin2theta);
         }
      }
   }
   if (eres.n >= 1){
      cout << "INFO: E:  " << E << " frac energy resolution:  " << eres.rms() << " bias:  " << eres.mean() << "\n";
      cout << "INFO: E:  " << E << " theta resolution:        " << tres.rms() << " bias:  " << tres.mean() << "\n";
      cout << "INFO: E:  " << E << " phi resolution:          " << pres.rms() << " bias:  " << pres.mean() << "\n";
      //cout << "INFO: E:  " << E << " theta pull:  " << tpull.rms() << " bias:  " << tpull.mean() << "\n";
      //cout << "INFO: E:  " << E << " phi pull:  " << ppull.rms() << " bias:  " << ppull.mean() << "\n";
      erms.push_back(eres.rms());
      eunc.push_back(eres.rms() / sqrt(cnt));
      trms.push_back(tres.rms());
      tunc.push_back(tres.rms() / sqrt(cnt));
      prms.push_back(pres.rms());
      punc.push_back(pres.rms() / sqrt(cnt));
   }
}

void mcarray::fill_histograms(double E){
   hereco      = new TH1F("hereco","",100, -30, 30.0);
   hfereco     = new TH1F("hfereco","",100, 0, 10.0);

   shower_fcn & sfcn = shower_fcn::instance();
   int n = 10;
   for (int i=0; i<n; i++){
      if (i%100 == 0) cout << i << "\n";
      generate(E, VERBOSE, 0.0);
      sfcn.gen_s_loge = 0.0;
      if (shower_fit(VERBOSE)){
         double e_fit = exp(sfcn.fit_s_loge);
         double e_gen = exp(sfcn.gen_s_loge);
         double e_frac = e_fit / e_gen; 
         //cout << "e_frac:  " << e_frac << "\n";
         if (e_frac > 9.0) e_frac = 9.0;  // make sure overflows in visible range
         hereco->Fill(log(e_fit) / log(10.0));
         hfereco->Fill(e_frac);
      }
   }
}
void mcarray::shower_histograms(){
   shower_fcn & sfcn = shower_fcn::instance();
   double log10_emin = 19.0;
   double log10_emax = 21.0;
   int nbins = 50;

   hshower_true = new TH1F("hshower_true", "", nbins, log10_emin-1, log10_emax+1); 
   hshower_true->Sumw2();
   hshower_reco = new TH1F("hshower_reco", "", nbins, log10_emin-1, log10_emax+1); 
   hshower_reco->Sumw2();
   hshower_near = new TH1F("hshower_near", "", nbins, log10_emin-1, log10_emax+1); 
   hshower_near->Sumw2();

   int n = 20;
   double w0 = nbins / (double) n;
   for (int i=0; i<n; i++){
      if (i%10 == 0) cout << i << "\n";
      double log10_e = log10_emin + (log10_emax - log10_emin)* rng.Uniform();
      double E = exp(log(10)*log10_e);
      double w = w0*flux(E)*E;
      hshower_true->Fill(log10_e, w);

      generate(E, QUIET, 0.0);
      if (shower_fit(QUIET)){
         double sig = sfcn.fit_s_loge / sfcn.unc_s_loge;
         double log10_efit = sfcn.fit_s_loge / log(10.0);
      
         if (sig > 10.0 && (log10_efit < 11.7)) {
            hshower_reco->Fill(log10_efit, w);            
            if (log10_efit > log10_e + log(1.1) / log(10.0)){
               hshower_near->Fill(log10_efit, w);
            }
         }
      }
   }
}


int main(int argc, char * argv[]){
   mcarray array;
   array.rng.SetSeed(0);
   shower_fcn & sfcn = shower_fcn::instance();

   char fname[200];

   if (argc < 4) { 
     cout << "USAGE:  <tag> d_n nruns \n";
     return 0;
   }

   const char * tag = argv[1];
   array.d_n = atoi(argv[2]);
   int nruns = atoi(argv[3]);
   cout << "INFO:  root file output with tag:  " << tag << "\n";
   cout << "INFO:  command line argument d_n = " << array.d_n << "\n";
   cout << "INFO:  command line argument nrums = " << nruns << "\n";


   array.d_size = 1.0;

   
   vector<double> _erms, _trms, _prms, _eunc, _tunc, _punc, _pen, _meannhits;
   int runs = 0;

   array.d_xs_gamma          = 0; 
   array.d_xs_mu             = 1E-5; 

   sprintf(fname, "resolution_%s.root", tag);
   cout << "Energy scan for d_xs_gamma = " << array.d_xs_gamma << "\n";
   cout << "Energy scan for d_xs_mu = " << array.d_xs_mu << "\n";
   for (double e=1E7; e<=1E12; e *= 2.0){
      double mhits = array.mean_hits(e);
      _meannhits.push_back(mhits);
      if (mhits > 5.0){
         array.resolutions(nruns, e, _erms, _eunc, _trms, _tunc, _prms, _punc);
         _pen.push_back(e*1E9);
         runs += 1;
      }
   }
   cout << _erms.size() << " " <<  _eunc.size() << " " << _pen.size() << " " << runs << endl;
 
   TGraphErrors *geres_err = new TGraphErrors(runs, &_pen[0], &_erms[0], NULL, &_eunc[0]);
   TGraphErrors *gtres_err = new TGraphErrors(runs, &_pen[0], &_trms[0], NULL, &_tunc[0]);
   TGraphErrors *gpres_err = new TGraphErrors(runs, &_pen[0], &_prms[0], NULL, &_punc[0]);
   TGraphErrors *gnhits_pen = new TGraphErrors(runs, &_pen[0], &_meannhits[0], NULL, NULL);
   
   TFile f(fname, "RECREATE");
   f.cd();
   geres_err->SetName("eres");
   gtres_err->SetName("tres");
   gpres_err->SetName("pres");
   gnhits_pen->SetName("meanNhits");
   geres_err->Write();
   gtres_err->Write();
   gpres_err->Write();
   gnhits_pen->Write();
   f.Close();

   return 0;
}


