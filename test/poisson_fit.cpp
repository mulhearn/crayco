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

#include "shower_dist.h"

using namespace std;

enum {QUIET=0, VERBOSE=1};

class calc_stats{
public:
   calc_stats(){tot_x = tot_x2 = n = max = 0; }
   void add(double x);
   double mean();
   double rms();
   double unc();

   double tot_x;
   double tot_x2;
   double n;
   double max;
};

void calc_stats::add(double x){
   if (x > max) max = x;
   tot_x  += x;
   tot_x2 += x*x;
   n+=1;
}

double calc_stats::mean(){
   double m = tot_x / n;
   // fix annoying double precision errors:
   if (m > max) m = max;
   return m;
}

double calc_stats::rms(){
   double m  = mean();
   double x2 = tot_x2 / ((double) n);
   return sqrt(fabs(x2 - m*m)); 
}

double calc_stats::unc(){
   return rms() / sqrt((double) n);
}


double sq(double x){ return x*x; }

double sqrtabs(double x){ return sqrt(fabs(x)); }


const double rM = 79.0;// m (EAS page 199).
const double sin2theta_max = sq(sin(1.22));

void cannonize_phi(double & phi){
   while(phi < 0.0)     phi += 2*M_PI;
   while(phi >= 2*M_PI) phi -= 2*M_PI;   
}

double delta_phi(double phi1, double phi2){
   double dphi = phi2 - phi1;
   while(dphi < -M_PI) dphi += 2*M_PI;
   while(dphi >= M_PI) dphi -= 2*M_PI;   
   return dphi;
}

double delta_sin2theta(double s2t1, double s2t2){
   double st1 = sqrt(fabs(s2t1));
   double st2 = sqrt(fabs(s2t2));
   if (st1 > 1.0) st1 = 1.0;
   if (st2 > 1.0) st2 = 1.0;
   double s1 = asin(st1);
   double s2 = asin(st2);
   return delta_phi(s1, s2);
}


double logn_gamma (double logE){
   // from Daniel's fit to Corsika...
   double a = -2.36;
   double b = 1.15;
   double logn = a + b*(logE);
   return logn;
}

double logn_mu(double logE){
  // old mu- only:
  //double a = -4.36;
  //double b = 0.93;
  // new mu- and mu+
  double a = -3.57;
  double b = 0.93;
  double logn = a + b*(logE);
  return logn;
}

double logE_gamma (double logn_gamma){
   // invert above:
   double a = -2.36;
   double b = 1.15;
   double logE = (logn_gamma - a) / b;
   return logE;
}

double logE_mu (double logn_gamma){
   // invert above:
   double a = -4.36;
   double b = 0.93;
   double logE = (logn_gamma - a) / b;
   return logE;
}

double shower_pdf(double x, double y, double s_sin2theta, double s_phi){
  // Debug with Gaussian:
  //const double gnorm = 1.0 / sqrt(2 * M_PI);
  //double dx = x/100*s_theta;
  //double dy = y/100*s_phi;
  //return gnorm * ( exp(-0.5*dx*dx)/(100*s_theta) * exp(-0.5*dy*dy)/(100*s_phi));

  //cout << s_theta << "\n";
  //cout << s_phi << "\n";

  // don't understand this yet, just copying from Daniel:
  double R = sqrt( x*x + y*y );
  double det_phi = std::atan2(y,x);      
  //double r = 0.5*R * sqrt(3.0-cos(2*(det_phi-s_phi)) + 2*cos(det_phi-s_phi)*cos(det_phi-s_phi)*cos(2*s_sin2theta));

  double cos2_dphi  = sq(cos(det_phi - s_phi));
  // simplified version of above:
  double r = R * sqrt(1.0 - cos2_dphi * s_sin2theta);

  // from EAS, page 199, eq 4.162
  // with overall factor of N_e removed
  const double s = 1.9;
  double norm = (1.0/(2*3.14159*rM*rM));
  double r1 = pow(r/rM,s-2);   
  double r2 = (double) pow( 1.0 + (r/rM),s-4.5);
  double gamma = tgammaf(4.5-s)/(tgammaf(s)*tgammaf(4.5-(2*s)));
  double pdf = norm*r1*r2*gamma;
  //cout << "r:  " << r << "\n";
  //cout << "pdf:  " << pdf << "\n";
  return pdf;
} 

double pull(double g, double x, double dx){ return (x-g) / dx; }

class shower_fcn {
public:
  // needed for TMinuit fcn interface
  static shower_fcn & instance(){ 
     static shower_fcn x;
     return x;
  }

  double count_hits(){
     double count_hits = 0;
     for (int i=0; i<d_x.size(); i++){
        if (d_h[i] > 0.0){ 
           count_hits+=1.0;
        }  
     }
     return count_hits;
  }

  // the actual fcn using data from the class shower_fcn:
  double calc_fcn(double * par){
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

  // detector data:
  vector<double> d_x;
  vector<double> d_y;
  vector<int>    d_h;

  //fixed parameters (for now)
  double d_xs_gamma;
  double d_xs_mu;
  double d_flat;
  
  // parameters used for generating shower data:
  double gen_s_loge;
  double gen_s_sin2theta;   
  double gen_s_phi;

  // fitted parameters:
  double fit_s_loge;
  double fit_s_sin2theta;
  double fit_s_phi;
  double unc_s_loge;
  double unc_s_sin2theta;
  double unc_s_phi;

private:
  shower_fcn(){}
};

// TMinuit fcn interface:
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  f = shower_fcn::instance().calc_fcn(par);
}


int do_fit(int verbose = VERBOSE){
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
   void noise_study(double E, const char * htag);  // fill histograms for noise study...
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
      if (do_fit(QUIET)){
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
      if (do_fit(VERBOSE)){
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

void mcarray::noise_study(double E, const char * htag){
   cout << "Energy of Generated Shower:  " << E << "\n";

   hereco      = new TH1F(htag,"",80, 9.0, 13.0);
   shower_fcn & sfcn = shower_fcn::instance();
   int n = 1000;
   for (int i=0; i<n; i++){
      if (i%100 == 0) cout << i << "\n";
      generate(E, QUIET, 0.0);
      //cout << "actual generated s_loge:  " << sfcn.gen_s_loge << "\n";
      sfcn.gen_s_loge = 10.0;
      //sfcn.d_flat     = 0.0; 
      if (do_fit(QUIET)){
         //cout << "fitted s_loge:  " << sfcn.fit_s_loge << "\n";
         //cout << "uncertainty s_loge:  " << sfcn.unc_s_loge << "\n";
         double e_fit = exp(sfcn.fit_s_loge);
         double sig = sfcn.fit_s_loge / sfcn.unc_s_loge;

         if (sig > 100.0) {
	   //cout << "logn significance:  " << sig << "\n";
	   //cout << "fitted energy:  " << e_fit << "\n";
            hereco->Fill(log(e_fit) / log(10.0));
         }
      }
   }
}

void mcarray::shower_histograms(){
   shower_fcn & sfcn = shower_fcn::instance();
   double log10_emin = 10.0;
   double log10_emax = 12.0;
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
      double w = w0*flux(E*1E9)*E;
      hshower_true->Fill(log10_e, w);

      generate(E, QUIET, 0.0);
      if (do_fit(QUIET)){
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

   //array.d_n        = 10; 
   //array.d_size     = 1.0;
   //array.d_xs_gamma       = 1E-9; 
   //1E-9 to 1E-7


   //array.emain      = 1E21;
   //array.escale     = 10.0;
   
   if (argc < 3) { return 0; }
   array.d_n = atoi(argv[1]);
   int nruns = atoi(argv[2]);
   cout << "INFO:  command line argument d_n = " << array.d_n << "\n";
   cout << "INFO:  command line argument nrums = " << nruns << "\n";

   array.d_size = 1.0;
   sprintf(fname, "resplots_%i.root", int(array.d_n));
   
   //array.print();
   //array.d_xs_gamma       = 1E-9; 
   //cout << "Energy scan for d_xs_gamma = " << array.d_xs_gamma << "\n";
   //array.mean_hits(1E7);
   //array.mean_hits(1E8);
   //array.mean_hits(1E9);
   //array.mean_hits(1E10);
   //array.mean_hits(1E11);
   //array.mean_hits(1E12);

   //array.d_xs_gamma       = 1E-8; 
   //cout << "Energy scan for d_xs_gamma = " << array.d_xs_gamma << "\n";
   //array.mean_hits(1E7);
   //array.mean_hits(1E8);
   //array.mean_hits(1E9);
   //array.mean_hits(1E10);
   //array.mean_hits(1E11);
   //array.mean_hits(1E12);

   vector<double> _erms7, _trms7, _prms7, _eunc7, _tunc7, _punc7;
   vector<double> _erms8, _trms8, _prms8, _eunc8, _tunc8, _punc8;
   vector<double> _erms9, _trms9, _prms9, _eunc9, _tunc9, _punc9;
   vector<double> _pen, _pen7, _pen8, _pen9; 
   vector<double> _meannhits7, _meannhits8, _meannhits9;
   int runs = 0;
   int runs7 = 0;
   int runs8 = 0;
   int runs9 = 0;
   
   array.d_xs_mu   = 0.0; 
   array.d_xs_gamma  = 1E-7; 
   cout << "Energy scan for d_xs = " << array.d_xs_gamma << "\n";
   cout << "Energy scan for d_xs_mu = " << array.d_xs_mu << "\n";
   for (double e=1E7; e<=1E12; e *= 2.0){
      _pen.push_back(e*1E9);
      runs += 1;
      double mhits = array.mean_hits(e);
      _meannhits7.push_back(mhits);
      if (mhits > 5.0){
         //array.resolutions(nruns, e, _erms7, _eunc7, _trms7, _tunc7, _prms7, _punc7);
         _pen7.push_back(e*1E9);
         runs7 += 1;
      }
   }
   //for (int i=0; i<_erms7.size(); i++){
   //  _eunc7.push_back(_erms7[i] / sqrt(runs7));
   //}
   cout << _erms7.size() << " " <<  _eunc7.size() << " " << _pen7.size() << " " << runs7 << endl;

   array.d_xs_gamma       = 5E-9; 
   array.d_xs_mu          = 5E-5; 
   cout << "Energy scan for d_xs_gamma = " << array.d_xs_gamma << "\n";
   cout << "Energy scan for d_xs_mu = " << array.d_xs_mu << "\n";
   for (double e=1E7; e<=1E12; e *= 2.0){
      double mhits = array.mean_hits(e);
      _meannhits8.push_back(mhits);
      if (mhits > 5.0){
         array.resolutions(nruns, e, _erms8, _eunc8, _trms8, _tunc8, _prms8, _punc8);
         _pen8.push_back(e*1E9);
         runs8 += 1;
      }
   }
   //for (int i=0; i<_erms8.size(); i++){
   //  _eunc8.push_back(_erms8[i] / sqrt(runs8));
   //}
   cout << _erms8.size() << " " <<  _eunc8.size() << " " << _pen8.size() << " " << runs8 << endl;

   array.d_xs_gamma          = 1E-9; 
   array.d_xs_mu             = 1E-5; 
   cout << "Energy scan for d_xs_gamma = " << array.d_xs_gamma << "\n";
   cout << "Energy scan for d_xs_mu = " << array.d_xs_mu << "\n";
   for (double e=1E7; e<=1E12; e *= 2.0){
      double mhits = array.mean_hits(e);
      _meannhits9.push_back(mhits);
      if (mhits > 5.0){
         array.resolutions(nruns, e, _erms9, _eunc9, _trms9, _tunc9, _prms9, _punc9);
         _pen9.push_back(e*1E9);
         runs9 += 1;
      }
   }
   //for (int i=0; i<_erms9.size(); i++){
   //  _eunc9.push_back(_erms9[i] / sqrt(runs9));
   //}
   cout << _erms9.size() << " " <<  _eunc9.size() << " " << _pen9.size() << " " << runs9 << endl;
 
   TGraphErrors *geres_err7 = new TGraphErrors(runs7, &_pen7[0], &_erms7[0], NULL, &_eunc7[0]);
   TGraphErrors *gtres_err7 = new TGraphErrors(runs7, &_pen7[0], &_trms7[0], NULL, &_eunc7[0]);
   TGraphErrors *gpres_err7 = new TGraphErrors(runs7, &_pen7[0], &_prms7[0], NULL, &_eunc7[0]);
   TGraphErrors *gnhits_pen7 = new TGraphErrors(runs, &_pen[0], &_meannhits7[0], NULL, NULL);
   TGraphErrors *geres_err8 = new TGraphErrors(runs8, &_pen8[0], &_erms8[0], NULL, &_eunc8[0]);
   TGraphErrors *gtres_err8 = new TGraphErrors(runs8, &_pen8[0], &_trms8[0], NULL, &_eunc8[0]);
   TGraphErrors *gpres_err8 = new TGraphErrors(runs8, &_pen8[0], &_prms8[0], NULL, &_eunc8[0]);
   TGraphErrors *gnhits_pen8 = new TGraphErrors(runs, &_pen[0], &_meannhits8[0], NULL, NULL);
   TGraphErrors *geres_err9 = new TGraphErrors(runs9, &_pen9[0], &_erms9[0], NULL, &_eunc9[0]);
   TGraphErrors *gtres_err9 = new TGraphErrors(runs9, &_pen9[0], &_trms9[0], NULL, &_eunc9[0]);
   TGraphErrors *gpres_err9 = new TGraphErrors(runs9, &_pen9[0], &_prms9[0], NULL, &_eunc9[0]);
   TGraphErrors *gnhits_pen9 = new TGraphErrors(runs, &_pen[0], &_meannhits9[0], NULL, NULL);
   
   TFile f(fname, "RECREATE");
   f.cd();
   geres_err7->SetName("eres7");
   gtres_err7->SetName("tres7");
   gpres_err7->SetName("pres7");
   gnhits_pen7->SetName("meanNhits7");
   geres_err7->Write();
   gtres_err7->Write();
   gpres_err7->Write();
   gnhits_pen7->Write();
   geres_err8->SetName("eres8");
   gtres_err8->SetName("tres8");
   gpres_err8->SetName("pres8");
   gnhits_pen8->SetName("meanNhits8");
   geres_err8->Write();
   gtres_err8->Write();
   gpres_err8->Write();
   gnhits_pen8->Write();
   geres_err9->SetName("eres9");
   gtres_err9->SetName("tres9");
   gpres_err9->SetName("pres9");
   gnhits_pen9->SetName("meanNhits9");
   geres_err9->Write();
   gtres_err9->Write();
   gpres_err9->Write();
   gnhits_pen9->Write();
   f.Close();


   //generate_pulls(rng);
   return 0;
}


int main_v2(){   
   int nruns = 20;
   mcarray array;
   array.rng.SetSeed(2014);
   shower_fcn & sfcn = shower_fcn::instance();

   array.d_flat     = 0.0;
   array.d_size     = 1.0;
   array.d_n = 1000;
   array.d_xs_gamma       = 1E-7; 
   array.d_xs_mu          = 0.0; 
   array.print();   

   TFile fout("noise.root", "RECREATE");
   fout.cd();

   array.d_flat     = 0.01;
   array.noise_study(1.0, "h_01_noise");
   array.hereco->Write();
   array.noise_study(1E12, "h_01_signal");
   array.hereco->Write();

   array.d_flat     = 0.05;
   array.noise_study(1.0, "h_05_noise");
   array.hereco->Write();
   array.noise_study(1E12, "h_05_signal");
   array.hereco->Write();

   array.d_flat     = 0.1;
   array.noise_study(1.0, "h_1_noise");
   array.hereco->Write();
   array.noise_study(1E12, "h_1_signal");
   array.hereco->Write();


   fout.Close();

   return 0;

   array.d_xs_gamma       = 1E-7; 
   cout << "Energy scan for d_xs_gamma = " << array.d_xs_gamma << "\n";
   for (double e=1E7; e<=1E12; e *= 2.0){
      double mhits = array.mean_hits(e);
      if (mhits > 10.0){
         //array.resolutions(nruns, e);
      }
   }

   array.d_xs_gamma       = 1E-8; 
   cout << "Energy scan for d_xs_gamma = " << array.d_xs_gamma << "\n";
   for (double e=1E7; e<=1E12; e *= 2.0){
      double mhits = array.mean_hits(e);
      if (mhits > 10.0){
         //array.resolutions(nruns, e);
      }
   }

   array.d_xs_gamma       = 1E-9; 
   cout << "Energy scan for d_xs_gamma = " << array.d_xs_gamma << "\n";
   for (double e=1E7; e<=1E12; e *= 2.0){
      double mhits = array.mean_hits(e);
      if (mhits > 10.0){
         //array.resolutions(nruns, e);
      }
   }




   return 0;
}
