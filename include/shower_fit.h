#ifndef __SHOWER_FIT_H__
#define __SHOWER_FIT_H__

#include <vector>

// shower_fit:  perform the shower fit to reconstruct shower parameters
// At the moment, you must set data into shower_fcn first...
enum {QUIET=0, VERBOSE=1};
int shower_fit(int verbose = VERBOSE);

// The shower poisson fit function, used with Minuit to fit shower parameters
class shower_fcn {
public:
  // needed for TMinuit fcn interface
  static shower_fcn & instance();
  double count_hits(int nmin=1);
  double count_hits_core(double radius);
  double calc_fcn(double * par);  // -2 log likelihood for shower with parameters <par>  
  double noise_only_fcn(); // -2 log likehood for noise only, depends only on d_flat

  // detector data:
  std::vector<double> d_x;
  std::vector<double> d_y;
  std::vector<int>    d_h;

  //fixed parameters (for now)
  double d_xs_gamma;
  double d_xs_mu;
  double d_flat;
  
  // parameters used for generating shower data:
  double gen_s_loge;
  double gen_s_sin2theta;   
  double gen_s_phi;
  double gen_s_x; // presently ignored in fit...
  double gen_s_y; // presently ignored in fit...

  // fitted parameters:
  double fit_s_loge;
  double fit_s_sin2theta;
  double fit_s_phi;
  double unc_s_loge;
  double unc_s_sin2theta;
  double unc_s_phi;

  // Minuit fit status variables
  double fmin;
  double fedm;
  double errdef;
  int npari;
  int nparx;
  int    istat;


private:
  shower_fcn(){}
};

// TMinuit fcn interface:
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);


#endif
