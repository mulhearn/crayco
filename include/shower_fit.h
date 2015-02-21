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
  double count_hits();
  double calc_fcn(double * par);

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
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);


#endif
