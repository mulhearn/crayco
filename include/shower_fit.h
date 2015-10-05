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

  int mode_cheat; // use generator values for starting point
  int mode_fix_theta_phi; // don't fit theta and phi
  int mode_fix_x_y; // don't fit the starting position

  double count_hits(int nmin=1);
  double count_hits_core(double radius);
  double calc_fcn(double * par);  // -2 log likelihood for shower with parameters <par>  
  double noise_only_fcn(); // -2 log likehood for noise only, depends only on d_flat

  // detector data:
  std::vector<double> d_x;
  std::vector<double> d_y;
  std::vector<int>    d_h;

  //fixed parameters:
  double d_xs_gamma;
  double d_xs_mu;
  double d_flat;
  double d_size;
  
  // parameters used for generating shower data:
  double gen_s_loge;	  
  double gen_s_sin2theta;   
  double gen_s_phi;	  
  double gen_s_x; 	  
  double gen_s_y;         

  // parameters used for generating shower data:
  double start_s_loge;
  double start_s_sin2theta;   
  double start_s_phi;
  double start_s_x; 
  double start_s_y; 
  void cheat(); // copy generator values to the fit starting values.
  void estimate(); // estimate or guess the starting values, for grid of size d_size in km

  // fitted parameters:
  double fit_s_loge;
  double fit_s_sin2theta;
  double fit_s_phi;
  double fit_s_x;
  double fit_s_y;
  double unc_s_loge;
  double unc_s_sin2theta;
  double unc_s_phi;
  double unc_s_x;
  double unc_s_y;

  // Minuit fit status variables
  double fmin;
  double fedm;
  double errdef;
  int npari;
  int nparx;
  int    istat;


private:
  shower_fcn();
};

// TMinuit fcn interface:
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);


#endif
