#ifndef __SHOWER_MC_H__
#define __SHOWER_MC_H__

#define SHOWER_MC_THETA_MAX 1.22  // max allowed theta for shower generation (used to calculate sin2theta_max)

class shower_mc {
 public:

  static const double sin2theta_max;

  shower_mc() {mode_manual_locations = 0; }

  TRandom rng;
  double d_n;     // number of detectors per square kilometer
  double d_size;  // array size, in kilometers
  double d_xs_gamma; // cross-section times gamma efficiency of detector
  double d_xs_mu;    // cross-section times muon efficiency of detector
  double d_flat;  // flat background rate (fakes)
    
  double tot_n() { return d_n * d_size * d_size; }

  void flat_in_circle(double r, double & x, double & y);

  // this is called at start of the the top-level generate routines,
  // but not if mode_manual_locations != 0.
  void generate_locations();
   
  // fill shower_fcn instance with MC data for the phone array
  void generate(double s_loge, double s_theta, double s_phi, int verbose=VERBOSE, double s_x=0.0, double s_y=0.0);    
  
  // throws a random phi, random theta, and generates MC data.
  // if fix_theta >= 0 then generates MC at that fixed theta value
  void generate(double E, int verbose=VERBOSE, double fix_theta=-1.0, double s_x=0.0, double s_y=0.0);    

  // probability for noise alone to produce k hits, set to exactly zero below min probability
  double prob_noise(int k, double min = 1E-50);
  
  // log of the probability for noise alone to produce k hits
  double log_prob_noise(int k);
  
  // generate an event with only noise.
  void noise(int verbose=VERBOSE);
  //
  // generate a noise event, k flat in [kmin,kmax] (the number of noise hits),
  // with member variable wgt set to prob(k) / pmax
  // NOTE:  call optimize_weighted_noise before first call with new parameters!
  double wgt;
  int kmin, kmax;
  double pmax;
  void weighted_noise(int verbose=VERBOSE);
  
  // find values of kmin and kmax that cover probabilities greater than pmin:
  // record highest proability pmax...
  void optimize_weighted_noise(double pmin, int verbose=VERBOSE);
    
  // this generates a shower with energy between such Emin and Emax:
  //  - uniformly distributed in log(E).
  //  - with variable wgt (above) set such that the integral across any energy region will
  //    will give the rate in Hz for a 1 km^2 grid.
  //  - (wgt = dr / dE  where r is the flux in Hz / km^2)  
  //  *** IF YOU USE THIS IN GEOMETRY DIFFERENT THAN 1km x 1km, you must weight it appropriately! ***
  double weighted_shower_energy(double Emin, double Emax, double ngen, int verbose=VERBOSE);

  // print
  void print();

  int mode_manual_locations;
 private:
  void init_shower_fcn();
  
};

#endif
