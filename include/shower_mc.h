#ifndef __SHOWER_MC_H__
#define __SHOWER_MC_H__

class shower_mc {
 public:
   TRandom rng;
   double d_n;     // number of detectors per square kilometer
   double d_size;  // array size, in kilometers
   double d_xs_gamma; // cross-section times gamma efficiency of detector
   double d_xs_mu;    // cross-section times muon efficiency of detector
   double d_flat;  // flat background rate (fakes)

   double tot_n() { return d_n * d_size * d_size; }

   // fill shower_fcn instance with MC data for the phone array
   void generate(double s_loge, double s_theta, double s_phi, int verbose=VERBOSE);    
   
   // throws a random phi, random theta, and generates MC data.
   // if fix_theta >= 0 then generates MC at that fixed theta value
   void generate(double E, int verbose=VERBOSE, double fix_theta=-1.0);

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

   // print
   void print();

 private:
   void init_shower_fcn();

};

#endif
