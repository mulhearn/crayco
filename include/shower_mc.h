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


   //
   // generate a noise event, flat in the number of phones with a hit,
   // with variable wgt set uppropriately.
   double wgt;
   void weighted_noise(int verbose=VERBOSE);

   // print
   void print();
};

#endif
