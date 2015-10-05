#ifndef __UTILS_H__
#define __UTILS_H__

// Convenience functions:
inline double sq(double x){ return x*x; }
inline double sqrtabs(double x){ return sqrt(fabs(x)); }

// Basic statistical analysis:
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

// calculate phi2 - phi1 in -pi to pi
double delta_phi(double phi1, double phi2);

// calculate minimum delta phi,
// allowing for a +/- pi ambiguity:
double delta_phi_ambig(double phi1, double phi2);

// given sin2theta1 and sin2theta2, calculate theta2 minus theta1 in -pi/2 to pi/2
double delta_sin2theta(double s2t1, double s2t2);

// calculate a pull:  (x - g) / dx
double pull(double g, double x, double dx);

#endif
