#include <cmath>
#include <math.h>

#include "utils.h"

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

double delta_phi(double phi1, double phi2){
   double dphi = phi2 - phi1;
   while(dphi < -M_PI) dphi += 2*M_PI;
   while(dphi >= M_PI) dphi -= 2*M_PI;   
   return dphi;
}

double delta_phi_ambig(double phi1, double phi2){
  double dphi = delta_phi(phi1,phi2);
  if (dphi < -M_PI/2.0) dphi += M_PI;
  if (dphi >= M_PI/2.0) dphi -= M_PI;   
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

double pull(double g, double x, double dx){ return (x-g) / dx; }
