#ifndef __SHOWER_MODEL_H__
#define __SHOWER_MODEL_H__

// All Energy is in units of eV

// log of number of total particles of gamma/mu given log of primary energy
double logn_gamma (double logE);
double logn_mu(double logE);

// log of primary energy given log of number of particles of a particular species
double logE_gamma (double logn);
double logE_mu (double logn);

// PDF for shower along ground at position x and y, given sin2theta and phi of shower.
// CHECK:  What value of age is assumed?
double shower_pdf(double x, double y, double s_sin2theta, double s_phi);

// Flux of showers in number / (km2 yr sr eV) 
// primary energy e in **eV**
double flux(double e);

#endif
