#include <math.h>

const int np=12;
float fluxd[np];

// eV
float energy[np] = {
   1e15,
   5e15,
   1e16,
   5e16,
   1e17,
   5e17,
   1e18,
   5e18,
   1e19,
   5e19,
   1e20,
    5e20
};

void fill_flux()
{
   //  E^{2.6} * flux, GeV^2.6/(GeV m^2 s sr)
   float Eflux_plot[np] = {
      1e4, // 1e15
      8e3,
      6e3, // 1e16
      3e3,
      2e3,// 1e17
      1e3,
      6e2,// 1e18
      2e2,
      1.5e2, // 1e19
      1e2,
      5e1, // 1e20
      5e0, // 5e20
   };

   for (int i=0;i<np;i++)
   {
      fluxd[i] = Eflux_plot[i] 
         * 1e6  //  m^2 /km^2
         * 3.155e7  //s/year
         / pow((energy[i]*1e-9),2.6); // remove E^2.6
   }

}

bool got_flux=false;

// e in eV
// flux in number / (km2 yr sr GeV)
double flux(double e)
{
   if (!got_flux)
   {
      fill_flux();
      got_flux=true;
   }
  
   if (e > 5E20) return 0.0;

   // look for our index
   int ilow=0;
   while ( energy[ilow+1] < e)
      ilow++;
   


   // linear interpolation of the logs
   double loge_low = log(energy[ilow]);
   double loge_hi = log(energy[ilow+1]);
   double logf_low = log(fluxd[ilow]);
   double logf_hi = log(fluxd[ilow+1]);
  
   double logf_int = logf_low + (log(e)-loge_low)*(logf_hi - logf_low)/(loge_hi - loge_low);
   return exp(logf_int);   
}
