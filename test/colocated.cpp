#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <math.h>
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TCanvas.h"
#include <iomanip>

#include "TMinuit.h"
#include <Math/ProbFuncMathCore.h>

#include "shower_model.h"
#include "utils.h"
#include "shower_fit.h"
#include "shower_mc.h"

using namespace std;

double mean_hits(shower_mc & mc, int nruns, double E){
   calc_stats hits;
   shower_fcn & sfcn = shower_fcn::instance();

   int npass = 0;

   for (int i=0; i<nruns; i++){
     mc.generate(E, QUIET, 0.0);  // fix theta=0
     //mc.generate(E, QUIET);
      double nhits = sfcn.count_hits(1);
      if (nhits >= 5) npass++;
      //cout << "INFO:  number of hits:  " << nhits << "\n";
      hits.add(nhits);
   }
   double eff = 100.0 * (double) npass / (double) nruns;
   
   //cout << "INFO: E:  " << E << " " << 9.0 + log(E)/log(10.0) << " phone hits mean:  " << hits.mean() << " +/- " << hits.unc() << " rms: " << hits.rms() << "\n";
   cout << "INFO: E:  " << E << " eV  phone hits mean:  " << hits.mean() << " +/- " << hits.unc() << " rms: " << hits.rms() << " eff: " << eff << "\n";
   cout << "INFO: eff/dense:  " << eff / mc.d_n << "\n";





   return hits.mean();
}

int main(int argc, char * argv[]){   
   shower_mc mc;
   mc.rng.SetSeed(0);
   shower_fcn & sfcn = shower_fcn::instance();

   mc.d_flat     = 0.0;
   mc.d_size     = 0.01;
   mc.d_n = 3;
   mc.d_xs_gamma          = 0.0; 
   mc.d_xs_mu             = 0.0; 

   char fname[200];
   if (argc < 6) { 
     cout << "usage:  mean_hits <nruns> <n> <xs_gamma> <xs_mu> <flat>\n";
     cout << "eg: mean_hits 10 1000 0.0 1E-5 0.02\n";
     return 0; 
   }
   int nruns     = atoi(argv[1]);
   mc.d_n        = atoi(argv[2]);
   mc.d_xs_gamma = atof(argv[3]);
   mc.d_xs_mu    = atof(argv[4]);
   mc.d_flat     = atof(argv[5]);

   mc.print();   
   mean_hits(mc, nruns, 1E18);
   mean_hits(mc, nruns, 5E18);
   mean_hits(mc, nruns, 1E19);
   mean_hits(mc, nruns, 5E19);
   mean_hits(mc, nruns, 1E20);
   mean_hits(mc, nruns, 5E20);

}
   

