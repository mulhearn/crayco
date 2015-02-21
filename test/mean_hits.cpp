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

   for (int i=0; i<nruns; i++){
      mc.generate(E, QUIET);
      double nhits = sfcn.count_hits();
      //cout << "INFO:  number of hits:  " << nhits << "\n";
      hits.add(nhits);
   }
   cout << "INFO: E:  " << E << " phone hits mean:  " << hits.mean() << " +/- " << hits.unc() << " rms: " << hits.rms() << "\n";
   return hits.mean();
}

int main(int argc, char * argv[]){   
   shower_mc mc;
   mc.rng.SetSeed(0);
   shower_fcn & sfcn = shower_fcn::instance();

   mc.d_flat     = 0.0;
   mc.d_size     = 1.0;
   mc.d_n = 1000;
   mc.d_xs_gamma          = 0.0; 
   mc.d_xs_mu             = 0.0; 

   char fname[200];
   if (argc < 5) { 
     cout << "usage:  mean_hits <nruns> <n> <xs_gamma> <xs_mu>\n";
     cout << "eg: mean_hits 10 1000 1E-9 1E-5\n";
     return 0; 
   }
   int nruns     = atoi(argv[1]);
   mc.d_n        = atoi(argv[2]);
   mc.d_xs_gamma = atof(argv[3]);
   mc.d_xs_mu    = atof(argv[4]);

   mc.print();   
   cout << "Energy scan:\n";
   mean_hits(mc, nruns, 1E7);
   mean_hits(mc, nruns, 1E8);
   mean_hits(mc, nruns, 1E9);
   mean_hits(mc, nruns, 1E10);
   mean_hits(mc, nruns, 1E11);
   mean_hits(mc, nruns, 1E12);
}
   

