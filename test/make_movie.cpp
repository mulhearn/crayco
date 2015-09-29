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
#include <fstream>

#include "TMinuit.h"
#include <Math/ProbFuncMathCore.h>

#include "shower_model.h"
#include "utils.h"
#include "shower_fit.h"
#include "shower_mc.h"

using namespace std;

void save_hits(shower_mc & mc, int nruns, double E){
   calc_stats hits;
   shower_fcn & sfcn = shower_fcn::instance();

   int npass = 0;

   for (int i=0; i<nruns; i++){
     mc.generate(E, QUIET, 0.0);  // fix theta=0
     char name[200];
     sprintf(name,"movie/frame%d.txt",i);
     cout << "INFO: outputting to " << name << "\n";
     ofstream output(name);
     int n = sfcn.d_x.size();
     for (int i=0; i<n; i++){
       output << sfcn.d_x[i] << " " << sfcn.d_y[i] << " " << sfcn.d_h[i] << "\n";
     }
     output.close();

   }
}

int main(int argc, char * argv[]){   
   shower_mc mc;
   mc.rng.SetSeed(0);
   shower_fcn & sfcn = shower_fcn::instance();

   mc.d_flat     = 0.0;
   mc.d_size     = 0.5;
   mc.d_n = 500;
   mc.d_xs_gamma          = 0.0; 
   mc.d_xs_mu             = 0.0; 

   char fname[200];
   if (argc < 6) { 
     cout << "usage:  make_movie <nruns> <n> <xs_gamma> <xs_mu> <flat>\n";
     cout << "eg: make_movie 10 1000 1E-9 1E-5 0.02\n";
     return 0; 
   }
   int nruns     = atoi(argv[1]);
   mc.d_n        = atoi(argv[2]);
   mc.d_xs_gamma = atof(argv[3]);
   mc.d_xs_mu    = atof(argv[4]);
   mc.d_flat     = atof(argv[5]);

   mc.print();   
   //cout << "Energy scan:\n";
   //cout << "ngamma:  " << exp(logn_gamma (log(1.07E11))) << "\n";
   //cout << "nmu:     " << exp(logn_mu (log(1.07E11))) << "\n";

   
   //save_hits(mc, nruns, 1E7);
   //save_hits(mc, nruns, 1E8);
   //save_hits(mc, nruns, 1E9);
   //save_hits(mc, nruns, 1E10);
   //save_hits(mc, nruns, 1E11);

   save_hits(mc, nruns, 1E12);
}
   

