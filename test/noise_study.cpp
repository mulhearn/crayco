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

TH1F * hereco;

void noise_study(shower_mc & mc, int nruns, double E, const char * htag){
   cout << "Energy of Generated Shower:  " << E << "\n";
   hereco      = new TH1F(htag,"",80, 9.0, 13.0);
   shower_fcn & sfcn = shower_fcn::instance();
   double tot_wgt = 0;
   int i=0;
   while (tot_wgt < nruns){
      if (i%100 == 0) cout << i << "\n";
      i++;

      if (E > 0.0) mc.generate(E, QUIET, 0.0);
      else         mc.weighted_noise(QUIET);
      tot_wgt += mc.wgt;

      //cout << "actual generated s_loge:  " << sfcn.gen_s_loge << "\n";
      sfcn.gen_s_loge = 10.0;
      //sfcn.d_flat     = 0.0; 
      if (shower_fit(QUIET)){
         //cout << "fitted s_loge:  " << sfcn.fit_s_loge << "\n";
         //cout << "uncertainty s_loge:  " << sfcn.unc_s_loge << "\n";
         double e_fit = exp(sfcn.fit_s_loge);
         double sig = sfcn.fit_s_loge / sfcn.unc_s_loge;

         if (sig > 100.0) {
	   //cout << "logn significance:  " << sig << "\n";
	   //cout << "fitted energy:  " << e_fit << "\n";
            hereco->Fill(log(e_fit) / log(10.0));
         }
      }
   }
}

int main(){   
   shower_fcn & sfcn = shower_fcn::instance();
   shower_mc mc;

   int nruns = 100;  // should take from command line...

   mc.rng.SetSeed(2014);

   mc.d_flat     = 0.0;
   mc.d_size     = 1.0;
   mc.d_n = 1000;
   mc.d_xs_gamma       = 1E-7; 
   mc.d_xs_mu          = 0.0; 
   mc.print();   

   TFile fout("noise.root", "RECREATE");
   fout.cd();

   mc.d_flat     = 0.02;
   noise_study(mc, nruns, 0, "h_02_noise");
   hereco->Write();
   noise_study(mc, nruns, 1E11, "h_02_signal");
   hereco->Write();

   mc.d_flat     = 0.2;
   noise_study(mc, nruns, 0, "h_2_noise");
   hereco->Write();
   noise_study(mc, nruns, 1E11, "h_2_signal");
   hereco->Write();

   fout.Close();

   return 0;
}
