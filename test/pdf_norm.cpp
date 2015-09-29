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

#include <Math/ProbFuncMathCore.h>

#include "shower_model.h"
#include "utils.h"
#include "shower_fit.h"
#include "shower_mc.h"

using namespace std;


int main(int argc, char * argv[]){   
  TRandom rng;

  cout << "INFO: testing...\n";

  cout << shower_pdf(0.1, 0.1, 0, 0) << "\n";

  double sum = 0;
  double count = 0;
  for(int i=0; i<100000; i++){    
    double x = (-500.0 + rng.Uniform() * 100000);
    double y = (-500.0 + rng.Uniform() * 100000);
    double p = shower_pdf(x, y, 0.0, 0.0);
    //cout << "x:   "  << x << "\n";
    //cout << "y:   "  << y << "\n";
    //cout << "p:   "  << p << "\n";
    sum += p;
    count += 1.0;
  }
  cout << "sum:        " << sum << "\n";
  cout << "integral:   " << sum * 100000 * 100000 / count << "\n";



}
   

