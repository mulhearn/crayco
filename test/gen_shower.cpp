#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <math.h>
#include "TFile.h"
#include "TH1.h"

const float max_theta = 1.22; // 70 degrees
//const float energy = 1.0e17;
const float rM = 79.0;// m (EAS page 199).

float H = 10e3; // 100 km up

//  log (E) = a + b log (N)
//  log(e) - a = b log N
//  (log(E)-a)/b = log N
// exp ( log (E) - a / b ) = N
// a = 1.28, b = 0.92 (proton, http://particle.astro.ru.nl/pub/icrc09-1044.pdf)

float nch (float E)
{
  float a = 1.28;
  float b = 0.92;
  float N = exp( (log(E) - a) / b );
  return N;
}


// from EAS, page 199, eq 4.162
float ldf(float N_e, float r, float s)
{
  float norm = (N_e/(2*3.14159*rM*rM));
  float r1 = pow(r/rM,s-2);
  float r2 = pow( 1+ (r/rM),s-4.5);
  float gamma = tgammaf(4.5-s)/(tgammaf(s)*tgammaf(4.5-(2*s)));

  return norm
    *r1
    *r2
    *gamma;

}

void read_det(char fname[1000],int &num,
	      std::vector<float> &x, std::vector<float> &y)
{
  std::ifstream fin(fname);
  char trash[100];
  fin >> trash;
  fin >> num;
  fin >> trash;
  int size;
  fin >> size;
  x.clear();
  y.clear();
  for (int i=0;i<num;i++)
    {
      int index; float xv,yv;
      fin >> index >> xv >> yv;
      x.push_back(xv);
      y.push_back(yv);
    }
}

int main(int argc, char *argv[])
{
  if (argc<3) { std::cout << "Usage: gen_shower det_file energy eff tag" << std::endl; exit(0);}

  char dfile[1000]; strcpy(dfile,argv[1]);
  float energy =atof(argv[2]);
  float eff = atof(argv[3]);
  char tag[1000]; strcpy(tag,argv[4]);

  float N = nch(energy);

  std::vector<float> x,y;
  int dsize;
  x.reserve(1000);
  y.reserve(1000);
  read_det(dfile,dsize,x,y);

  int i=0;

  //  TFile *tfo = new TFile(Form("phi_templates_%s.root",tag),"RECREATE");

  for (float theta = 0; theta <= max_theta; theta += max_theta/5.0)
    for (float age = 1.0;age <= 1.2; age += 0.2)
      for (float phi = 0; phi <= 3.14159/2; phi += 3.14159/8.0)
	{
	  for (int ix=0;ix<5;ix++)
	    {
	  
	  //	float theta = max_theta*rand()/(RAND_MAX*1.0);
	  //	float phi = 2*3.14159*rand()/(RAND_MAX*1.0);
	  //      float age = 0.4 + 1.2*(rand()/(RAND_MAX*1.0));

	  std::cout << " Theta = " << theta << " N = " << N << " phi = " << phi << " age = " << age << "\t";
	  
	  char fout_name[1000];
	  sprintf(fout_name,"shower%d_theta%1.2f_phi%1.1f_age%1.1f_%s_eff%s_energy%s_%s",i,theta,phi,age,tag,argv[3],argv[2],dfile);
	  std::ofstream fout(fout_name);
	  
	  int iout=0;

	  TH1F *h = new TH1F(Form("h_%d",i),Form("h_theta%1.2f_phi%1.1f",theta,phi),100,0,3.14159);
	  h->Sumw2();
	  for (int is=0;is<dsize;is++)
	    {
	      float R = sqrt( x[is]*x[is] + y[is]*y[is] );
	      //	  float r = R * cos(theta);
	      
	      float det_phi = std::atan2(y[is],x[is]);
	      
	      float r = 0.5*R * sqrt(3.0-cos(2*(det_phi-phi)) + 2*cos(det_phi-phi)*cos(det_phi-phi)*cos(2*theta));
	      
	      
	      float Ht = H*tan(theta);
	      float D = fabs(Ht - R);
	      float l = sqrt(H*H+D*D);
	      float time = l/(3.0e8);
	      
	      float density = ldf( N, r, age );

	      float size = 0.5e-4; // m^2;
	      float prob = size*density*eff;
	      float trand = rand()/(1.0*RAND_MAX); // between 0-1
	      if (trand < prob) 
		{
		  fout << x[is] << "\t" << y[is] << "\t" << density << "\t" << time <<std::endl;
		  h->Fill(fabs(det_phi));
		  //std::cout << " detector phi = " << det_phi << std::endl;
		  iout++;
		}
	    }
	  fout.close();
	  std::cout << " done with shower " << i << " " << fout_name << " hits = " << iout << std::endl;
	  i++;
	    }
	}

}

