#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include "global.h"
#include "fermi_bonn.h"

#include <iomanip>
#include <string>
#include <math.h>
#include <iostream>
void fermi_bonn() {



Double_t x[] = {
0.,
0.0522648,    
0.11324,      
0.156794,     
0.200348,     
0.243902,     
0.278746,     
0.3223,       
0.357143,     
0.400697,     
0.43554,      
0.479094,     
0.522648,     
0.574913,     
0.618467,     
0.670732,     
0.722997,     
0.766551,     
0.818815,     
0.879791,     
0.932056,     
0.984321,     
1.0453,       
1.10627,      
1.15854,      
1.22822,      
1.2892,       
1.35017,      
1.41115,      
1.48084,      
1.55052,      
1.62021,      
1.6899,       
1.76829,      
1.84669,      
1.92509,      
2.00348,      
2.08188,      
2.16028,      
2.23868,      
2.32578,      
2.40418,      
2.48258,      
2.56969,      
2.64808,      
2.73519,      
2.81359,      
2.89199,      
2.97909,      
3.05749,      
3.1446,       
3.223,        
3.30139,      
3.3885,       
3.4669,       
3.5453,       
3.62369,      
3.70209,      
3.77178,      
3.85889,      
3.91115,      
3.98955,      
4.04181,      
4.12021,      
4.17247,      
4.25087,      
4.32056,      
4.39895,      
4.44251,      
4.52091,      
4.5993,       
4.68641,      
4.76481,      
4.84321,      
4.9216		
	};

Double_t y[] = {
    0.,
    10.6082,
  8.37678,
   5.87802,
   4.12463,
   3.01046,
   2.11245,
 1.48231,
   1.04015,
   0.729874,
  0.512155,
   0.359381,
   0.262303,
   0.184059,
   0.13434,
   0.0980512,
   0.071565,
   0.0502175,
   0.0366524,
   0.0278256,
   0.0203092,
   0.0148231,
 0.0112534,
  0.00821353,
  0.00623551,
  0.00473384,
 0.00359381,
  0.00272833,
  0.00215443,
  0.00170125,
  0.0013434,
  0.00106082,
 0.000837678,
  0.000688029,
  0.000587802,
  0.000482793,
  0.000412463,
  0.000352378,
  0.000301046,
  0.000257191,
  0.000228546,
  0.000203092,
  0.000180472,
  0.000160372,
  0.00013701,
  0.00012175,
  0.00010819,
  9.61404e-05,
  8.54327e-05,
  7.59175e-05,
 6.48583e-05,
5.76347e-05,
  5.12155e-05,
 4.37548e-05,
 4.04425e-05,
 3.45511e-05,
  3.07029e-05,
  2.72833e-05,
  2.33089e-05,
  2.15443e-05,
  1.84059e-05,
  1.70125e-05,
  1.45343e-05,
  1.3434e-05,
  1.1477e-05,
  1.06082e-05,
  9.42668e-06,
  8.05347e-06,
  7.4438e-06,
  6.35944e-06,
 5.65115e-06,
  4.82793e-06,
  4.12463e-06,
  3.52378e-06,
 3.01046e-06
};
Int_t i;

for (i=0; i<75; i++) {
y[i] = y[i]*x[i]*x[i]*x[i];
x[i] = x[i]*197.;
};



ROOT::Math::Interpolator inter(74, ROOT::Math::Interpolation::kCSPLINE);



inter.SetData(74, x, y);

TH1* h_interpol = new TH1D("Paris interpol","Paris interpol",1000,0.,940.);

for (i=1; i<1000; i++) {
h_interpol->SetBinContent(i,inter.Eval(h_interpol->GetBinCenter(i)));
};



Double_t PI = 3.1415;
Double_t px,py,pz,theta,phi;
Double_t momentum; 

Double_t q[1];
Double_t probSum[1];




probSum[0]= ((float) rand() / (RAND_MAX));

h_interpol->GetQuantiles(1,q,probSum);
momentum =q[0]/1000.;

theta = PI*((double) rand() / (RAND_MAX));
phi = 2.*PI*((double) rand() / (RAND_MAX));
px_fermi = momentum*sin(theta)*cos(phi);
py_fermi = momentum*sin(theta)*sin(phi);
pz_fermi = momentum*cos(theta);



h_interpol->Delete();



};











