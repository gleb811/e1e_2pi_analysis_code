#include <TMath.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include "cuts_empty.h"
#include "beta_func_empty.h"
#include "global.h"
#include <iostream>



using namespace std;


 void beta_func_empty() {
 
 Float_t beta_nom_pip,beta_nom_pim,beta_nom_p,delta_t_pip,delta_t_pim,delta_t_p;
 Float_t m_proton, m_pip;
 Float_t p_fid_a_1, p_fid_b_1;
 Float_t pip_fid_a_1,pip_fid_b_1;
 Float_t th_min_1,par1_1,par2_1, pim_fid_a_1,pim_fid_b_1;
 cuts_empty particle_ID_empty;
 
 m_proton =  0.938272;
 m_pip = 0.13957;
   
beta_nom_pip = P_PIp/sqrt(m_pip*m_pip+P_PIp*P_PIp);
beta_nom_pim = P_PIm/sqrt(m_pip*m_pip+P_PIm*P_PIm);
beta_nom_p = P_P/sqrt(m_proton*m_proton+P_P*P_P);
delta_t_pip = PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.;
delta_t_pim = PIm_dist*(1./beta_nom_pim-1/beta_PIm)/30.;
delta_t_p = P_dist*(1./beta_nom_p-1/beta_P)/30.; 
 
 
pip_fid_a_1 = 24.*(1-exp(-1.*0.08*(th_PIp-9.)));
  pip_fid_b_1 = -1.*25.*(1-exp(-1.*0.1*(th_PIp-10.)));
  
if ((ph_PIp >= 330) && (ph_PIp <=360)){
//    if ((ph_PIp > pip_fid_b_1+360) && (ph_PIp < pip_fid_a_1+360)){
    
//time_pip[0][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.); 

// if ((PdHit_PIp==36)&&(delta_t_pip > 0.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.39)*30./PIp_dist);
 
// if ((PdHit_PIp==41)&&(delta_t_pip < -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 3.8)*30./PIp_dist);
//if ((PdHit_PIp==41)&&(delta_t_pip > -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.2)*30./PIp_dist);
 
// if ((PdHit_PIp==42)&&(delta_t_pip < 1.21)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.23)*30./PIp_dist);
if ((PdHit_PIp==42)&&(delta_t_pip > 1.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.)*30./PIp_dist);

if ((PdHit_PIp==45)&&(delta_t_pip < -4.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 4.8)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < 1.)&&(delta_t_pip > -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.45)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < -2.)&&(delta_t_pip > -4.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 3.)*30./PIp_dist);

//if ((PdHit_PIp==46)&&(delta_t_pip > 1.03)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.04)*30./PIp_dist);
//if ((PdHit_PIp==46)&&(delta_t_pip < 1.03)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.02)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip < 1.)&&(delta_t_pip > -1.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.26)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 6.)*30./PIp_dist);

//   };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIp >= 0) && (ph_PIp <=30)){
// if ((ph_PIp > pip_fid_b_1) && (ph_PIp < pip_fid_a_1)){ 

//time_pip[0][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);

if ((PdHit_PIp==42)&&(delta_t_pip > 1.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.)*30./PIp_dist);

if ((PdHit_PIp==45)&&(delta_t_pip < -4.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 4.8)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < 1.)&&(delta_t_pip > -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.45)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < -2.)&&(delta_t_pip > -4.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 3.)*30./PIp_dist);

//if ((PdHit_PIp==46)&&(delta_t_pip > 1.03)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.04)*30./PIp_dist);
//if ((PdHit_PIp==46)&&(delta_t_pip < 1.03)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.02)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip < 1.)&&(delta_t_pip > -1.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.26)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 6.)*30./PIp_dist);



  // };//end of the fiducial cut for second part of sector1
   }; //end of the second part of sector1
  
  
  if ((ph_PIp >= 30) && (ph_PIp <=90)) {
//  if ((ph_PIp > pip_fid_b_1+60) && (ph_PIp < pip_fid_a_1+60)){

//time_pip[1][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);



if ((PdHit_PIp==40)&&(delta_t_pip > 6.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 7.)*30./PIp_dist);
if ((PdHit_PIp==40)&&(delta_t_pip -2.5)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 3.)*30./PIp_dist);

if ((PdHit_PIp==41)&&(delta_t_pip < -3.6)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip +4.1)*30./PIp_dist);
if ((PdHit_PIp==41)&&(delta_t_pip < -2.)&&(delta_t_pip > -3.6)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.1)*30./PIp_dist);




if ((PdHit_PIp==42)&&(delta_t_pip < -4.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.8)*30./PIp_dist);



if ((PdHit_PIp==43)&&(delta_t_pip < -6.)&&(delta_t_pip > -7.7)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 6.95)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip < -7.7)&&(delta_t_pip > -10.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 8.3)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip < -11.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 11.5)*30./PIp_dist);


if ((PdHit_PIp==44)&&(delta_t_pip > -5.5)&&(delta_t_pip < -3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.)*30./PIp_dist);
if ((PdHit_PIp==44)&&(delta_t_pip > -8.)&&(delta_t_pip < -5.5)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 6.7)*30./PIp_dist);

if ((PdHit_PIp==45)&&(delta_t_pip > -2.)&&(delta_t_pip < 2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.8)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip > -7.)&&(delta_t_pip < -4.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 5.2)*30./PIp_dist);

if ((PdHit_PIp==46)&&(delta_t_pip < 1.)&&(delta_t_pip > -3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.9)*30./PIp_dist);
if ((PdHit_PIp==46)&&(delta_t_pip < -7.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 8.)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip < 3.)&&(delta_t_pip > 0.7)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 1.6)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < 0.7)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.12)*30./PIp_dist);


// };//end of the fiducial cut for sector2
};//end of the sector2

  
if ((ph_PIp >= 90) && (ph_PIp <=150)) {
//if ((ph_PIp > pip_fid_b_1+120) && (ph_PIp < pip_fid_a_1+120)){

//time_pip[2][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);

//if (delta_t_pip < 2.055) beta_new_vs_p_pip[2][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1.);


if ((PdHit_PIp==25)&&(delta_t_pip < -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.5)*30./PIp_dist);


if ((PdHit_PIp==40)&&(delta_t_pip > 2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip-4.)*30./PIp_dist);


if ((PdHit_PIp==41)&&(delta_t_pip > 0.37)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.86)*30./PIp_dist);


if ((PdHit_PIp==42)&&(delta_t_pip < 3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 1.9)*30./PIp_dist);
if ((PdHit_PIp==42)&&(delta_t_pip > 3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 5.36)*30./PIp_dist);


if ((PdHit_PIp==44)&&(delta_t_pip < 9.)&&(delta_t_pip > 7.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 8.1)*30./PIp_dist);
if ((PdHit_PIp==44)&&(delta_t_pip > 9.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 10.5)*30./PIp_dist);

if ((PdHit_PIp==46)&&(delta_t_pip < 2.0)&&(delta_t_pip > -2.0)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.5)*30./PIp_dist);
if ((PdHit_PIp==46)&&(delta_t_pip > 2.0)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 3.25)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip < -1.2)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 2.)*30./PIp_dist);


// };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIp >= 150) && (ph_PIp <=210)){
// if ((ph_PIp > pip_fid_b_1+180) && (ph_PIp < pip_fid_a_1+180)){

//time_pip[3][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);


if ((PdHit_PIp==39)&&(delta_t_pip < -3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.3)*30./PIp_dist);


if ((PdHit_PIp==41)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.6)*30./PIp_dist);

if ((PdHit_PIp==42)&&(delta_t_pip < -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip +3.)*30./PIp_dist);


if ((PdHit_PIp==43)&&(delta_t_pip > 0.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.9)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip < 0.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.)*30./PIp_dist);


if ((PdHit_PIp==46)&&(delta_t_pip > -0.85)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.2)*30./PIp_dist);
if ((PdHit_PIp==46)&&(delta_t_pip < -0.85)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.3)*30./PIp_dist);

if ((PdHit_PIp==48)&&(delta_t_pip > -0.1)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.45)*30./PIp_dist);
if ((PdHit_PIp==48)&&(delta_t_pip < -0.1)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 1.1)*30./PIp_dist);

// };//end of the fiducial cut for sector4
  }; //end of the sector4
  
if ((ph_PIp >= 210) && (ph_PIp <=270)) {
//if ((ph_PIp > pip_fid_b_1+240) && (ph_PIp < pip_fid_a_1+240)){

//time_pip[4][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);


if ((PdHit_PIp==40)&&(delta_t_pip > -3.3)&&(delta_t_pip < -0.8)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.6)*30./PIp_dist);
if ((PdHit_PIp==40)&&(delta_t_pip < -3.3)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.85)*30./PIp_dist);


if ((PdHit_PIp==42)&&(delta_t_pip < 0.65)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.33)*30./PIp_dist);
if ((PdHit_PIp==42)&&(delta_t_pip > 0.65)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.0)*30./PIp_dist);


if ((PdHit_PIp==44)&&(delta_t_pip < 0.2)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.6)*30./PIp_dist);
if ((PdHit_PIp==44)&&(delta_t_pip > 0.2)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.73)*30./PIp_dist);

if ((PdHit_PIp==46)&&(delta_t_pip > 1.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 2.)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip < 3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 2.0)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip > 3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 4.0)*30./PIp_dist);


// };//end of the fiducial cut for sector5
}; //end of the sector5
  
  
 if ((ph_PIp >= 270) && (ph_PIp <=330)){
// if ((ph_PIp > pip_fid_b_1+300) && (ph_PIp < pip_fid_a_1+300)){

//time_pip[5][PdHit_PIp-1]->Fill(P_PIp,PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.,1.);

if ((PdHit_PIp==31)&&(delta_t_pip < 0.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.6)*30./PIp_dist);
if ((PdHit_PIp==31)&&(delta_t_pip > 0.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.47)*30./PIp_dist);


if ((PdHit_PIp==40)&&(delta_t_pip < -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.)*30./PIp_dist);

if ((PdHit_PIp==44)&&(delta_t_pip < -1.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.5)*30./PIp_dist);


if ((PdHit_PIp==45)&&(delta_t_pip < -2.8)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.7)*30./PIp_dist);


if ((PdHit_PIp==47)&&(delta_t_pip > -2.)&&(delta_t_pip < 0.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.84)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -2.)&&(delta_t_pip > -4.0)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 2.7)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -4.0)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 5.3)*30./PIp_dist);

if ((PdHit_PIp==48)&&(delta_t_pip > 1.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip -2.78)*30./PIp_dist);


// };//end of the fiducial cut for sector6
  }; //end of the sector6





th_min_1=(11.09+8./(0.472*P_PIm+0.117));
  par1_1=0.705+1.1*P_PIm;
  par2_1=-63.2-29.3*P_PIm;       
   pim_fid_a_1=30.5*pow((sin((th_PIm-th_min_1)*0.01745)),(par1_1+par2_1/th_PIm+1530./th_PIm/th_PIm))-1;   pim_fid_b_1=-30.5*pow((sin((th_PIm-th_min_1)*0.01745)),(par1_1+par2_1/th_PIm+1530./th_PIm/th_PIm))+1;
   









// pi-

if ((ph_PIm >= 330) && (ph_PIm <=360)){
//    if ((ph_PIm > pim_fid_b_1+360) && (ph_PIm < pim_fid_a_1+360)){
    
//time_pim[0][PdHit_PIm-1]->Fill(P_PIm,delta_t_pim,1.); 

// if ((PdHit_PIm==36)&&(delta_t_pim > 0.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.39)*30./PIm_dist);
 
// if ((PdHit_PIm==41)&&(delta_t_pim < -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 3.8)*30./PIm_dist);
//if ((PdHit_PIm==41)&&(delta_t_pim > -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.2)*30./PIm_dist);
 
// if ((PdHit_PIm==42)&&(delta_t_pim < 1.21)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 0.23)*30./PIm_dist);
if ((PdHit_PIm==42)&&(delta_t_pim > 1.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.)*30./PIm_dist);

if ((PdHit_PIm==45)&&(delta_t_pim < -4.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 4.8)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < 1.)&&(delta_t_pim > -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.45)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < -2.)&&(delta_t_pim > -4.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 3.)*30./PIm_dist);

//if ((PdHit_PIm==46)&&(delta_t_pim > 1.03)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.04)*30./PIm_dist);
//if ((PdHit_PIm==46)&&(delta_t_pim < 1.03)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 0.02)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim < 1.)&&(delta_t_pim > -1.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.26)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 6.)*30./PIm_dist);

//   };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIm >= 0) && (ph_PIm <=30)){
// if ((ph_PIm > pim_fid_b_1) && (ph_PIm < pim_fid_a_1)){ 

//time_pim[0][PdHit_PIm-1]->Fill(P_PIm,delta_t_pim,1.);

if ((PdHit_PIm==42)&&(delta_t_pim > 1.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.)*30./PIm_dist);

if ((PdHit_PIm==45)&&(delta_t_pim < -4.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 4.8)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < 1.)&&(delta_t_pim > -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.45)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < -2.)&&(delta_t_pim > -4.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 3.)*30./PIm_dist);

//if ((PdHit_PIm==46)&&(delta_t_pim > 1.03)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.04)*30./PIm_dist);
//if ((PdHit_PIm==46)&&(delta_t_pim < 1.03)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 0.02)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim < 1.)&&(delta_t_pim > -1.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.26)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 6.)*30./PIm_dist);



  // };//end of the fiducial cut for second part of sector1
   }; //end of the second part of sector1
  
  
  if ((ph_PIm >= 30) && (ph_PIm <=90)) {
//  if ((ph_PIm > pim_fid_b_1+60) && (ph_PIm < pim_fid_a_1+60)){

//time_pim[1][PdHit_PIm-1]->Fill(P_PIm,delta_t_pim,1.);



if ((PdHit_PIm==40)&&(delta_t_pim > 6.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 7.)*30./PIm_dist);
if ((PdHit_PIm==40)&&(delta_t_pim -2.5)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 3.)*30./PIm_dist);

if ((PdHit_PIm==41)&&(delta_t_pim < -3.6)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim +4.1)*30./PIm_dist);
if ((PdHit_PIm==41)&&(delta_t_pim < -2.)&&(delta_t_pim > -3.6)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.1)*30./PIm_dist);




if ((PdHit_PIm==42)&&(delta_t_pim < -4.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.8)*30./PIm_dist);



if ((PdHit_PIm==43)&&(delta_t_pim < -6.)&&(delta_t_pim > -7.7)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 6.95)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim < -7.7)&&(delta_t_pim > -10.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 8.3)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim < -11.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 11.5)*30./PIm_dist);


if ((PdHit_PIm==44)&&(delta_t_pim > -5.5)&&(delta_t_pim < -3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.)*30./PIm_dist);
if ((PdHit_PIm==44)&&(delta_t_pim > -8.)&&(delta_t_pim < -5.5)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 6.7)*30./PIm_dist);

if ((PdHit_PIm==45)&&(delta_t_pim > -2.)&&(delta_t_pim < 2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.8)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim > -7.)&&(delta_t_pim < -4.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 5.2)*30./PIm_dist);

if ((PdHit_PIm==46)&&(delta_t_pim < 1.)&&(delta_t_pim > -3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.9)*30./PIm_dist);
if ((PdHit_PIm==46)&&(delta_t_pim < -7.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 8.)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim < 3.)&&(delta_t_pim > 0.7)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 1.6)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < 0.7)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.12)*30./PIm_dist);


// };//end of the fiducial cut for sector2
};//end of the sector2

  
if ((ph_PIm >= 90) && (ph_PIm <=150)) {
//if ((ph_PIm > pim_fid_b_1+120) && (ph_PIm < pim_fid_a_1+120)){

//time_pim[2][PdHit_PIm-1]->Fill(P_PIm,delta_t_pim,1.);

//if (delta_t_pim < 2.055) beta_new_vs_p_pim[2][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1.);


if ((PdHit_PIm==25)&&(delta_t_pim < -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.5)*30./PIm_dist);


if ((PdHit_PIm==40)&&(delta_t_pim > 2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim-4.)*30./PIm_dist);


if ((PdHit_PIm==41)&&(delta_t_pim > 0.37)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.86)*30./PIm_dist);


if ((PdHit_PIm==42)&&(delta_t_pim < 3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 1.9)*30./PIm_dist);
if ((PdHit_PIm==42)&&(delta_t_pim > 3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 5.36)*30./PIm_dist);


if ((PdHit_PIm==44)&&(delta_t_pim < 9.)&&(delta_t_pim > 7.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 8.1)*30./PIm_dist);
if ((PdHit_PIm==44)&&(delta_t_pim > 9.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 10.5)*30./PIm_dist);

if ((PdHit_PIm==46)&&(delta_t_pim < 2.0)&&(delta_t_pim > -2.0)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.5)*30./PIm_dist);
if ((PdHit_PIm==46)&&(delta_t_pim > 2.0)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 3.25)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim < -1.2)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 2.)*30./PIm_dist);


// };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIm >= 150) && (ph_PIm <=210)){
// if ((ph_PIm > pim_fid_b_1+180) && (ph_PIm < pim_fid_a_1+180)){

//time_pim[3][PdHit_PIm-1]->Fill(P_PIm,delta_t_pim,1.);


if ((PdHit_PIm==39)&&(delta_t_pim < -3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.3)*30./PIm_dist);


if ((PdHit_PIm==41)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.6)*30./PIm_dist);

if ((PdHit_PIm==42)&&(delta_t_pim < -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim +3.)*30./PIm_dist);


if ((PdHit_PIm==43)&&(delta_t_pim > 0.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.9)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim < 0.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.)*30./PIm_dist);


if ((PdHit_PIm==46)&&(delta_t_pim > -0.85)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.2)*30./PIm_dist);
if ((PdHit_PIm==46)&&(delta_t_pim < -0.85)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.3)*30./PIm_dist);

if ((PdHit_PIm==48)&&(delta_t_pim > -0.1)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 0.45)*30./PIm_dist);
if ((PdHit_PIm==48)&&(delta_t_pim < -0.1)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 1.1)*30./PIm_dist);

// };//end of the fiducial cut for sector4
  }; //end of the sector4
  
if ((ph_PIm >= 210) && (ph_PIm <=270)) {
//if ((ph_PIm > pim_fid_b_1+240) && (ph_PIm < pim_fid_a_1+240)){

//time_pim[4][PdHit_PIm-1]->Fill(P_PIm,delta_t_pim,1.);


if ((PdHit_PIm==40)&&(delta_t_pim > -3.3)&&(delta_t_pim < -0.8)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.6)*30./PIm_dist);
if ((PdHit_PIm==40)&&(delta_t_pim < -3.3)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.85)*30./PIm_dist);


if ((PdHit_PIm==42)&&(delta_t_pim < 0.65)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.33)*30./PIm_dist);
if ((PdHit_PIm==42)&&(delta_t_pim > 0.65)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.0)*30./PIm_dist);


if ((PdHit_PIm==44)&&(delta_t_pim < 0.2)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.6)*30./PIm_dist);
if ((PdHit_PIm==44)&&(delta_t_pim > 0.2)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.73)*30./PIm_dist);

if ((PdHit_PIm==46)&&(delta_t_pim > 1.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 2.)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim < 3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 2.0)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim > 3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 4.0)*30./PIm_dist);


// };//end of the fiducial cut for sector5
}; //end of the sector5
  
  
 if ((ph_PIm >= 270) && (ph_PIm <=330)){
// if ((ph_PIm > pim_fid_b_1+300) && (ph_PIm < pim_fid_a_1+300)){

//time_pim[5][PdHit_PIm-1]->Fill(P_PIm,PIm_dist*(1./beta_nom_pim-1/beta_PIm)/30.,1.);

if ((PdHit_PIm==31)&&(delta_t_pim < 0.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.6)*30./PIm_dist);
if ((PdHit_PIm==31)&&(delta_t_pim > 0.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.47)*30./PIm_dist);


if ((PdHit_PIm==40)&&(delta_t_pim < -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.)*30./PIm_dist);

if ((PdHit_PIm==44)&&(delta_t_pim < -1.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.5)*30./PIm_dist);


if ((PdHit_PIm==45)&&(delta_t_pim < -2.8)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.7)*30./PIm_dist);


if ((PdHit_PIm==47)&&(delta_t_pim > -2.)&&(delta_t_pim < 0.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.84)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -2.)&&(delta_t_pim > -4.0)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 2.7)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -4.0)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 5.3)*30./PIm_dist);

if ((PdHit_PIm==48)&&(delta_t_pim > 1.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim -2.78)*30./PIm_dist);


// };//end of the fiducial cut for sector6
  }; //end of the sector6



  
 
 
 };
 
 
