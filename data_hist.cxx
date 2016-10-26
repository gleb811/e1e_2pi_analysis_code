#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include "cuts_data.h"
#include "data_hist.h"
#include "global.h"
#include <iostream>



using namespace std;


 void data_hist() {
 Float_t beta_nom_pip,beta_nom_pim,beta_nom_p;
 beta_nom_pip = P_PIp/sqrt(m_pip*m_pip+P_PIp*P_PIp);
beta_nom_pim = P_PIm/sqrt(m_pip*m_pip+P_PIm*P_PIm);
beta_nom_p = P_P/sqrt(m_proton*m_proton+P_P*P_P);
 
 /*if (sector ==1) {
if (pmt_hit == -1) ph_el_left[0][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[0][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[0][segment]->Fill(nphe,1.);
 };
if (sector ==2)  {
if (pmt_hit == -1) ph_el_left[1][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[1][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[1][segment]->Fill(nphe,1.);

};

if (sector ==3) {
 if (pmt_hit == -1) ph_el_left[2][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[2][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[2][segment]->Fill(nphe,1.);
 
 };
if (sector ==4)  {
if (pmt_hit == -1) ph_el_left[3][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[3][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[3][segment]->Fill(nphe,1.);

};
if (sector ==5) {
if (pmt_hit == -1) ph_el_left[4][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[4][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[4][segment]->Fill(nphe,1.);

 };
if (sector ==6) {
if (pmt_hit == -1) ph_el_left[5][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[5][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[5][segment]->Fill(nphe,1.);

 };*/



/*
if (sector ==1) {
h_cc_nphe_total_s1->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s1->Fill(theta_cc,ph_cc,1.);
 };
if (sector ==2)  {
h_cc_nphe_total_s2->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s2->Fill(theta_cc,ph_cc,1.);

};

if (sector ==3) {
h_cc_nphe_total_s3->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s3->Fill(theta_cc,ph_cc,1.); 
 };
if (sector ==4)  {
h_cc_nphe_total_s4->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s4->Fill(theta_cc,ph_cc,1.);
};
if (sector ==5) {
h_cc_nphe_total_s5->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s5->Fill(theta_cc,ph_cc,1.);
 };
if (sector ==6) {
h_cc_nphe_total_s6->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s6->Fill(theta_cc,ph_cc,1.);
 };
*/


//cout << P_dist << "  "<< beta_nom_p << " "<< beta_P<<"\n";
W_2pi_selection->Fill(W,Q2,1.);

//if ((Q2>0.45)&&(Q2<0.5)){
if ((W>1.3)&&(W<1.9)){

if ((ph_P >= 330)&& (ph_P <= 360)) {
time_p[0][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
beta_vs_p_p [0][PdHit_P-1]->Fill(P_P,beta_P,1.);
if (W>1.3) th_vs_p_p_1[0]->Fill(P_P,th_P,1.); 
 };//end of the first part of the first sector
 
 if ((ph_P >= 0) && (ph_P <= 30)) {
time_p[0][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.); 
beta_vs_p_p [0][PdHit_P-1]->Fill(P_P,beta_P,1.);
if (W>1.3) th_vs_p_p_1[0]->Fill(P_P,th_P,1.);
 };//end of the second part of the first sector
  
 if ((ph_P >= 30) && (ph_P <=90)){ 
time_p[1][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
beta_vs_p_p [1][PdHit_P-1]->Fill(P_P,beta_P,1.);
if (W>1.3) th_vs_p_p_1[1]->Fill(P_P,th_P,1.);
 };//end of the sector2
 
 if ((ph_P >=90) && (ph_P <=150)) {
time_p[2][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
beta_vs_p_p [2][PdHit_P-1]->Fill(P_P,beta_P,1.);
if (W>1.3) th_vs_p_p_1[2]->Fill(P_P,th_P,1.);
 };//end of the sector3
 
 if ((ph_P >= 150) && (ph_P <= 210)) {
time_p[3][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
beta_vs_p_p [3][PdHit_P-1]->Fill(P_P,beta_P,1.);
if (W>1.3) th_vs_p_p_1[3]->Fill(P_P,th_P,1.);
};//end of the sector4

if ((ph_P >= 210) && (ph_P <=270)) {
time_p[4][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
beta_vs_p_p [4][PdHit_P-1]->Fill(P_P,beta_P,1.);
if (W>1.3) th_vs_p_p_1[4]->Fill(P_P,th_P,1.);
 };//end of the sector5
 
 if ((ph_P >= 270) && (ph_P <=330)) {
 
time_p[5][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
beta_vs_p_p [5][PdHit_P-1]->Fill(P_P,beta_P,1.);
if (W>1.3) th_vs_p_p_1[5]->Fill(P_P,th_P,1.);
 };//end of the sector6


//-------------------------------------
if ((ph_PIp >= 330)&& (ph_PIp <= 360)) {
time_pip[0][PdHit_PIp-1]->Fill(P_PIp,PIp_dist/30.*(1./beta_nom_pip-1./beta_PIp),1.);
th_vs_p_pip_1[0]->Fill(P_PIp,th_PIp,1.);
beta_vs_p_pip [0][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1.);

 };//end of the first part of the first sector
 
 if ((ph_PIp >= 0) && (ph_PIp <= 30)) {
time_pip[0][PdHit_PIp-1]->Fill(P_PIp,PIp_dist/30.*(1./beta_nom_pip-1./beta_PIp),1.); 
th_vs_p_pip_1[0]->Fill(P_PIp,th_PIp,1.);
beta_vs_p_pip [0][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1.);

 };//end of the second part of the first sector
  
 if ((ph_PIp >= 30) && (ph_PIp <=90)){ 
time_pip[1][PdHit_PIp-1]->Fill(P_PIp,PIp_dist/30.*(1./beta_nom_pip-1./beta_PIp),1.);
th_vs_p_pip_1[1]->Fill(P_PIp,th_PIp,1.);

beta_vs_p_pip [1][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1.);
 };//end of the sector2
 
 if ((ph_PIp >=90) && (ph_PIp <=150)) {
time_pip[2][PdHit_PIp-1]->Fill(P_PIp,PIp_dist/30.*(1./beta_nom_pip-1./beta_PIp),1.);
th_vs_p_pip_1[2]->Fill(P_PIp,th_PIp,1.);
beta_vs_p_pip [2][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1.);
 };//end of the sector3
 
 if ((ph_PIp >= 150) && (ph_PIp <= 210)) {
 time_pip[3][PdHit_PIp-1]->Fill(P_PIp,PIp_dist/30.*(1./beta_nom_pip-1./beta_PIp),1.);
th_vs_p_pip_1[3]->Fill(P_PIp,th_PIp,1.);
beta_vs_p_pip [3][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1.);
};//end of the sector4

if ((ph_PIp >= 210) && (ph_PIp <=270)) {
time_pip[4][PdHit_PIp-1]->Fill(P_PIp,PIp_dist/30.*(1./beta_nom_pip-1./beta_PIp),1.);
th_vs_p_pip_1[4]->Fill(P_PIp,th_PIp,1.);
beta_vs_p_pip [4][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1.);
 };//end of the sector5
 
 if ((ph_PIp >= 270) && (ph_PIp <=330)) {
 
 time_pip[5][PdHit_PIp-1]->Fill(P_PIp,PIp_dist/30.*(1./beta_nom_pip-1./beta_PIp),1.);
th_vs_p_pip_1[5]->Fill(P_PIp,th_PIp,1.);
beta_vs_p_pip [5][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1.);
 };//end of the sector6

//-----------------------------------
if ((ph_PIm >= 330)&& (ph_PIm <= 360)) {
time_pim[0][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
th_vs_p_pim_1[0]->Fill(P_PIm,th_PIm,1.);
beta_vs_p_pim [0][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1.); 
 };//end of the first part of the first sector
 
 if ((ph_PIm >= 0) && (ph_PIm <= 30)) {
time_pim[0][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.); 
th_vs_p_pim_1[0]->Fill(P_PIm,th_PIm,1.);
beta_vs_p_pim [0][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1.);
 };//end of the second part of the first sector
  
 if ((ph_PIm >= 30) && (ph_PIm <=90)){ 
time_pim[1][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
th_vs_p_pim_1[1]->Fill(P_PIm,th_PIm,1.);
beta_vs_p_pim [1][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1.);
 };//end of the sector2
 
 if ((ph_PIm >=90) && (ph_PIm <=150)) {
time_pim[2][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
th_vs_p_pim_1[2]->Fill(P_PIm,th_PIm,1.);
beta_vs_p_pim [2][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1.);
 };//end of the sector3
 
 if ((ph_PIm >= 150) && (ph_PIm <= 210)) {
time_pim[3][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
th_vs_p_pim_1[3]->Fill(P_PIm,th_PIm,1.);
beta_vs_p_pim [3][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1.);
};//end of the sector4

if ((ph_PIm >= 210) && (ph_PIm <=270)) {
time_pim[4][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
th_vs_p_pim_1[4]->Fill(P_PIm,th_PIm,1.);
beta_vs_p_pim [4][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1.);
 };//end of the sector5
 
 if ((ph_PIm >= 270) && (ph_PIm <=330)) {
 
time_pim[5][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
th_vs_p_pim_1[5]->Fill(P_PIm,th_PIm,1.);
beta_vs_p_pim [5][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1.);
 };//end of the sector6



}; //end of 1.35<W<1.4
//}; // end of 0.45<Q2<0.5

/*
if (sector ==1) {
 hist_z_el_1->Fill(z_EL,1.);
//th_vs_p_e_1[0]->Fill(P_EL,th_EL,1.);
 };
if (sector ==2)  {
//th_vs_p_e_1[1]->Fill(P_EL,th_EL,1.);
hist_z_el_2->Fill(z_EL,1.);
};

if (sector ==3) {
 hist_z_el_3->Fill(z_EL,1.);
//th_vs_p_e_1[2]->Fill(P_EL,th_EL,1.);
 };
if (sector ==4)  {
//th_vs_p_e_1[3]->Fill(P_EL,th_EL,1.);
hist_z_el_4->Fill(z_EL,1.);
};
if (sector ==5) {
//th_vs_p_e_1[4]->Fill(P_EL,th_EL,1.);
 hist_z_el_5->Fill(z_EL,1.);
 };
if (sector ==6) {
//th_vs_p_e_1[5]->Fill(P_EL,th_EL,1.);
 hist_z_el_6->Fill(z_EL,1.);
 };
 
*/
 };
