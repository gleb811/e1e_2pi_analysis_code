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
#include "global.h"
#include <iostream>




using namespace std;

//ostringstream qqq1;



 bool cuts_empty::Electron_cuts_empty(){
   
   
   

   
   bool cuts_empty;
   Float_t th_min,par1,par2,fid_a,fid_b,a,b;
   Short_t i;
   Float_t ph_el_arr[3][6][18] = {{{1000.,1000.,23.,24.,25.,23.,25.,27.,25.,34.,26.,30.,28.,27.,34.,31.,28.,1000.},
                             {1000.,1000.,30.,25.,25.,22.,30.,21.,22.,24.,23.,25.,31.,27.,27.,29.,29.,1000.},
                             {1000.,1000.,23.,25.,24.,23.,23.,23.,31.,28.,27.,28.,31.,32.,32.,32.,30.,1000.},
                             {1000.,1000.,30.,20.,24.,1.,25.,26.,24.,27.,25.,27.,22.,20.,32.,39.,30.,1000.},
                             {1000.,1000.,31.,30.,25.,28.,27.,23.,22.,26.,31.,26.,47.,25.,35.,42.,36.,1000.},
                             {1000.,1000.,28.,28.,22.,28.,27.,26.,26.,28.,26.,29.,28.,27.,31.,36.,34.,1000.}},
			     {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
	                     {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
			     {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}},
			     {{1000.,1000.,22.,20.,24.,25.,22.,24.,26.,26.,32.,26.,23.,30.,28.,25.,25.,1000.},
                             {1000.,1000.,28.,20.,22.,30.,25.,23.,30.,33.,23.,24.,27.,35.,30.,33.,36.,1000.},
			     {1000.,1000.,23.,26.,28.,23.,23.,25.,24.,27.,25.,27.,26.,30.,33.,29.,36.,1000.},
			     {1000.,1000.,1.,1.,1.,1.,20.,20.,22.,26.,1.,27.,20.,23.,27.,25.,25.,1000.},
			     {1000.,1000.,27.,26.,25.,51.,38.,32.,29.,30.,31.,36.,35.,27.,38.,41.,31.,1000.},
			     {1000.,1000.,29.,21.,26.,33.,29.,26.,29.,32.,36.,33.,26.,32.,50.,35.,36.,1000.}}};


 cuts_empty = false; 
   
   
   if ((dc_z_EL > -2.3)&&(dc_z_EL < 2.)) {
     if (P_EL > 0.461) {
//      if ((sc_y < (0.52*sc_x-47.)) && (sc_y > (-0.52*sc_x+47.))) {  
   th_min=(9.5+17./(P_EL+0.2));
  par1=0.85+1.1*P_EL;
  par2=-62.8-30.*P_EL;       
   fid_a=37.3*(0.85+1.1*P_EL)*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1525./th_EL/th_EL));
     
   fid_b=-37.3*(0.85+1.1*P_EL)*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1525./th_EL/th_EL));   
   a = fid_a;
   b = fid_b; 
//   P_EL = floor((P_EL*1000) + 0.5)/1000;


       switch (sector) {
case 1 : 


 
 if (((ECT/P_EL) > (-0.03606*P_EL*P_EL+0.1357*P_EL+0.08835)) && ((ECT/P_EL) < (0.008336*P_EL*P_EL-0.0301*P_EL)+0.38)) {
 
// if (nphe > 25.) {


if ((ph_EL >= 330) && (ph_EL <= 360)){


if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+360) && (ph_EL < fid_a+360)){

// if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_1_empty->Fill(z_EL,1.);
 //~fid cuts in Cherenkov plane
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_1_empty->Fill(dc_z_EL-0.2,1.); 
 if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons
//if (avrg_nphe_sector1->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector1_after->Fill(nphe,1.);
//if ((theta_cc >1.64786+2.16008*segment-0.00433375*segment*segment)&&(theta_cc <11.3593+1.38551*segment+0.0451710*segment*segment)){

if (norm_nphe_s1->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./50.+1.)) > 0.8) {
//if (nphe > ph_el_arr[pmt_hit+1][0][segment]){

if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
   cuts_empty = true;
   };
}; 
  // };
   
   };
   };
   };
   };
  
   
   }; //fiducial
   }; //second part of sector 1
   
  if ((ph_EL >= 0) && (ph_EL <= 30)) {


if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b) && (ph_EL < fid_a)){

 
//if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_1_empty->Fill(z_EL,1.);
//~fid cuts in Cherenkov plane
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_1_empty->Fill(dc_z_EL,1.);
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
//cut on average number of photoelectrons
//if (avrg_nphe_sector1->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector1_after->Fill(nphe,1.);
//if ((theta_cc >1.64786+2.16008*segment-0.00433375*segment*segment)&&(theta_cc <11.3593+1.38551*segment+0.0451710*segment*segment)){

if (norm_nphe_s1->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./50.+1.)) > 0.8) {

//if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
   cuts_empty = true; 
   };
  };
//  }; 
   };
   };
   };
   };

 }; //fiducial
   }; //first part of sector 1


// }; //nphe cut
//??? 

     
 
 }; // ectot vs p cut
 
 break;
 
case 2 : 
 
 if (((ECT/P_EL) > (-0.03116*P_EL*P_EL+0.1228*P_EL+0.1099)) && ((ECT/P_EL) < (0.001805*P_EL*P_EL-0.01003*P_EL)+0.419)) {
  
// if (nphe > 25.) {

 if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+60) && (ph_EL < fid_a+60)){

// if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_2_empty->Fill(z_EL,1.);
 
 //~fid cuts in Cherenkov plane
 if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_2_empty->Fill(dc_z_EL,1.);
  if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons
//if (avrg_nphe_sector2->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector2_after->Fill(nphe,1.);
//if ((theta_cc >1.95003+2.00182*segment+0.00390572*segment*segment)&&(theta_cc <11.1869+1.37368*segment+0.047233*segment*segment)){
if (norm_nphe_s2->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./50.+1.)) > 0.8) {
//if (nphe > ph_el_arr[pmt_hit+1][1][segment]){
if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
   cuts_empty = true; 
   };
   };
 //  };
   
   };
   };
   };
   }; 
   
   }; //fiducial 

 //}; //nphe cut
 
 
  
 
 }; // ectot vs p cut

 break; 
 
case 3 : 
 
if (((ECT/P_EL) > (-0.024*P_EL*P_EL+0.1014*P_EL+0.09593)) && ((ECT/P_EL) < (0.003918*P_EL*P_EL-0.01219*P_EL)+0.3616)) { 
  
// if (nphe > 25.) {

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+120) && (ph_EL < fid_a+120)){
//if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_3_empty->Fill(z_EL,1.);
 //~fid cuts in Cherenkov plane
 
  if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_3_empty->Fill(dc_z_EL-0.1,1.);
//cut on average number of photoelectrons
//if (avrg_nphe_sector3->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector3_after->Fill(nphe,1.);
//if ((theta_cc >1.9837+2.0442*segment+0.0022649*segment*segment)&&(theta_cc <11.2718+1.35251*segment+0.047262*segment*segment)){
if (norm_nphe_s3->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./50.+1.)) > 0.8) {
//if (nphe > ph_el_arr[pmt_hit+1][2][segment]){

  if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
   cuts_empty = true; 
   };
  };
 //  };
   
   };
   };
   };
   }; 
   
   }; //fiducial  

// }; //nphe cut
 
 
 
 }; // ectot vs p cut

 break;  
 
case 4 : 
 
 if (((ECT/P_EL) > (-0.03557*P_EL*P_EL+0.1322*P_EL+0.0951)) && ((ECT/P_EL) < (-0.001583*P_EL*P_EL-0.008766*P_EL)+0.3724)) {
  
// if (nphe > 25.) {

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+180) && (ph_EL < fid_a+180)){
//if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_4_empty->Fill(z_EL,1.);
//~fid cuts in Cherenkov plane
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_4_empty->Fill(dc_z_EL-0.15,1.);
  if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
//cut on average number of photoelectrons
//if (avrg_nphe_sector4->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector4_after->Fill(nphe,1.);
//if ((theta_cc >1.36283+2.27581*segment-0.0137087*segment*segment)&&(theta_cc <11.7047+1.27222*segment+0.0517623*segment*segment)){
//if (nphe > ph_el_arr[pmt_hit+1][3][segment]){
if (norm_nphe_s4->GetBinContent(int((theta_cc+5.+1.)*200./60.),int((ph_cc+25.)*200./50.+1.)) > 0.8) {
   if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
   cuts_empty = true; 
   };
 }; 
   //};
   
   };
   };
   };
   };
   
   
     }; //fiducial 
     
        

 //}; //nphe cut
 
  
 }; // ectot vs p cut

 break;  
 
case 5 : 
 
 if (((ECT/P_EL) > (-0.03453*P_EL*P_EL+0.132*P_EL+0.08782)) && ((ECT/P_EL) < (0.004611*P_EL*P_EL-0.01966*P_EL)+0.3846)) {
//  if ((sc_x < 287.) || (sc_x > 307.)) {
// if ((sc_y >  (-9.*sc_x+2060.)) || (sc_y < (-9.*sc_x+1920.))) { 
 //if (nphe > 25.) {

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+240) && (ph_EL < fid_a+240)){
//if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_5_empty->Fill(z_EL,1.);
//~fid cuts in Cherenkov plane
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_5_empty->Fill(dc_z_EL,1.);
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
//cut on average number of photoelectrons
//if (avrg_nphe_sector5->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector5_after->Fill(nphe,1.);
//if ((theta_cc >2.17772+1.95055*segment+0.00993131*segment*segment)&&(theta_cc <11.9184+1.34684*segment+0.0471248*segment*segment)){
if (norm_nphe_s5->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./50.+1.)) > 0.8) {
//if (nphe > ph_el_arr[pmt_hit+1][4][segment]){
   if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
   if ((th_EL > (9.5 + 17./(P_EL+0.72 + 0.2) +16.3)) || (th_EL < (9.5 + 17./(P_EL+0.72 + 0.2) +14.7))) {
   cuts_empty = true; 
   };
   };
};
 //  };
   
   };
   };
   };
   };
   
     }; //fiducial     

 //}; //nphe cuts
 
//  }; // sc plane cuts
// }; // sc plane cuts
 }; // ectot vs p cut

 break;   
 
case 6 : 

 
// if (((ECT/P_EL) > (-0.0504534*P_EL*P_EL+0.176674*P_EL+0.0627479)) && ((ECT/P_EL) < (0.0112167*P_EL*P_EL-0.0273361*P_EL)+0.382789)) {
 if (((ECT/P_EL) > (-0.033*P_EL*P_EL+0.1336*P_EL+0.08992)) && ((ECT/P_EL) < (0.002539*P_EL*P_EL-0.005779*P_EL)+0.3669)) {  
 //if (nphe > 25.) {

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+300) && (ph_EL < fid_a+300)){
//if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_6_empty->Fill(z_EL,1.);
//~fid cuts in Cherenkov plane
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_6_empty->Fill(dc_z_EL,1.);
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons
//if (avrg_nphe_sector6->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector6_after->Fill(nphe,1.);
//if ((theta_cc >2.18831+2.01682*segment+0.00288938*segment*segment)&&(theta_cc <11.5751+1.25302*segment+0.0536114*segment*segment)){

//if (nphe > ph_el_arr[pmt_hit+1][5][segment]){
if (norm_nphe_s6->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./50.+1.)) > 0.8) {
if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
   cuts_empty = true; 
   };
};
  // };
   
   };
   };
   };
   };
   
      }; //fiducial   

 //}; //nphe cut
     
 
 }; // ectot vs p cut
 break;      
 
   
   }; // end of switch
   }; // end of z_EL cut
   };  // end of calorimeter threshold cut
  
   

   
   
   return cuts_empty;
   
   };
   
   
   
   
   
   
   
   
   ////////////////////////////////////////////
      bool cuts_empty::Proton_cuts_empty(){
       
   bool cuts_empty;
   Float_t m_p,p_fid_a,p_fid_b;
   m_p=0.938272;   
   p_fid_a = 24.*(1-exp(-1.*0.08*(th_P-9.)));
   p_fid_b = -1.*25.*(1-exp(-1.*0.1*(th_P-10.)));
   cuts_empty = false; 
   


if (fabs(dc_z_EL - dc_z_P) < 4.) {
//if ((n_P == 1)&& (beta_P < 0.9496*P_P/sqrt(m_p*m_p+0.9497*P_P*P_P-0.06649) + 0.04136) && (beta_P > 1.045*P_P/sqrt(m_p*m_p+0.896*P_P*P_P - 0.2) - 0.139)){

 if (( th_P > pow((P_P-0.04-0.304992),(0.0758186))*91.5643-48.2057 + 3.5 ) || ( th_P < pow((P_P-0.07-0.304992),(0.0758186))*91.5643-48.2057 + 1. ) || (P_P < 0.04+0.304992)) {    
if (( th_P > pow((P_P-0.04-0.374992),(0.0758186))*91.5643-48.2057 + 0.5 ) || ( th_P < pow((P_P-0.04-0.410992),(0.0758186))*91.5643-48.2057 + 0.5 ) || (P_P < 0.04+0.374992)) { 

if ((PdHit_P !=48)&& (((beta_P < ((P_P/(sqrt(P_P*P_P+0.938*0.938))+0.02)*(1.2+0.92*P_P)/(1.+P_P))) && (beta_P > ((P_P/(sqrt(P_P*P_P+0.938*0.938))-0.05)/((1.+P_P)/(0.9+1.06*P_P)))))||(PdHit_P > 40))){ 

//&& (beta_P < 0.9675*P_P/sqrt(m_p*m_p+0.9386*P_P*P_P-0.1723) + 0.0063) && (beta_P > 0.9408*P_P/sqrt(m_p*m_p+0.7455*P_P*P_P - 0.2544) - 0.1126)

if ((ph_P >= 330)&& (ph_P <= 360)) {

//&&( PdHit_P != 42)&&(PdHit_P != 45)&&(PdHit_P != 47)&&(PdHit_P !=48)
// &&( PdHit_P != 42)&&(PdHit_P != 45)&&(PdHit_P != 47)&&(PdHit_P !=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)) {

  if ((ph_P > p_fid_b+360) && (ph_P < p_fid_a+360)){

cuts_empty = true;
};//end of fiducial cut for the first part of the first sector
  
 //};//end of W-cut
 
 };//end of the first part of the first sector
 
 if ((ph_P >= 0) && (ph_P <= 30)) {
 
 //&&( PdHit_P != 42)&&(PdHit_P != 45)&&(PdHit_P != 47)&&(PdHit_P !=48)
 //&&( PdHit_P != 42)&&(PdHit_P != 45)&&(PdHit_P != 47)&&(PdHit_P !=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)) {
if ((ph_P > p_fid_b) && (ph_P < p_fid_a)){


cuts_empty = true;
};//end of fiducial cut for the second part of the first sector
//};//end of W-cut
  };//end of the second part of the first sector
 
 
 
 if ((ph_P >= 30) && (ph_P <=90)){ 
 //&&(PdHit_P!=24)&&(PdHit_P!=41)&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=44)&&(PdHit_P!=45)&&(PdHit_P!=46)&&(PdHit_P!=48)
 //&&(PdHit_P!=45)&&(PdHit_P!=46)&&(PdHit_P!=47)&&(PdHit_P!=41)&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=44)&&(PdHit_P!=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)) {
 if ((ph_P > p_fid_b+60) && (ph_P < p_fid_a+60)){

if (( th_P > 26.5087*P_P*P_P*P_P -116.557*P_P*P_P+ 175.167*P_P-61.7717+2.7 ) || ( th_P < 26.5087*P_P*P_P*P_P -116.557*P_P*P_P+ 175.167*P_P-61.7717-2.8 )) {
cuts_empty = true;
};
};//end of the fiducial cut for sector2
//};//end of W-cut
 
 };//end of the sector2
 
 if ((ph_P >=90) && (ph_P <=150)) {
 
 //&&(PdHit_P!=25)&&(PdHit_P!=40)&&(PdHit_P!=42)&&(PdHit_P!=44)&&(PdHit_P!=46)&&(PdHit_P!=48)
 //&&(PdHit_P!=40)&&(PdHit_P!=41)&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=46)&&(PdHit_P!=47)&&(PdHit_P!=44)&&(PdHit_P!=25)&&(PdHit_P!=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)){

if ((ph_P > p_fid_b+120) && (ph_P < p_fid_a+120)){

cuts_empty = true;
};//end of the fiducial cut for sector3
//}; //end of W-cut
  };//end of the sector3
 
 if ((ph_P >= 150) && (ph_P <= 210)) {
 
//&&(PdHit_P!=42)&&(PdHit_P!=39)&&(PdHit_P!=48) 
 //&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=46)&&(PdHit_P!=39)&&(PdHit_P!=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)){
if ((ph_P > p_fid_b+180) && (ph_P < p_fid_a+180)){

cuts_empty = true;
};//end of the fiducial cut for sector4
//}; //end of W-cut

};//end of the sector4

if ((ph_P >= 210) && (ph_P <=270)) {

//&&(PdHit_P!=40)&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=17)&&(PdHit_P!=48)
//&&(PdHit_P!=40)&&(PdHit_P!=42)&&(PdHit_P!=44)&&(PdHit_P!=46)&&(PdHit_P!=47)&&(PdHit_P!=17)&&(PdHit_P!=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)) {
 if ((ph_P > p_fid_b+240) && (ph_P < p_fid_a+240)){

if (( th_P > 31.2482*(P_P+0.045)*(P_P+0.045)*(P_P+0.045) -135.817*(P_P+0.045)*(P_P+0.045)+ 198.038*(P_P+0.045)-66.968+1.8) || ( th_P < 31.2482*(P_P-0.01)*(P_P-0.01)*(P_P-0.01) -135.817*(P_P-0.01)*(P_P-0.01)+ 198.038*(P_P-0.01)-66.968-2.5)) {
 cuts_empty = true;
 };
 };//end of the fiducial cut for sector5
// };//end of W-cut
 
 };//end of the sector5
 
 if ((ph_P >= 270) && (ph_P <=330)) {
 
 //&&(PdHit_P!=40)&&(PdHit_P!=45)&&(PdHit_P!=31)&&(PdHit_P!=48)&&(PdHit_P!=47)
 //&&(PdHit_P!=40)&&(PdHit_P!=43)&&(PdHit_P!=44)&&(PdHit_P!=45)&&(PdHit_P!=31)&&(PdHit_P!=48)&&(PdHit_P!=47)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)){
if ((ph_P > p_fid_b+300) && (ph_P < p_fid_a+300)){

cuts_empty = true;
 };//end of the fiducial cut for sector6
// };//end of W-cut
 };//end of the sector6
 
 
 
 };
 };
 };
 };
 

    return cuts_empty;
    };
    
    
    
    /////////////////////
    
    
   bool cuts_empty::PIp_cuts_empty(){
       
   bool cuts_empty;
   Float_t m_pip,pip_fid_a,pip_fid_b,beta_nom_pip;
   m_pip = 0.13957;
  pip_fid_a = 24.*(1-exp(-1.*0.08*(th_PIp-9.)));
  pip_fid_b = -1.*25.*(1-exp(-1.*0.1*(th_PIp-10.)));
  beta_nom_pip = P_PIp/sqrt(m_pip*m_pip+P_PIp*P_PIp);
   cuts_empty = false; 
   
      if (fabs(dc_z_EL - dc_z_PIp) < 4.) {
   
     //cout << th_PIp<< " rgdgdf "<<ph_PIp<<" iiiiii "<<beta_PIp<<" riuthy "<<P_PIp<<"\n";
//   if ((n_PIp == 1)&& (PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30. < 0.0001769/(P_PIp*P_PIp*P_PIp*P_PIp+0.0001471)+0.8465)&&(PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30. > -0.0002121/(P_PIp*P_PIp*P_PIp*P_PIp+5.685e-05)-0.8411)) {
   
 
   //(beta_PIp < 0.8*P_PIp/sqrt(m_pip*m_pip+0.9108*P_PIp*P_PIp-0.001768) + 0.2) && (beta_PIp > 1.054*P_PIp/sqrt(m_pip*m_pip+0.7001*P_PIp*P_PIp - 0.006497) - 0.2999)
   //&& (beta_PIp < 0.8*P_PIp/sqrt(m_pip*m_pip+0.91*P_PIp*P_PIp-0.0034) + 0.2) && (beta_PIp > 1.05*P_PIp/sqrt(m_pip*m_pip+0.7*P_PIp*P_PIp - 0.0056) - 0.297)
   
   //&& (beta_PIp < 0.8*P_PIp/sqrt(m_pip*m_pip+0.91*P_PIp*P_PIp-0.0034) + 0.2) && (beta_PIp > 1.05*P_PIp/sqrt(m_pip*m_pip+0.7*P_PIp*P_PIp - 0.0056) - 0.297)
   
   if  ((PdHit_PIp !=48)&& (((beta_PIp < ((206.-P_PIp-0.02)*(pow(((200.-P_PIp)/(200.+P_PIp)),0.7)*(P_PIp+0.05))/sqrt((P_PIp+0.05)*(P_PIp+0.05)+0.019)/(200.+P_PIp+0.02)+0.019)) && (beta_PIp > (((((1.+5.*1.07*(P_PIp-0.07))*(P_PIp-0.07))))/sqrt((P_PIp-0.07)*(P_PIp-0.07)+0.138*0.138)/(1+5.*(P_PIp-0.07))-0.1)) && (PdHit_PIp <= 40))||((PdHit_PIp > 40)&&((beta_PIp > ((P_PIp/(sqrt(P_PIp*P_PIp+0.938*0.938))+0.02)*(1.2+0.92*P_PIp)/(1.+P_PIp))))))) {
      
      
   if ((ph_PIp >= 330) && (ph_PIp <=360)){
   
   //&&(PdHit_PIp != 42)&&(PdHit_PIp != 45)&&(PdHit_PIp != 46)&&(PdHit_PIp != 47)&&(PdHit_PIp !=48)
   //&&(PdHit_PIp != 42)&&(PdHit_PIp != 45)&&(PdHit_PIp != 47)&&(PdHit_PIp !=48)
   
   if ((ph_PIp > pip_fid_b+360) && (ph_PIp < pip_fid_a+360)){

  if (( th_PIp > (pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1) || (th_PIp < (pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1) || (P_PIp < 0.0575818)) {

if (( th_PIp > (304.23*(P_PIp+0.1)*(P_PIp+0.1)*(P_PIp+0.1) -255.798*(P_PIp+0.1)*(P_PIp+0.1)+497.462*(P_PIp+0.1) +38.0385)*exp(-1.85*(P_PIp+0.1)) +11. -1.5) || ( th_PIp <  (304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) +8. -11.)) { 


   cuts_empty = true;
  };
}; 
   };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIp >= 0) && (ph_PIp <=30)){
  
  //&&(PdHit_PIp != 42)&&(PdHit_PIp != 45)&&(PdHit_PIp != 46)&&(PdHit_PIp != 47)&&(PdHit_PIp !=48)
  //&&( PdHit_PIp != 42)&&(PdHit_PIp != 45)&&(PdHit_PIp != 47)&&(PdHit_PIp !=48)

 if ((ph_PIp > pip_fid_b) && (ph_PIp < pip_fid_a)){


if (( th_PIp > (pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1) || (th_PIp < (pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1) || (P_PIp < 0.0575818)) {

if (( th_PIp > (304.23*(P_PIp+0.1)*(P_PIp+0.1)*(P_PIp+0.1) -255.798*(P_PIp+0.1)*(P_PIp+0.1)+497.462*(P_PIp+0.1) +38.0385)*exp(-1.85*(P_PIp+0.1)) +11. -1.5) || ( th_PIp <  (304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) +8. -11.)) {


cuts_empty = true;
   };//end of the fiducial cut for second part of sector1
   }; //end of the second part of sector1
  };
};
  
  if ((ph_PIp >= 30) && (ph_PIp <=90)) {
  
 //&&(PdHit_PIp!=24)&&(PdHit_PIp!=45)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=41)&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=44)&&(PdHit_PIp!=48) //&&(PdHit_PIp!=45)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=41)&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=44)&&(PdHit_PIp!=48)

if ((ph_PIp > pip_fid_b+60) && (ph_PIp < pip_fid_a+60)){

if (( th_PIp > pow((P_PIp-0.415068),(0.226449))*48.7564 + 2.79478-1.) || ( th_PIp <pow((P_PIp-0.449975),( 0.315164 ))*36.608 +  9.74262-1.) || (P_PIp < 0.415068)) {
if (( th_PIp > (387.289*P_PIp*P_PIp*P_PIp -758.466*P_PIp*P_PIp+ 842.881*P_PIp-299.953+15.)*exp(-2*P_PIp)) || ( th_PIp < (387.289*P_PIp*P_PIp*P_PIp -758.466*P_PIp*P_PIp+ 842.881*P_PIp-299.953-15.)*exp(-2*P_PIp))) {


cuts_empty = true;

};
};
 };//end of the fiducial cut for sector2
};//end of the sector2
  
  
if ((ph_PIp >= 90) && (ph_PIp <=150)) {

//&&(PdHit_PIp!=40)&&(PdHit_PIp!=41)&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=44)&&(PdHit_PIp!=25)&&(PdHit_PIp!=48)
//&&(PdHit_PIp!=40)&&(PdHit_PIp!=41)&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=44)&&(PdHit_PIp!=25)&&(PdHit_PIp!=48)

if ((ph_PIp > pip_fid_b+120) && (ph_PIp < pip_fid_a+120)){

if (( th_PIp > pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374) || (th_PIp < pow((P_PIp -0.454898),( 0.289291))* 35.7267+6.65908+1.5) || (P_PIp < 0.416536)) {
 if (( th_PIp > ((10000*P_PIp*P_PIp*P_PIp-3607.41*P_PIp*P_PIp+ 1725.72*P_PIp-10.6776)*exp(-4.7*P_PIp)) || (th_PIp < (10000*P_PIp*P_PIp*P_PIp-4505.62*P_PIp*P_PIp+  2056.24  *P_PIp -77.4077 +5.)*exp(-4.7*P_PIp)))) {
  if (( th_PIp > pow((P_PIp-0.097536),(0.188376))*67.4593-21.4374+260.5-206.) || (th_PIp < pow((P_PIp-0.14),(0.188373))*67.4593-21.4374+260.5-210.) || (P_PIp < 0.097536)) {
 
 

cuts_empty = true;

};
};
};

 };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIp >= 150) && (ph_PIp <=210)){
 
 //&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=46)&&(PdHit_PIp!=39)&&(PdHit_PIp!=48)
 //&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=46)&&(PdHit_PIp!=39)&&(PdHit_PIp!=48)
 
 if ((ph_PIp > pip_fid_b+180) && (ph_PIp < pip_fid_a+180)){

if ((th_PIp > pow((P_PIp-0.412699),(0.214407))*52.0544 -0.0995427 - 2.3) || (th_PIp < pow((P_PIp-0.452908),(0.102883))*84.0374  -40.301+ 1.5) || (P_PIp < 0.412699)) {
if ((th_PIp > (1600*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1068.36*(P_PIp+0.03)*(P_PIp+0.03)+ 775.016*(P_PIp+0.03)-1.13034)*exp(-2.75*(P_PIp+0.03))) || (th_PIp < (pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.45*P_PIp)-7.)) { 

cuts_empty = true;

};
};

 };//end of the fiducial cut for sector4
 
 }; //end of the sector4
  
if ((ph_PIp >= 210) && (ph_PIp <=270)) {

//&&(PdHit_PIp!=40)&&(PdHit_PIp!=42)&&(PdHit_PIp!=44)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=17)&&(PdHit_PIp!=48)
//&&(PdHit_PIp!=40)&&(PdHit_PIp!=42)&&(PdHit_PIp!=44)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=17)&&(PdHit_PIp!=48)

if ((ph_PIp > pip_fid_b+240) && (ph_PIp < pip_fid_a+240)){

if (( th_PIp > pow((P_PIp-0.304992),(0.0758186))*91.5643-48.2057 -1.) || ( th_PIp < pow((P_PIp -0.36848),( 0.0864219))*70.4769  -34.9998+ 1.5) || (P_PIp < 0.30499)) {
if (( th_PIp > (525.498*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1284.98*(P_PIp+0.03)*(P_PIp+0.03)+1460.67*(P_PIp+0.03)-499.999)*exp(-1.94*(P_PIp+0.03))) || ( th_PIp < (525.498*(P_PIp-0.02)*(P_PIp-0.02)*(P_PIp-0.02) -1284.98*(P_PIp-0.02)*(P_PIp-0.02)+1460.67*(P_PIp-0.02)-499.999)*exp(-1.94*(P_PIp-0.02)) - 4.7)) {
if (( th_PIp > (304.23*(P_PIp)*(P_PIp)*(P_PIp) -255.798*(P_PIp)*(P_PIp)+497.462*(P_PIp) +38.0385)*exp(-1.85*(P_PIp))) || ( th_PIp < (304.23*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -255.798*(P_PIp+0.03)*(P_PIp+0.03)+497.462*(P_PIp+0.03) +38.0385)*exp(-1.85*(P_PIp+0.03)) -11.)) {



cuts_empty = true;

};
};
};




 };//end of the fiducial cut for sector5

}; //end of the sector5
  
  
 if ((ph_PIp >= 270) && (ph_PIp <=330)){
 
 //&&(PdHit_PIp!=40)&&(PdHit_PIp!=44)&&(PdHit_PIp!=45)&&(PdHit_PIp!=31)&&(PdHit_PIp!=48)&&(PdHit_PIp!=47)
 //&&(PdHit_PIp!=40)&&(PdHit_PIp!=43)&&(PdHit_PIp!=44)&&(PdHit_PIp!=45)&&(PdHit_PIp!=31)&&(PdHit_PIp!=48)&&(PdHit_PIp!=47)

 if ((ph_PIp > pip_fid_b+300) && (ph_PIp < pip_fid_a+300)){

if (( th_PIp > pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374-1.) || ( th_PIp < pow((P_PIp -0.454098),(0.0912936))*58.2946-20.4843+1.5) || (P_PIp < 0.415)) {
if (( th_PIp > pow((P_PIp-0.097536),(0.188376))*67.4593-21.4374+260.5-210.) || ( th_PIp < pow((P_PIp-0.14),(0.188373))*67.4593-21.4374+260.5-214.)|| (P_PIp < 0.1) ) {


cuts_empty = true;

};
};


 };//end of the fiducial cut for sector6
 
 }; //end of the sector6
  
  };
  
  
   };
    return cuts_empty;
   };
   
   ////////////////////
      
   bool cuts_empty::PIm_cuts_empty(){
       
   bool cuts_empty;
   Float_t m_pim,th_min,par1,par2,pim_fid_a,pim_fid_b,beta_nom_pim;
   
    
     m_pim = 0.13957; 
  th_min=(11.09+8./(0.472*P_PIm+0.117));
  par1=0.705+1.1*P_PIm;
  par2=-63.2-33.3*P_PIm;       
   pim_fid_a=30.5*pow((sin((th_PIm-th_min)*0.01745)),(par1+par2/th_PIm+1530./th_PIm/th_PIm))-1;
     
   pim_fid_b=-30.5*pow((sin((th_PIm-th_min)*0.01745)),(par1+par2/th_PIm+1530./th_PIm/th_PIm))+1; 
   beta_nom_pim = P_PIm/sqrt(m_pim*m_pim+P_PIm*P_PIm);
   cuts_empty = false; 
   
   if (fabs(dc_z_EL - dc_z_PIm) < 4.) {
 
   if  ((PdHit_PIm !=48)&& (((beta_PIm < ((206.-P_PIm-0.02)*(pow(((200.-P_PIm)/(200.+P_PIm)),0.7)*(P_PIm+0.05))/sqrt((P_PIm+0.05)*(P_PIm+0.05)+0.019)/(200.+P_PIm+0.02)+0.019)) && (beta_PIm > (((((1.+5.*1.07*(P_PIm-0.07))*(P_PIm-0.07))))/sqrt((P_PIm-0.07)*(P_PIm-0.07)+0.138*0.138)/(1+5.*(P_PIm-0.07))-0.1)) && (PdHit_PIm <= 40))||((PdHit_PIm > 40)&&((beta_PIm > ((P_PIm/(sqrt(P_PIm*P_PIm+0.938*0.938))+0.02)*(1.2+0.92*P_PIm)/(1.+P_PIm))))))) {
   
   
 if ((ph_PIm >= 330) && (ph_PIm <=360)){
 
 //&&(PdHit_PIm!=42)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
 //&&(PdHit_PIm!=42)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)



  if ((ph_PIm > pim_fid_b+360) && (ph_PIm < pim_fid_a+360)){
 
  if (( th_PIm >  (11.09+8./(0.472*(P_PIm+0.25)+0.117))+85.) || (th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+75.)) {
 if (( th_PIm >  (11.09+8./(0.472*(P_PIm+0.25)+0.117))+101.) || (th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+96.5)) {
 
   cuts_empty = true;
   };
   };
   
   };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIm >= 0) && (ph_PIm <=30)){
  
  //&&(PdHit_PIm!=42)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
  //&&(PdHit_PIm!=42)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)

 
 if ((ph_PIm > pim_fid_b) && (ph_PIm < pim_fid_a)){


 if (( th_PIm >  (11.09+8./(0.472*(P_PIm+0.25)+0.117))+85.) || (th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+75.)) {
 if (( th_PIm >  (11.09+8./(0.472*(P_PIm+0.25)+0.117))+101.) || (th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+96.5)) {
 
cuts_empty = true;

};
};

   };//end of the fiducial cut for second part of sector1
   };//end of the second part of sector1
  
    if ((ph_PIm >= 30) && (ph_PIm <=90)) {
    
 //  &&(PdHit_PIm!=16)&&(PdHit_PIm!=24)&&(PdHit_PIm!=41)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48) //&&(PdHit_PIm!=16)&&(PdHit_PIm!=24)&&(PdHit_PIm!=41)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)


if ((ph_PIm > pim_fid_b+60) && (ph_PIm < pim_fid_a+60)){

cuts_empty = true;
 };//end of the fiducial cut for sector2
};//end of the sector2
    
if ((ph_PIm >= 90) && (ph_PIm <=150)) {

//&&(PdHit_PIm!=25)&&(PdHit_PIm!=40)&&(PdHit_PIm!=41)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
//&&(PdHit_PIm!=25)&&(PdHit_PIm!=40)&&(PdHit_PIm!=41)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)

 
 if ((ph_PIm > pim_fid_b+120) && (ph_PIm < pim_fid_a+120)){

if ((th_PIm > 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)+2.) || (th_PIm < 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)-2.5)) {
cuts_empty = true;
};
 };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIm >= 150) && (ph_PIm <=210)){
 
 //&&(PdHit_PIm!=39)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
 //&&(PdHit_PIm!=39)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
 
 
 if ((ph_PIm > pim_fid_b+180) && (ph_PIm < pim_fid_a+180)){


if (( th_PIm > 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.15)+1.81169e-07)+86.) || (th_PIm < 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.15)+1.81169e-07)+77.)) {
if (( th_PIm > 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.35)+1.81169e-07)+66.) || (th_PIm < 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.31)+1.81169e-07)+62.)) {
cuts_empty = true;
};
};
 };//end of the fiducial cut for sector4
 }; //end of the sector4
  
if ((ph_PIm >= 210) && (ph_PIm <=270)) {

//&&(PdHit_PIm!=17)&&(PdHit_PIm!=40)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
//&&(PdHit_PIm!=17)&&(PdHit_PIm!=40)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)

if ((ph_PIm > pim_fid_b+240) && (ph_PIm < pim_fid_a+240)){


if ((th_PIm > 36.152+3.69909e-05/(5.40783e-06*(P_PIm)+1.81169e-07)+4.) || (th_PIm < 36.152+3.69909e-05/(5.40783e-06*(P_PIm)+1.81169e-07)-1.5)) {

cuts_empty = true;
};
 };//end of the fiducial cut for sector5
}; //end of the sector5
    
 if ((ph_PIm >= 270) && (ph_PIm <=330)){
 
 //&&(PdHit_PIm!=31)&&(PdHit_PIm!=40)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48) //&&(PdHit_PIm!=31)&&(PdHit_PIm!=40)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)

 
 if ((ph_PIm > pim_fid_b+300) && (ph_PIm < pim_fid_a+300)){

cuts_empty = true;
 };//end of the fiducial cut for sector6
 }; //end of the sector6
  
   
  
  };
   };
    return cuts_empty;
   }; 
   
   
