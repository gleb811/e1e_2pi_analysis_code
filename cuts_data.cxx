#include <TMath.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include "cuts_data.h"
#include "global.h"
#include <iostream>



using namespace std;

ostringstream qqq1;


   bool cuts_data::Electron_cuts_data(){
   
   
   

   
   bool cuts_data;
   Float_t th_min,par1,par2,fid_a,fid_b,a,b;
   Short_t i;
   Float_t ph_el_arr[3][6][18] = {{{1000.,1000.,26.,24.,25.,28.,24.,30.,30.,34.,26.,30.,32.,27.,36.,31.,28.,1000.},
                             {1000.,1000.,31.,20.,25.,22.,23.,21.,22.,24.,21.,31.,31.,31.,21.,29.,29.,1000.},
                             {1000.,1000.,27.,32.,24.,23.,19.,23.,31.,28.,27.,28.,32.,34.,27.,37.,37.,1000.},
                             {1000.,1000.,28.,20.,24.,19.,28.,26.,24.,22.,25.,29.,22.,23.,34.,43.,27.,1000.},
                             {1000.,1000.,31.,33.,25.,28.,29.,23.,22.,28.,35.,29.,51.,29.,30.,44.,36.,1000.},
                             {1000.,1000.,32.,28.,26.,29.,30.,30.,28.,46.,27.,39.,32.,29.,32.,33.,36.,1000.}},
			     {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
	                     {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
			     {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}},
			     {{1000.,1000.,18.,0.,24.,22.,22.,27.,26.,28.,32.,26.,23.,30.,28.,29.,25.,1000.},
                             {1000.,1000.,31.,20.,22.,33.,27.,26.,35.,35.,23.,24.,30.,45.,40.,40.,46.,1000.},
			     {1000.,1000.,23.,26.,28.,23.,23.,28.,24.,29.,25.,27.,29.,35.,37.,29.,36.,1000.},
			     {1000.,1000.,0.,0.,0.,0.,25.,30.,26.,27.,0.,27.,20.,28.,32.,27.,34.,1000.},
			     {1000.,1000.,27.,26.,27.,51.,40.,35.,31.,33.,33.,36.,40.,37.,35.,44.,31.,1000.},
			     {1000.,1000.,31.,21.,28.,33.,31.,28.,25.,28.,38.,36.,29.,35.,50.,35.,36.,1000.}}};
			     
			     
			     
/* Float_t th_vs_seg_cc_arr[2][6][18] = {{{0., 14.9943, 16.1178, 17.944, 19.7594, 21.7965, 24.0849, 25.7348, 27.287, 29.7387, 32.118, 34.554, 37.6616, 39.9555, 41.8899, 43.9917, 44.9401,0.},
                             {0., 14.4619, 15.8294, 17.3894, 19.5803, 21.6848, 23.4877, 25.6689, 27.2811, 29.7402, 32.1088, 34.3413, 37.7965, 39.9164, 41.93, 43.8781, 44.9396,0.},
                             {0., 14.6396, 15.8791, 17.681, 19.5967, 21.4913, 23.4704, 25.7742, 27.3183, 29.4058, 32.1693, 34.3859, 37.7187, 39.9258, 41.8853, 43.7546, 44.8999,0.},
                             {0., 15.0222, 16.1285, 17.8297, 19.7349, 21.93, 24.0715, 25.7447, 27.4221, 29.2554, 31.7778, 34.6932, 37.5713, 39.6431, 41.9178, 43.9135, 44.9405,0.},
                             {0., 15.1934, 17.3943, 18.4212, 19.8698, 21.9927, 24.0866, 25.9066, 28.4957, 30.6615, 33.2479, 35.3405, 37.4675, 38.6754, 41.6704, 44.0829, 44.9731,0.},
                             {0., 14.7309, 16.1419, 17.8516, 19.5722, 21.5672, 23.7528, 25.4896, 27.317, 29.8693, 32.1174, 34.5841, 37.7443, 40.01, 42.0101, 44.3091, 45.0681,0.}},
			     
			     
			     {{0., 9.70842, 11.6354, 13.4193, 15.059, 16.8693, 18.3642, 20.481, 22.5206, 23.9727, 25.4384, 27.9494, 29.6408, 33.1416, 35.9524, 38.0618, 40.6736,0.},
                             {0., 10.0549, 11.7277, 13.3703, 14.522, 16.1756, 18.3656, 20.0143, 22.3117, 23.8772, 25.3239, 27.825, 29.4949, 32.6421, 35.7132, 38.164, 40.3459,0.},
		             {0., 10.0805, 11.7872, 13.4138, 14.878, 16.6669, 18.492, 20.1909, 22.4518, 24.0873, 25.337, 27.7718, 29.7123, 33.1764, 36.0684, 38.3921, 40.5504,0.},
	                     {0., 9.89467, 11.8022, 13.4676, 14.8875, 16.5666, 18.234, 20.494, 22.3551, 24.1051, 25.5014, 27.4897, 29.1652, 32.9201, 35.3268, 37.8445, 40.393,0.},
		             {0., 10.2692, 11.1835, 13.4246, 14.8611, 16.5345, 18.2405, 20.3903, 22.0994, 23.9946, 26.2842, 28.3204, 31.0017, 32.6091, 36.7202, 38.2352, 40.8262,0.},
			     {0., 10.084, 11.641, 13.3442, 14.8961, 16.6584, 18.0017, 20.4061, 22.2266, 23.8465, 25.3002, 27.8278, 29.4927, 32.6927, 35.9115, 38.1519, 40.59,0.}}};*/
			     
			     
 Float_t th_vs_seg_cc_arr[2][6][18] = {{{0, 14.1133, 15.3753, 17.1888, 18.9788, 20.9753, 23.1317, 24.8555, 26.4884, 28.7626, 31.0185, 33.4327, 36.3182, 38.8195, 40.9003, 43.0034, 44.229,0.},
                             {0., 13.7274, 15.1581, 16.7182, 18.7325, 20.7666, 22.6319, 24.7182, 26.4459, 28.7519, 30.9952, 33.2372, 36.4083, 38.7015, 40.8938, 42.9257, 44.174,0.},
                             {0., 13.8798, 15.2084, 16.9716, 18.8154, 20.6872, 22.6378, 24.8429, 26.5011, 28.5067, 31.0552, 33.2685, 36.3783, 38.8005, 40.9158, 42.8609, 44.175,0.},
                             {0., 14.1676, 15.4142, 17.1058, 18.9349, 21.0361, 23.0975, 24.8628, 26.5717, 28.3813, 30.7322, 33.4418, 36.1686, 38.5218, 40.8193, 42.902, 44.1826,0.},
                             {0., 14.4373, 16.3017, 17.6014, 18.9795, 20.9873, 23.0632, 24.9872, 27.4241, 29.5447, 32.082, 34.165, 36.3894, 37.6569, 40.8454, 43.1083, 44.282,0.},
                             {0., 13.9564, 15.3918, 17.1004, 18.7928, 20.749, 22.7943, 24.6424, 26.4686, 28.8655, 30.9812, 33.4581, 36.369, 38.7905, 40.9936, 43.2829, 44.3217,0.}},
			     
			     
			     {{0., 10.5894, 12.3756, 14.1748, 15.8384, 17.6905, 19.3173, 21.3626, 23.3202, 24.9398, 26.5138, 29.0621, 30.9914, 34.278, 36.9419, 39.0501, 41.3847,0.},
                             {0., 10.7894, 12.3926, 14.0418, 15.373, 17.0938, 19.2219, 20.9716, 23.1494, 24.8626, 26.403, 28.9296, 30.8862, 33.8615, 36.7493, 39.1163, 41.1115,0.},
		             {0., 10.8404, 12.4489, 14.1227, 15.6575, 17.471, 19.3248, 21.1219, 23.2713, 24.9794, 26.4234, 28.8907, 31.0595, 34.3022, 37.0379, 39.2859, 41.2753,0.},
	                     {0., 10.7492, 12.5137, 14.1913, 15.6846, 17.4605, 19.2089, 21.3768, 23.2076, 24.9714, 26.5435, 28.7212, 30.5629, 34.0419, 36.4253, 38.856, 41.1509,0.},
		             {0., 10.9971, 12.3757, 14.2543, 15.852, 17.6056, 19.2968, 21.3097, 23.1726, 25.1176, 27.4525, 29.5022, 32.0802, 33.6269, 37.5453, 39.2098, 41.5173,0.},
			     {0., 10.8585, 12.3911, 14.0955, 15.6754, 17.4765, 18.9602, 21.2534, 23.075, 24.8503, 26.4364, 28.9539, 30.868, 33.9122, 36.9279, 39.1781, 41.3363,0.}}};			     
			     			     
			     



  // Float_t ph_el_r[6][18] = {{1000.,1000.,18.,20.,24.,22.,22.,27.,26.,28.,32.,26.,23.,30.,28.,29.,25.,1000.},
  //                           {1000.,1000.,31.,20.,22.,33.,27.,26.,35.,35.,23.,24.,30.,45.,40.,40.,46.,1000.},
//			     {1000.,1000.,23.,26.,28.,23.,23.,28.,24.,29.,25.,27.,29.,35.,37.,29.,36.,1000.},
	//		     {1000.,1000.,20.,21.,22.,24.,25.,34.,26.,27.,24.,27.,20.,28.,32.,27.,34.,1000.},
		//	     {1000.,1000.,27.,26.,27.,51.,40.,35.,31.,33.,33.,46.,40.,37.,48.,44.,31.,1000.},
			///     {1000.,1000.,31.,21.,28.,33.,31.,28.,32.,38.,38.,36.,29.,35.,90.,35.,36.,1000.}};


//for (i=0;i<18;i++){
//cout << ph_el[0][5][i] << "\n";
//};
 cuts_data = false; 
   
   
// if ((LiveTime > 0.89) && (LiveTime <0.925) && (inclusive > 48000) &&(inclusive < 60000) && (elastic > 16400) && (elastic < 20500)&&(hist_twopi_cut->GetBinContent(two_pions_block+1) > 120.)&&(hist_twopi_cut->GetBinContent(two_pions_block+1) < 400.)){



   if ((LiveTime > 0.89) && (LiveTime <0.925) && (inclusive > 48000) &&(inclusive < 60000) && (elastic > 16400) && (elastic < 20500)){
     if ((dc_z_EL > -2.3)&&(dc_z_EL < 2.)) {
     if (P_EL > 0.461) {
//     if ((sc_y < (0.52*sc_x-47.)) && (sc_y > (-0.52*sc_x+47.))) {
//     if ((z_EL > -2.)&&(z_EL < 2.)) {
     
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
  hist_ectot_sector1->Fill(P_EL,ECT/P_EL,1.); 
  
W_incl->Fill(W,Q2,1.);

if ((ph_EL >= 330) && (ph_EL <= 360)){

ph_vs_th_1 -> Fill(th_EL,ph_EL-360,1.);

if ((P_EL < 1.75999) && (P_EL > 0.4)){
ph_vs_th_1pe[int((P_EL*100-40)/8)]->Fill(th_EL,ph_EL-360,1);
};

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+360) && (ph_EL < fid_a+360)){
if ((P_EL < 1.75999) && (P_EL > 0.4)){
ph_vs_th_1pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-360,1.);
};

 W_incl_fid->Fill(W,Q2,1.);
 if (W>1.3) th_vs_p_e_2[0]->Fill(P_EL,th_EL,1.);
 nphe_sector1->Fill(nphe,1.);
 if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_1->Fill(dc_z_EL,1.);
 //~fid cuts in Cherenkov plane
 if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons

//if ((theta_cc >1.64786+2.16008*segment-0.00433375*segment*segment)&&(theta_cc <11.3593+1.38551*segment+0.0451710*segment*segment)){

//if ((theta_cc >th_vs_seg_cc_arr[1][0][segment])&&(theta_cc <th_vs_seg_cc_arr[0][0][segment])){

if (norm_nphe_s1->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./60.+1.)) > 0.8) {


nphe_sector1_after->Fill(nphe,1.);
//cout << pmt_hit << "     " << segment << "\n";

if (pmt_hit == -1) ph_el_left[0][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[0][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[0][segment]->Fill(nphe,1.);

if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
//if (nphe > 50){

//th_cc_vs_seg[0]->Fill(segment+1,theta_cc,1.);
//if (pmt_hit == -1) ph_el_left[0][segment]->Fill(nphe,1.);
//if (pmt_hit == 0) ph_el_both[0][segment]->Fill(nphe,1.);
//if (pmt_hit == 1) ph_el_right[0][segment]->Fill(nphe,1.);
   Wbefore[0]->Fill(W,1.);
   cuts_data = true;
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
ph_vs_th_1 -> Fill(th_EL,ph_EL,1.); 
if ((P_EL < 1.75999) &&(P_EL > 0.4)){
ph_vs_th_1pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL,1.);

qqq1.str("");
qqq1 << 400+80*int((P_EL*100-40)/8) << " < P_{el} <" << 400+80*(int((P_EL*100-40)/8)+1) << " MeV";
ph_vs_th_1pe[int((P_EL*100-40)/8)]->SetTitle(qqq1.str().c_str());
ph_vs_th_1pe[int((P_EL*100-40)/8)]->SetTitleSize(0.09);
};
if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b) && (ph_EL < fid_a)){
if ((P_EL < 1.75999) &&(P_EL > 0.4)){
ph_vs_th_1pe_fid[int((P_EL*100-40)/8)] -> Fill(th_EL,ph_EL,1.);
};

 W_incl_fid->Fill(W,Q2,1.);
if (W>1.3) th_vs_p_e_2[0]->Fill(P_EL,th_EL,1.); 
nphe_sector1->Fill(nphe,1.); 
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_1->Fill(dc_z_EL+0.2,1.);
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons


//if ((theta_cc >1.64786+2.16008*segment-0.00433375*segment*segment)&&(theta_cc <11.3593+1.38551*segment+0.0451710*segment*segment)){

//if ((theta_cc >th_vs_seg_cc_arr[1][0][segment])&&(theta_cc <th_vs_seg_cc_arr[0][0][segment])){

if (norm_nphe_s1->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./60.+1.)) > 0.8) {


nphe_sector1_after->Fill(nphe,1.);

if (pmt_hit == -1) ph_el_left[0][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[0][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[0][segment]->Fill(nphe,1.);

if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
//if (nphe > 50){
//th_cc_vs_seg[0]->Fill(segment+1,theta_cc,1.);
//if (pmt_hit == -1) ph_el_left[0][segment]->Fill(nphe,1.);
//if (pmt_hit == 0) ph_el_both[0][segment]->Fill(nphe,1.);
//if (pmt_hit == 1) ph_el_right[0][segment]->Fill(nphe,1.);
   Wbefore[0]->Fill(W,1.);
   cuts_data = true; 
  };
};
// }; 
 };
   };
   };
   };
 }; //fiducial
   }; //first part of sector 1


// }; //nphe cut

     
 
 }; // ectot vs p cut
 
 break;
 
case 2 : 

 
 
 if (((ECT/P_EL) > (-0.03116*P_EL*P_EL+0.1228*P_EL+0.1099)) && ((ECT/P_EL) < (0.001805*P_EL*P_EL-0.01003*P_EL)+0.419)) {
 hist_ectot_sector2->Fill(P_EL,ECT/P_EL,1.);  

W_incl->Fill(W,Q2,1.);

 ph_vs_th_2 -> Fill(th_EL,ph_EL-60,1.);

if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_2pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-60,1.);
qqq1.str("");
qqq1 << 400+80*int((P_EL*100-40)/8) << " < P_{el} <" << 400+80*(int((P_EL*100-40)/8)+1) << " MeV";
ph_vs_th_2pe[int((P_EL*100-40)/8)]->SetTitle(qqq1.str().c_str());
ph_vs_th_2pe[int((P_EL*100-40)/8)]->SetTitleSize(0.09);
}; 


 if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+60) && (ph_EL < fid_a+60)){
if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_2pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-60,1.);
};

 W_incl_fid->Fill(W,Q2,1.);
 if (W>1.3) th_vs_p_e_2[1]->Fill(P_EL,th_EL,1.);
 //cout << pmt_hit << "     " << segment << "\n";nphe_sector2->Fill(nphe,1.);
   if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_2->Fill(dc_z_EL+0.15,1.);
 
 //~fid cuts in Cherenkov plane
  if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons

//if ((theta_cc >1.95003+2.00182*segment+0.00390572*segment*segment)&&(theta_cc <11.1869+1.37368*segment+0.047233*segment*segment)){

//if ((theta_cc >th_vs_seg_cc_arr[1][1][segment])&&(theta_cc <th_vs_seg_cc_arr[0][1][segment])){

if (norm_nphe_s2->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./60.+1.)) > 0.8) {
nphe_sector2_after->Fill(nphe,1.);
//cout << pmt_hit << "     " << segment << "\n";

if (pmt_hit == -1) ph_el_left[1][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[1][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[1][segment]->Fill(nphe,1.);

if (nphe > ph_el_arr[pmt_hit+1][1][segment]){
//if (nphe > 50.){
//th_cc_vs_seg[1]->Fill(segment+1,theta_cc,1.);
//if (pmt_hit == -1) ph_el_left[1][segment]->Fill(nphe,1.);
//if (pmt_hit == 0) ph_el_both[1][segment]->Fill(nphe,1.);
//if (pmt_hit == 1) ph_el_right[1][segment]->Fill(nphe,1.);
   Wbefore[1]->Fill(W,1.);
   cuts_data = true; 
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
  hist_ectot_sector3->Fill(P_EL,ECT/P_EL,1.);  

W_incl->Fill(W,Q2,1.);

ph_vs_th_3 -> Fill(th_EL,ph_EL-120,1.);

if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_3pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-120,1.);
qqq1.str("");
qqq1 << 400+80*int((P_EL*100-40)/8) << " < P_{el} <" << 400+80*(int((P_EL*100-40)/8)+1) << " MeV";
ph_vs_th_3pe[int((P_EL*100-40)/8)]->SetTitle(qqq1.str().c_str());
ph_vs_th_3pe[int((P_EL*100-40)/8)]->SetTitleSize(0.09);
};

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+120) && (ph_EL < fid_a+120)){
if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_3pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-120,1.);
}; 

 W_incl_fid->Fill(W,Q2,1.);
 if (W>1.3) th_vs_p_e_2[2]->Fill(P_EL,th_EL,1.);
 nphe_sector3->Fill(nphe,1.);
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_3->Fill(dc_z_EL-0.05,1.);
 //~fid cuts in Cherenkov plane
 
  if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons

//nphe_sector3_after->Fill(nphe,1.);
//if ((theta_cc >1.9837+2.0442*segment+0.0022649*segment*segment)&&(theta_cc <11.2718+1.35251*segment+0.047262*segment*segment)){

//if ((theta_cc >th_vs_seg_cc_arr[1][2][segment])&&(theta_cc <th_vs_seg_cc_arr[0][2][segment])){

if (norm_nphe_s3->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./60.+1.)) > 0.8) {
nphe_sector3_after->Fill(nphe,1.);

if (pmt_hit == -1) ph_el_left[2][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[2][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[2][segment]->Fill(nphe,1.);

if (nphe > ph_el_arr[pmt_hit+1][2][segment]){
//if (nphe > 50.){
//th_cc_vs_seg[2]->Fill(segment+1,theta_cc,1.);
//if (pmt_hit == -1) ph_el_left[2][segment]->Fill(nphe,1.);
//if (pmt_hit == 0) ph_el_both[2][segment]->Fill(nphe,1.);
//if (pmt_hit == 1) ph_el_right[2][segment]->Fill(nphe,1.);
   Wbefore[2]->Fill(W,1.);
   cuts_data = true; 
};
};
//   };
   
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
 hist_ectot_sector4->Fill(P_EL,ECT/P_EL,1.);  

W_incl->Fill(W,Q2,1.);

ph_vs_th_4 -> Fill(th_EL,ph_EL-180,1);


if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_4pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-180,1.);
qqq1.str("");
qqq1 << 400+80*int((P_EL*100-40)/8) << " < P_{el} <" << 400+80*(int((P_EL*100-40)/8)+1) << " MeV";
ph_vs_th_4pe[int((P_EL*100-40)/8)]->SetTitle(qqq1.str().c_str());
ph_vs_th_4pe[int((P_EL*100-40)/8)]->SetTitleSize(0.09);
};

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+180) && (ph_EL < fid_a+180)){
if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_4pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-180,1.);
};


 W_incl_fid->Fill(W,Q2,1.);
if (W>1.3) th_vs_p_e_2[3]->Fill(P_EL,th_EL,1.); 
nphe_sector4->Fill(nphe,1.);
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_4->Fill(dc_z_EL-0.15-0.1,1.);
//~fid cuts in Cherenkov plane
  if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons


//if ((theta_cc >1.36283+2.27581*segment-0.0137087*segment*segment)&&(theta_cc <11.7047+1.27222*segment+0.0517623*segment*segment)){

//if ((theta_cc >th_vs_seg_cc_arr[1][3][segment])&&(theta_cc <th_vs_seg_cc_arr[0][3][segment])){

if (norm_nphe_s4->GetBinContent(int((theta_cc+5.+1.)*200./60.),int((ph_cc+25.)*200./60.+1.)) > 0.8) {
nphe_sector4_after->Fill(nphe,1.);

if (pmt_hit == -1) ph_el_left[3][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[3][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[3][segment]->Fill(nphe,1.);

if (nphe > ph_el_arr[pmt_hit+1][3][segment]){
//if (nphe > 50.){
//th_cc_vs_seg[3]->Fill(segment+1,theta_cc,1.);
//if (pmt_hit == -1) ph_el_left[3][segment]->Fill(nphe,1.);
//if (pmt_hit == 0) ph_el_both[3][segment]->Fill(nphe,1.);
//if (pmt_hit == 1) ph_el_right[3][segment]->Fill(nphe,1.);
   Wbefore[3]->Fill(W,1.);
   cuts_data = true; 
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
 
case 5 : 

 
 
 if (((ECT/P_EL) > (-0.03453*P_EL*P_EL+0.132*P_EL+0.08782)) && ((ECT/P_EL) < (0.004611*P_EL*P_EL-0.01966*P_EL)+0.3846)) {
 hist_ectot_sector5->Fill(P_EL,ECT/P_EL,1.);   
// if ((sc_x < 287.) || (sc_x > 307.)) {
// if ((sc_y >  (-9.*sc_x+2060.)) || (sc_y < (-9.*sc_x+1920.))) {

W_incl->Fill(W,Q2,1.);

 ph_vs_th_5 -> Fill(th_EL,ph_EL-240,1);
 
 
if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_5pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-240,1.);
qqq1.str("");
qqq1 << 400+80*int((P_EL*100-40)/8) << " < P_{el} <" << 400+80*(int((P_EL*100-40)/8)+1) << " MeV";
ph_vs_th_5pe[int((P_EL*100-40)/8)]->SetTitle(qqq1.str().c_str());
ph_vs_th_5pe[int((P_EL*100-40)/8)]->SetTitleSize(0.09);
};

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+240) && (ph_EL < fid_a+240)){

//if ((th_EL > (9.5+17./((P_EL+0.32)+0.2))+25.)||(th_EL < (9.5+17./((P_EL+0.3)+0.2))+21.8)){
if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_5pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-240,1.);
};

 W_incl_fid->Fill(W,Q2,1.);

nphe_sector5->Fill(nphe,1.);
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_5->Fill(dc_z_EL-0.2,1.);
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons


//if ((theta_cc >2.17772+1.95055*segment+0.00993131*segment*segment)&&(theta_cc <11.9184+1.34684*segment+0.0471248*segment*segment)){

//if ((theta_cc >th_vs_seg_cc_arr[1][4][segment])&&(theta_cc <th_vs_seg_cc_arr[0][4][segment])){

if (norm_nphe_s5->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./60.+1.)) > 0.8) {
nphe_sector5_after->Fill(nphe,1.);

if (pmt_hit == -1) ph_el_left[4][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[4][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[4][segment]->Fill(nphe,1.);

//cout << pmt_hit << "   " << segment << "  " << ph_el_arr[pmt_hit+1][4][segment] << "\n";
if (nphe > ph_el_arr[pmt_hit+1][4][segment]){
//if (nphe > 50.){
//th_cc_vs_seg[4]->Fill(segment+1,theta_cc,1.);
//if (pmt_hit == -1) ph_el_left[4][segment]->Fill(nphe,1.);
//if (pmt_hit == 0) ph_el_both[4][segment]->Fill(nphe,1.);
//if (pmt_hit == 1) ph_el_right[4][segment]->Fill(nphe,1.);
   Wbefore[4]->Fill(W,1.);
   
   sc_sector5->Fill(sc_y,sc_x,1.);
   if ((th_EL > (9.5 + 17./(P_EL+0.72 + 0.2) +16.3)) || (th_EL < (9.5 + 17./(P_EL+0.72 + 0.2) +14.7))) { 
   if (W>1.3) th_vs_p_e_2[4]->Fill(P_EL,th_EL,1.);
   cuts_data = true; 
   };
};
};
//   };
   
   };
   };
   };
//   };
//   };//th_vs_p
     }; //fiducial     

 //}; //nphe cut
 
 
// }; // sc plane cuts
// }; // sc plane cuts
 
 }; // ectot vs p cut

 break;   
 
case 6 : 


 
 if (((ECT/P_EL) > (-0.033*P_EL*P_EL+0.1336*P_EL+0.08992)) && ((ECT/P_EL) < (0.002539*P_EL*P_EL-0.005779*P_EL)+0.3669)) {
  hist_ectot_sector6->Fill(P_EL,ECT/P_EL,1.); 

W_incl->Fill(W,Q2,1.);

ph_vs_th_6 -> Fill(th_EL,ph_EL-300,1); 



if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_6pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-300,1.);
qqq1.str("");
qqq1 << 400+80*int((P_EL*100-40)/8) << " < P_{el} <" << 400+80*(int((P_EL*100-40)/8)+1) << " MeV";
ph_vs_th_6pe[int((P_EL*100-40)/8)]->SetTitle(qqq1.str().c_str());
ph_vs_th_6pe[int((P_EL*100-40)/8)]->SetTitleSize(0.09);
};
if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+300) && (ph_EL < fid_a+300)){
if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
ph_vs_th_6pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-300,1.);
};


 W_incl_fid->Fill(W,Q2,1.);
if (W>1.3) th_vs_p_e_2[5]->Fill(P_EL,th_EL,1.); 
nphe_sector6->Fill(nphe,1.);
if ((W>1.65)&&(W<1.85)&&(Q2>0.4)&&(Q2<0.7)) hist_z_el_6->Fill(dc_z_EL,1.);
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons


//if ((theta_cc >2.18831+2.01682*segment+0.00288938*segment*segment)&&(theta_cc <11.5751+1.25302*segment+0.0536114*segment*segment)){

//if ((theta_cc >th_vs_seg_cc_arr[1][5][segment])&&(theta_cc <th_vs_seg_cc_arr[0][5][segment])){

if (norm_nphe_s6->GetBinContent(int((theta_cc+5.)*200./60.+1.),int((ph_cc+25.)*200./60.+1.)) > 0.8) {
nphe_sector6_after->Fill(nphe,1.);

if (pmt_hit == -1) ph_el_left[5][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[5][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[5][segment]->Fill(nphe,1.);

if (nphe > ph_el_arr[pmt_hit+1][5][segment]){
//if (nphe > 50.){
//th_cc_vs_seg[5]->Fill(segment+1,theta_cc,1.);
//if (pmt_hit == -1) ph_el_left[5][segment]->Fill(nphe,1.);
//if (pmt_hit == 0) ph_el_both[5][segment]->Fill(nphe,1.);
//if (pmt_hit == 1) ph_el_right[5][segment]->Fill(nphe,1.);

   Wbefore[5]->Fill(W,1.);
   cuts_data = true; 
};
};
//  };
   
//   };
   };
  };
   };
   
      }; //fiducial   

 //}; //nphe cut
 
   
 
 }; // ectot vs p cut
 break;      
 
   };
   }; // end of switch
    }; //end of z_EL cut
   }; //end of sc plane cuts
   };  // end of calorimeter threshold cut
   };  // end of quality check cuts
   

   
   
   return cuts_data;
   
   };




///////////////////////////////////////////////////////////////////////////////////////////////

   bool cuts_data::Proton_cuts_data(){
       
   bool cuts_data;
   Float_t m_p,p_fid_a,p_fid_b;
   m_p=0.938272;   
   p_fid_a = 24.*(1-exp(-1.*0.08*(th_P-9.)));
   p_fid_b = -1.*25.*(1-exp(-1.*0.1*(th_P-10.)));
   cuts_data = false; 

if ((W>1.3)&&(W<1.9)) {

if (fabs(dc_z_EL - dc_z_P) < 4.) {

//if ((PdHit_P !=48)&& (((beta_P < ((P_P/(sqrt(P_P*P_P+0.938*0.938))+0.02)*(1.2+0.92*P_P)/(1.+P_P))) && (beta_P > ((P_P/(sqrt(P_P*P_P+0.938*0.938))-0.05)/((1.+P_P)/(0.9+1.06*P_P)))))||(PdHit_P > 40))){
     
 if (( th_P > pow((P_P-0.04-0.304992),(0.0758186))*91.5643-48.2057 + 3.5 ) || ( th_P < pow((P_P-0.07-0.304992),(0.0758186))*91.5643-48.2057 + 1. ) || (P_P < 0.04+0.304992)) {    
if (( th_P > pow((P_P-0.04-0.374992),(0.0758186))*91.5643-48.2057 + 0.5 ) || ( th_P < pow((P_P-0.04-0.410992),(0.0758186))*91.5643-48.2057 + 0.5 ) || (P_P < 0.04+0.374992)) { 
 
if ((PdHit_P !=48)&& (((beta_P < ((P_P/(sqrt(P_P*P_P+0.938*0.938))+0.02)*(1.2+0.92*P_P)/(1.+P_P))) && (beta_P > ((P_P/(sqrt(P_P*P_P+0.938*0.938))-0.05)/((1.+P_P)/(0.9+1.06*P_P))))))){

if ((ph_P >= 330)&& (ph_P <= 360)) {

if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)&&(n_PIp == 1)) ph_th_p_1[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-360,1);

  if ((ph_P > p_fid_b+360) && (ph_P < p_fid_a+360)){
//if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)) ph_th_p_1[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-360,1);
ph_vs_th_p_1_w -> Fill(th_P,ph_P-360,1);

if ((W > 1.3)&&(n_PIp == 1)) th_vs_p_p_2[0]->Fill(P_P,th_P,1.);
padd_rate[0]->Fill(PdHit_P,1.);
padd_prot_mass[0]->Fill(PdHit_P,P_P/beta_P*sqrt(1.-beta_P*beta_P),1.);
cuts_data = true;
};//end of fiducial cut for the first part of the first sector
  
 //};//end of W-cut
 ph_vs_th_p_1 -> Fill(th_P,ph_P-360,1);
 
 };//end of the first part of the first sector
 
 if ((ph_P >= 0) && (ph_P <= 30)) {

if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)&&(n_PIp == 1)) ph_th_p_1[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P,1);
 
if ((ph_P > p_fid_b) && (ph_P < p_fid_a)){
//if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)) ph_th_p_1[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P,1);
ph_vs_th_p_1_w -> Fill(th_P,ph_P,1);
padd_prot_mass[0]->Fill(PdHit_P,P_P/beta_P*sqrt(1.-beta_P*beta_P),1.);
if ((W > 1.3)&&(n_PIp == 1)) th_vs_p_p_2[0]->Fill(P_P,th_P,1.);
padd_rate[0]->Fill(PdHit_P,1.);

cuts_data = true;
};//end of fiducial cut for the second part of the first sector
//};//end of W-cut
 ph_vs_th_p_1 -> Fill(th_P,ph_P,1);
 };//end of the second part of the first sector
 
 
 
 if ((ph_P >= 30) && (ph_P <=90)){ 

if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)&&(n_PIp == 1)) ph_th_p_2[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-60,1);
 if ((ph_P > p_fid_b+60) && (ph_P < p_fid_a+60)){
 
//  if ((th_P >26.5087*P_P*P_P*P_P -116.557*P_P*P_P+ 175.167*P_P-61.7717+2.7 ) || (th_P <26.5087*P_P*P_P*P_P -116.557*P_P*P_P+ 175.167*P_P-61.7717-2.8)){
//if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)) ph_th_p_2[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-60,1);
ph_vs_th_p_2_w -> Fill(th_P,ph_P-60,1);


if (( th_P > 26.5087*P_P*P_P*P_P -116.557*P_P*P_P+ 175.167*P_P-61.7717+2.7 ) || ( th_P < 26.5087*P_P*P_P*P_P -116.557*P_P*P_P+ 175.167*P_P-61.7717-2.8 )) {
if ((W > 1.3)&&(n_PIp == 1)) th_vs_p_p_2[1]->Fill(P_P,th_P,1.);
padd_rate[1]->Fill(PdHit_P,1.);
padd_prot_mass[1]->Fill(PdHit_P,P_P/beta_P*sqrt(1.-beta_P*beta_P),1.);
cuts_data = true;
};
//};//th_vs_p
};//end of the fiducial cut for sector2
//};//end of W-cut
 
ph_vs_th_p_2 -> Fill(th_P,ph_P-60,1);
 };//end of the sector2
 
 if ((ph_P >=90) && (ph_P <=150)) {
 
 if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)&&(n_PIp == 1)) ph_th_p_3[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-120,1);
if ((ph_P > p_fid_b+120) && (ph_P < p_fid_a+120)){
//if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)) ph_th_p_3[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-120,1);
ph_vs_th_p_3_w -> Fill(th_P,ph_P-120,1);

if ((W > 1.3)&&(n_PIp == 1)) th_vs_p_p_2[2]->Fill(P_P,th_P,1.);
padd_rate[2]->Fill(PdHit_P,1.);
padd_prot_mass[2]->Fill(PdHit_P,P_P/beta_P*sqrt(1.-beta_P*beta_P),1.);
cuts_data = true;
};//end of the fiducial cut for sector3
//}; //end of W-cut
 ph_vs_th_p_3 -> Fill(th_P,ph_P-120,1);
 };//end of the sector3
 
 if ((ph_P >= 150) && (ph_P <= 210)) {
 
if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)&&(n_PIp == 1)) ph_th_p_4[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-180,1);
if ((ph_P > p_fid_b+180) && (ph_P < p_fid_a+180)){
//if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)) ph_th_p_4[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-180,1);
ph_vs_th_p_4_w -> Fill(th_P,ph_P-180,1);

if ((W > 1.3)&&(n_PIp == 1)) th_vs_p_p_2[3]->Fill(P_P,th_P,1.);
padd_rate[3]->Fill(PdHit_P,1.);
padd_prot_mass[3]->Fill(PdHit_P,P_P/beta_P*sqrt(1.-beta_P*beta_P),1.);
cuts_data = true;
};//end of the fiducial cut for sector4
//}; //end of W-cut
ph_vs_th_p_4 -> Fill(th_P,ph_P-180,1);
};//end of the sector4

if ((ph_P >= 210) && (ph_P <=270)) {

if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)&&(n_PIp == 1)) ph_th_p_5[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-240,1);
 if ((ph_P > p_fid_b+240) && (ph_P < p_fid_a+240)){
 
//if ((P_P<=0.321436)||(th_P > pow((P_P -0.321436),(0.0704348))*88.0419-46.9342) || (th_P < 31.2482*(P_P-0.01)*(P_P-0.01)*(P_P-0.01) -135.817*(P_P-0.01)*(P_P-0.01)+ 198.038*(P_P-0.01)-66.968-2.5)||((th_P > 31.2482*(P_P+0.045)*(P_P+0.045)*(P_P+0.045) -135.817*(P_P+0.045)*(P_P+0.045)+ 198.038*(P_P+0.045)-66.968+1.8)&&(th_P < pow((P_P-0.371051),(  0.0649747))*87.0943  -49.9895 -1.))){ 

//if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)) ph_th_p_5[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-240,1);
 ph_vs_th_p_5_w -> Fill(th_P,ph_P-240,1);



if (( th_P > 31.2482*(P_P+0.045)*(P_P+0.045)*(P_P+0.045) -135.817*(P_P+0.045)*(P_P+0.045)+ 198.038*(P_P+0.045)-66.968+1.8) || ( th_P < 31.2482*(P_P-0.01)*(P_P-0.01)*(P_P-0.01) -135.817*(P_P-0.01)*(P_P-0.01)+ 198.038*(P_P-0.01)-66.968-2.5)) {
if ((W > 1.3)&&(n_PIp == 1)) th_vs_p_p_2[4]->Fill(P_P,th_P,1.);
padd_rate[4]->Fill(PdHit_P,1.);
padd_prot_mass[4]->Fill(PdHit_P,P_P/beta_P*sqrt(1.-beta_P*beta_P),1.);
 cuts_data = true;
 };
// };//th_vs_p
 };//end of the fiducial cut for sector5
// };//end of W-cut
 ph_vs_th_p_5 -> Fill(th_P,ph_P-240,1);
 };//end of the sector5
 
 if ((ph_P >= 270) && (ph_P <=330)) {
 
 if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)&&(n_PIp == 1)) ph_th_p_6[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-300,1);
if ((ph_P > p_fid_b+300) && (ph_P < p_fid_a+300)){
//if ((P_P > 0.2) && (P_P < 1.2)&&(W > 1.3)) ph_th_p_6[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-300,1);
ph_vs_th_p_6_w -> Fill(th_P,ph_P-300,1);

if ((W > 1.3)&&(n_PIp == 1)) th_vs_p_p_2[5]->Fill(P_P,th_P,1.);
padd_rate[5]->Fill(PdHit_P,1.);
padd_prot_mass[5]->Fill(PdHit_P,P_P/beta_P*sqrt(1.-beta_P*beta_P),1.);
cuts_data = true;
 };//end of the fiducial cut for sector6
// };//end of W-cut
 ph_vs_th_p_6 -> Fill(th_P,ph_P-300,1);
 };//end of the sector6
 
 
 
 };
   };
   };
   };
   };
    return cuts_data;
    };
    
    
/////////////////////////////
    
    bool cuts_data::PIp_cuts_data(){
       
   bool cuts_data;
   Float_t m_pip,pip_fid_a,pip_fid_b, beta_nom_pip;
   m_pip = 0.13957;
   pip_fid_a = 24.*(1-exp(-1.*0.08*(th_PIp-9.)));
   pip_fid_b = -1.*25.*(1-exp(-1.*0.1*(th_PIp-10.)));

  beta_nom_pip = P_PIp/sqrt(m_pip*m_pip+P_PIp*P_PIp);
   cuts_data = false; 
   
   if (fabs(dc_z_EL - dc_z_PIp) < 4.) {
   
   //btvsp_test->Fill(P_P,beta_P,1);
   //cout << th_PIp<< " rgdgdf "<<ph_PIp<<" iiiiii "<<beta_PIp<<" riuthy "<<P_PIp<<"\n";
   if  ((PdHit_PIp !=48)&& (((beta_PIp < ((206.-P_PIp-0.02)*(pow(((200.-P_PIp)/(200.+P_PIp)),0.7)*(P_PIp+0.05))/sqrt((P_PIp+0.05)*(P_PIp+0.05)+0.019)/(200.+P_PIp+0.02)+0.019)) && (beta_PIp > (((((1.+5.*1.07*(P_PIp-0.07))*(P_PIp-0.07))))/sqrt((P_PIp-0.07)*(P_PIp-0.07)+0.138*0.138)/(1+5.*(P_PIp-0.07))-0.1))))) {
   
//     if  ((PdHit_PIp !=48)&&((((beta_PIp < ((206.-P_PIp-0.02)*(pow((200.-P_PIp/(200.+P_PIp)),0.7)*(P_PIp+0.05))/sqrt((P_PIp+0.05)*(P_PIp+0.05)+0.019)/(200.+P_PIp+0.02)+0.02)) && (beta_PIp > (((((1.+5.*1.07*(P_PIp-0.07))*(P_PIp-0.07))))/sqrt((P_PIp-0.07)*(P_PIp-0.07)+0.138*0.138)/(1+5.*(P_PIp-0.07))-0.1)))&&(PdHit_PIp <= 40)))) { 
  //cout << "n_PIp = "<< n_PIp << "\n"; 
  
 //  if ((PdHit_PIp !=48)&&(beta_PIp < 0.8*P_PIp/sqrt(m_pip*m_pip+0.9108*P_PIp*P_PIp-0.001768) + 0.2) && (beta_PIp > 1.054*P_PIp/sqrt(m_pip*m_pip+0.7001*P_PIp*P_PIp - 0.006497) - 0.2999)){

   
 //  if ((PdHit_PIp !=48) && (beta_PIp < 0.8*P_PIp/sqrt(m_pip*m_pip+0.91*P_PIp*P_PIp-0.0034) + 0.2) && (beta_PIp > 1.05*P_PIp/sqrt(m_pip*m_pip+0.7*P_PIp*P_PIp - 0.0056) - 0.297)){
   
   //&& (beta_PIp < 0.8*P_PIp/sqrt(m_pip*m_pip+0.91*P_PIp*P_PIp-0.0034) + 0.2) && (beta_PIp > 1.05*P_PIp/sqrt(m_pip*m_pip+0.7*P_PIp*P_PIp - 0.0056) - 0.297)
   
 //  if ((PdHit_PIp !=48) && ( beta_PIp <  0.2392*P_PIp/sqrt(m_pip*m_pip+0.1164*P_PIp*P_PIp-0.0177) + 0.3346+0.1)&&( beta_PIp >  0.2351*P_PIp/sqrt(m_pip*m_pip+0.0264*P_PIp*P_PIp-0.0191) - 0.4943-0.1)){ 
      
   if ((ph_PIp >= 330) && (ph_PIp <=360)){
   
//   if ((th_PIp > (304.23*(P_PIp+0.1)*(P_PIp+0.1)*(P_PIp+0.1) -255.798*(P_PIp+0.1)*(P_PIp+0.1)+497.462*(P_PIp+0.1) +38.0385)*exp(-1.85*(P_PIp+0.1)) +9.5) || (th_PIp < pow((P_PIp -0.454098),(0.0912936))*58.2946-20.4843+0.9)||((th_PIp < (pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1)&&((th_PIp > pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374)||(P_PIp<=0.416536)))||((th_PIp < (304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) -3)&&(th_PIp > (pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1))){ 
   
     if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_1[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-360,1);
 ph_vs_th_pip_1 -> Fill(th_PIp,ph_PIp-360,1);
   if ((ph_PIp > pip_fid_b+360) && (ph_PIp < pip_fid_a+360)){
 
  
   
  if (( th_PIp > (pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1) || (th_PIp < (pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1) || (P_PIp < 0.0575818)) {

if (( th_PIp > (304.23*(P_PIp+0.1)*(P_PIp+0.1)*(P_PIp+0.1) -255.798*(P_PIp+0.1)*(P_PIp+0.1)+497.462*(P_PIp+0.1) +38.0385)*exp(-1.85*(P_PIp+0.1)) +11. -1.5) || ( th_PIp <  (304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) +8. -11.)) { 

 th_vs_p_pip_2[0]->Fill(P_PIp,th_PIp,1.);
 padd_prot_mass[0]->Fill(PdHit_PIp,P_PIp/beta_PIp*sqrt(1.-beta_PIp*beta_PIp),1.);
   cuts_data = true;
   };
   };
//   };//th_vs_p
   };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIp >= 0) && (ph_PIp <=30)){
  
//  if ((th_PIp > (304.23*(P_PIp+0.1)*(P_PIp+0.1)*(P_PIp+0.1) -255.798*(P_PIp+0.1)*(P_PIp+0.1)+497.462*(P_PIp+0.1) +38.0385)*exp(-1.85*(P_PIp+0.1)) +9.5) || (th_PIp < pow((P_PIp -0.454098),(0.0912936))*58.2946-20.4843+0.9)||((th_PIp < (pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1)&&((th_PIp > pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374)||(P_PIp<=0.416536)))||((th_PIp < (304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) -3)&&(th_PIp > (pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1))){ 
  
if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_1[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp,1);

 ph_vs_th_pip_1 -> Fill(th_PIp,ph_PIp,1);
 if ((ph_PIp > pip_fid_b) && (ph_PIp < pip_fid_a)){

if (( th_PIp > (pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1) || (th_PIp < (pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1) || (P_PIp < 0.0575818)) {

if (( th_PIp > (304.23*(P_PIp+0.1)*(P_PIp+0.1)*(P_PIp+0.1) -255.798*(P_PIp+0.1)*(P_PIp+0.1)+497.462*(P_PIp+0.1) +38.0385)*exp(-1.85*(P_PIp+0.1)) +11. -1.5) || ( th_PIp <  (304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) +8. -11.)) {

th_vs_p_pip_2[0]->Fill(P_PIp,th_PIp,1.);
padd_prot_mass[0]->Fill(PdHit_PIp,P_PIp/beta_PIp*sqrt(1.-beta_PIp*beta_PIp),1.);
cuts_data = true;
};
};
//};//th_vs_p
   };//end of the fiducial cut for second part of sector1
   }; //end of the second part of sector1
  
  
  if ((ph_PIp >= 30) && (ph_PIp <=90)) {
  
 if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_2[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-60,1);
ph_vs_th_pip_2 -> Fill(th_PIp,ph_PIp-60,1);

if ((ph_PIp > pip_fid_b+60) && (ph_PIp < pip_fid_a+60)){

//if ((P_PIp<=0.415068)||(th_PIp >  pow((P_PIp-0.415068),(0.226449))*48.7564 + 2.79478-1.) || (th_PIp < (387.289*P_PIp*P_PIp*P_PIp -758.466*P_PIp*P_PIp+ 842.881*P_PIp-299.953-15.)*exp(-2*P_PIp))||((th_PIp > (387.289*P_PIp*P_PIp*P_PIp -758.466*P_PIp*P_PIp+ 842.881*P_PIp-299.953+15.)*exp(-2*P_PIp) )&&(th_PIp < pow((P_PIp-0.449975),( 0.315164 ))*36.608 +  9.74262-1.))){ 

if (( th_PIp > pow((P_PIp-0.415068),(0.226449))*48.7564 + 2.79478-1.) || ( th_PIp <pow((P_PIp-0.449975),( 0.315164 ))*36.608 +  9.74262-1.) || (P_PIp < 0.415068)) {
if (( th_PIp > (387.289*P_PIp*P_PIp*P_PIp -758.466*P_PIp*P_PIp+ 842.881*P_PIp-299.953+15.)*exp(-2*P_PIp)) || ( th_PIp < (387.289*P_PIp*P_PIp*P_PIp -758.466*P_PIp*P_PIp+ 842.881*P_PIp-299.953-15.)*exp(-2*P_PIp))) {

th_vs_p_pip_2[1]->Fill(P_PIp,th_PIp,1.);
padd_prot_mass[1]->Fill(PdHit_PIp,P_PIp/beta_PIp*sqrt(1.-beta_PIp*beta_PIp),1.);
cuts_data = true;

};
};
//};//th_vs_p
 };//end of the fiducial cut for sector2
};//end of the sector2
  
  
if ((ph_PIp >= 90) && (ph_PIp <=150)) {

if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_3[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-120,1);
ph_vs_th_pip_3 -> Fill(th_PIp,ph_PIp-120,1);

if ((ph_PIp > pip_fid_b+120) && (ph_PIp < pip_fid_a+120)){

//if ((th_PIp < pow((P_PIp -0.454898),( 0.289291))* 35.7267+6.65908+1.5) || (th_PIp > (10000*P_PIp*P_PIp*P_PIp-3607.41*P_PIp*P_PIp+ 1725.72*P_PIp-10.6776)*exp(-4.7*P_PIp))||((th_PIp < (10000*P_PIp*P_PIp*P_PIp-4505.62*P_PIp*P_PIp+  2056.24  *P_PIp -77.4077 +5.)*exp(-4.7*P_PIp))&&((th_PIp > pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374)||(P_PIp <= 0.416536)))){ 

if (( th_PIp > pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374) || (th_PIp < pow((P_PIp -0.454898),( 0.289291))* 35.7267+6.65908+1.5) || (P_PIp < 0.416536)) {
 if (( th_PIp > ((10000*P_PIp*P_PIp*P_PIp-3607.41*P_PIp*P_PIp+ 1725.72*P_PIp-10.6776)*exp(-4.7*P_PIp)) || (th_PIp < (10000*P_PIp*P_PIp*P_PIp-4505.62*P_PIp*P_PIp+  2056.24  *P_PIp -77.4077 +5.)*exp(-4.7*P_PIp)))) {
  if (( th_PIp > pow((P_PIp-0.097536),(0.188376))*67.4593-21.4374+260.5-206.) || (th_PIp < pow((P_PIp-0.14),(0.188373))*67.4593-21.4374+260.5-210.) || (P_PIp < 0.097536)) {
 

th_vs_p_pip_2[2]->Fill(P_PIp,th_PIp,1.);
padd_prot_mass[2]->Fill(PdHit_PIp,P_PIp/beta_PIp*sqrt(1.-beta_PIp*beta_PIp),1.);
cuts_data = true;

};
};
};
//};//th_vs_p
 };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIp >= 150) && (ph_PIp <=210)){
 
 
if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_4[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-180,1);
 ph_vs_th_pip_4 -> Fill(th_PIp,ph_PIp-180,1);
 
 if ((ph_PIp > pip_fid_b+180) && (ph_PIp < pip_fid_a+180)){

// if ((th_PIp <pow((P_PIp-0.452908),(0.102883))*84.0374  -40.301+ 1.3) ||(th_PIp> (304.23*(P_PIp+0.18)*(P_PIp+0.18)*(P_PIp+0.18) -255.798*(P_PIp+0.18)*(P_PIp+0.18)+497.462*(P_PIp+0.18) +38.0385)*exp(-1.85*(P_PIp+0.18)) +5.)|| ((th_PIp< (304.23*(P_PIp+0.18)*(P_PIp+0.18)*(P_PIp+0.18) -255.798*(P_PIp+0.18)*(P_PIp+0.18)+497.462*(P_PIp+0.18) +38.0385)*exp(-1.85*(P_PIp+0.18)) - 2.)&&(th_PIp >(1600*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1068.36*(P_PIp+0.03)*(P_PIp+0.03)+ 775.016*(P_PIp+0.03)-1.13034)*exp(-2.75*(P_PIp+0.03))))||((th_PIp < (pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.45*P_PIp)-7.)&&((th_PIp >pow((P_PIp-0.412699),(0.214407))*52.0544 -0.0995427-2.1 )||(P_PIp <=0.412699)))){ 

if ((th_PIp > pow((P_PIp-0.412699),(0.214407))*52.0544 -0.0995427 - 2.3) || (th_PIp < pow((P_PIp-0.452908),(0.102883))*84.0374  -40.301+ 1.5) || (P_PIp < 0.412699)) {
if ((th_PIp > (1600*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1068.36*(P_PIp+0.03)*(P_PIp+0.03)+ 775.016*(P_PIp+0.03)-1.13034)*exp(-2.75*(P_PIp+0.03))) || (th_PIp < (pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.45*P_PIp)-7.)) { 

th_vs_p_pip_2[3]->Fill(P_PIp,th_PIp,1.);
padd_prot_mass[3]->Fill(PdHit_PIp,P_PIp/beta_PIp*sqrt(1.-beta_PIp*beta_PIp),1.);
cuts_data = true;

};
};
//};//th_vs_p
 };//end of the fiducial cut for sector4
 
 }; //end of the sector4
  
if ((ph_PIp >= 210) && (ph_PIp <=270)) {

if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_5[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-240,1);
ph_vs_th_pip_5 -> Fill(th_PIp,ph_PIp-240,1);
if ((ph_PIp > pip_fid_b+240) && (ph_PIp < pip_fid_a+240)){

//if ((th_PIp < (525.498*(P_PIp-0.02)*(P_PIp-0.02)*(P_PIp-0.02) -1284.98*(P_PIp-0.02)*(P_PIp-0.02)+1460.67*(P_PIp-0.02)-499.999)*exp(-1.94*(P_PIp-0.02)) -4.7) || (th_PIp > (304.23*(P_PIp)*(P_PIp)*(P_PIp) -255.798*(P_PIp)*(P_PIp)+497.462*(P_PIp) +38.0385)*exp(-1.85*(P_PIp)) ) ||((th_PIp < (304.23*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -255.798*(P_PIp+0.03)*(P_PIp+0.03)+497.462*(P_PIp+0.03) +38.0385)*exp(-1.85*(P_PIp+0.03)) -11.)&&((th_PIp > pow((P_PIp-0.304992),(0.0758186))*91.5643-48.2057 - 1.)||(P_PIp <=0.304992))) || ((th_PIp < pow((P_PIp -0.36848),( 0.0864219))*70.4769  -34.9998+ 1.5)&&(th_PIp > (525.498*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1284.98*(P_PIp+0.03)*(P_PIp+0.03)+1460.67*(P_PIp+0.03)-499.999)*exp(-1.94*(P_PIp+0.03))))){

if (( th_PIp > pow((P_PIp-0.304992),(0.0758186))*91.5643-48.2057 -1.) || ( th_PIp < pow((P_PIp -0.36848),( 0.0864219))*70.4769  -34.9998+ 1.5) || (P_PIp < 0.30499)) {
if (( th_PIp > (525.498*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1284.98*(P_PIp+0.03)*(P_PIp+0.03)+1460.67*(P_PIp+0.03)-499.999)*exp(-1.94*(P_PIp+0.03))) || ( th_PIp < (525.498*(P_PIp-0.02)*(P_PIp-0.02)*(P_PIp-0.02) -1284.98*(P_PIp-0.02)*(P_PIp-0.02)+1460.67*(P_PIp-0.02)-499.999)*exp(-1.94*(P_PIp-0.02)) - 4.7)) {
if (( th_PIp > (304.23*(P_PIp)*(P_PIp)*(P_PIp) -255.798*(P_PIp)*(P_PIp)+497.462*(P_PIp) +38.0385)*exp(-1.85*(P_PIp))) || ( th_PIp < (304.23*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -255.798*(P_PIp+0.03)*(P_PIp+0.03)+497.462*(P_PIp+0.03) +38.0385)*exp(-1.85*(P_PIp+0.03)) -11.)) {


th_vs_p_pip_2[4]->Fill(P_PIp,th_PIp,1.);
padd_prot_mass[4]->Fill(PdHit_PIp,P_PIp/beta_PIp*sqrt(1.-beta_PIp*beta_PIp),1.);
cuts_data = true;

};
};
};
//};//th_vs_p
 };//end of the fiducial cut for sector5

}; //end of the sector5
  
  
 if ((ph_PIp >= 270) && (ph_PIp <=330)){
 
 

if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_6[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-300,1);
 ph_vs_th_pip_6 -> Fill(th_PIp,ph_PIp-300,1);
 if ((ph_PIp > pip_fid_b+300) && (ph_PIp < pip_fid_a+300)){

//if ((th_PIp < pow((P_PIp -0.454098),(0.0912936))*58.2946-20.4843+1.5) || (P_PIp<=0.05+0.0942469)||(th_PIp > pow((P_PIp-0.05-0.0942469),( 0.0582707))*114.358-50 -0.5)||((th_PIp < pow((P_PIp-0.05-0.126994),( 0.0706829))* 110.073-50+2.) &&((th_PIp > pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374-1.)||(P_PIp<=0.416536)))){

if (( th_PIp > pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374-1.) || ( th_PIp < pow((P_PIp -0.454098),(0.0912936))*58.2946-20.4843+1.5) || (P_PIp < 0.415)) {
if (( th_PIp > pow((P_PIp-0.097536),(0.188376))*67.4593-21.4374+260.5-210.) || ( th_PIp < pow((P_PIp-0.14),(0.188373))*67.4593-21.4374+260.5-214.)|| (P_PIp < 0.1) ) {

 th_vs_p_pip_2[5]->Fill(P_PIp,th_PIp,1.);
 padd_prot_mass[5]->Fill(PdHit_PIp,P_PIp/beta_PIp*sqrt(1.-beta_PIp*beta_PIp),1.);
cuts_data = true;

};
};
//};//th_vs_p
 };//end of the fiducial cut for sector6
 
 }; //end of the sector6
  };
//  if (th_PIp < 20.)cuts_data = false;
  
  
   };
    return cuts_data;
   };
   
   
   /////////////////////////////////////

  
    bool cuts_data::PIm_cuts_data(){
       
   bool cuts_data;
   Float_t m_pim,th_min,par1,par2,pim_fid_a,pim_fid_b,beta_nom_pim;
   
    
   m_pim = 0.13957; 
  th_min=(11.09+8./(0.472*P_PIm+0.117));
  par1=0.705+1.1*P_PIm;
  par2=-63.2-33.3*P_PIm;       
   pim_fid_a=30.5*pow((sin((th_PIm-th_min)*0.01745)),(par1+par2/th_PIm+1530./th_PIm/th_PIm))-1;
     
   pim_fid_b=-30.5*pow((sin((th_PIm-th_min)*0.01745)),(par1+par2/th_PIm+1530./th_PIm/th_PIm))+1; 
   beta_nom_pim = P_PIm/sqrt(m_pim*m_pim+P_PIm*P_PIm);
   cuts_data = false; 
   
   if (fabs(dc_z_EL - dc_z_PIm) < 4.) {
 
   if  ((PdHit_PIm !=48)&& (((beta_PIm < ((206.-P_PIm-0.02)*(pow(((200.-P_PIm)/(200.+P_PIm)),0.7)*(P_PIm+0.05))/sqrt((P_PIm+0.05)*(P_PIm+0.05)+0.019)/(200.+P_PIm+0.02)+0.019)) && (beta_PIm > (((((1.+5.*1.07*(P_PIm-0.07))*(P_PIm-0.07))))/sqrt((P_PIm-0.07)*(P_PIm-0.07)+0.138*0.138)/(1+5.*(P_PIm-0.07))-0.1))))) {
   
//   if  ((PdHit_PIm !=48)&&((((beta_PIm < ((206.-P_PIm-0.02)*(pow((200.-P_PIm/(200.+P_PIm)),0.7)*(P_PIm+0.05))/sqrt((P_PIm+0.05)*(P_PIm+0.05)+0.019)/(200.+P_PIm+0.02)+0.02)) && (beta_PIm > (((((1.+5.*1.07*(P_PIm-0.07))*(P_PIm-0.07))))/sqrt((P_PIm-0.07)*(P_PIm-0.07)+0.138*0.138)/(1+5.*(P_PIm-0.07))-0.1)))&&(PdHit_PIm <= 40)))) {   
   
   //(beta_PIm < 0.8*P_PIm/sqrt(m_pim*m_pim+1.2*P_PIm*P_PIm+0.008756) + 0.2963) && (beta_PIm > 0.7686*P_PIm/sqrt(m_pim*m_pim+0.4188*P_PIm*P_PIm - 0.01057)-0.2111)
   //&& (beta_PIm < 0.8*P_PIm/sqrt(m_pim*m_pim+1.2*P_PIm*P_PIm+0.008756) + 0.02963) && (beta_PIm > 0.7686*P_PIm/sqrt(m_pim*m_pim+0.4188*P_PIm*P_PIm - 0.01057)-0.2111)
   
   //&& (beta_PIm < 0.2338*P_PIm/sqrt(m_pim*m_pim+0.0518*P_PIm*P_PIm-0.0187) + 0.0001728) && (beta_PIm > 0.9659*P_PIm/sqrt(m_pim*m_pim+0.9729*P_PIm*P_PIm + 0.008634)-0.0003043)
   //&& (beta_PIm < 0.1717*P_PIm/sqrt(m_pim*m_pim+0.028*P_PIm*P_PIm-0.02) + 0.00023) && (beta_PIm > 0.112*P_PIm/sqrt(m_pim*m_pim+0.015*P_PIm*P_PIm - 0.02))
   
   
 if ((ph_PIm >= 330) && (ph_PIm <=360)){
 
 
 if ((P_PIm > 0) && (P_PIm < 2.)) ph_vs_th_pim[0][int(P_PIm/0.2)]->Fill(th_PIm,ph_PIm-360,1);
 ph_th_pim_all_p[0] -> Fill(th_PIm,ph_PIm-360,1);
 

  if ((ph_PIm > pim_fid_b+360) && (ph_PIm < pim_fid_a+360)){

// if ((th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+75.) || (th_PIm >(11.09+8./(0.472*(P_PIm+0.25)+0.117))+85. )){ 

 
 
 if (( th_PIm >  (11.09+8./(0.472*(P_PIm+0.25)+0.117))+85.) || (th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+75.)) {
 if (( th_PIm >  (11.09+8./(0.472*(P_PIm+0.25)+0.117))+101.) || (th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+96.5)) {
 
  th_vs_p_pim_2[0]->Fill(P_PIm,th_PIm,1.);
   cuts_data = true;
   
   };
   };
   
   
   
//   };//th_vs_p
   };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIm >= 0) && (ph_PIm <=30)){
  
 
  
if ((P_PIm > 0) && (P_PIm < 2.)) ph_vs_th_pim[0][int(P_PIm/0.2)]->Fill(th_PIm,ph_PIm,1);
 ph_th_pim_all_p[0] -> Fill(th_PIm,ph_PIm,1);

 
 if ((ph_PIm > pim_fid_b) && (ph_PIm < pim_fid_a)){
// if ((th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+75.) || (th_PIm >(11.09+8./(0.472*(P_PIm+0.25)+0.117))+85. )){ 

 if (( th_PIm >  (11.09+8./(0.472*(P_PIm+0.25)+0.117))+85.) || (th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+75.)) {
 if (( th_PIm >  (11.09+8./(0.472*(P_PIm+0.25)+0.117))+101.) || (th_PIm < (11.09+8./(0.472*(P_PIm+0.25)+0.117))+96.5)) {
 

 th_vs_p_pim_2[0]->Fill(P_PIm,th_PIm,1.);
cuts_data = true;

   };
   };

//};//th_vs_p
   };//end of the fiducial cut for second part of sector1
   };//end of the second part of sector1
  
    if ((ph_PIm >= 30) && (ph_PIm <=90)) {
    
    
if ((P_PIm > 0) && (P_PIm < 2.)) ph_vs_th_pim[1][int(P_PIm/0.2)]->Fill(th_PIm,ph_PIm-60,1);
ph_th_pim_all_p[1] -> Fill(th_PIm,ph_PIm-60,1);


if ((ph_PIm > pim_fid_b+60) && (ph_PIm < pim_fid_a+60)){
//if ((th_PIm < 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)-2.5) || (th_PIm > 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)+2.)){ 

if ((th_PIm > 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)+2.) || (th_PIm < 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)-2.5)) {

 th_vs_p_pim_2[1]->Fill(P_PIm,th_PIm,1.);
cuts_data = true;

};

//};//th_vs_p
 };//end of the fiducial cut for sector2
};//end of the sector2
    
if ((ph_PIm >= 90) && (ph_PIm <=150)) {


if ((P_PIm > 0) && (P_PIm < 2.)) ph_vs_th_pim[2][int(P_PIm/0.2)]->Fill(th_PIm,ph_PIm-120,1);
ph_th_pim_all_p[2] -> Fill(th_PIm,ph_PIm-120,1);

 
 if ((ph_PIm > pim_fid_b+120) && (ph_PIm < pim_fid_a+120)){

if ((th_PIm > 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.1)+1.81169e-07)+80.) || (th_PIm <  36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.08)+1.81169e-07)+71.)){ 


  th_vs_p_pim_2[2]->Fill(P_PIm,th_PIm,1.);
cuts_data = true;


};//th_vs_p
 };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIm >= 150) && (ph_PIm <=210)){
 
 
if ((P_PIm > 0) && (P_PIm < 2.)) ph_vs_th_pim[3][int(P_PIm/0.2)]->Fill(th_PIm,ph_PIm-180,1);
 ph_th_pim_all_p[3] -> Fill(th_PIm,ph_PIm-180,1);

 
 if ((ph_PIm > pim_fid_b+180) && (ph_PIm < pim_fid_a+180)){



if (( th_PIm > 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.15)+1.81169e-07)+86.) || (th_PIm < 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.15)+1.81169e-07)+77.)) {
if (( th_PIm > 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.35)+1.81169e-07)+66.) || (th_PIm < 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.31)+1.81169e-07)+62.)) {

 th_vs_p_pim_2[3]->Fill(P_PIm,th_PIm,1.);
cuts_data = true;
};
};//th_vs_p


 };//end of the fiducial cut for sector4
 }; //end of the sector4
  
if ((ph_PIm >= 210) && (ph_PIm <=270)) {


if ((P_PIm > 0) && (P_PIm < 2.)) ph_vs_th_pim[4][int(P_PIm/0.2)]->Fill(th_PIm,ph_PIm-240,1);
ph_th_pim_all_p[4] -> Fill(th_PIm,ph_PIm-240,1);


if ((ph_PIm > pim_fid_b+240) && (ph_PIm < pim_fid_a+240)){



if ((th_PIm > 36.152+3.69909e-05/(5.40783e-06*(P_PIm)+1.81169e-07)+4.) || (th_PIm < 36.152+3.69909e-05/(5.40783e-06*(P_PIm)+1.81169e-07)-1.5)) {

th_vs_p_pim_2[4]->Fill(P_PIm,th_PIm,1.);
cuts_data = true;

};//th_vs_p
 };//end of the fiducial cut for sector5
}; //end of the sector5
    
 if ((ph_PIm >= 270) && (ph_PIm <=330)){
 

 
if ((P_PIm > 0) && (P_PIm < 2.)) ph_vs_th_pim[5][int(P_PIm/0.2)]->Fill(th_PIm,ph_PIm-300,1);
 ph_th_pim_all_p[5] -> Fill(th_PIm,ph_PIm-300,1);
 
 
 if ((ph_PIm > pim_fid_b+300) && (ph_PIm < pim_fid_a+300)){
 

 th_vs_p_pim_2[5]->Fill(P_PIm,th_PIm,1.);
cuts_data = true;

//};//th_vs_p
 };//end of the fiducial cut for sector6
 }; //end of the sector6
  
 
  
 
  
  };
  
   };
    return cuts_data;
   };
  
    
   
   
