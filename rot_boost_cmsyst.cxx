#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include "cuts_data.h"
#include "rot_boost_cmsyst.h"
#include "global.h"
#include <iostream>



using namespace std;


 void rot_boost_cmsyst() {


 Float_t m_proton,m_pip,beta;
 Float_t a_gamma, b_gamma, a_beta,b_beta;
 TVector3 Vect3_gamma, Vect3_beta,V3_anti_z(0,0,-1);
TLorentzVector P4_gamma;
 TRotation rot;
 m_proton = 0.938272;
 m_pip = 0.13957;
 
 
 
inv_m_pip_pim = sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg));
inv_m_pip_p = sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg));
inv_m_pim_p = sqrt((P4_PP_reg+P4_PIm_reg)*(P4_PP_reg+P4_PIm_reg));

 ///ROTATION;
P4_gamma = P4_EL - P4_ELP_reg;

/*
th_gamma = acos(P4_gamma[2]/sqrt(P4_gamma[0]*P4_gamma[0]+P4_gamma[1]*P4_gamma[1]+P4_gamma[2]*P4_gamma[2]));

if (P4_gamma[0] != 0.) {
phi_gamma = atan(P4_gamma[1]/P4_gamma[0]);
}
 else {
if(P4_gamma[1] > 0.) phi_gamma = M_PI/2;
if(P4_gamma[1] < 0.) phi_gamma = 3*M_PI/2;
};

if ((P4_gamma[0] < 0.) && (P4_gamma[1] > 0)) phi_gamma = phi_gamma+M_PI;
if (( P4_gamma[0]< 0.) && ( P4_gamma[1]< 0)) phi_gamma = phi_gamma+M_PI;
if ((P4_gamma[0] > 0.) && ( P4_gamma[1]< 0)) phi_gamma = phi_gamma+2*M_PI;

phi_gamma = M_PI;

P4_PP_rot[0] = P4_PP_reg[0]*cos(ph_EL*M_PI/180.) + P4_PP_reg[1]*sin(ph_EL*M_PI/180.);
P4_PP_rot[1] = -P4_PP_reg[0]*sin(ph_EL*M_PI/180.) + P4_PP_reg[1]*cos(ph_EL*M_PI/180.);
P4_PP_rot[2] = P4_PP_reg[2];
P4_PP_rot[3] = P4_PP_reg[3];

P4_PP_rot_1[0] = -P4_PP_rot[0]*cos(th_gamma)*cos(phi_gamma)-P4_PP_rot[1]*sin(phi_gamma)- P4_PP_rot[2]*sin(th_gamma)*cos(phi_gamma);
P4_PP_rot_1[1] = P4_PP_rot[0]*cos(th_gamma)*sin(phi_gamma) - P4_PP_rot[1]*cos(phi_gamma)+P4_PP_rot[2]*sin(th_gamma)*sin(phi_gamma);
P4_PP_rot_1[2] = -P4_PP_rot[0]*sin(th_gamma)+ P4_PP_rot[2]*cos(th_gamma);
P4_PP_rot_1[3] = P4_PP_rot[3];

*/
 TVector3 uz = P4_gamma.Vect().Unit();
 TVector3 ux = (P4_EL.Vect().Cross(P4_ELP_reg.Vect())).Unit();
 ux.Rotate(3.*M_PI/2,uz);
 rot.SetZAxis(uz,ux).Invert();
 P4_PP_reg.Transform(rot);
 P4_PIm_reg.Transform(rot);
 P4_PIp_reg.Transform(rot);
 P4_gamma.Transform(rot);
 
 
 //BOOST
 // E_gamma = (W*W+Q2-m_proton*m_proton)/(2*m_proton);
// E_p_gamma_lab = (W*W+Q2+m_proton*m_proton)/(2*m_proton);
// P_p_gamma_lab = sqrt(E_gamma*E_gamma -Q2);
//beta = P_p_gamma_lab/E_p_gamma_lab;
 beta = sqrt(P4_gamma[3]*P4_gamma[3]+Q2)/(P4_gamma[3]+m_proton);
//gamma = 1./sqrt(1-beta*beta);
 /*
 P4_PP_rot_2_boost[0] = P4_PP_reg[0];
 P4_PP_rot_2_boost[1] = P4_PP_reg[1];
 P4_PP_rot_2_boost[2] = gamma*(P4_PP_reg[2]-beta*P4_PP_reg[3]);
 P4_PP_rot_2_boost[3] = gamma*(P4_PP_reg[3]-beta*P4_PP_reg[2]);

*/

//cout << P4_PIm_reg[0]<<" "<<beta<<" "<<P4_gamma[3]*P4_gamma[3]+Q2<<" "<<P4_gamma.Vect().Mag2()<<"\n";
P4_PP_reg.Boost(0,0,-beta);
P4_PIm_reg.Boost(0,0,-beta);
P4_PIp_reg.Boost(0,0,-beta);
P4_gamma.Boost(0,0,-beta);

//cout<<"1"<< P4_PP_rot_2_boost[0]<< " " << P4_PP_rot_2_boost[1] << " " <<P4_PP_rot_2_boost[2]<< " "<< P4_PP_rot_2_boost[3]<< "\n";

//cout<<"2"<< P4_PP_reg[0]<< " " << P4_PP_reg[1] << " " <<P4_PP_reg[2]<< " "<< P4_PP_reg[3]<< "\n";

theta_PIm_cm =(180./M_PI)*P4_PIm_reg.Theta();
theta_PIp_cm =(180./M_PI)*P4_PIp_reg.Theta();
theta_P_cm =(180./M_PI)*P4_PP_reg.Theta(); //acos((P4_PIm_reg[0]*P4_gamma[0]+P4_PIm_reg[1]*P4_gamma[1]+P4_PIm_reg[2]*P4_gamma[2])/(sqrt(P4_PIm_reg[0]*P4_PIm_reg[0]+P4_PIm_reg[1]*P4_PIm_reg[1]+P4_PIm_reg[2]*P4_PIm_reg[2])*sqrt(P4_gamma[0]*P4_gamma[0]+P4_gamma[1]*P4_gamma[1]+P4_gamma[2]*P4_gamma[2])));

//phi_PIm_cm = atan(P4_PIm_reg[1]/P4_PIm_reg[0]);
//if((P4_PIm_reg[0]>0)&&(P4_PIm_reg[1]<0)) phi_PIm_cm = phi_PIm_cm+2*M_PI;
//if((P4_PIm_reg[0]<0)&&(P4_PIm_reg[1]<0)) phi_PIm_cm = phi_PIm_cm+M_PI;
//if((P4_PIm_reg[0]<0)&&(P4_PIm_reg[1]>0)) phi_PIm_cm = phi_PIm_cm+M_PI;
//if((P4_PIm_reg[0]==0)&&(P4_PIm_reg[1]>0)) phi_PIm_cm = M_PI/2;
//if((P4_PIm_reg[0]==0)&&(P4_PIm_reg[1]<0)) phi_PIm_cm = 3*M_PI/2;


if (P4_PIm_reg.Phi()>0) phi_PIm_cm = (180./M_PI)*P4_PIm_reg.Phi();
if (P4_PIm_reg.Phi()<0) phi_PIm_cm = (180./M_PI)*(P4_PIm_reg.Phi()+2*M_PI);

if (P4_PIp_reg.Phi()>0) phi_PIp_cm = (180./M_PI)*P4_PIp_reg.Phi();
if (P4_PIp_reg.Phi()<0) phi_PIp_cm = (180./M_PI)*(P4_PIp_reg.Phi()+2*M_PI);

if (P4_PP_reg.Phi()>0) phi_P_cm = (180./M_PI)*P4_PP_reg.Phi();
if (P4_PP_reg.Phi()<0) phi_P_cm = (180./M_PI)*(P4_PP_reg.Phi()+2*M_PI);

///1
a_gamma = sqrt(1./(1-pow((P4_PIm_reg.Vect().Unit() * V3_anti_z),2)));
b_gamma = -(P4_PIm_reg.Vect().Unit() * V3_anti_z)*a_gamma;
Vect3_gamma = a_gamma*V3_anti_z +b_gamma*P4_PIm_reg.Vect().Unit();

a_beta = sqrt(1./(1-pow((P4_PIm_reg.Vect().Unit() * P4_PIp_reg.Vect().Unit()),2)));
b_beta = -(P4_PIm_reg.Vect().Unit() * P4_PIp_reg.Vect().Unit())*a_beta;
Vect3_beta = a_beta*P4_PIp_reg.Vect().Unit() + b_beta*P4_PIm_reg.Vect().Unit();

alpha_PPIp_piPIm = (180./M_PI)*acos(Vect3_gamma * Vect3_beta);
if (Vect3_gamma.Cross(Vect3_beta) * P4_PIm_reg.Vect() < 0) alpha_PPIp_piPIm = 360. - alpha_PPIp_piPIm;

///2
a_gamma = sqrt(1./(1-pow((P4_PP_reg.Vect().Unit() * V3_anti_z),2)));
b_gamma = -(P4_PP_reg.Vect().Unit() * V3_anti_z)*a_gamma;
Vect3_gamma = a_gamma*V3_anti_z +b_gamma*P4_PP_reg.Vect().Unit();

a_beta = sqrt(1./(1-pow((P4_PP_reg.Vect().Unit() * P4_PIp_reg.Vect().Unit()),2)));
b_beta = -(P4_PP_reg.Vect().Unit() * P4_PIp_reg.Vect().Unit())*a_beta;
Vect3_beta = a_beta*P4_PIp_reg.Vect().Unit() + b_beta*P4_PP_reg.Vect().Unit();

alpha_PIpPIm_pipf = (180./M_PI)*acos(Vect3_gamma * Vect3_beta);

if (Vect3_gamma.Cross(Vect3_beta) * P4_PP_reg.Vect() < 0) alpha_PIpPIm_pipf = 360. - alpha_PIpPIm_pipf;

///3
a_gamma = sqrt(1./(1-pow((P4_PIp_reg.Vect().Unit() * V3_anti_z),2)));
b_gamma = -(P4_PIp_reg.Vect().Unit() * V3_anti_z)*a_gamma;
Vect3_gamma = a_gamma*V3_anti_z +b_gamma*P4_PIp_reg.Vect().Unit();

a_beta = sqrt(1./(1-pow((P4_PIp_reg.Vect().Unit() * P4_PIm_reg.Vect().Unit()),2)));
b_beta = -(P4_PIp_reg.Vect().Unit() * P4_PIm_reg.Vect().Unit())*a_beta;
Vect3_beta = a_beta*P4_PIm_reg.Vect().Unit() + b_beta*P4_PIp_reg.Vect().Unit();

alpha_PPIm_piPIp = (180./M_PI)*acos(Vect3_gamma * Vect3_beta);

if (Vect3_gamma.Cross(Vect3_beta) * P4_PIp_reg.Vect() < 0) alpha_PPIm_piPIp = 360. - alpha_PPIm_piPIp;
return;
//return (phi_P_cm, phi_PIp_cm, phi_PIm_cm, theta_PIm_cm,theta_PIp_cm, theta_P_cm, alpha_PPIp_piPIm,alpha_PIpPIm_pipf,alpha_PPIm_piPIp); 
 
 };
