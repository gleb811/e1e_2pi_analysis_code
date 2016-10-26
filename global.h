#ifndef GLOBAL_H
#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include <math.h>
#include <TLorentzVector.h>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TText.h"
#include "TStyle.h"
#include "TGObject.h"
#include "TObject.h"
#include "TSystem.h"
#include "TMinuit.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#define GLOBAL_H

#endif

extern Int_t two_pions_block;
extern Int_t npart,segment,sector,indtype,n_P,n_PIp,n_PIm;
extern Short_t pmt_hit;
extern Int_t PdHit_EL,PdHit_PIp,PdHit_P,PdHit_PIm;

extern Float_t m_proton,m_pip,beta;

extern Float_t LiveTime,inclusive,elastic,P_EL,th_EL,z_EL,dc_z_EL,dc_z_PIm,dc_z_PIp,dc_z_P,ph_EL,ECT,nphe,theta_cc,ph_cc;
extern Float_t ph_P,th_P,P_P,beta_P,beta_PIp,beta_PIm;

extern Float_t P_time,PIp_time,PIm_time;
extern Float_t P_dist,PIp_dist,PIm_dist;
extern Float_t P_PIp,ph_PIp,th_PIp,Nphe_PIp;
extern Float_t P_PIm,ph_PIm,th_PIm;
extern Float_t px_fermi,py_fermi,pz_fermi;
extern TLorentzVector  P4_inprot_miss;



extern Float_t sc_x,sc_y,sc_x_p,sc_y_p,sc_x_pip,sc_y_pip,sc_x_pim,sc_y_pim;

extern Double_t theta_PIm_cm,theta_PIp_cm,theta_P_cm,phi_PIm_cm,phi_P_cm,phi_PIp_cm,alpha_PPIp_piPIm, alpha_PIpPIm_pipf,alpha_PPIm_piPIp; 
extern TLorentzVector P4_EL,P4_ELP_reg,P4_PP_reg,P4_PIp_reg,P4_PIm_reg;
extern Double_t inv_m_pip_pim,inv_m_pip_p,inv_m_pim_p;

extern TH1F *padd_rate[6];
extern TH2F *padd_prot_mass[6];

extern TH2F *sc_sector5;
extern TH1F *Wbefore[6],*Wafter[6];
extern TH2F *hist_ectot_sector1_sim,*hist_ectot_sector2_sim,*hist_ectot_sector3_sim,*hist_ectot_sector4_sim,*hist_ectot_sector5_sim,*hist_ectot_sector6_sim;
extern TH2F *hist_ectot_sector1,*hist_ectot_sector2,*hist_ectot_sector3,*hist_ectot_sector4,*hist_ectot_sector5,*hist_ectot_sector6;
extern TH2F *W_incl,*W_2pi,*W_incl_fid,*W_2pi_selection;

extern TH2F *eout_vs_ein_before,*eout_vs_ein_after;

extern TH2F *W_2pi_old, *W_2pi_sim,*W_2pi_old_sim;
extern TH1F *hist_PIm_miss,*hist_PIm_miss_en,*hist_PIp_miss,*hist_PIp_miss_d,*hist_PIp_miss_d_bef, *hist_PIp_miss_en,*hist_P_miss,*hist_P_miss_en, *h_miss_mass_0;
extern TH1D *h_pim_mis_radcorr[22],*h_pim_mis_radcorr_sim[22],*h_pim_mis_radcorr_sim_bckg[22],*h_pr_mom_test[22],*h_pr_mom_sim_test[22],*h_pip_mom_test[22],*h_pip_mom_sim_test[22];
extern TH1D*h_pip_mis_radcorr[22],*h_pip_mis_radcorr_sim[22],*h_pip_mis_radcorr_sim_bckg[22];
extern TH1D*h_pr_mis_radcorr[22],*h_pr_mis_radcorr_sim[22],*h_pr_mis_radcorr_sim_bckg[22];
extern TH1D*h_0_mis_radcorr[22],*h_0_mis_radcorr_sim[22],*h_0_mis_radcorr_sim_bckg[22];
extern TH2F *h_thp_vs_thp[22];
extern TH1F *h_miss_mass_0_d, *h_miss_mass_0_d_sim;
extern TH1F *h_miss_mass_0_d_mmcut, *h_miss_mass_0_d_sim_mmcut;
extern TH1F *h_inprot_miss,*h_inprot_miss_en;
extern TH1F *h_inprot_miss_sim,*h_inprot_miss_en_sim;

extern TH1F *hist_PIm_miss_all_reg_1, *hist_PIm_miss_all_reg_2;
extern TH1F *hist_PIm_miss_all_reg_1_sim, *hist_PIm_miss_all_reg_2_sim;

extern TH1F *hist_miss_en_0,*h_miss_mom_0,*h_miss_mom_0_sim,*h_miss_mom_0_d,*h_miss_mom_0_d_sim; 
extern TH1F *h_miss_mom_0_d_mmcut,*h_miss_mom_0_d_sim_mmcut; 
extern TH1F *h_miss_mom_0_nocut, *h_miss_mom_0_cut_on0, *h_miss_mom_0_cut_onpim, *h_miss_mom_0_nocut_sim, *h_miss_mom_0_cut_on0_sim, *h_miss_mom_0_cut_onpim_sim;

extern TH1F *h_PIp_miss_d_sim,*h_PIp_miss_d_bef_sim,*h_PIm_miss_sim;
extern TH1F *h_PIp_miss_sim, *h_miss_mass_0_sim, *h_miss_en_0_sim, *h_PIm_miss_en_sim, *h_PIp_miss_en_sim, *hist_P_miss_en_sim;
extern TH1F *hist_w_hadr_all_reg, *hist_w_el_all_reg, *hist_w_sim_new_1dim,*hist_w_sim_old_1dim ;
extern Float_t W,Q2;
extern TH2F *ph_vs_th_1,*ph_vs_th_2,*ph_vs_th_3,*ph_vs_th_4,*ph_vs_th_5,*ph_vs_th_6;
extern TH2F *ph_vs_th_1pe[17],*ph_vs_th_2pe[17],*ph_vs_th_3pe[17],*ph_vs_th_4pe[17],*ph_vs_th_5pe[17],*ph_vs_th_6pe[17];
extern TH2F *h_mm_0_vs_npart, *h_mm_pim_vs_npart, *h_mm_pip_vs_npart;
extern TH2F *h_mm_0_vs_npart_sim, *h_mm_pim_vs_npart_sim, *h_mm_pip_vs_npart_sim;

extern TH2F *ph_vs_th_el_sim[6][7];
extern TH2F *ph_vs_th_p_sim[6];
extern TH2F *ph_vs_th_pip_sim[6];
extern TH2F *ph_vs_th_pim_sim[6][5];

extern TH2F  *h_cc_nphe_total_s1,*h_cc_nphe_total_s2,*h_cc_nphe_total_s3,*h_cc_nphe_total_s4,*h_cc_nphe_total_s5,*h_cc_nphe_total_s6;

extern TH2F  *h_cc_nphe_final_s1,*h_cc_nphe_final_s2,*h_cc_nphe_final_s3,*h_cc_nphe_final_s4,*h_cc_nphe_final_s5,*h_cc_nphe_final_s6;


extern TH2F *ph_vs_th_1pe_fid[17],*ph_vs_th_2pe_fid[17],*ph_vs_th_3pe_fid[17],*ph_vs_th_4pe_fid[17],*ph_vs_th_5pe_fid[17],*ph_vs_th_6pe_fid[17];

extern TH2F *hist_sector1,*hist_sector2,*hist_sector3,*hist_sector4,*hist_sector5,*hist_sector6;
extern TH2F *hist_nphe_sector1,*hist_nphe_sector2,*hist_nphe_sector3,*hist_nphe_sector4,*hist_nphe_sector5,*hist_nphe_sector6;
extern TH1F *nphe_sector1,*nphe_sector2,*nphe_sector3,*nphe_sector4,*nphe_sector5,*nphe_sector6;
extern TH1F *nphe_sector1_after,*nphe_sector2_after,*nphe_sector3_after,*nphe_sector4_after,*nphe_sector5_after,*nphe_sector6_after;
extern TH2F  *ph_vs_th_p_1,*ph_vs_th_p_2,*ph_vs_th_p_3,*ph_vs_th_p_4,*ph_vs_th_p_5,*ph_vs_th_p_6;
extern TH2F  *ph_vs_th_p_1_w,*ph_vs_th_p_2_w,*ph_vs_th_p_3_w,*ph_vs_th_p_4_w,*ph_vs_th_p_5_w,*ph_vs_th_p_6_w;
extern TH2F  *ph_th_p_1[5],*ph_th_p_2[5],*ph_th_p_3[5],*ph_th_p_4[5],*ph_th_p_5[5],*ph_th_p_6[5];
extern TH1F *hist_z_el_1, *hist_z_el_2, *hist_z_el_3, *hist_z_el_4, *hist_z_el_5, *hist_z_el_6;
extern TH1F *hist_z_el_1_empty, *hist_z_el_2_empty, *hist_z_el_3_empty, *hist_z_el_4_empty, *hist_z_el_5_empty, *hist_z_el_6_empty;
extern TH1F *hist_z_el_1_sim_1, *hist_z_el_2_sim_1, *hist_z_el_3_sim_1, *hist_z_el_4_sim_1, *hist_z_el_5_sim_1, *hist_z_el_6_sim_1;
extern TH1F *hist_z_el_1_sim_2, *hist_z_el_2_sim_2, *hist_z_el_3_sim_2, *hist_z_el_4_sim_2, *hist_z_el_5_sim_2, *hist_z_el_6_sim_2;

extern TH1F *h_mixed_prod,*h_mixed_prod_2;
extern TH1F *h_mixed_prod_sim,*h_mixed_prod_2_sim;

extern TH2F  *ph_th_pip_1[6],*ph_th_pip_2[6],*ph_th_pip_3[6],*ph_th_pip_4[6],*ph_th_pip_5[6],*ph_th_pip_6[6];

extern TH2F  *ph_vs_th_pip_1,*ph_vs_th_pip_2,*ph_vs_th_pip_3,*ph_vs_th_pip_4,*ph_vs_th_pip_5,*ph_vs_th_pip_6;

extern TH2F  *norm_nphe_s1,*norm_nphe_s2,*norm_nphe_s3,*norm_nphe_s4,*norm_nphe_s5,*norm_nphe_s6;

extern TH2F  *avrg_nphe_sector1,*avrg_nphe_sector2,*avrg_nphe_sector3,*avrg_nphe_sector4,*avrg_nphe_sector5,*avrg_nphe_sector6;
extern TFile *outFile;
extern TH1F *hist_ltime,*hist_ltime_1d,*hist_n_incl,*hist_n_incl_1d,*hist_n_twopions,*hist_n_twopions_1d,*hist_n_elast,*hist_n_elast_1d,*hist_twopi_cut;
 
extern TH2F  *W_2pi_fid_p;
extern TH2F  *h_delta_w_vs_w_old_data;
extern TH2F  *h_delta_w_vs_w_old_data_mmcut;
extern TH2F  *h_delta_w_vs_w_old_gen;
extern TH2F  *h_delta_w_vs_w_old_rec;

 extern TH2F *beta_vs_p_p[6][48];
 extern TH2F *beta_vs_p_pim[6][48];
  extern TH2F *beta_vs_p_pip[6][48];
  
  extern TH2F *beta_vs_p_p_sim[6][48];
 extern TH2F *beta_vs_p_pim_sim[6][48];
  extern TH2F *beta_vs_p_pip_sim[6][48]; 
  
  
  extern TH2F *time_pip[6][48];
  extern TH2F *time_pim[6][48];
  extern TH2F *time_p[6][48];
  
  extern TH2F *time_pip_sim[6][48];
  extern TH2F *time_pim_sim[6][48];
  extern TH2F *time_p_sim[6][48];  
  
  
  extern TH2F *ph_vs_th_pim[6][10];

  //photoelectrons
  extern TH1F *ph_el_left[6][20];
  extern TH1F *ph_el_both[6][20];
  extern TH1F *ph_el_right[6][20];
  
  extern TH2F *ph_th_pim_all_p[6];

extern TH1F *h_inv_m_pip_pim,*h_inv_m_pip_p,*h_inv_m_pim_p; 
extern TH1F *h_inv_m_pip_pim_bin, *h_inv_m_pip_p_bin, *h_inv_m_pim_p_bin, *h_theta_PIm_cm, *h_theta_PIp_cm, *h_theta_P_cm, *h_phi_PIm_cm, *h_phi_PIp_cm, *h_phi_P_cm, *h_alpha_PIpPIm_pipf, *h_alpha_PPIp_piPIm, *h_alpha_PPIm_piPIp;


extern THnSparseD *h_5dim_0_1[12][21];
extern THnSparseD *h_5dim_0_2[12][21];
extern THnSparseD *h_5dim_0_3[12][21];

extern THnSparseD *h_5dim_pim_1[12][21];
extern THnSparseD *h_5dim_pim_2[12][21];
extern THnSparseD *h_5dim_pim_3[12][21];

extern THnSparseD *h_5dim_pip_1[12][21];
extern THnSparseD *h_5dim_pip_2[12][21];
extern THnSparseD *h_5dim_pip_3[12][21];

extern THnSparseD *h_5dim_pr_1[12][21];
extern THnSparseD *h_5dim_pr_2[12][21];
extern THnSparseD *h_5dim_pr_3[12][21];



extern THnSparseD *h_5dim_0_1_empty[12][21];
extern THnSparseD *h_5dim_0_2_empty[12][21];
extern THnSparseD *h_5dim_0_3_empty[12][21];

extern THnSparseD *h_5dim_pim_1_empty[12][21];
extern THnSparseD *h_5dim_pim_2_empty[12][21];
extern THnSparseD *h_5dim_pim_3_empty[12][21];

extern THnSparseD *h_5dim_pip_1_empty[12][21];
extern THnSparseD *h_5dim_pip_2_empty[12][21];
extern THnSparseD *h_5dim_pip_3_empty[12][21];

extern THnSparseD *h_5dim_pr_1_empty[12][21];
extern THnSparseD *h_5dim_pr_2_empty[12][21];
extern THnSparseD *h_5dim_pr_3_empty[12][21];





extern THnSparseD *h_5dim_0_1_sim_1[12][21];
extern THnSparseD *h_5dim_0_2_sim_1[12][21];
extern THnSparseD *h_5dim_0_3_sim_1[12][21];

extern THnSparseD *h_5dim_pim_1_sim_1[12][21];
extern THnSparseD *h_5dim_pim_2_sim_1[12][21];
extern THnSparseD *h_5dim_pim_3_sim_1[12][21];

extern THnSparseD *h_5dim_pip_1_sim_1[12][21];
extern THnSparseD *h_5dim_pip_2_sim_1[12][21];
extern THnSparseD *h_5dim_pip_3_sim_1[12][21];

extern THnSparseD *h_5dim_pr_1_sim_1[12][21];
extern THnSparseD *h_5dim_pr_2_sim_1[12][21];
extern THnSparseD *h_5dim_pr_3_sim_1[12][21];






extern THnSparseD *h_5dim_1_sim_2[12][21];
extern THnSparseD *h_5dim_2_sim_2[12][21];
extern THnSparseD *h_5dim_3_sim_2[12][21];

extern THnSparseD *h_5dim_excl_1[9];
extern THnSparseD *h_5dim_excl_2[9];
extern THnSparseD *h_5dim_excl_3[9];

extern THnSparseD *h_5dim_excl_1_empty[9];
extern THnSparseD *h_5dim_excl_2_empty[9];
extern THnSparseD *h_5dim_excl_3_empty[9];

extern THnSparseD *h_5dim_1_sim_excl_1[9];
extern THnSparseD *h_5dim_2_sim_excl_1[9];
extern THnSparseD *h_5dim_3_sim_excl_1[9];

extern THnSparseD *h_5dim_1_sim_excl_2[9];
extern THnSparseD *h_5dim_2_sim_excl_2[9];
extern THnSparseD *h_5dim_3_sim_excl_2[9];


extern TH1D *h1prj_inv_m_pip_p[12][21]; 
extern TH1D *h1prj_inv_m_pip_pim[12][21]; 
extern TH1D *h1prj_th_P[12][21]; 
extern TH1D *h1prj_phi_P[12][21]; 
extern TH1D *h1prj_alpha_PIpPIm_pipf[12][21];

extern TH1D *h2prj_inv_m_pip_p[12][21];
extern TH1D *h2prj_inv_m_pip_pim[12][21];
extern TH1D *h2prj_th_PIm[12][21];
extern TH1D *h2prj_phi_PIm[12][21];
extern TH1D *h2prj_alpha_PPIp_piPIm[12][21];

extern TH1D *h3prj_inv_m_pim_p[12][21];
extern TH1D *h3prj_inv_m_pip_pim[12][21];
extern TH1D *h3prj_th_PIp[12][21];
extern TH1D *h3prj_phi_PIp[12][21];
extern TH1D *h3prj_alpha_PPIm_piPIp[12][21];



extern TH1D *h1prj_inv_m_pip_p_sim[12][21]; 
extern TH1D *h1prj_inv_m_pip_pim_sim[12][21]; 
extern TH1D *h1prj_th_P_sim[12][21]; 
extern TH1D *h1prj_phi_P_sim[12][21]; 
extern TH1D *h1prj_alpha_PIpPIm_pipf_sim[12][21];

extern TH1D *h2prj_inv_m_pip_p_sim[12][21];
extern TH1D *h2prj_inv_m_pip_pim_sim[12][21];
extern TH1D *h2prj_th_PIm_sim[12][21];
extern TH1D *h2prj_phi_PIm_sim[12][21];
extern TH1D *h2prj_alpha_PPIp_piPIm_sim[12][21];

extern TH1D *h3prj_inv_m_pim_p_sim[12][21];
extern TH1D *h3prj_inv_m_pip_pim_sim[12][21];
extern TH1D *h3prj_th_PIp_sim[12][21];
extern TH1D *h3prj_phi_PIp_sim[12][21];
extern TH1D *h3prj_alpha_PPIm_piPIp_sim[12][21];


extern TH1D *h1prj_inv_m_pip_p_sim_2[12][21]; 
extern TH1D *h1prj_inv_m_pip_pim_sim_2[12][21]; 
extern TH1D *h1prj_th_P_sim_2[12][21]; 
extern TH1D *h1prj_phi_P_sim_2[12][21]; 
extern TH1D *h1prj_alpha_PIpPIm_pipf_sim_2[12][21];

extern TH1D *h2prj_inv_m_pip_p_sim_2[12][21];
extern TH1D *h2prj_inv_m_pip_pim_sim_2[12][21];
extern TH1D *h2prj_th_PIm_sim_2[12][21];
extern TH1D *h2prj_phi_PIm_sim_2[12][21];
extern TH1D *h2prj_alpha_PPIp_piPIm_sim_2[12][21];

extern TH1D *h3prj_inv_m_pim_p_sim_2[12][21];
extern TH1D *h3prj_inv_m_pip_pim_sim_2[12][21];
extern TH1D *h3prj_th_PIp_sim_2[12][21];
extern TH1D *h3prj_phi_PIp_sim_2[12][21];
extern TH1D *h3prj_alpha_PPIm_piPIp_sim_2[12][21];





extern TH1F *h_pim_mm_q2_w[12][21];
extern TH1F *h_pim_mm_q2_w_sim[12][21];

extern TH1F *h_0_mm_q2_w[21];
extern TH1F *h_0_mm_q2_w_sim[21];

extern TH1F *h_0_d_mm_q2_w[21];
extern TH1F *h_0_d_mm_q2_w_sim[21];

extern TH1F *h_0_mmom_q2_w[21];
extern TH1F *h_0_mmom_q2_w_sim[21];

extern TH1F *h_0_d_mmom_q2_w[21];
extern TH1F *h_0_d_mmom_q2_w_sim[21];





extern TH1D *h_cos_th;

extern TH1D *h_w_int[12];
extern TH2F *th_cc_vs_seg[6];
extern TH2F *th_cc_vs_seg_sim[6];
extern TH2F *th_vs_p_e_1[6], *th_vs_p_e_2[6], *th_vs_p_p_1[6],*th_vs_p_p_2[6],*th_vs_p_pip_1[6],*th_vs_p_pip_2[6],*th_vs_p_pim_1[6],*th_vs_p_pim_2[6];

extern TH2F *th_vs_p_e_1_sim[6], *th_vs_p_e_2_sim[6], *th_vs_p_p_1_sim[6],*th_vs_p_p_2_sim[6],*th_vs_p_pip_1_sim[6],*th_vs_p_pip_2_sim[6],*th_vs_p_pim_1_sim[6],*th_vs_p_pim_2_sim[6];

int global();
