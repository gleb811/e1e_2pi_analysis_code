#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include <math.h>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TGObject.h"
#include "TObject.h"
#include "TSystem.h"
#include "TMinuit.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include "global.h"
#include <stdio.h>
#include <dlfcn.h>
#include <sstream>
#include <TLorentzVector.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring> 
#include <cstdlib>

 using namespace std; 
 
Int_t npart,segment,sector,indtype,n_P,n_PIp,n_PIm;
Short_t pmt_hit;
Int_t PdHit_EL,PdHit_PIp,PdHit_P,PdHit_PIm; 
 
 Float_t m_proton,m_pip,beta;
 
 Float_t LiveTime,inclusive,elastic,P_EL,th_EL,z_EL,dc_z_EL,dc_z_PIm,dc_z_PIp,dc_z_P,ph_EL,ECT,nphe,theta_cc,ph_cc;
 

Float_t W,Q2;
Float_t ph_P,th_P,P_P,beta_P,beta_PIp,beta_PIm;
//Float_t beta_P_time,beta_PIp_time,beta_PIm_time;
Float_t P_time,PIp_time,PIm_time;
Float_t P_dist,PIp_dist,PIm_dist;
Float_t P_PIp,ph_PIp,th_PIp,Nphe_PIp;
Float_t P_PIm,ph_PIm,th_PIm;
 Float_t px_fermi,py_fermi,pz_fermi;
TLorentzVector  P4_inprot_miss;

Float_t sc_x,sc_y,sc_x_p,sc_y_p,sc_x_pip,sc_y_pip,sc_x_pim,sc_y_pim;
Double_t theta_PIm_cm,theta_PIp_cm,theta_P_cm,phi_PIm_cm,phi_P_cm,phi_PIp_cm,alpha_PPIp_piPIm, alpha_PIpPIm_pipf,alpha_PPIm_piPIp; 
TLorentzVector P4_EL,P4_ELP_reg,P4_PP_reg,P4_PIp_reg,P4_PIm_reg;
Double_t inv_m_pip_pim,inv_m_pip_p,inv_m_pim_p;
//booking electron histograms

THnSparse* whatever;

/*TH2F *hist_sector1 = new TH2F("cc_sector1","cc_sector1",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector2 = new TH2F("cc_sector2","cc_sector2",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector3 = new TH2F("cc_sector3","cc_sector3",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector4 = new TH2F("cc_sector4","cc_sector4",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector5 = new TH2F("cc_sector5","cc_sector5",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector6 = new TH2F("cc_sector6","cc_sector6",200, -5., 55.,200., -25., 25.);*/
TH2F *sc_sector5 = new TH2F("sc_sector5","sc_sector5",100, -180., 180., 100, 50., 400.);

TH1F *nphe_sector1 = new TH1F("nphe_sector1","nphe_sector1",501, -1., 500.);
TH1F *nphe_sector1_after = new TH1F("nphe_sector1_after","nphe_sector1_after",501, -1., 500.);
TH1F *nphe_sector2 = new TH1F("nphe_sector2","nphe_sector2",501, -1., 500.);
TH1F *nphe_sector2_after = new TH1F("nphe_sector2_after","nphe_sector2_after",501, -1., 500.);
TH1F *nphe_sector3 = new TH1F("nphe_sector3","nphe_sector3",501, -1., 500.);
TH1F *nphe_sector3_after = new TH1F("nphe_sector3_after","nphe_sector3_after",501, -1., 500.);
TH1F *nphe_sector4 = new TH1F("nphe_sector4","nphe_sector4",501, -1., 500.);
TH1F *nphe_sector4_after = new TH1F("nphe_sector4_after","nphe_sector4_after",501, -1., 500.);
TH1F *nphe_sector5 = new TH1F("nphe_sector5","nphe_sector5",501, -1., 500.);
TH1F *nphe_sector5_after = new TH1F("nphe_sector5_after","nphe_sector5_after",501, -1., 500.);
TH1F *nphe_sector6 = new TH1F("nphe_sector6","nphe_sector6",501, -1., 500.);
TH1F *nphe_sector6_after = new TH1F("nphe_sector6_after","nphe_sector6_after",501, -1., 500.);


/*TH2F *hist_nphe_sector1 = new TH2F("cc_nphe_sector1","cc_nphe_sector1",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector2 = new TH2F("cc_nphe_sector2","cc_nphe_sector2",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector3 = new TH2F("cc_nphe_sector3","cc_nphe_sector3",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector4 = new TH2F("cc_nphe_sector4","cc_nphe_sector4",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector5 = new TH2F("cc_nphe_sector5","cc_nphe_sector5",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector6 = new TH2F("cc_nphe_sector6","cc_nphe_sector6",200, -5., 55.,200., -25., 25.);*/

TH2F  *avrg_nphe_sector1,*avrg_nphe_sector2,*avrg_nphe_sector3,*avrg_nphe_sector4,*avrg_nphe_sector5,*avrg_nphe_sector6;

TH2F  *norm_nphe_s1,*norm_nphe_s2,*norm_nphe_s3,*norm_nphe_s4,*norm_nphe_s5,*norm_nphe_s6;





TH2F *h_cc_nphe_total_s1 = new TH2F("h_cc_nphe_total_s1","h_cc_nphe_total_s1",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s2 = new TH2F("h_cc_nphe_total_s2","h_cc_nphe_total_s2",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s3 = new TH2F("h_cc_nphe_total_s3","h_cc_nphe_total_s3",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s4 = new TH2F("h_cc_nphe_total_s4","h_cc_nphe_total_s4",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s5 = new TH2F("h_cc_nphe_total_s5","h_cc_nphe_total_s5",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s6 = new TH2F("h_cc_nphe_total_s6","h_cc_nphe_total_s6",200, -5., 55.,200., -25., 25.);


TH2F *h_cc_nphe_final_s1 = new TH2F("h_cc_nphe_final_s1","h_cc_nphe_final_s1",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s2 = new TH2F("h_cc_nphe_final_s2","h_cc_nphe_final_s2",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s3 = new TH2F("h_cc_nphe_final_s3","h_cc_nphe_final_s3",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s4 = new TH2F("h_cc_nphe_final_s4","h_cc_nphe_final_s4",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s5 = new TH2F("h_cc_nphe_final_s5","h_cc_nphe_final_s5",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s6 = new TH2F("h_cc_nphe_final_s6","h_cc_nphe_final_s6",200, -5., 55.,200., -25., 25.);




TH2F *hist_ectot_sector1 = new TH2F("ectot_sector1","ectot_sector1",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector2 = new TH2F("ectot_sector2","ectot_sector2",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector3 = new TH2F("ectot_sector3","ectot_sector3",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector4 = new TH2F("ectot_sector4","ectot_sector4",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector5 = new TH2F("ectot_sector5","ectot_sector5",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector6 = new TH2F("ectot_sector6","ectot_sector6",220, -0.005, 2.195,80, -0.005, 0.795);

TH2F *hist_ectot_sector1_sim = new TH2F("ectot_sector1_sim","ectot_sector1_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector2_sim = new TH2F("ectot_sector2_sim","ectot_sector2_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector3_sim = new TH2F("ectot_sector3_sim","ectot_sector3_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector4_sim = new TH2F("ectot_sector4_sim","ectot_sector4_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector5_sim = new TH2F("ectot_sector5_sim","ectot_sector5_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector6_sim = new TH2F("ectot_sector6_sim","ectot_sector6_sim",440, -0.005, 2.195,160, -0.005, 0.795);


TH2F *W_incl = new TH2F("W_incl","W_incl",200, 0.8, 1.9,200., 0.2, 1.9);
TH2F *W_2pi = new TH2F("W_2pi_new","W_2pi_new",210., 1.2, 1.9,240., 0.2, 1.4);
TH2F *W_2pi_old = new TH2F("W_2pi_old","W_2pi_old",210., 1.2, 1.9,240., 0.2, 1.4);

TH2F *W_2pi_sim= new TH2F("W_2pi_new_sim","W_2pi_new_sim",210., 1.2, 1.9,240., 0.2, 1.4);
TH2F *W_2pi_old_sim = new TH2F("W_2pi_old_sim","W_2pi_old_sim",210., 1.2, 1.9,240., 0.2, 1.4);


TH2F *W_incl_fid = new TH2F("W_incl_fid","W_incl_fid",200, 0.8, 1.9,200., 0.2, 1.9);
TH2F *W_2pi_selection = new TH2F("W_2pi_selection","W_2pi_selection",28, 1.2, 1.9,24., 0.2, 1.4);
TH2F *W_2pi_fid_p = new TH2F("W_2pi_fid_p","W_2pi_fid_p",200, 1.2, 1.9,200., 0.2, 1.4);

TH1F *h_inprot_miss = new TH1F("h_inprot_miss","h_inprot_miss",800, -2., 2.);
TH1F *h_inprot_miss_en = new TH1F("h_inprot_miss_en","h_inprot_miss_en",800, -2., 2.);

TH1F *h_inprot_miss_sim = new TH1F("h_inprot_miss_sim","h_inprot_miss_sim",800, -2., 2.);
TH1F *h_inprot_miss_en_sim = new TH1F("h_inprot_miss_en_sim","h_inprot_miss_en_sim",800, -2., 2.);



TH1F *hist_PIm_miss = new TH1F("hist_PIm_miss","hist_PIm_miss",2000, -2., 2.);
TH1F *h_PIm_miss_sim = new TH1F("h_PIm_miss_sim","h_PIm_miss_sim",2000, -2., 2.);

TH1D *h_pim_mis_radcorr[22];
TH1D *h_pim_mis_radcorr_sim[22];
TH1D *h_pim_mis_radcorr_sim_bckg[22];
TH1D *h_pip_mis_radcorr[22];
TH1D *h_pip_mis_radcorr_sim[22];
TH1D *h_pip_mis_radcorr_sim_bckg[22];
TH1D *h_pr_mis_radcorr[22];
TH1D *h_pr_mis_radcorr_sim[22];
TH1D *h_pr_mis_radcorr_sim_bckg[22];
TH1D *h_0_mis_radcorr[22];
TH1D *h_0_mis_radcorr_sim[22];
TH1D *h_0_mis_radcorr_sim_bckg[22];
TH1D *h_pr_mom_test[22];
TH1D *h_pr_mom_sim_test[22];
TH1D *h_pip_mom_test[22];
TH1D *h_pip_mom_sim_test[22];

TH2F *h_thp_vs_thp[22];

TH1F *hist_w_hadr_all_reg = new TH1F("hist_w_hadr_all_reg","hist_w_hadr_all_reg",210, 1.2, 1.9);
TH1F *hist_w_el_all_reg = new TH1F("hist_w_el_all_reg","hist_w_el_all_reg",210, 1.2, 1.9);
TH1F *hist_w_sim_old_1dim = new TH1F("hist_w_sim_old_1dim","hist_w_sim_old_1dim",210, 1.2, 1.9);
TH1F *hist_w_sim_new_1dim = new TH1F("hist_w_sim_new_1dim","hist_w_sim_new_1dim",210, 1.2, 1.9);


TH1F *hist_PIm_miss_all_reg_1 = new TH1F("hist_PIm_miss_all_reg_1","hist_PIm_miss_all_reg_1",400, -0.2, 0.2);
TH1F *hist_PIm_miss_all_reg_2 = new TH1F("hist_PIm_miss_all_reg_2","hist_PIm_miss_all_reg_2",400, -0.2, 0.2);

TH1F *hist_PIm_miss_all_reg_1_sim = new TH1F("hist_PIm_miss_all_reg_1_sim","hist_PIm_miss_all_reg_1_sim",400, -0.2, 0.2);
TH1F *hist_PIm_miss_all_reg_2_sim = new TH1F("hist_PIm_miss_all_reg_2_sim","hist_PIm_miss_all_reg_2_sim",400, -0.2, 0.2);


TH1F *h_PIp_miss_sim = new TH1F("h_PIp_miss_sim","h_PIp_miss_sim",500, -0.7, 0.3);
TH1F *h_miss_mass_0_sim = new TH1F("h_miss_mass_0_sim","h_miss_mass_0_sim",400, -0.2, 0.2);
TH1F *h_miss_en_0_sim = new TH1F("h_miss_en_0_sim","h_miss_en_0_sim",400, -2, 2);


TH1F *h_miss_mass_0 = new TH1F("h_miss_mass_0","h_miss_mass_0",400, -0.2, 0.2);
TH1F *h_PIm_miss_en_sim = new TH1F("h_PIm_miss_en_sim","h_PIm_miss_en_sim",400, -2, 2);
TH1F *h_PIp_miss_en_sim = new TH1F("h_PIp_miss_en_sim","h_PIp_miss_en_sim",400, -2, 2);

TH1F *hist_PIm_miss_en = new TH1F("hist_PIm_miss_en","hist_PIm_miss_en",400, -2, 2);

TH1F *h_miss_mass_0_d = new TH1F("h_miss_mass_0_d","h_miss_mass_0_d",2000, -0., 2.);
TH1F *h_miss_mass_0_d_sim = new TH1F("h_miss_mass_0_d_sim","h_miss_mass_0_d_sim",2000, -0., 2.);


TH1F *h_miss_mass_0_d_mmcut = new TH1F("h_miss_mass_0_d_mmcut","h_miss_mass_0_d_mmcut",2000, -0., 2.);
TH1F *h_miss_mass_0_d_sim_mmcut = new TH1F("h_miss_mass_0_d_sim_mmcut","h_miss_mass_0_d_sim_mmcut",2000, -0., 2.);

TH1F *h_miss_mom_0 = new TH1F("h_miss_mom_0","h_miss_mom_0",400, -0.1, 0.9);
TH1F *h_miss_mom_0_sim = new TH1F("h_miss_mom_0_sim","h_miss_mom_0_sim",400, -0.1, 0.9);


TH1F *h_miss_mom_0_nocut = new TH1F("h_miss_mom_0_nocut","h_miss_mom_0_nocut",400, -0.1, 0.9);
TH1F *h_miss_mom_0_cut_on0 = new TH1F("h_miss_mom_0_cut_on0","h_miss_mom_0_cut_on0",400, -0.1, 0.9);
TH1F *h_miss_mom_0_cut_onpim = new TH1F("h_miss_mom_0_cut_onpim","h_miss_mom_0_cut_onpim",400, -0.1, 0.9);

TH1F *h_miss_mom_0_nocut_sim = new TH1F("h_miss_mom_0_nocut_sim","h_miss_mom_0_nocut_sim",400, -0.1, 0.9);
TH1F *h_miss_mom_0_cut_on0_sim = new TH1F("h_miss_mom_0_cut_on0_sim","h_miss_mom_0_cut_on0_sim",400, -0.1, 0.9);
TH1F *h_miss_mom_0_cut_onpim_sim = new TH1F("h_miss_mom_0_cut_onpim_sim","h_miss_mom_0_cut_onpim_sim",400, -0.1, 0.9);


TH1F *h_miss_mom_0_d = new TH1F("h_miss_mom_0_d","h_miss_mom_0_d",400, -0.1, 0.9);
TH1F *h_miss_mom_0_d_sim = new TH1F("h_miss_mom_0_d_sim","h_miss_mom_0_d_sim",400, -0.1, 0.9);

TH1F *h_miss_mom_0_d_mmcut = new TH1F("h_miss_mom_0_d_mmcut","h_miss_mom_0_d_mmcut",400, -0.1, 0.9);
TH1F *h_miss_mom_0_d_sim_mmcut = new TH1F("h_miss_mom_0_d_sim_mmcut","h_miss_mom_0_d_sim_mmcut",400, -0.1, 0.9);
//extern 


TH1F *hist_PIp_miss = new TH1F("hist_PIp_miss","hist_PIp_miss",500, -0.7, 0.3);
TH1F *hist_PIp_miss_d = new TH1F("hist_PIp_miss_d","hist_PIp_miss_d",500, -1., 3.);
TH1F *hist_PIp_miss_d_bef = new TH1F("hist_PIp_miss_d_bef","hist_PIp_miss_d_bef",500, -1., 3.);
TH1F *h_PIp_miss_d_sim = new TH1F("h_PIp_miss_d_sim","h_PIp_miss_d_sim",500, -1., 3.);
TH1F *h_PIp_miss_d_bef_sim = new TH1F("h_PIp_miss_d_bef_sim","h_PIp_miss_d_bef_sim",500, -1., 3.);

TH1F *hist_PIp_miss_en = new TH1F("hist_PIp_miss_en","hist_PIp_miss_en",400, -2, 2);
TH1F *hist_P_miss = new TH1F("hist_P_miss","hist_P_miss",400, -0.5, 1.5);
TH1F *hist_P_miss_en = new TH1F("hist_P_miss_en","hist_P_miss_en",400, -2, 2);
TH1F *hist_P_miss_en_sim = new TH1F("hist_P_miss_en_sim","hist_P_miss_en_sim",400, -2, 2);
TH1F *hist_miss_en_0 = new TH1F("hist_miss_en_0","hist_miss_en_0",400, -2, 2);
 

TH2F *h_mm_0_vs_npart = new TH2F("mm_0_vs_npart","mm_0_vs_npart",200, 0,10,300, -0.4, 0.2);
TH2F *h_mm_pim_vs_npart = new TH2F("mm_pim_vs_npart","mm_pim_vs_npart",200, 0,10,300, -0.4, 0.2);
TH2F *h_mm_pip_vs_npart = new TH2F("mm_pip_vs_npart","mm_pip_vs_npart",200, 0,10,300, -0.4, 0.2);

TH2F *h_mm_0_vs_npart_sim = new TH2F("mm_0_vs_npart_sim","mm_0_vs_npart_sim",200, 0,10,300, -0.4, 0.2);
TH2F *h_mm_pim_vs_npart_sim = new TH2F("mm_pim_vs_npart_sim","mm_pim_vs_npart_sim",200, 0,10,300, -0.4, 0.2);
TH2F *h_mm_pip_vs_npart_sim = new TH2F("mm_pip_vs_npart_sim","mm_pip_vs_npart_sim",200, 0,10,300, -0.4, 0.2);

TH1F *h_mixed_prod = new TH1F("h_mixed_prod","h_mixed_prod",400, -0.5, 0.5);
TH1F *h_mixed_prod_2 = new TH1F("h_mixed_prod_2","h_mixed_prod_2",400, -0.5, 0.5);

TH1F *h_mixed_prod_sim = new TH1F("h_mixed_prod_sim","h_mixed_prod_sim",400, -0.5, 0.5);
TH1F *h_mixed_prod_2_sim = new TH1F("h_mixed_prod_2_sim","h_mixed_prod_2_sim",400, -0.5, 0.5);



TH1F *hist_z_el_1 = new TH1F("hist_z_el_1","hist_z_el_1",400, -10.,10.);
TH1F *hist_z_el_2 = new TH1F("hist_z_el_2","hist_z_el_2",400, -10.,10.);
TH1F *hist_z_el_3 = new TH1F("hist_z_el_3","hist_z_el_3",400, -10.,10.);
TH1F *hist_z_el_4 = new TH1F("hist_z_el_4","hist_z_el_4",400, -10.,10.);
TH1F *hist_z_el_5 = new TH1F("hist_z_el_5","hist_z_el_5",400, -10.,10.);
TH1F *hist_z_el_6 = new TH1F("hist_z_el_6","hist_z_el_6",400, -10.,10.);


TH1F *hist_z_el_1_empty = new TH1F("hist_z_el_1_empty","hist_z_el_1_empty",400, -10.,10.);
TH1F *hist_z_el_2_empty = new TH1F("hist_z_el_2_empty","hist_z_el_2_empty",400, -10.,10.);
TH1F *hist_z_el_3_empty = new TH1F("hist_z_el_3_empty","hist_z_el_3_empty",400, -10.,10.);
TH1F *hist_z_el_4_empty = new TH1F("hist_z_el_4_empty","hist_z_el_4_empty",400, -10.,10.);
TH1F *hist_z_el_5_empty = new TH1F("hist_z_el_5_empty","hist_z_el_5_empty",400, -10.,10.);
TH1F *hist_z_el_6_empty = new TH1F("hist_z_el_6_empty","hist_z_el_6_empty",400, -10.,10.);


TH1F *hist_z_el_1_sim_1 = new TH1F("hist_z_el_1_sim_1","hist_z_el_1_sim_1",400, -10.,10.);
TH1F *hist_z_el_2_sim_1 = new TH1F("hist_z_el_2_sim_1","hist_z_el_2_sim_1",400, -10.,10.);
TH1F *hist_z_el_3_sim_1 = new TH1F("hist_z_el_3_sim_1","hist_z_el_3_sim_1",400, -10.,10.);
TH1F *hist_z_el_4_sim_1 = new TH1F("hist_z_el_4_sim_1","hist_z_el_4_sim_1",400, -10.,10.);
TH1F *hist_z_el_5_sim_1 = new TH1F("hist_z_el_5_sim_1","hist_z_el_5_sim_1",400, -10.,10.);
TH1F *hist_z_el_6_sim_1 = new TH1F("hist_z_el_6_sim_1","hist_z_el_6_sim_1",400, -10.,10.);

TH1F *hist_z_el_1_sim_2 = new TH1F("hist_z_el_1_sim_2","hist_z_el_1_sim_2",400, -10.,10.);
TH1F *hist_z_el_2_sim_2 = new TH1F("hist_z_el_2_sim_2","hist_z_el_2_sim_2",400, -10.,10.);
TH1F *hist_z_el_3_sim_2 = new TH1F("hist_z_el_3_sim_2","hist_z_el_3_sim_2",400, -10.,10.);
TH1F *hist_z_el_4_sim_2 = new TH1F("hist_z_el_4_sim_2","hist_z_el_4_sim_2",400, -10.,10.);
TH1F *hist_z_el_5_sim_2 = new TH1F("hist_z_el_5_sim_2","hist_z_el_5_sim_2",400, -10.,10.);
TH1F *hist_z_el_6_sim_2 = new TH1F("hist_z_el_6_sim_2","hist_z_el_6_sim_2",400, -10.,10.);




//extern TH1F *inv_m_pip_pim_bin, *inv_m_pip_p_bin, *inv_m_pim_p_bin, *theta_PIm_cm, *theta_PIp_cm, *theta_P_cm, *phi_PIm_cm, *phi_PIp_cm, *phi_P_cm, *alpha_PIpPIm_pipf, *alpha_PPIp_piPIm, *alpha_PPIm_piPIp;

TH1F *h_inv_m_pip_pim_bin =  new TH1F("inv_m_pip_pim_bin","inv_m_pip_pim_bin",12.,0.228,0.69);
TH1F *h_inv_m_pip_p_bin =  new TH1F("inv_m_pip_p_bin","inv_m_pip_p_bin",12.,1.,1.485);
TH1F *h_inv_m_pim_p_bin =  new TH1F("inv_m_pim_p_bin","inv_m_pim_p_bin",12.,1.,1.485);
TH1F *h_theta_PIm_cm = new TH1F("theta_PIm_cm","theta_PIm_cm",10.,0.,180.);
TH1F *h_theta_PIp_cm = new TH1F("theta_PIp_cm","theta_PIp_cm",10.,0.,180.);
TH1F *h_theta_P_cm = new TH1F("theta_P_cm","theta_P_cm",10.,0.,180.);
TH1F *h_phi_PIm_cm = new TH1F("phi_PIm_cm","phi_PIm_cm",8.,0.,360.);
TH1F *h_phi_PIp_cm = new TH1F("phi_PIp_cm","phi_PIp_cm",8.,0.,360.);
TH1F *h_phi_P_cm = new TH1F("phi_P_cm","phi_P_cm",8.,0.,360.);
TH1F *h_alpha_PIpPIm_pipf = new TH1F("alpha_PIpPIm_pipf","alpha_PIpPIm_pipf",8.,0.,360.);
TH1F *h_alpha_PPIp_piPIm = new TH1F("alpha_PPIp_piPIm","alpha_PPIp_piPIm",8.,0.,360.);
TH1F *h_alpha_PPIm_piPIp = new TH1F("alpha_PPIm_piPIp","alpha_PPIm_piPIp",8.,0.,360.);

static const Int_t ndims = 5;
  Int_t bins[5];
  
  
  Double_t W_bin [21];
  
   Double_t xmin[5] = {(0.938272 + 0.13957), (0.13957 + 0.13957),0.,0.,0.};
  
   Double_t xmax_1[21];
   Double_t xmax_2[21];
   
   Double_t xmax[5];






TH1F *h_inv_m_pip_pim = new TH1F("inv_m_pip_pim","inv_m_pip_pim",210, 0.,2.1);
TH1F *h_inv_m_pip_p = new TH1F("inv_m_pip_p","inv_m_pip_p",210, 0.,2.1);
TH1F *h_inv_m_pim_p = new TH1F("inv_m_pim_p","inv_m_pim_p",210, 0.,2.1);


TH2F *ph_vs_th_1 = new TH2F("ph_vs_th_1","ph_vs_th_1",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_2 = new TH2F("ph_vs_th_2","ph_vs_th_2",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_3 = new TH2F("ph_vs_th_3","ph_vs_th_3",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_4 = new TH2F("ph_vs_th_4","ph_vs_th_4",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_5 = new TH2F("ph_vs_th_5","ph_vs_th_5",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_6 = new TH2F("ph_vs_th_6","ph_vs_th_6",200, 0., 60.,150., -35., 35.);

TH1D *h_cos_th;

ostringstream qqq;
Short_t i,j;
TH2F *ph_vs_th_1pe[17],*ph_vs_th_2pe[17],*ph_vs_th_3pe[17],*ph_vs_th_4pe[17],*ph_vs_th_5pe[17],*ph_vs_th_6pe[17];
TH2F *ph_vs_th_1pe_fid[17],*ph_vs_th_2pe_fid[17],*ph_vs_th_3pe_fid[17],*ph_vs_th_4pe_fid[17],*ph_vs_th_5pe_fid[17],*ph_vs_th_6pe_fid[17];
TH2F *beta_vs_p_p[6][48];
TH2F *beta_vs_p_pip[6][48];
TH2F *beta_vs_p_pim[6][48]; 

TH2F *beta_vs_p_p_sim[6][48];
TH2F *beta_vs_p_pip_sim[6][48];
TH2F *beta_vs_p_pim_sim[6][48]; 


TH2F *ph_vs_th_pim[6][10];
TH2F *ph_th_pim_all_p[6]; 

TH2F *time_pip[6][48];
TH2F *time_pim[6][48];
TH2F *time_p[6][48];

TH2F *time_pip_sim[6][48];
TH2F *time_pim_sim[6][48];
TH2F *time_p_sim[6][48];

TH2F *ph_vs_th_el_sim[6][7];
TH2F *ph_vs_th_p_sim[6];
TH2F *ph_vs_th_pip_sim[6];
TH2F *ph_vs_th_pim_sim[6][5];

TH1F *Wbefore[6];
TH1F *Wafter[6];

//photoelectrons
TH1F *ph_el_left[6][20];
TH1F *ph_el_both[6][20];
TH1F *ph_el_right[6][20];


THnSparseD *h_5dim_0_1[12][21];
THnSparseD *h_5dim_0_2[12][21];
THnSparseD *h_5dim_0_3[12][21];

THnSparseD *h_5dim_pim_1[12][21];
THnSparseD *h_5dim_pim_2[12][21];
THnSparseD *h_5dim_pim_3[12][21];

THnSparseD *h_5dim_pip_1[12][21];
THnSparseD *h_5dim_pip_2[12][21];
THnSparseD *h_5dim_pip_3[12][21];

THnSparseD *h_5dim_pr_1[12][21];
THnSparseD *h_5dim_pr_2[12][21];
THnSparseD *h_5dim_pr_3[12][21];


THnSparseD *h_5dim_0_1_empty[12][21];
THnSparseD *h_5dim_0_2_empty[12][21];
THnSparseD *h_5dim_0_3_empty[12][21];

THnSparseD *h_5dim_pim_1_empty[12][21];
THnSparseD *h_5dim_pim_2_empty[12][21];
THnSparseD *h_5dim_pim_3_empty[12][21];

THnSparseD *h_5dim_pip_1_empty[12][21];
THnSparseD *h_5dim_pip_2_empty[12][21];
THnSparseD *h_5dim_pip_3_empty[12][21];

THnSparseD *h_5dim_pr_1_empty[12][21];
THnSparseD *h_5dim_pr_2_empty[12][21];
THnSparseD *h_5dim_pr_3_empty[12][21];

THnSparseD *h_5dim_0_1_sim_1[12][21];
THnSparseD *h_5dim_0_2_sim_1[12][21];
THnSparseD *h_5dim_0_3_sim_1[12][21];

THnSparseD *h_5dim_pim_1_sim_1[12][21];
THnSparseD *h_5dim_pim_2_sim_1[12][21];
THnSparseD *h_5dim_pim_3_sim_1[12][21];

THnSparseD *h_5dim_pip_1_sim_1[12][21];
THnSparseD *h_5dim_pip_2_sim_1[12][21];
THnSparseD *h_5dim_pip_3_sim_1[12][21];

THnSparseD *h_5dim_pr_1_sim_1[12][21];
THnSparseD *h_5dim_pr_2_sim_1[12][21];
THnSparseD *h_5dim_pr_3_sim_1[12][21];


THnSparseD *h_5dim_1_sim_2[12][21];
THnSparseD *h_5dim_2_sim_2[12][21];
THnSparseD *h_5dim_3_sim_2[12][21];

THnSparseD *h_5dim_excl_1[9];
THnSparseD *h_5dim_excl_2[9];
THnSparseD *h_5dim_excl_3[9];

THnSparseD *h_5dim_excl_1_empty[9];
THnSparseD *h_5dim_excl_2_empty[9];
THnSparseD *h_5dim_excl_3_empty[9];

THnSparseD *h_5dim_1_sim_excl_1[9];
THnSparseD *h_5dim_2_sim_excl_1[9];
THnSparseD *h_5dim_3_sim_excl_1[9];

THnSparseD *h_5dim_1_sim_excl_2[9];
THnSparseD *h_5dim_2_sim_excl_2[9];
THnSparseD *h_5dim_3_sim_excl_2[9];


TH1D *h1prj_inv_m_pip_p[12][21]; 
TH1D *h1prj_inv_m_pip_pim[12][21]; 
TH1D *h1prj_th_P[12][21]; 
TH1D *h1prj_phi_P[12][21]; 
TH1D *h1prj_alpha_PIpPIm_pipf[12][21];

TH1D *h2prj_inv_m_pip_p[12][21]; 
TH1D *h2prj_inv_m_pip_pim[12][21]; 
TH1D *h2prj_th_PIm[12][21]; 
TH1D *h2prj_phi_PIm[12][21]; 
TH1D *h2prj_alpha_PPIp_piPIm[12][21];

TH1D *h3prj_inv_m_pim_p[12][21]; 
TH1D *h3prj_inv_m_pip_pim[12][21]; 
TH1D *h3prj_th_PIp[12][21]; 
TH1D *h3prj_phi_PIp[12][21]; 
TH1D *h3prj_alpha_PPIm_piPIp[12][21];  


TH1D *h1prj_inv_m_pip_p_sim[12][21]; 
TH1D *h1prj_inv_m_pip_pim_sim[12][21]; 
TH1D *h1prj_th_P_sim[12][21]; 
TH1D *h1prj_phi_P_sim[12][21]; 
TH1D *h1prj_alpha_PIpPIm_pipf_sim[12][21];

TH1D *h2prj_inv_m_pip_p_sim[12][21]; 
TH1D *h2prj_inv_m_pip_pim_sim[12][21]; 
TH1D *h2prj_th_PIm_sim[12][21]; 
TH1D *h2prj_phi_PIm_sim[12][21]; 
TH1D *h2prj_alpha_PPIp_piPIm_sim[12][21];

TH1D *h3prj_inv_m_pim_p_sim[12][21]; 
TH1D *h3prj_inv_m_pip_pim_sim[12][21]; 
TH1D *h3prj_th_PIp_sim[12][21]; 
TH1D *h3prj_phi_PIp_sim[12][21]; 
TH1D *h3prj_alpha_PPIm_piPIp_sim[12][21]; 

TH1D *h1prj_inv_m_pip_p_sim_2[12][21]; 
TH1D *h1prj_inv_m_pip_pim_sim_2[12][21]; 
TH1D *h1prj_th_P_sim_2[12][21]; 
TH1D *h1prj_phi_P_sim_2[12][21]; 
TH1D *h1prj_alpha_PIpPIm_pipf_sim_2[12][21];

TH1D *h2prj_inv_m_pip_p_sim_2[12][21]; 
TH1D *h2prj_inv_m_pip_pim_sim_2[12][21]; 
TH1D *h2prj_th_PIm_sim_2[12][21]; 
TH1D *h2prj_phi_PIm_sim_2[12][21]; 
TH1D *h2prj_alpha_PPIp_piPIm_sim_2[12][21];

TH1D *h3prj_inv_m_pim_p_sim_2[12][21]; 
TH1D *h3prj_inv_m_pip_pim_sim_2[12][21]; 
TH1D *h3prj_th_PIp_sim_2[12][21]; 
TH1D *h3prj_phi_PIp_sim_2[12][21]; 
TH1D *h3prj_alpha_PPIm_piPIp_sim_2[12][21]; 




TH1F *h_pim_mm_q2_w[12][21];
TH1F *h_0_mm_q2_w[21];
TH1F *h_0_d_mm_q2_w[21];

TH1F *h_pim_mm_q2_w_sim[12][21];
TH1F *h_0_mm_q2_w_sim[21];
TH1F *h_0_d_mm_q2_w_sim[21];

TH1F *h_0_mmom_q2_w[21];
TH1F *h_0_mmom_q2_w_sim[21];
TH1F *h_0_d_mmom_q2_w[21];
TH1F *h_0_d_mmom_q2_w_sim[21];

TH1D *h_w_int[12];
TH2F *th_cc_vs_seg[6];
TH2F *th_cc_vs_seg_sim[6];
TH2F *th_vs_p_e_1[6], *th_vs_p_e_2[6], *th_vs_p_p_1[6],*th_vs_p_p_2[6],*th_vs_p_pip_1[6],*th_vs_p_pip_2[6],*th_vs_p_pim_1[6],*th_vs_p_pim_2[6];
TH2F *th_vs_p_e_1_sim[6], *th_vs_p_e_2_sim[6], *th_vs_p_p_1_sim[6],*th_vs_p_p_2_sim[6],*th_vs_p_pip_1_sim[6],*th_vs_p_pip_2_sim[6],*th_vs_p_pim_1_sim[6],*th_vs_p_pim_2_sim[6];
// booking proton histograms

TH2F *eout_vs_ein_before = new TH2F("eout_vs_ein_before","eout_vs_ein_before",100, 0., 0.5,100, 0., 0.5);
TH2F *eout_vs_ein_after = new TH2F("eout_vs_ein_after","eout_vs_ein_after",100, 0., 0.5,100, 0., 0.5);

TH2F *ph_vs_th_p_1 = new TH2F("ph_vs_th_p_1","ph_vs_th_p_1",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_2 = new TH2F("ph_vs_th_p_2","ph_vs_th_p_2",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_3 = new TH2F("ph_vs_th_p_3","ph_vs_th_p_3",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_4 = new TH2F("ph_vs_th_p_4","ph_vs_th_p_4",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_5 = new TH2F("ph_vs_th_p_5","ph_vs_th_p_5",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_6 = new TH2F("ph_vs_th_p_6","ph_vs_th_p_6",220, 0., 110.,150., -35., 35.);

TH2F *ph_vs_th_p_1_w = new TH2F("ph_vs_th_p_1_w","ph_vs_th_p_1_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_2_w = new TH2F("ph_vs_th_p_2_w","ph_vs_th_p_2_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_3_w = new TH2F("ph_vs_th_p_3_w","ph_vs_th_p_3_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_4_w = new TH2F("ph_vs_th_p_4_w","ph_vs_th_p_4_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_5_w = new TH2F("ph_vs_th_p_5_w","ph_vs_th_p_5_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_6_w = new TH2F("ph_vs_th_p_6_w","ph_vs_th_p_6_w",220, 0., 110.,150., -35., 35.);


TH2F *ph_vs_th_pip_1 = new TH2F("ph_vs_th_pip_1","ph_vs_th_pip_1",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_2 = new TH2F("ph_vs_th_pip_2","ph_vs_th_pip_2",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_3 = new TH2F("ph_vs_th_pip_3","ph_vs_th_pip_3",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_4 = new TH2F("ph_vs_th_pip_4","ph_vs_th_pip_4",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_5 = new TH2F("ph_vs_th_pip_5","ph_vs_th_pip_5",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_6 = new TH2F("ph_vs_th_pip_6","ph_vs_th_pip_6",220, 0., 110.,150., -35., 35.);

TH2F  *h_delta_w_vs_w_old_data = new TH2F("h_delta_w_vs_w_old_data","h_delta_w_vs_w_old_data",180, 1.3, 1.9,300., -0.3,0.3);
TH2F  *h_delta_w_vs_w_old_data_mmcut = new TH2F("h_delta_w_vs_w_old_data_mmcut","h_delta_w_vs_w_old_data_mmcut",180, 1.3, 1.9,300., -0.3,0.3);
TH2F  *h_delta_w_vs_w_old_gen = new TH2F("h_delta_w_vs_w_old_gen","h_delta_w_vs_w_old_gen",180, 1.3, 1.9,300., -0.3,0.3);
TH2F  *h_delta_w_vs_w_old_rec = new TH2F("h_delta_w_vs_w_old_rec","h_delta_w_vs_w_old_rec",180, 1.3, 1.9,300., -0.3,0.3);

TH2F  *ph_th_p_1[5],*ph_th_p_2[5],*ph_th_p_3[5],*ph_th_p_4[5],*ph_th_p_5[5],*ph_th_p_6[5];

TH2F  *ph_th_pip_1[6],*ph_th_pip_2[6],*ph_th_pip_3[6],*ph_th_pip_4[6],*ph_th_pip_5[6],*ph_th_pip_6[6];

TH1F *padd_rate[6];
TH2F *padd_prot_mass[6];
//output file

TFile *outFile;

// bookking quality check histograms

TH1F *hist_ltime = new TH1F("ltime","ltime",35000, -.5, 34999.5);
TH1F *hist_ltime_1d = new TH1F("ltime_1d","ltime_1d",1100, 0., 1.1);
TH1F *hist_n_incl = new TH1F("n_incl","n_incl",35000, -.5, 34999.5);
TH1F *hist_n_incl_1d = new TH1F("n_incl_1d","n_incl_1d",1000, 0., 300000.);
TH1F *hist_n_twopions = new TH1F("n_twopions","n_twopions",35000, -.5, 34999.5);
TH1F *hist_n_twopions_1d = new TH1F("n_twopions_1d","n_twopions_1d",5000, -0.5, 4999.5);
TH1F *hist_n_elast = new TH1F("n_elast","n_elast",35000, -.5, 34999.5);
TH1F *hist_n_elast_1d = new TH1F("n_elast_1d","n_elast_1d",1000, 0., 80000.);
TH1F *hist_twopi_cut;
Int_t two_pions_block;


int global() {

for (i=0; i<6; i++) {
qqq.str("");
qqq << "padd_rate_" << i;
padd_rate[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),48, .5, 48.5);
qqq.str("");
qqq << "padd_prot_mass_" << i;
padd_prot_mass[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),48, .5, 48.5,100,0.,1.5);
};
// loops for electron histograms

for (i=0; i<6; i++) {
qqq.str("");
qqq << "W_before_sector" << i;
Wbefore[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),200, 0.8, 1.2);
qqq.str("");
qqq.str("");
qqq << "W_after_sector" << i << "]";
Wafter[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),200, 0.8, 1.2);
qqq.str("");
};

for (i=0; i<17; i++) {
qqq << "ph_vs_th_1pe[" << i << "]";
ph_vs_th_1pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<22; i++) {

qqq << "pim_mis_radcorr_" << i;
h_pim_mis_radcorr[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, -0.05, 0.07);
qqq.str("");

qqq << "pim_mis_radcorr_sim_" << i;
h_pim_mis_radcorr_sim[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, -0.05, 0.07);
qqq.str("");

qqq << "pim_mis_radcorr_sim_bckg_" << i;
h_pim_mis_radcorr_sim_bckg[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, -0.05, 0.07);
qqq.str("");


qqq << "pip_mis_radcorr_" << i;
h_pip_mis_radcorr[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, -0.05, 0.07);
qqq.str("");

qqq << "pip_mis_radcorr_sim_" << i;
h_pip_mis_radcorr_sim[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, -0.05, 0.07);
qqq.str("");

qqq << "pip_mis_radcorr_sim_bckg_" << i;
h_pip_mis_radcorr_sim_bckg[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, -0.05, 0.07);
qqq.str("");


qqq << "pr_mis_radcorr_" << i;
h_pr_mis_radcorr[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, 0.8, 0.95);
qqq.str("");

qqq << "pr_mis_radcorr_sim_" << i;
h_pr_mis_radcorr_sim[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, 0.8, 0.95);
qqq.str("");

qqq << "pr_mis_radcorr_sim_bckg_" << i;
h_pr_mis_radcorr_sim_bckg[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, 0.8, 0.95);
qqq.str("");

qqq << "0_mis_radcorr_" << i;
h_0_mis_radcorr[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, -0.06, 0.06);
qqq.str("");

qqq << "0_mis_radcorr_sim_" << i;
h_0_mis_radcorr_sim[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, -0.06, 0.06);
qqq.str("");

qqq << "0_mis_radcorr_sim_bckg_" << i;
h_0_mis_radcorr_sim_bckg[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),200, -0.06, 0.06);
qqq.str("");



qqq << "h_thp_vs_thp_" << i;
h_thp_vs_thp[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),30, 0., 180.,30, 0., 180.);
qqq.str("");

};

for (i=0; i<22; i++) {

qqq << "pr_mom_test_" << i;
h_pr_mom_test[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),100, 0., 2.);
qqq.str("");

qqq << "pip_mom_test_" << i;
h_pip_mom_test[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),100, 0., 2.);
qqq.str("");


};

for (i=0; i<22; i++) {

qqq << "pr_mom_sim_test_" << i;
h_pr_mom_sim_test[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),100, 0., 2.);
qqq.str("");

qqq << "pip_mom_sim_test_" << i;
h_pip_mom_sim_test[i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),100, 0., 2.);
qqq.str("");

};



for (i=0; i<17; i++) {
qqq << "ph_vs_th_2pe[" <<i << "]";
ph_vs_th_2pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_3pe[" << i << "]";
ph_vs_th_3pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_4pe[" << i << "]";
ph_vs_th_4pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_5pe[" << i << "]";
ph_vs_th_5pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_6pe[" << i << "]";
ph_vs_th_6pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_1pe_fid[" << i << "]";
ph_vs_th_1pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_2pe_fid[" << i << "]";
ph_vs_th_2pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_3pe_fid[" << i << "]";
ph_vs_th_3pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_4pe_fid[" << i << "]";
ph_vs_th_4pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_5pe_fid[" << i << "]";
ph_vs_th_5pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_6pe_fid[" << i << "]";
ph_vs_th_6pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};

//loops for proton histograms


for (i=0; i<5; i++) {
qqq << "ph_th_p_1[" << i << "]";
ph_th_p_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_2[" << i << "]";
ph_th_p_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_3[" << i << "]";
ph_th_p_3[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_4[" << i << "]";
ph_th_p_4[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_5[" << i << "]";
ph_th_p_5[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_6[" << i << "]";
ph_th_p_6[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

//pip fid

for (i=0; i<6; i++) {
qqq << "ph_th_pip_1[" << i << "]";
ph_th_pip_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_2[" << i << "]";
ph_th_pip_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_3[" << i << "]";
ph_th_pip_3[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_4[" << i << "]";
ph_th_pip_4[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_5[" << i << "]";
ph_th_pip_5[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_6[" << i << "]";
ph_th_pip_6[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};




for (i=0; i<6; i++) {
for(j=0; j<7; j++){
qqq << "ph_vs_th_el_sim_" << i+1 << "_" << j+1;
ph_vs_th_el_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),240, 0., 60.,210., -35., 35.);
qqq.str("");
};
//for(j=0; j<5; j++){
qqq << "ph_vs_th_p_sim_" << i+1 << "_" << j+1;
ph_vs_th_p_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),330, 0., 110.,210., -35., 35.);
qqq.str("");
//};
//for(j=0; j<6; j++){
qqq << "ph_vs_th_pip_sim_" << i+1 << "_" << j+1;
ph_vs_th_pip_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),330, 0., 110.,210., -35., 35.);
qqq.str("");
//};
for(j=0; j<5; j++){
qqq << "ph_vs_th_pim_sim_" << i+1 << "_" << j+1;
ph_vs_th_pim_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),360, 0., 120.,210., -35., 35.);
qqq.str("");

};
};




//beta_vs_p_pos[48][6]
for (i=0; i<6; i++) {
for(j=0; j<48; j++){
qqq << "beta_vs_p_p_" << i+1 << "_" << j+1;
beta_vs_p_p[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,300, 0., 1.5);
qqq.str("");
qqq << "beta_vs_p_pip_" << i+1 << "_" << j+1;
beta_vs_p_pip[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,300, 0., 1.5);
qqq.str("");
qqq << "beta_vs_p_pim_" << i+1 << "_" << j+1;
beta_vs_p_pim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,300, 0., 1.5);
qqq.str("");


qqq << "time_p_" << i+1 << "_" << j+1;
time_p[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");
qqq << "time_pip_" << i+1 << "_" << j+1;
time_pip[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");
qqq << "time_pim_" << i+1 << "_" << j+1;
time_pim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");


qqq << "beta_vs_p_p_sim_" << i+1 << "_" << j+1;
beta_vs_p_p_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,240, 0., 1.2);
qqq.str("");
qqq << "beta_vs_p_pip_sim_" << i+1 << "_" << j+1;
beta_vs_p_pip_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,240, 0., 1.2);
qqq.str("");
qqq << "beta_vs_p_pim_sim_" << i+1 << "_" << j+1;
beta_vs_p_pim_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,240, 0., 1.2);
qqq.str("");


qqq << "time_p_sim_" << i+1 << "_" << j+1;
time_p_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");
qqq << "time_pip_sim_" << i+1 << "_" << j+1;
time_pip_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");
qqq << "time_pim_sim_" << i+1 << "_" << j+1;
time_pim_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");







};
};



//TH2F *ph_vs_th_pim[6][10]
//TH2F *ph_th_pim_all_p[6];
for (i=0; i<6; i++) {
qqq << "ph_th_pim_all_p[" << i+1;
ph_th_pim_all_p[i]=new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0., 120.,150., -35., 35.);
qqq.str("");
for(j=0; j<10; j++){
qqq << "ph_vs_th_pim[" << i+1 << "][" << j+1 <<"]" ;
ph_vs_th_pim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0., 120.,150., -35., 35.);
qqq.str("");
};
};


/*for (i=0; i<6; i++) {
for(j=0; j<57; j++){
qqq << "y_sc_vs_x_sc_p" << i+1 << "_" << j+1;
y_sc_vs_x_sc_p[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),600, -300.,300.,500, 0.,500.);
qqq.str("");
qqq << "y_sc_vs_x_sc_pim" << i+1 << "_" << j+1;
y_sc_vs_x_sc_pim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),600, -300.,300.,500, 0.,500.);
qqq.str("");
qqq << "y_sc_vs_x_sc_pip" << i+1 << "_" << j+1;
y_sc_vs_x_sc_pip[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),600, -300.,300.,500, 0.,500.);
qqq.str("");
qqq << "y_sc_vs_x_sc_el" << i+1 << "_" << j+1;
y_sc_vs_x_sc_el[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),600, -300.,300.,500, 0.,500.);
qqq.str("");

};
};*/

 //Int_t bins[5] = {12, 12, 10, 8, 8};
  
  
 /* Double_t W_bin [20];
  for (i=0; i<20; i++) W_bin[i] = 1.325+0.025*i;
  
   Double_t xmin[5] = {1., 0.288,0.,0.,0.};
   Double_t xmax_1[5];
   Double_t xmax_1[5];
   for (i=0; i<20; i++) {
   xmax_1[i] = W_bin[i] - 0.13957 - 0.05;
   xmax_2[i] = W_bin[i] - 0.938272 - 0.05;
   };
   Double_t xmax[5] = {1.485,0.69, 180.,360.,360.};*/
   
   



xmax[2] = 180.;
xmax[3] = 360.;
xmax[4] = 360.;

/*bins[2]=10;
bins[3]=8;
bins[4]=8;*/

bins[2]=10;
bins[3]=5;
bins[4]=8;


for (i=0; i<9; i++) {


W_bin[i] = 1.425+0.025*i;//right edge


xmax[0] = xmax_1[i]; 
xmax[1] = xmax_2[i];

bins[0]=12;
bins[1]=12;


xmax_1[i] =  (1.4125+0.025*i - 0.13957)+((1.4125+0.025*i - 0.13957)-(0.938272 + 0.13957))/(bins[0]-1);
xmax_2[i] = (1.4125+0.025*i - 0.938272)+((1.4125+0.025*i - 0.938272)-(0.13957 + 0.13957))/(bins[1]-1);

qqq << "h_5dim_excl_1"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_2"<<"_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_3"<<"_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_3[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_excl_1_empty" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_1_empty[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_2_empty" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_2_empty[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_3_empty" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_3_empty[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_1_sim_excl_1" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_1_sim_excl_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_2_sim_excl_1"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_2_sim_excl_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_3_sim_excl_1" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_3_sim_excl_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_1_sim_excl_2"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_1_sim_excl_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_2_sim_excl_2"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_2_sim_excl_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_3_sim_excl_2"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_3_sim_excl_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


};






for (i=0; i<21; i++) {





if ((i>=0)&&(i<=1)){
bins[0]=8;
bins[1]=8;
bins[2]=6;
bins[3]=5;
bins[4]=5;
};
if ((i>=2)&&(i<=3)){
bins[0]=10;
bins[1]=10;
bins[2]=8;
bins[3]=5;
bins[4]=6;
};

if ((i>=4)&&(i<=6)){
bins[0]=12;
bins[1]=12;
bins[2]=10;
bins[3]=5;
bins[4]=8;
};

if ((i>=7)&&(i<=14)){
bins[0]=12;
bins[1]=12;
bins[2]=10;
bins[3]=8;
bins[4]=8;
};

if ((i>=15)&&(i<=22)){
bins[0]=12;
bins[1]=12;
bins[2]=10;
bins[3]=8;
bins[4]=8;
};

W_bin[i] = 1.325+0.025*i;//right edge
xmax_1[i] =  (1.3125+0.025*i - 0.13957)+((1.3125+0.025*i - 0.13957)-(0.938272 + 0.13957))/(bins[0]-1);
xmax_2[i] = (1.3125+0.025*i - 0.938272)+((1.3125+0.025*i - 0.938272)-(0.13957 + 0.13957))/(bins[1]-1);

xmax[0] = xmax_1[i]; 
xmax[1] = xmax_2[i];

//cout << i<<" bins0= " << bins[0]<<  " bins1= " << bins[1]<< " bins2= " << bins[2]<< " bins3= " << bins[3]<< " bins4= " << bins[4]<< "\n";

qqq << "h_0_mm"<<"_w_" << 1000*(1.3125+0.025*i);
h_0_mm_q2_w[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.2, 0.2);
qqq.str("");

qqq << "h_0_mm_sim""_w_" << 1000*(1.3125+0.025*i);
h_0_mm_q2_w_sim[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.2, 0.2);
qqq.str("");


qqq << "h_0_d_mm"<<"_w_" << 1000*(1.3125+0.025*i);
h_0_d_mm_q2_w[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -1., 3.);
qqq.str("");

qqq << "h_0_d_mm_sim"<< "_w_" << 1000*(1.3125+0.025*i);
h_0_d_mm_q2_w_sim[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -1., 3.);
qqq.str("");

///////mom
qqq << "h_0_mmom"<<"_w_" << 1000*(1.3125+0.025*i);
h_0_mmom_q2_w[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),500, -0.1, 0.9);
qqq.str("");

qqq << "h_0_mmom_sim""_w_" << 1000*(1.3125+0.025*i);
h_0_mmom_q2_w_sim[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),500, -0.1, 0.9);
qqq.str("");


qqq << "h_0_d_mmom"<<"_w_" << 1000*(1.3125+0.025*i);
h_0_d_mmom_q2_w[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),600, -0.1, 1.7);
qqq.str("");

qqq << "h_0_d_mmom_sim_"<< "_w_" << 1000*(1.3125+0.025*i);
h_0_d_mmom_q2_w_sim[i]= new TH1F(qqq.str().c_str(),qqq.str().c_str(),600, -0.1, 1.7);
qqq.str("");






for(j=0; j<12; j++){

qqq << "h_pim_mm_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_pim_mm_q2_w[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),2000, -2., 2.);
qqq.str("");

qqq << "h_pim_mm_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_pim_mm_q2_w_sim[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),2000, -2., 2.);
qqq.str("");


qqq << "h_5dim_0_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_0_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_0_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_0_2[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_0_3_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_0_3[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pim_2[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_3_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pim_3[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pip_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pip_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pip_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pip_2[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pip_3_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pip_3[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pr_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pr_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pr_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pr_2[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pr_3_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pr_3[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");




qqq << "h_5dim_0_1_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_0_1_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_0_2_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_0_2_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_0_3_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_0_3_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_pim_1_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pim_1_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_2_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pim_2_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_3_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pim_3_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pip_1_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pip_1_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pip_2_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pip_2_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pip_3_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pip_3_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pr_1_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pr_1_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pr_2_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pr_2_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pr_3_empty_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pr_3_empty[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");










qqq << "h_5dim_0_1_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_0_1_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_0_2_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_0_2_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_0_3_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_0_3_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_1_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pim_1_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_2_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pim_2_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_3_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pim_3_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pip_1_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pip_1_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pip_2_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pip_2_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pip_3_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pip_3_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pr_1_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pr_1_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pr_2_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pr_2_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pr_3_sim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_pr_3_sim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");



qqq << "h_5dim_1_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_1_sim_2[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_2_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_2_sim_2[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_3_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.3125+0.025*i);
h_5dim_3_sim_2[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");



/*
qqq << "inv_m_pip_p_1_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i;
h1prj_inv_m_pip_p[j][i] =  new TH1D(qqq.str().c_str(),qqq.str().c_str(),bins[0],xmin[0],xmax[0]);
qqq.str("");
 
qqq << "inv_m_pip_pim_1_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i; 
h1prj_inv_m_pip_pim[j][i] =  new TH1D(qqq.str().c_str(),qqq.str().c_str(),bins[1],xmin[1],xmax[1]); 
qqq.str("");

qqq << "theta_P_cm_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i; 
h1prj_th_P[j][i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),10.,0.,180.); 
qqq.str("");

qqq << "phi_P_cm_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i; 
h1prj_phi_P[j][i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),8.,0.,360.); 
qqq.str("");

qqq << "alpha_PIpPIm_pipf_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i; 
h1prj_alpha_PIpPIm_pipf[j][i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),8.,0.,360.);
qqq.str("");

//////
qqq << "inv_m_pip_p_2_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i; 
h2prj_inv_m_pip_p[j][i] =  new TH1D(qqq.str().c_str(),qqq.str().c_str(),bins[0],xmin[0],xmax[0]);  
qqq.str("");

qqq << "inv_m_pip_pim_2_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i; 
h2prj_inv_m_pip_pim[j][i] =  new TH1D(qqq.str().c_str(),qqq.str().c_str(),bins[1],xmin[1],xmax[1]);
qqq.str("");
  
qqq << "theta_PIm_cm_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i;   
h2prj_th_PIm[j][i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),10.,0.,180.); 
qqq.str("");

qqq << "phi_PIm_cm_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i; 
h2prj_phi_PIm[j][i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),8.,0.,360.); 
qqq.str("");

qqq << "alpha_PPIp_piPIm_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i; 
h2prj_alpha_PPIp_piPIm[j][i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),8.,0.,360.);
qqq.str("");

///////////////

qqq << "inv_m_pim_p_3_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i;
h3prj_inv_m_pim_p[j][i] =  new TH1D(qqq.str().c_str(),qqq.str().c_str(),bins[0],xmin[0],xmax[0]);
qqq.str("");

qqq << "inv_m_pip_pim_3_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i; 
h3prj_inv_m_pip_pim[j][i] =  new TH1D(qqq.str().c_str(),qqq.str().c_str(),bins[1],xmin[1],xmax[1]); 
qqq.str("");
  
qqq << "theta_PIp_cm_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i;     
h3prj_th_PIp[j][i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),10.,0.,180.); 
qqq.str("");

qqq << "phi_PIp_cm_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i;
h3prj_phi_PIp[j][i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),8.,0.,360.); 
qqq.str("");

qqq << "alpha_PPIm_piPIp_"<<"q2_" << 0.425+0.05*j << "_w_" << 1.3125+0.025*i;
h3prj_alpha_PPIm_piPIp[j][i] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),8.,0.,360.);
qqq.str("");*/

};
};


for(j=0; j<12; j++){
qqq << "w_int_" << j+1;
h_w_int[j] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str("");
};

for (i=0; i<6; i++) {
for(j=0; j<20; j++){
qqq << "sector_" << i+1 << "_left_" << j+1 ;
ph_el_left[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),500, 0.,500);
qqq.str("");
};
};

for (i=0; i<6; i++) {
for(j=0; j<20; j++){
qqq << "sector_" << i+1 << "_both_" << j+1 ;
ph_el_both[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),500, 0.,500);
qqq.str("");
};
};

for (i=0; i<6; i++) {
for(j=0; j<20; j++){
qqq << "sector_" << i+1 << "_right_" << j+1 ;
ph_el_right[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),500, 0.,500);
qqq.str("");
};
};

for (i=0; i<6; i++) {
qqq << "th_cc_vs_seg_" << i+1;
th_cc_vs_seg[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),20, -0.5,19.5,200,0.,50.);
qqq.str("");

qqq << "th_cc_vs_seg_sim_" << i+1;
th_cc_vs_seg_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),20, -0.5,19.5,200,0.,50.);
qqq.str("");


};

for (i=0; i<6; i++) {
qqq << "th_vs_p_e_1_" << i+1;
th_vs_p_e_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),440, 0.,2.2,180,0.,60.);
qqq.str("");
qqq << "th_vs_p_e_2_" << i+1;
th_vs_p_e_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),440, 0.,2.2,180,0.,60.);
qqq.str("");
qqq << "th_vs_p_p_1_" << i+1;
th_vs_p_p_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_p_2_" << i+1;
th_vs_p_p_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pip_1_" << i+1;
th_vs_p_pip_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pip_2_" << i+1;
th_vs_p_pip_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pim_1_" << i+1;
th_vs_p_pim_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pim_2_" << i+1;
th_vs_p_pim_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");

qqq << "th_vs_p_e_1_sim_" << i+1;
th_vs_p_e_1_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),440, 0.,2.2,180,0.,60.);
qqq.str("");
qqq << "th_vs_p_e_2_sim_" << i+1;
th_vs_p_e_2_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),440, 0.,2.2,180,0.,60.);
qqq.str("");
qqq << "th_vs_p_p_1_sim_" << i+1;
th_vs_p_p_1_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_p_2_sim_" << i+1;
th_vs_p_p_2_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pip_1_sim_" << i+1;
th_vs_p_pip_1_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pip_2_sim_" << i+1;
th_vs_p_pip_2_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pim_1_sim_" << i+1;
th_vs_p_pim_1_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pim_2_sim_" << i+1;
th_vs_p_pim_2_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");



};


};
