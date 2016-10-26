#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TH3.h"
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
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGButtonGroup.h>
#include <RQ_OBJECT.h>
#include <TGNumberEntry.h>
#include <TGProgressBar.h>
#include <TGLabel.h>
#include <stdio.h>
#include <dlfcn.h>
#include "MyMainFrame.h"
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <TGFileDialog.h>
#include <GuiTypes.h>
#include <TGDoubleSlider.h>
#include <TGComboBox.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TError.h> 
#include <auto_ptr.h>
#ifndef __CINT__
#include <cstdlib>
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include <RooLandau.h> 
#include <RooNumConvPdf.h>
#include <RooDataHist.h>
#include "RooBinning.h"
#include <sys/types.h>
#include <wait.h>
#include <unistd.h>
#include <cstring>
#include <getopt.h>
#include <cstdlib>
#include "global.h"
#include "beta_func_data.h"
#include "beta_func_empty.h"
#include "rot_boost_cmsyst.h"
#include "fermi_bonn.h"
#include "data_hist.h"
#include "output.h"



 using namespace std; 
 

#define _USE_MATH_DEFINES


    void MyMainFrame::MainFrame(UChar_t flag, Float_t E_beam, Short_t nfiles, Short_t nfiles_empty, Short_t nfiles_sim, string inp_files[],string inp_files_empty[], string inp_files_sim[], string outfile_in) { 


	inpfile_inp = inp_files[0];
	outfile_inp = outfile_in;
	
	E0 = E_beam;

        n_files = nfiles;
	n_files_empty = nfiles_empty;
	n_files_sim = nfiles_sim;

//sozdaem massiv strok (novie)
        file = new string[nfiles];
	file_empty = new string[nfiles_empty];
	file_sim = new string[nfiles_sim];

//prisvaevaem novim massivam strok massivi strok iz input faila 
        file = inp_files;
	file_empty = inp_files_empty;
	file_sim = inp_files_sim;
	data_sim = flag;
	




   

     gROOT->SetBatch(true);
gROOT->ProcessLine( "gErrorIgnoreLevel = kWarning; " );
    ostringstream adc_num;
    ostringstream tdc_num;
    ostringstream ref_tdc;
    ostringstream qqq;
    ostringstream qqq1;
    
	 
    Long64_t j;
 Double_t fract_integ[3][6][18];   
 Double_t * adc_offset;
 adc_offset = new Double_t [12];
 Double_t * adc_cut;
 adc_cut = new Double_t [12]; 
 Short_t m, ti;
 Long64_t i,nstart,nstop,n_incl,n_elast, k_long;
 
 
TH1I *hist_adc_off[12]; 

 TFile *finp;
 Int_t block_total = 0;
 Int_t block_last = 0;
 Float_t Qfull = 0.;
  Float_t Qfull_empty = 0.;
  Float_t Qfull_sim = 0.; 
 
  Int_t block_curr = 0;
 Int_t block_prev = 0;
  Short_t m_old=0;
  
  
  
  
   Float_t P_EL_old,th_EL_old,Q2_old,P_EL_new,th_EL_new,W_new,Q2_new;
 bool selection;
 bool selection_pim_miss, selection_0_miss,selection_pip_miss,selection_proton_miss;
 bool selection_pim_miss_sim, selection_0_miss_sim,selection_pip_miss_sim,selection_proton_miss_sim;
 bool selection_pim_miss_empt, selection_0_miss_empt,selection_pip_miss_empt,selection_proton_miss_empt; 
 

 
  
	UChar_t pdhit;
	
	
	Float_t * p; 
        p = new Float_t [20];
	
        Float_t m_deutron,NpheCC_EL,Nphe_pip,Nphe_pim,ECtot_EL;
	Float_t p_miss_1,p_miss_2,p_miss_3,p_miss_4,p_miss_sqr;
	Float_t beta_nom_pip,beta_nom_pim,beta_nom_p,delta_t_pip,delta_t_pim,delta_t_p;
	Float_t p_fid_a_1, p_fid_b_1;
	Float_t  x_EL,y_EL;
	Float_t th_PIm_miss,ph_PIm_miss;
	Float_t  ECin_EL,ECout_EL;
	Int_t block, block_tot;
        Long64_t gpart,k,last_i,last_k;
	Float_t  q_l,Qdiff,Qcurr,Qprev,Qtotal,deltaQ;
	Int_t sc_part_local;
	Float_t sc_pd_local,sc_sect_local;
	Float_t sc_z,fid_a,fid_b,a,b;
	Float_t nx,ny,nz,par1,par2,th_min;
	Float_t sx,sy,sz,px,py,pz;
	Float_t delta_mom_p; 
	Float_t E_gamma,E_p_gamma_lab,P_p_gamma_lab,beta1,gamma;
	Float_t pip_fid_a_1,pip_fid_b_1;
	Float_t th_min_1,par1_1,par2_1, pim_fid_a_1,pim_fid_b_1;
	Float_t W_old;
	Float_t pf_x,pf_y,pf_z,pxel_new,pyel_new,pzel_new;
   	Double_t integ, err, err1, old_bin_cont, new_bin_cont;
	Float_t n_twopions = 0.;
	Float_t n_twopions_old = 0.;	
	bool two_pions_flag;
	//Int_t two_pions_block;
	
	Double_t Var1[5],Var2[5],Var3[5]; 
	Double_t Var_1[5],Var_2[5],Var_3[5]; 
	
 	bool cut_fiduch;
	bool bool_el_id_data, bool_proton_id_data, bool_pip_id_data,bool_pim_id_data;
	bool bool_el_id_sim, bool_proton_id_sim, bool_pip_id_sim,bool_pim_id_sim;
	bool bool_el_id_empt, bool_proton_id_empt, bool_pip_id_empt,bool_pim_id_empt;

	
	m_proton = 0.938272;
	//m_proton = 0.93957;
	m_deutron = 1.875612;
	m_pip = 0.13957;
	  
  
	TLorentzVector P4_D,P4_P,P4_PIm_miss,P4_PIp_miss,P4_PIp_miss_d,P4_P_miss,P4_PIm_miss_0;
	TLorentzVector P4_PP_rot,P4_PP_rot_1,P4_PP_rot_2,P4_PP_rot_3,P4_PP_rot_2_boost;
	TLorentzVector  P4_PP_cor, P4_PIp_cor, P4_PIm_cor;
	TLorentzVector  P4_miss_0,P4_miss_0_d,P4_PIm_miss_d;
	TLorentzVector  P4_miss_0_en_comp,P4_PIm_miss_en_comp;
	TLorentzVector  P4_ELP_for_miss_en_comp;
	
	TVector3 V3_dir_gamma,P3_PP_rot,uz,ux;
		
	Float_t th_gamma, phi_gamma;
		 
	P4_EL.SetXYZT(0,0,2.039,2.039);
	P4_P.SetXYZT(0,0,0,m_proton);
	P4_D.SetXYZT(0,0,0,m_deutron);
	 


global();

     //   TFile *twopicutfile = new TFile("two_pi_rate.root","READ");
     //   hist_twopi_cut = (TH1F*)twopicutfile->Get("n_twopions");
	
	
        TFile *nphefile = new TFile("new_ratio.root","READ");
	
	
        norm_nphe_s1 = (TH2F*)nphefile->Get("h301");
	norm_nphe_s2 = (TH2F*)nphefile->Get("h302");
	norm_nphe_s3 = (TH2F*)nphefile->Get("h303");
	norm_nphe_s4 = (TH2F*)nphefile->Get("h304");
	norm_nphe_s5 = (TH2F*)nphefile->Get("h305");
	norm_nphe_s6 = (TH2F*)nphefile->Get("h306");	

        //cout << " bin content1 = " << norm_nphe_s1->GetBinContent(100,100) << "\n";
	
      
	

      // cout << " bin content2 = " << avrg_nphe_sector1->GetBinContent(100,100) << "\n";

ifstream input("phel_integr_fract.txt");
if(input.is_open()){
i=0;
    while(!input.eof()){
          string line1,t_str, e_str,r_str,fr_str;
	   Int_t t,e,r;
	   Double_t fr;
           getline(input,line1); //read number
	   if (line1.length() != 0){ 
              t_str= line1.substr(0,line1.find(","));
            t = atof(t_str.c_str());
		   
	    e_str = line1.substr(t_str.length()+1,line1.substr(t_str.length()+1).find(","));
            e = atof(e_str.c_str());
	    	    
	    r_str = line1.substr(t_str.length()+e_str.length()+2, line1.substr(t_str.length()+e_str.length()+2).find(","));
            r = atof(r_str.c_str());
	    
	    fr_str = line1.substr(t_str.length()+e_str.length()+r_str.length()+3);
	    fr = atof(fr_str.c_str());
	    	    	    
	    //cout << t<< "   " << e << "   " << r << "   " << fr <<" \n";
	    fract_integ[t][e][r] = fr;
	    i=i+1;
	    	    };
	    
    };
};

input.close();
//for(k=0; k<3; k++){
//for(i=0; i<6; i++){
//for(j=0; j<18; j++){

//cout << k << "," << i << "," <<j << "," << fract_integ[k][i][j]<< "\n";
//};
//};
//};



  for (m=1; m<=n_files; m++) {

//sozdaem fail s imenem, vzyatim iz masiva strok  
  finp = new TFile(file[m-1].c_str()); 
  
  
  
 cout << "Processing file " << m << "\n"; 
 
 
 //cout << " bin content3 = " << norm_nphe_s1->GetBinContent(100,100) << "\n";
//berem derevo iz faila    
  TTree *t21 = (TTree*)finp->Get("t21");
  
  TBranch *br_indtype = t21->GetBranch("indtype");
    TBranch *br_npart = t21->GetBranch("npart");
    TBranch *br_p_el = t21->GetBranch("P_EL");
    TBranch *br_block = t21->GetBranch("block");
    TBranch *br_deltaQ = t21->GetBranch("deltaQ");
    TBranch *br_LiveTime = t21->GetBranch("LiveTime");
    TBranch *br_n_incl = t21->GetBranch("n_incl");
    TBranch *br_n_elast = t21->GetBranch("n_elast");
    TBranch *br_ph_EL = t21->GetBranch("ph_EL");
    TBranch *br_th_EL = t21->GetBranch("th_EL");
    TBranch *br_W = t21->GetBranch("W");
    TBranch *br_Q2 = t21->GetBranch("Q2");
    TBranch *br_nphe_el = t21->GetBranch("NpheCC_EL");
    TBranch *br_ectot_el = t21->GetBranch("ECtot_EL");
    TBranch *br_ecin_el = t21->GetBranch("ECin_EL");
    TBranch *br_ecout_el = t21->GetBranch("ECout_EL");
    TBranch *br_x_EL = t21->GetBranch("x_EL");
    TBranch *br_y_EL = t21->GetBranch("y_EL");
    TBranch *br_z_EL = t21->GetBranch("z_EL");
    TBranch *br_dc_z_EL = t21->GetBranch("dc_z_EL");
    TBranch *br_dc_z_PIm = t21->GetBranch("dc_z_PIm");
    TBranch *br_dc_z_PIp = t21->GetBranch("dc_z_PIp");
    TBranch *br_dc_z_P = t21->GetBranch("dc_z_P");    
    TBranch *br_z_P = t21->GetBranch("z_P");
    TBranch *br_z_PIp = t21->GetBranch("z_PIp");
    TBranch *br_z_PIm = t21->GetBranch("z_PIm");
    TBranch *br_PdHit_EL = t21->GetBranch("PdHit_EL");
    TBranch *br_PdHit_PIp = t21->GetBranch("PdHit_PIp");
    TBranch *br_PdHit_PIm = t21->GetBranch("PdHit_PIm");
    TBranch *br_PdHit_P = t21->GetBranch("PdHit_P");
    TBranch *br_sc_x = t21->GetBranch("sc_x");
    TBranch *br_sc_y = t21->GetBranch("sc_y");
    TBranch *br_sc_z = t21->GetBranch("sc_z");
    TBranch *br_pmt_hit = t21->GetBranch("pmt_hit");
    TBranch *br_segment = t21->GetBranch("segment");
    TBranch *br_theta_cc = t21->GetBranch("theta_cc");
    TBranch *br_ph_cc = t21->GetBranch("ph_cc");
    TBranch *br_sector = t21->GetBranch("sector");
    TBranch *br_n_PIp = t21->GetBranch("n_PIp");
    TBranch *br_n_PIm = t21->GetBranch("n_PIm");
    TBranch *br_n_P = t21->GetBranch("n_P");
    TBranch *br_p_pip = t21->GetBranch("P_PIp");
    TBranch *br_p_pim = t21->GetBranch("P_PIm");
    TBranch *br_p_p = t21->GetBranch("P_P");
    TBranch *br_th_PIp = t21->GetBranch("th_PIp");
    TBranch *br_th_PIm = t21->GetBranch("th_PIm");
    TBranch *br_th_P = t21->GetBranch("th_P");
    TBranch *br_ph_PIp = t21->GetBranch("ph_PIp");
    TBranch *br_ph_PIm = t21->GetBranch("ph_PIm");
    TBranch *br_ph_P = t21->GetBranch("ph_P");
    TBranch *br_beta_PIm = t21->GetBranch("beta_PIm");
    TBranch *br_beta_PIp = t21->GetBranch("beta_PIp");
    TBranch *br_beta_P = t21->GetBranch("beta_P");
    TBranch *br_beta_PIm_time = t21->GetBranch("beta_PIm_time");
    TBranch *br_beta_PIp_time = t21->GetBranch("beta_PIp_time");
    TBranch *br_beta_P_time = t21->GetBranch("beta_P_time");    
    TBranch *br_nphe_pip = t21->GetBranch("NpheCC_PIp");
    TBranch *br_PIp_time = t21->GetBranch("PIp_time");
    TBranch *br_PIp_dist = t21->GetBranch("PIp_dist");
    TBranch *br_P_time = t21->GetBranch("P_time");
    TBranch *br_P_dist = t21->GetBranch("P_dist");
    TBranch *br_PIm_time = t21->GetBranch("PIm_time");
    TBranch *br_PIm_dist = t21->GetBranch("PIm_dist");
 
    
  Bool_t adc_cut_switch,tdc_cut_switch;
  
  Qdiff = 0.;
  Qcurr = 0.;
  Qprev = 0.;
  Qtotal = 0.;
  k = 0;
  block = 0;
  last_k = 0;
  nstart = 0;
  nstop = 0;
  n_incl = 0;
  n_elast = 0;
  

//tsikl po sobitiyam. ih odinakovoe kol-vo v kazdoi peremennoi  
  for (i=0; i<br_sector->GetEntries(); i++) { 
  
  Qprev = Qcurr;

br_segment->GetEntry(i);
br_pmt_hit->GetEntry(i);
 br_deltaQ->GetEntry(i);
  br_n_incl->GetEntry(i);
  br_n_elast->GetEntry(i);  
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_W->GetEntry(i);
  br_Q2->GetEntry(i); 
  br_npart->GetEntry(i);
  br_n_PIm->GetEntry(i);
  br_n_PIp->GetEntry(i);
  br_n_P->GetEntry(i);
  br_sector->GetEntry(i);
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_ph_cc->GetEntry(i);
  br_theta_cc->GetEntry(i);  
  br_p_el->GetEntry(i);  
  br_p_p->GetEntry(i); 
  br_p_pim->GetEntry(i);
  br_p_pip->GetEntry(i); 
  br_ectot_el->GetEntry(i); 
  br_ecin_el->GetEntry(i); 
  br_ecout_el->GetEntry(i);
  br_nphe_el->GetEntry(i);
  br_nphe_pip->GetEntry(i);
  br_th_EL->GetEntry(i);
  br_ph_EL->GetEntry(i); 
  br_th_P->GetEntry(i);
  br_ph_P->GetEntry(i);
  br_th_PIp->GetEntry(i);
  br_ph_PIp->GetEntry(i);
  br_th_PIm->GetEntry(i);
  br_ph_PIm->GetEntry(i);
  br_beta_PIp->GetEntry(i);
  br_beta_PIm->GetEntry(i);
  br_beta_P->GetEntry(i);
   br_beta_PIp_time->GetEntry(i);
  br_beta_PIm_time->GetEntry(i);
  br_beta_P_time->GetEntry(i); 
  br_PdHit_EL->GetEntry(i);
  br_PdHit_PIp->GetEntry(i);
  br_PdHit_PIm->GetEntry(i);
  br_PdHit_P->GetEntry(i);
  br_sc_x->GetEntry(i);
  br_sc_y->GetEntry(i); 
  br_sc_z->GetEntry(i); 
  br_z_EL->GetEntry(i);
  br_dc_z_EL->GetEntry(i);
  br_dc_z_PIm->GetEntry(i);  
  br_dc_z_PIp->GetEntry(i);   
  br_dc_z_P->GetEntry(i);     
  br_W->GetEntry(i);
  br_Q2->GetEntry(i);
  br_indtype->GetEntry(i);
  
  
  br_PIp_time->GetEntry(i);
  br_PIm_time->GetEntry(i);
  br_P_time->GetEntry(i); 
  br_PIp_dist->GetEntry(i);
  br_PIm_dist->GetEntry(i);
  br_P_dist->GetEntry(i);
  

  
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  P_P = br_p_p->GetLeaf("P_P")->GetValue();
  P_PIp = br_p_pip->GetLeaf("P_PIp")->GetValue();
  P_PIm = br_p_pim->GetLeaf("P_PIm")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();
  th_P = br_th_P->GetLeaf("th_P")->GetValue();
  ph_P = br_ph_P->GetLeaf("ph_P")->GetValue();
  th_PIp = br_th_PIp->GetLeaf("th_PIp")->GetValue();
  ph_PIp = br_ph_PIp->GetLeaf("ph_PIp")->GetValue();
  npart = br_npart->GetLeaf("npart")->GetValue();
  n_P = br_n_P->GetLeaf("n_P")->GetValue();
  n_PIp = br_n_PIp->GetLeaf("n_PIp")->GetValue();
  Nphe_pip = br_nphe_pip->GetLeaf("NpheCC_PIp")->GetValue();
  PdHit_EL = br_PdHit_EL->GetLeaf("PdHit_EL")->GetValue();
  PdHit_PIp = br_PdHit_PIp->GetLeaf("PdHit_PIp")->GetValue();
  PdHit_PIm = br_PdHit_PIm->GetLeaf("PdHit_PIm")->GetValue();
  PdHit_P = br_PdHit_P->GetLeaf("PdHit_P")->GetValue();
  sc_x = br_sc_x->GetLeaf("sc_x")->GetValue();
  sc_y = br_sc_y->GetLeaf("sc_y")->GetValue();
  sc_z = br_sc_z->GetLeaf("sc_z")->GetValue();  
  beta_PIp = br_beta_PIp_time->GetLeaf("beta_PIp_time")->GetValue();
  beta_PIm = br_beta_PIm_time->GetLeaf("beta_PIm_time")->GetValue();
  beta_P = br_beta_P_time->GetLeaf("beta_P_time")->GetValue();
  n_PIm = br_n_PIm->GetLeaf("n_PIm")->GetValue();
  ph_PIm = br_ph_PIm->GetLeaf("ph_PIm")->GetValue();
  th_PIm = br_th_PIm->GetLeaf("th_PIm")->GetValue();
  z_EL = br_z_EL->GetLeaf("z_EL")->GetValue();
  dc_z_EL = br_dc_z_EL->GetLeaf("dc_z_EL")->GetValue();
  dc_z_PIm = br_dc_z_PIm->GetLeaf("dc_z_PIm")->GetValue();
  dc_z_PIp = br_dc_z_PIp->GetLeaf("dc_z_PIp")->GetValue();
  dc_z_P = br_dc_z_P->GetLeaf("dc_z_P")->GetValue();  
  W = br_W->GetLeaf("W")->GetValue();
  Q2 = br_Q2->GetLeaf("Q2")->GetValue();
  segment = br_segment->GetLeaf("segment")->GetValue();
  pmt_hit = br_pmt_hit->GetLeaf("pmt_hit")->GetValue();

  PIp_time = br_PIp_time->GetLeaf("PIp_time")->GetValue();
  PIm_time = br_PIm_time->GetLeaf("PIm_time")->GetValue();
  P_time = br_P_time->GetLeaf("P_time")->GetValue();
  PIp_dist = br_PIp_dist->GetLeaf("PIp_dist")->GetValue();
  PIm_dist = br_PIm_dist->GetLeaf("PIm_dist")->GetValue();
  P_dist = br_P_dist->GetLeaf("P_dist")->GetValue();
  indtype = br_indtype->GetLeaf("indtype")->GetValue();
  
  
  
//p_fid_a = 24.*(1-exp(-1.*0.08*(th_P-9.)));
//p_fid_b = -1.*25.*(1-exp(-1.*0.1*(th_P-10.)));
  //pip_fid_a = 24.*(1-exp(-1.*0.08*(th_PIp-9.)));
  // pip_fid_b = -1.*25.*(1-exp(-1.*0.1*(th_PIp-10.))); 
  
 
  if (br_ectot_el->GetLeaf("ECtot_EL")->GetValue() > (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue())) { 
          ECT  = br_ectot_el->GetLeaf("ECtot_EL")->GetValue();  
         } else {
          ECT  = (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue());    
          };
	
	  if ((br_ecin_el->GetLeaf("ECin_EL")->GetValue()) > 0.) {
          if ((br_ecout_el->GetLeaf("ECout_EL")->GetValue()) > 0.) {
	  eout_vs_ein_before->Fill(br_ecin_el->GetLeaf("ECin_EL")->GetValue(),br_ecout_el->GetLeaf("ECout_EL")->GetValue(),1.);
	  };
	  };
	  
	  //quatity check   
	     
	     block_curr = br_block->GetLeaf("block")->GetValue(); 

    if ((m_old != 0) && (m != m_old))block_total=block_total+block_prev;
    
    m_old=m;


   LiveTime=br_LiveTime->GetLeaf("LiveTime")->GetValue();
   
   inclusive = br_n_incl->GetLeaf("n_incl")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue();
   elastic =  br_n_elast->GetLeaf("n_elast")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue(); 
   
   two_pions_flag = false; 
   if (block_curr != block_prev) {

   hist_ltime->Fill(block_total+block_curr,LiveTime);
   
   hist_ltime_1d->Fill(LiveTime);
   hist_n_incl->Fill(block_total+block_curr,inclusive);
hist_n_incl_1d->Fill(inclusive);  
 hist_n_elast->Fill(block_total+block_curr,elastic);
hist_n_elast_1d->Fill(elastic); 
 
 
//if ((LiveTime > 0.89) && (LiveTime <0.925) && (inclusive > 48000) &&(inclusive < 60000) && (elastic > 16400) && (elastic < 20500)){

   
   two_pions_block = block_total+block_prev;
   two_pions_flag = true; 
 
    hist_n_twopions_1d->Fill(n_twopions/br_deltaQ->GetLeaf("deltaQ")->GetValue());
   hist_n_twopions->Fill(two_pions_block,(n_twopions_old/br_deltaQ->GetLeaf("deltaQ")->GetValue()));
//    cout << "N twopions = " << n_twopions/br_deltaQ->GetLeaf("deltaQ")->GetValue() << "  N twopions from hist = " << hist_twopi_cut->GetBinContent(two_pions_block+2) << "\n";
    n_twopions = 0.;

   
// if ((LiveTime > 0.89) && (LiveTime <0.925) && (inclusive > 48000) &&(inclusive < 60000) && (elastic > 16400) && (elastic < 20500)&&(hist_twopi_cut->GetBinContent(two_pions_block+2) > 120.)&&(hist_twopi_cut->GetBinContent(two_pions_block+2) < 400.)){
   if ((LiveTime > 0.89) && (LiveTime <0.925) && (inclusive > 48000) &&(inclusive < 60000) && (elastic > 16400) && (elastic < 20500)){
   
 Qfull = Qfull + br_deltaQ->GetLeaf("deltaQ")->GetValue();
  };   
   block_prev = block_curr;
   };   
	  




	  
//cout<<segment<<"\n";	  

//MOM CORR
/*
P_EL_new = corrfunc.correct_pel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
th_EL_new = corrfunc.correct_thel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
P_EL_old = P_EL;
P_EL = P_EL_new;
th_EL_old = th_EL;
th_EL = th_EL_new;



W_new = pow(float((E0+m_proton-P_EL)),2);
W_new = W_new -pow(P_EL*sin(th_EL*M_PI/180.)*cos(ph_EL*M_PI/180.),2);
W_new = W_new -pow(P_EL*sin(th_EL*M_PI/180.)*sin(ph_EL*M_PI/180.),2);
W_new = W_new -pow(float(E0 - P_EL*cos(th_EL*M_PI/180.)),2);
W_new = sqrt(W_new);
W = W_new;  

Q2_new = pow(float(E0-P_EL),2);
Q2_new = Q2_new -pow(P_EL*sin(th_EL*M_PI/180.)*cos(ph_EL*M_PI/180.),2);
Q2_new = Q2_new -pow(P_EL*sin(th_EL*M_PI/180.)*sin(ph_EL*M_PI/180.),2);
Q2_new = Q2_new -pow(float(E0 - P_EL*cos(th_EL*M_PI/180.)),2);
Q2 = -Q2_new;

*/   
 
sector =  int(br_sector->GetLeaf("sector")->GetValue());
nphe = br_nphe_el->GetLeaf("NpheCC_EL")->GetValue();
theta_cc = br_theta_cc->GetLeaf("theta_cc")->GetValue();
ph_cc = br_ph_cc->GetLeaf("ph_cc")->GetValue();  
  

//Missing mass pim
/*p_miss_1 = -P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.)-P_P*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.)-P_PIp*cos(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.);
p_miss_2 = -P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.)-P_P*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.)-P_PIp*sin(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.);
p_miss_3 = E0-P_EL*cos(th_EL*M_PI/180.)-P_P*cos(th_P*M_PI/180.)-P_PIp*cos(th_PIp*M_PI/180.);
p_miss_4 = E0+m_proton-P_EL-sqrt(m_proton*m_proton+P_P*P_P)-sqrt(m_pip*m_pip+P_PIp*P_PIp);

p_miss_sqr = p_miss_4*p_miss_4-p_miss_3*p_miss_3-p_miss_2*p_miss_2-p_miss_1*p_miss_1;*/

P4_ELP_reg.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),P_EL);

P4_PP_reg.SetXYZT(P_P*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*cos(th_P*M_PI/180.),sqrt(m_proton*m_proton+P_P*P_P));

P4_PIp_reg.SetXYZT(P_PIp*cos(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*sin(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*cos(th_PIp*M_PI/180.),sqrt(m_pip*m_pip+P_PIp*P_PIp));


P4_PIm_reg.SetXYZT(P_PIm*cos(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*sin(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*cos(th_PIm*M_PI/180.),sqrt(m_pip*m_pip+P_PIm*P_PIm));




P4_PIm_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;
//P4_PIm_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;
P4_PIp_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;

P4_PIp_miss_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_PIm_miss_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;

P4_P_miss = P4_EL + P4_P - P4_ELP_reg - P4_PIp_reg - P4_PIm_reg;
P4_miss_0 = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;

P4_miss_0_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;

//W_old = W;
//W = sqrt((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2());
//W=sqrt((P4_P+P4_EL-P4_ELP_reg).Mag2());


///ENERGY LOSS

/*
delta_mom_p = corrfunc.correct_energy_theta_pf(P_P, th_P);

//cout << "theta = " << th_P << " P_P = " << P_P << " delta_p = " << delta_mom_p << "\n";

P4_PP_cor.SetXYZT((delta_mom_p+P_P)*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),(delta_mom_p+P_P)*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),(delta_mom_p+P_P)*cos(th_P*M_PI/180.),sqrt((delta_mom_p+P_P)*(delta_mom_p+P_P)+m_proton*m_proton));

P4_PIm_cor = P4_EL + P4_P - P4_ELP_reg - P4_PP_cor - P4_PIp_reg;
P4_PIp_cor = P4_EL + P4_P - P4_ELP_reg - P4_PP_cor - P4_PIm_reg;
*/






th_PIm_miss = (180./M_PI)*acos(P4_PIm_miss[2]/sqrt(P4_PIm_miss[0]*P4_PIm_miss[0]+P4_PIm_miss[1]*P4_PIm_miss[1]+P4_PIm_miss[2]*P4_PIm_miss[2]));

/*
if(P4_PIm_miss[0] == 0.) {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;

};


if ((P4_PIm_miss[1]>0)&&(P4_PIm_miss[0]>0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
};


if ((P4_PIm_miss[1]>0)&&(P4_PIm_miss[0]<0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 180. - (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};

if ((P4_PIm_miss[1]<0)&&(P4_PIm_miss[0]<0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 180. + (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};

if ((P4_PIm_miss[1]<0)&&(P4_PIm_miss[0]>0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 360. - (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};*/


if (P4_PIm_miss[0] != 0.) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
}
 else {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;
};
				   

if ((P4_PIm_miss[0] < 0.) && (P4_PIm_miss[1] > 0)) ph_PIm_miss = ph_PIm_miss+180.;
if (( P4_PIm_miss[0]< 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+180.;
if ((P4_PIm_miss[0] > 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+360.;

//E_miss = E0-m_proton-P_EL-sqrt(m_proton^2+P_P^2)-sqrt(m_pip^2+P_PIp^2)-sqrt(m_pip^2+P_PIm^2);

//cout << "theta_cc = " << theta_cc << " ph_cc = " << ph_cc << " bin content = " << norm_nphe_s1->GetBinContent(100,100) << "\n";


//if (PdHit_PIp==46) cout << "1 "<< beta_PIp << "\n";
//beta_func_data();
//if (PdHit_PIp==46) cout << "2 "<< beta_PIp << "\n";


bool_el_id_data=particle_ID_data.Electron_cuts_data(); 




selection = false;
selection_pim_miss = false;
selection_pip_miss = false;
selection_proton_miss = false;
selection_0_miss = false;

if (bool_el_id_data) {

beta_func_data();

data_hist();

if ((br_ecin_el->GetLeaf("ECin_EL")->GetValue()) > 0.) {
if ((br_ecout_el->GetLeaf("ECout_EL")->GetValue()) > 0.) {
eout_vs_ein_after->Fill(br_ecin_el->GetLeaf("ECin_EL")->GetValue(),br_ecout_el->GetLeaf("ECout_EL")->GetValue(),1.);
};
};



bool_proton_id_data=particle_ID_data.Proton_cuts_data();
bool_pip_id_data=particle_ID_data.PIp_cuts_data();
bool_pim_id_data=particle_ID_data.PIm_cuts_data();

//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart >=3)){
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){


if ((bool_pip_id_data)&&(bool_pim_id_data)&&(!bool_proton_id_data)){
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) < W - m_pip + 0.05) {
if ((abs(dc_z_EL - dc_z_PIp) < 4.) && (abs(dc_z_EL - dc_z_PIm) < 4.) && (abs(dc_z_PIp - dc_z_PIm) < 4.)) {
hist_P_miss-> Fill(P4_P_miss.Mag2(),1.);
if (P4_P_miss.Mag2() > m_proton*m_proton - 0.04 ){
if (P4_P_miss.Mag2() < m_proton*m_proton + 0.04 ){
if ((W > 1.58) && (W < 1.63) && (Q2 > 0.5) && (Q2 < 0.55)) hist_P_miss_en->Fill(P4_P_miss[3],1.);
if (P4_P_miss[3] > m_proton  - 0.05) {
selection_proton_miss = true;
P4_PP_reg = P4_P_miss;
};
};
};
};
};
};
};
};
};
};
};

};




//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=4)){
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if ((bool_proton_id_data)&&(bool_pip_id_data)&&(bool_pim_id_data)){

hist_PIm_miss_all_reg_1-> Fill(P4_PIm_miss.Mag2(),1.);



if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {







if ((abs(dc_z_EL - dc_z_PIp) < 4.) && (abs(dc_z_EL - dc_z_PIm) < 4.) && (abs(dc_z_PIp - dc_z_PIm) < 4.) && (abs(dc_z_PIp - dc_z_P) < 4.) && (abs(dc_z_PIm - dc_z_P) < 4.) && (abs(dc_z_EL - dc_z_P) < 4.)) {


if ((P4_miss_0.Mag2()>-0.04)&&(P4_miss_0.Mag2()<0.04)){
h_miss_mom_0->Fill((P4_miss_0.Vect()).Mag(),1.);


h_miss_mom_0_cut_on0->Fill((P4_miss_0.Vect()).Mag(),1.);

if ((W > 1.58) && (W < 1.63) && (Q2 > 0.5) && (Q2 < 0.55)) hist_miss_en_0->Fill(P4_miss_0[3],1.);

if ((P4_miss_0[3] >  - 0.05)) {

selection_0_miss = true;
if ((Q2>0.4)&&(Q2<0.45)) h_pr_mom_test[int((W-1.3)/0.025)]->Fill((P4_PP_reg.Vect()).Mag());
if ((Q2>0.4)&&(Q2<0.45)) h_pip_mom_test[int((W-1.3)/0.025)]->Fill((P4_PIp_reg.Vect()).Mag());

};
};
};
};
};
};
};
};
};
};

};


//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=3)){
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if ((bool_pip_id_data)&&(bool_proton_id_data)&&(!bool_pim_id_data)){
//if ((bool_pip_id_data)&&(bool_proton_id_data)){

if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) < W - m_pip + 0.05) {

h_mm_pim_vs_npart->Fill(npart,P4_PIm_miss.Mag2(),1.);



if ((abs(dc_z_EL - dc_z_PIp) < 4.) && (abs(dc_z_EL - dc_z_P) < 4.) && (abs(dc_z_PIp - dc_z_P) < 4.)) {

h_pim_mm_q2_w[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]-> Fill(P4_PIm_miss.Mag2(),1.);
hist_PIm_miss-> Fill(P4_PIm_miss.Mag2(),1.);

if ((P4_PIm_miss.Mag2()> -0.04)&&(P4_PIm_miss.Mag2()< 0.06)){

if ((W > 1.58) && (W < 1.63) && (Q2 > 0.5) && (Q2 < 0.55)) hist_PIm_miss_en->Fill(P4_PIm_miss[3],1.);

if((P4_PIm_miss[3] > m_pip  - 0.05)) {
selection_pim_miss = true;
P4_PIm_reg = P4_PIm_miss;
if ((Q2>0.4)&&(Q2<0.45)) h_pr_mom_test[int((W-1.3)/0.025)]->Fill((P4_PP_reg.Vect()).Mag());
if ((Q2>0.4)&&(Q2<0.45)) h_pip_mom_test[int((W-1.3)/0.025)]->Fill((P4_PIp_reg.Vect()).Mag());
};
};

};
};
};
};
};
};

};
};


if ((bool_pim_id_data)&&(bool_proton_id_data)&&(!bool_pip_id_data)){

if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {



if ((abs(dc_z_EL - dc_z_P) < 4.) && (abs(dc_z_EL - dc_z_PIm) < 4.) && (abs(dc_z_P - dc_z_PIm) < 4.)) {

if ((P4_PIp_miss.Mag2()> -0.04)&&(P4_PIp_miss.Mag2()< 0.06)){

if ((W > 1.58) && (W < 1.63) && (Q2 > 0.5) && (Q2 < 0.55)) hist_PIp_miss_en->Fill(P4_PIp_miss[3],1.);

if((P4_PIp_miss[3] > m_pip  - 0.05)) {
selection_pip_miss = true;
P4_PIp_reg = P4_PIp_miss;
};
};

};
};
};
};
};
};
};
};





};


};//konets ifa el cutov


if ((selection_pim_miss)||(selection_0_miss)||(selection_pip_miss)||(selection_proton_miss)) {
//if ((selection_pim_miss)||(selection_0_miss)||(selection_proton_miss)) {
//if ((selection_pim_miss)) {
//if ((selection_pim_miss)) {
//data_hist();
W_2pi_fid_p->Fill(W,Q2,1.);



rot_boost_cmsyst();

//if ((selection_0_miss)) h_delta_w_vs_w_old_data->Fill(W, W-W_old,1.);

if((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){

Var_1[0] = inv_m_pip_p;
Var_1[1] = inv_m_pip_pim;
Var_1[2] = theta_P_cm;
Var_1[3] = phi_P_cm;
Var_1[4] = alpha_PIpPIm_pipf;

Var_2[0] = inv_m_pip_p;
Var_2[1] = inv_m_pip_pim;
Var_2[2] = theta_PIm_cm;
Var_2[3] = phi_PIm_cm;
Var_2[4] = alpha_PPIp_piPIm;

Var_3[0] = inv_m_pim_p;
Var_3[1] = inv_m_pip_pim;
Var_3[2] = theta_PIp_cm;
Var_3[3] = phi_PIp_cm;
Var_3[4] = alpha_PPIm_piPIp;


if (selection_0_miss) {
h_5dim_0_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_0_2[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_0_3[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 

if ((Q2>0.45)&&(Q2<0.5)) {

h_0_mis_radcorr[int((W-1.3)/0.025)]->Fill(P4_miss_0.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);


};

};
if (selection_pim_miss){
h_5dim_pim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pim_2[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pim_3[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 






if ((Q2>0.45)&&(Q2<0.5)) {

if (n_PIm == 0)h_pim_mis_radcorr[int((W-1.3)/0.025)]->Fill(P4_PIm_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);

h_thp_vs_thp[int((W-1.3)/0.025)]->Fill(theta_P_cm,th_P,1.);

};


};
if (selection_pip_miss){
h_5dim_pip_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pip_2[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pip_3[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 

if ((Q2>0.45)&&(Q2<0.5)) {

if (n_PIp == 0)h_pip_mis_radcorr[int((W-1.3)/0.025)]->Fill(P4_PIp_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);


};

};
if (selection_proton_miss){
h_5dim_pr_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pr_2[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pr_3[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 

if ((Q2>0.45)&&(Q2<0.5)) {

if (n_P == 0)h_pr_mis_radcorr[int((W-1.3)/0.025)]->Fill(P4_P_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);


};

};

n_twopions = n_twopions + fract_integ[pmt_hit+1][sector-1][segment];




h_inv_m_pip_pim->Fill(sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)),1.);
h_inv_m_pip_p->Fill(sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)),1.);
h_inv_m_pim_p->Fill(sqrt((P4_PP_reg+P4_PIm_reg)*(P4_PP_reg+P4_PIm_reg)),1.);


};// end W & Q2 cuts



 


};//end selection??




    };// konets tsicla po sobitiyam (po i)
  
  
  t21->Delete();

 finp->Close();  
 
     };//konets tsicla po failam (po m)
    
    
 /////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////////// 
 
 block_curr = 0;
 block_prev = 0;
 m_old=0;
  
 for (m=1; m<=n_files_empty; m++) {
  
  finp = new TFile(file_empty[m-1].c_str()); 
  
  
  
 cout << "Processing file with empty target " << m << "\n"; 
  
  
  // cout << " bin content3 = " << avrg_nphe_sector1->GetBinContent(100,100) << "\n";
  
  TTree *t21 = (TTree*)finp->Get("t21");
  
   TBranch *br_indtype = t21->GetBranch("indtype");
    TBranch *br_npart = t21->GetBranch("npart");
    TBranch *br_p_el = t21->GetBranch("P_EL");
    TBranch *br_block = t21->GetBranch("block");
    TBranch *br_deltaQ = t21->GetBranch("deltaQ");
    TBranch *br_LiveTime = t21->GetBranch("LiveTime");
    TBranch *br_n_incl = t21->GetBranch("n_incl");
    TBranch *br_n_elast = t21->GetBranch("n_elast");
    TBranch *br_ph_EL = t21->GetBranch("ph_EL");
    TBranch *br_th_EL = t21->GetBranch("th_EL");
    TBranch *br_W = t21->GetBranch("W");
    TBranch *br_Q2 = t21->GetBranch("Q2");
    TBranch *br_nphe_el = t21->GetBranch("NpheCC_EL");
    TBranch *br_ectot_el = t21->GetBranch("ECtot_EL");
    TBranch *br_ecin_el = t21->GetBranch("ECin_EL");
    TBranch *br_ecout_el = t21->GetBranch("ECout_EL");
    TBranch *br_x_EL = t21->GetBranch("x_EL");
    TBranch *br_y_EL = t21->GetBranch("y_EL");
    TBranch *br_z_EL = t21->GetBranch("z_EL");
    TBranch *br_dc_z_EL = t21->GetBranch("dc_z_EL"); 
    TBranch *br_dc_z_PIm = t21->GetBranch("dc_z_PIm");
    TBranch *br_dc_z_PIp = t21->GetBranch("dc_z_PIp");
    TBranch *br_dc_z_P = t21->GetBranch("dc_z_P"); 
    TBranch *br_z_P = t21->GetBranch("z_P");
    TBranch *br_z_PIp = t21->GetBranch("z_PIp");
    TBranch *br_z_PIm = t21->GetBranch("z_PIm");
    TBranch *br_PdHit_EL = t21->GetBranch("PdHit_EL");
    TBranch *br_PdHit_PIp = t21->GetBranch("PdHit_PIp");
    TBranch *br_PdHit_PIm = t21->GetBranch("PdHit_PIm");
    TBranch *br_PdHit_P = t21->GetBranch("PdHit_P");
    TBranch *br_sc_x = t21->GetBranch("sc_x");
    TBranch *br_sc_y = t21->GetBranch("sc_y");
    TBranch *br_sc_z = t21->GetBranch("sc_z");
    TBranch *br_pmt_hit = t21->GetBranch("pmt_hit");
    TBranch *br_segment = t21->GetBranch("segment");
    TBranch *br_theta_cc = t21->GetBranch("theta_cc");
    TBranch *br_ph_cc = t21->GetBranch("ph_cc");
    TBranch *br_sector = t21->GetBranch("sector");
    TBranch *br_n_PIp = t21->GetBranch("n_PIp");
    TBranch *br_n_PIm = t21->GetBranch("n_PIm");
    TBranch *br_n_P = t21->GetBranch("n_P");
    TBranch *br_p_pip = t21->GetBranch("P_PIp");
    TBranch *br_p_pim = t21->GetBranch("P_PIm");
    TBranch *br_p_p = t21->GetBranch("P_P");
    TBranch *br_th_PIp = t21->GetBranch("th_PIp");
    TBranch *br_th_PIm = t21->GetBranch("th_PIm");
    TBranch *br_th_P = t21->GetBranch("th_P");
    TBranch *br_ph_PIp = t21->GetBranch("ph_PIp");
    TBranch *br_ph_PIm = t21->GetBranch("ph_PIm");
    TBranch *br_ph_P = t21->GetBranch("ph_P");
    TBranch *br_beta_PIm = t21->GetBranch("beta_PIm");
    TBranch *br_beta_PIp = t21->GetBranch("beta_PIp");
    TBranch *br_beta_P = t21->GetBranch("beta_P");
    TBranch *br_beta_PIm_time = t21->GetBranch("beta_PIm_time");
    TBranch *br_beta_PIp_time = t21->GetBranch("beta_PIp_time");
    TBranch *br_beta_P_time = t21->GetBranch("beta_P_time");   
    TBranch *br_nphe_pip = t21->GetBranch("NpheCC_PIp");
    TBranch *br_PIp_time = t21->GetBranch("PIp_time");
    TBranch *br_PIp_dist = t21->GetBranch("PIp_dist");
    TBranch *br_P_time = t21->GetBranch("P_time");
    TBranch *br_P_dist = t21->GetBranch("P_dist");
    TBranch *br_PIm_time = t21->GetBranch("PIm_time");
    TBranch *br_PIm_dist = t21->GetBranch("PIm_dist");  
    
   
   

  Bool_t adc_cut_switch,tdc_cut_switch;
  
  Qdiff = 0.;
  Qcurr = 0.;
  Qprev = 0.;
  Qtotal = 0.;
  k = 0;
  block = 0;
  last_k = 0;
  nstart = 0;
  nstop = 0;
  n_incl = 0;
  n_elast = 0;
  
  
  for (i=0; i<br_sector->GetEntries(); i++) { 
  
  Qprev = Qcurr;

  br_segment->GetEntry(i);
  br_pmt_hit->GetEntry(i);
  br_deltaQ->GetEntry(i);
  br_n_incl->GetEntry(i);
  br_n_elast->GetEntry(i);  
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_W->GetEntry(i);
  br_Q2->GetEntry(i); 
  br_npart->GetEntry(i);
  br_n_PIm->GetEntry(i);
  br_n_PIp->GetEntry(i);
  br_n_P->GetEntry(i);
  br_sector->GetEntry(i);
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_ph_cc->GetEntry(i);
  br_theta_cc->GetEntry(i);  
  br_p_el->GetEntry(i);  
  br_p_p->GetEntry(i); 
  br_p_pim->GetEntry(i);
  br_p_pip->GetEntry(i); 
  br_ectot_el->GetEntry(i); 
  br_ecin_el->GetEntry(i); 
  br_ecout_el->GetEntry(i);
  br_nphe_el->GetEntry(i);
  br_nphe_pip->GetEntry(i);
  br_th_EL->GetEntry(i);
  br_ph_EL->GetEntry(i); 
  br_th_P->GetEntry(i);
  br_ph_P->GetEntry(i);
  br_th_PIp->GetEntry(i);
  br_ph_PIp->GetEntry(i);
  br_th_PIm->GetEntry(i);
  br_ph_PIm->GetEntry(i);
  br_beta_PIp->GetEntry(i);
  br_beta_PIm->GetEntry(i);
  br_beta_P->GetEntry(i);
  br_beta_PIp_time->GetEntry(i);
  br_beta_PIm_time->GetEntry(i);
  br_beta_P_time->GetEntry(i);  
  br_PdHit_EL->GetEntry(i);
  br_PdHit_PIp->GetEntry(i);
  br_PdHit_PIm->GetEntry(i);
  br_PdHit_P->GetEntry(i);
  br_sc_x->GetEntry(i);
  br_sc_y->GetEntry(i);
  br_sc_z->GetEntry(i);  
  br_z_EL->GetEntry(i);
  br_dc_z_EL->GetEntry(i);  
  br_dc_z_PIm->GetEntry(i);   
  br_dc_z_PIp->GetEntry(i); 
  br_dc_z_P->GetEntry(i);        
  br_W->GetEntry(i);
  br_Q2->GetEntry(i);
  br_indtype->GetEntry(i);
  br_PIp_time->GetEntry(i);
  br_PIm_time->GetEntry(i);
  br_P_time->GetEntry(i); 
  br_PIp_dist->GetEntry(i);
  br_PIm_dist->GetEntry(i);
  br_P_dist->GetEntry(i);
  
  
  
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  P_P = br_p_p->GetLeaf("P_P")->GetValue();
  P_PIp = br_p_pip->GetLeaf("P_PIp")->GetValue();
  P_PIm = br_p_pim->GetLeaf("P_PIm")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();
  th_P = br_th_P->GetLeaf("th_P")->GetValue();
  ph_P = br_ph_P->GetLeaf("ph_P")->GetValue();
  th_PIp = br_th_PIp->GetLeaf("th_PIp")->GetValue();
  ph_PIp = br_ph_PIp->GetLeaf("ph_PIp")->GetValue();
  npart = br_npart->GetLeaf("npart")->GetValue();
  n_P = br_n_P->GetLeaf("n_P")->GetValue();
  n_PIp = br_n_PIp->GetLeaf("n_PIp")->GetValue();
  Nphe_pip = br_nphe_pip->GetLeaf("NpheCC_PIp")->GetValue();
  PdHit_EL = br_PdHit_EL->GetLeaf("PdHit_EL")->GetValue();
  PdHit_PIp = br_PdHit_PIp->GetLeaf("PdHit_PIp")->GetValue();
  PdHit_PIm = br_PdHit_PIm->GetLeaf("PdHit_PIm")->GetValue();
  PdHit_P = br_PdHit_P->GetLeaf("PdHit_P")->GetValue();
  sc_x = br_sc_x->GetLeaf("sc_x")->GetValue();
  sc_y = br_sc_y->GetLeaf("sc_y")->GetValue();  
  sc_z = br_sc_z->GetLeaf("sc_z")->GetValue();    
  beta_PIp = br_beta_PIp_time->GetLeaf("beta_PIp_time")->GetValue();
  beta_PIm = br_beta_PIm_time->GetLeaf("beta_PIm_time")->GetValue();
  beta_P = br_beta_P_time->GetLeaf("beta_P_time")->GetValue();
  n_PIm = br_n_PIm->GetLeaf("n_PIm")->GetValue();
  ph_PIm = br_ph_PIm->GetLeaf("ph_PIm")->GetValue();
  th_PIm = br_th_PIm->GetLeaf("th_PIm")->GetValue();
  z_EL = br_z_EL->GetLeaf("z_EL")->GetValue();
  dc_z_EL = br_dc_z_EL->GetLeaf("dc_z_EL")->GetValue(); 
  dc_z_PIm = br_dc_z_PIm->GetLeaf("dc_z_PIm")->GetValue();    
  dc_z_PIp = br_dc_z_PIp->GetLeaf("dc_z_PIp")->GetValue();  
  dc_z_P = br_dc_z_P->GetLeaf("dc_z_P")->GetValue();     
  W = br_W->GetLeaf("W")->GetValue();
  Q2 = br_Q2->GetLeaf("Q2")->GetValue();
  segment = br_segment->GetLeaf("segment")->GetValue();
  pmt_hit = br_pmt_hit->GetLeaf("pmt_hit")->GetValue();

  PIp_time = br_PIp_time->GetLeaf("PIp_time")->GetValue();
  PIm_time = br_PIm_time->GetLeaf("PIm_time")->GetValue();
  P_time = br_P_time->GetLeaf("P_time")->GetValue();
  PIp_dist = br_PIp_dist->GetLeaf("PIp_dist")->GetValue();
  PIm_dist = br_PIm_dist->GetLeaf("PIm_dist")->GetValue();
  P_dist = br_P_dist->GetLeaf("P_dist")->GetValue();
  indtype = br_indtype->GetLeaf("indtype")->GetValue();
  
  
  if (br_ectot_el->GetLeaf("ECtot_EL")->GetValue() > (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue())) { 
          ECT  = br_ectot_el->GetLeaf("ECtot_EL")->GetValue();  
         } else {
          ECT  = (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue());    
          };
	  
	  //quatity check   
	     
	     block_curr = br_block->GetLeaf("block")->GetValue(); 

    if ((m_old != 0) && (m != m_old))block_total=block_total+block_prev;
    
    m_old=m;


   LiveTime=br_LiveTime->GetLeaf("LiveTime")->GetValue();
   
   inclusive = br_n_incl->GetLeaf("n_incl")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue();
   elastic =  br_n_elast->GetLeaf("n_elast")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue(); 

   if (block_curr != block_prev) {


 
 Qfull_empty = Qfull_empty + br_deltaQ->GetLeaf("deltaQ")->GetValue();
 
   
  
   
   block_prev = block_curr;
   };   
	  
  
 
sector =  int(br_sector->GetLeaf("sector")->GetValue());
nphe = br_nphe_el->GetLeaf("NpheCC_EL")->GetValue();
theta_cc = br_theta_cc->GetLeaf("theta_cc")->GetValue();
ph_cc = br_ph_cc->GetLeaf("ph_cc")->GetValue();  
  


P4_ELP_reg.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),P_EL);

P4_PP_reg.SetXYZT(P_P*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*cos(th_P*M_PI/180.),sqrt(m_proton*m_proton+P_P*P_P));

P4_PIp_reg.SetXYZT(P_PIp*cos(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*sin(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*cos(th_PIp*M_PI/180.),sqrt(m_pip*m_pip+P_PIp*P_PIp));


P4_PIm_reg.SetXYZT(P_PIm*cos(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*sin(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*cos(th_PIm*M_PI/180.),sqrt(m_pip*m_pip+P_PIm*P_PIm));


P4_PIm_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;
P4_PIp_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_PIp_miss_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_P_miss = P4_EL + P4_P - P4_ELP_reg - P4_PIp_reg - P4_PIm_reg;
P4_miss_0 = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;

//W_old = W;
//W = sqrt((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2());

th_PIm_miss = (180./M_PI)*acos(P4_PIm_miss[2]/sqrt(P4_PIm_miss[0]*P4_PIm_miss[0]+P4_PIm_miss[1]*P4_PIm_miss[1]+P4_PIm_miss[2]*P4_PIm_miss[2]));


if (P4_PIm_miss[0] != 0.) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
}
 else {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;
};
				   

if ((P4_PIm_miss[0] < 0.) && (P4_PIm_miss[1] > 0)) ph_PIm_miss = ph_PIm_miss+180.;
if (( P4_PIm_miss[0]< 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+180.;
if ((P4_PIm_miss[0] > 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+360.;

//beta_func_empty();

bool_el_id_empt=particle_ID_empty.Electron_cuts_empty(); 
bool_proton_id_empt=particle_ID_empty.Proton_cuts_empty();
bool_pip_id_empt=particle_ID_empty.PIp_cuts_empty();
bool_pim_id_empt=particle_ID_empty.PIm_cuts_empty();


selection = false;
selection_pim_miss_empt = false;
selection_pip_miss_empt = false;
selection_proton_miss_empt = false;
selection_0_miss_empt = false;

if (bool_el_id_empt) {

beta_func_empty();

//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart >=3)){
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if ((bool_pip_id_empt)&&(bool_pim_id_empt)&&(!bool_proton_id_empt)){

if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) < W - m_pip + 0.05) {
if (P4_P_miss[3] > m_proton  - 0.05) {
if ((abs(dc_z_EL - dc_z_PIp) < 4.) && (abs(dc_z_EL - dc_z_PIm) < 4.) && (abs(dc_z_PIp - dc_z_PIm) < 4.)) {
if (P4_P_miss.Mag2() >m_proton*m_proton - 0.04 ){
if (P4_P_miss.Mag2() <m_proton*m_proton + 0.04 ){
selection_proton_miss_empt = true;
P4_PP_reg = P4_P_miss;
};
};
};
};
};
};
};
};
};
};

};
};

//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=4)){
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if ((bool_pip_id_empt)&&(bool_proton_id_empt)&&(bool_pim_id_empt)){

/*
if (sector ==1)  hist_z_el_1_empty->Fill(z_EL-0.15,1.); 
if (sector ==2)  hist_z_el_2_empty->Fill(z_EL-0.15,1.);
if (sector ==3)  hist_z_el_3_empty->Fill(z_EL-0.15,1.);
if (sector ==4)  hist_z_el_4_empty->Fill(z_EL-0.15,1.);
if (sector ==5)  hist_z_el_5_empty->Fill(z_EL-0.15,1.);
if (sector ==6)  hist_z_el_6_empty->Fill(z_EL-0.15,1.);
*/


if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {

if ((P4_miss_0[3] >  - 0.05)) {
if ((abs(dc_z_EL - dc_z_PIp) < 4.) && (abs(dc_z_EL - dc_z_PIm) < 4.) && (abs(dc_z_PIp - dc_z_PIm) < 4.) && (abs(dc_z_PIp - dc_z_P) < 4.) && (abs(dc_z_PIm - dc_z_P) < 4.) && (abs(dc_z_EL - dc_z_P) < 4.)) {
if ((P4_miss_0.Mag2()>-0.04)&&(P4_miss_0.Mag2()<0.04)){
selection_0_miss_empt = true;
};
};
};
};
};
};
};
};
};
};
};


//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart >=3)){
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if ((bool_pip_id_empt)&&(bool_proton_id_empt)&&(!bool_pim_id_empt)){
//if ((bool_pip_id_empt)&&(bool_proton_id_empt)){
if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) < W - m_pip + 0.05) {


if (P4_PIm_miss[3] > m_pip  - 0.05) {
if ((abs(dc_z_EL - dc_z_PIp) < 4.) && (abs(dc_z_EL - dc_z_P) < 4.) && (abs(dc_z_PIp - dc_z_P) < 4.)) {
if ((P4_PIm_miss.Mag2()> -0.04)&&(P4_PIm_miss.Mag2()< 0.06)){
selection_pim_miss_empt = true;
P4_PIm_reg = P4_PIm_miss;
};
};
}
};
};
};
};
};
};
};

if ((bool_pim_id_empt)&&(bool_proton_id_empt)&&(!bool_pip_id_empt)){
if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {



if (P4_PIp_miss[3] > m_pip  - 0.05) {
if ((abs(dc_z_EL - dc_z_P) < 4.) && (abs(dc_z_EL - dc_z_PIm) < 4.) && (abs(dc_z_P - dc_z_PIm) < 4.)) {
if ((P4_PIp_miss.Mag2()> -0.04)&&(P4_PIp_miss.Mag2()< 0.06)){
selection_pip_miss_empt = true;
P4_PIp_reg = P4_PIp_miss;
};
};
};
};
};
};
};
};
};
};

};



if ((selection_pim_miss_empt)||(selection_0_miss_empt)||(selection_pip_miss_empt)||(selection_proton_miss_empt)) {
//if ((selection_pim_miss_empt)||(selection_0_miss_empt)||(selection_proton_miss_empt)) {
//if ((selection_pim_miss_empt)) {
//if ((selection_pim_miss_empt)) {
rot_boost_cmsyst();

if((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){

Var_1[0] = inv_m_pip_p;
Var_1[1] = inv_m_pip_pim;
Var_1[2] = theta_P_cm;
Var_1[3] = phi_P_cm;
Var_1[4] = alpha_PIpPIm_pipf;

Var_2[0] = inv_m_pip_p;
Var_2[1] = inv_m_pip_pim;
Var_2[2] = theta_PIm_cm;
Var_2[3] = phi_PIm_cm;
Var_2[4] = alpha_PPIp_piPIm;

Var_3[0] = inv_m_pim_p;
Var_3[1] = inv_m_pip_pim;
Var_3[2] = theta_PIp_cm;
Var_3[3] = phi_PIp_cm;
Var_3[4] = alpha_PPIm_piPIp;


if (selection_0_miss_empt) {
h_5dim_0_1_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_0_2_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_0_3_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 
};
if (selection_pim_miss_empt) {
h_5dim_pim_1_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pim_2_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pim_3_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 
};
if (selection_pip_miss_empt) {
h_5dim_pip_1_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pip_2_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pip_3_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 
};
if (selection_proton_miss_empt) {
h_5dim_pr_1_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pr_2_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pr_3_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 
};
/*if ((selection_0_miss_empt)) {
if ((W > 1.4)&&(W < 1.625)&&(Q2 > 0.4)&&(Q2 < 0.9)) {

h_5dim_excl_1_empty[int((W-1.4)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_excl_2_empty[int((W-1.4)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_excl_3_empty[int((W-1.4)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]);
};
};*/


};
};







 }; //konec ifa electronnih cutov
    };// konets tsicla po sobitiyam (po i)
  
  
  t21->Delete();

 finp->Close();  
 
     };//konets tsicla po failam (po m)  
  
  
  
  
  
  
//////////////////////////////////////////-sim  
 block_curr = 0;
 block_prev = 0;
 m_old=0;
 
for (m=1; m<=n_files_sim; m++) {
  
  finp = new TFile(file_sim[m-1].c_str()); 
  
  
  
 cout << "Processing file with simulation " << m << "\n"; 
  
  
  // cout << " bin content3 = " << avrg_nphe_sector1->GetBinContent(100,100) << "\n";
    
  TTree *t21 = (TTree*)finp->Get("t21");
  
    TBranch *br_indtype = t21->GetBranch("indtype");
    TBranch *br_npart = t21->GetBranch("npart");
    TBranch *br_p_el = t21->GetBranch("P_EL");
    TBranch *br_block = t21->GetBranch("block");
    TBranch *br_deltaQ = t21->GetBranch("deltaQ");
    TBranch *br_LiveTime = t21->GetBranch("LiveTime");
    TBranch *br_n_incl = t21->GetBranch("n_incl");
    TBranch *br_n_elast = t21->GetBranch("n_elast");
    TBranch *br_ph_EL = t21->GetBranch("ph_EL");
    TBranch *br_th_EL = t21->GetBranch("th_EL");
    TBranch *br_W = t21->GetBranch("W");
    TBranch *br_Q2 = t21->GetBranch("Q2");
    TBranch *br_nphe_el = t21->GetBranch("NpheCC_EL");
    TBranch *br_ectot_el = t21->GetBranch("ECtot_EL");
    TBranch *br_ecin_el = t21->GetBranch("ECin_EL");
    TBranch *br_ecout_el = t21->GetBranch("ECout_EL");
    TBranch *br_x_EL = t21->GetBranch("x_EL");
    TBranch *br_y_EL = t21->GetBranch("y_EL");
    TBranch *br_z_EL = t21->GetBranch("z_EL");
    TBranch *br_dc_z_EL = t21->GetBranch("dc_z_EL");
    TBranch *br_dc_z_PIm = t21->GetBranch("dc_z_PIm");
    TBranch *br_dc_z_PIp = t21->GetBranch("dc_z_PIp");
    TBranch *br_dc_z_P = t21->GetBranch("dc_z_P"); 
    TBranch *br_z_P = t21->GetBranch("z_P");
    TBranch *br_z_PIp = t21->GetBranch("z_PIp");
    TBranch *br_z_PIm = t21->GetBranch("z_PIm");
    TBranch *br_PdHit_EL = t21->GetBranch("PdHit_EL");
    TBranch *br_PdHit_PIp = t21->GetBranch("PdHit_PIp");
    TBranch *br_PdHit_PIm = t21->GetBranch("PdHit_PIm");
    TBranch *br_PdHit_P = t21->GetBranch("PdHit_P");
    TBranch *br_sc_x = t21->GetBranch("sc_x");
    TBranch *br_sc_y = t21->GetBranch("sc_y");
    TBranch *br_sc_z = t21->GetBranch("sc_z");
    TBranch *br_pmt_hit = t21->GetBranch("pmt_hit");
    TBranch *br_segment = t21->GetBranch("segment");
    TBranch *br_theta_cc = t21->GetBranch("theta_cc");
    TBranch *br_ph_cc = t21->GetBranch("ph_cc");
    TBranch *br_sector = t21->GetBranch("sector");
    TBranch *br_n_PIp = t21->GetBranch("n_PIp");
    TBranch *br_n_PIm = t21->GetBranch("n_PIm");
    TBranch *br_n_P = t21->GetBranch("n_P");
    TBranch *br_p_pip = t21->GetBranch("P_PIp");
    TBranch *br_p_pim = t21->GetBranch("P_PIm");
    TBranch *br_p_p = t21->GetBranch("P_P");
    TBranch *br_th_PIp = t21->GetBranch("th_PIp");
    TBranch *br_th_PIm = t21->GetBranch("th_PIm");
    TBranch *br_th_P = t21->GetBranch("th_P");
    TBranch *br_ph_PIp = t21->GetBranch("ph_PIp");
    TBranch *br_ph_PIm = t21->GetBranch("ph_PIm");
    TBranch *br_ph_P = t21->GetBranch("ph_P");
    TBranch *br_beta_PIm = t21->GetBranch("beta_PIm");
    TBranch *br_beta_PIp = t21->GetBranch("beta_PIp");
    TBranch *br_beta_P = t21->GetBranch("beta_P");
    TBranch *br_beta_PIm_time = t21->GetBranch("beta_PIm_time");
    TBranch *br_beta_PIp_time = t21->GetBranch("beta_PIp_time");
    TBranch *br_beta_P_time = t21->GetBranch("beta_P_time");    
    TBranch *br_nphe_pip = t21->GetBranch("NpheCC_PIp");
    TBranch *br_PIp_time = t21->GetBranch("PIp_time");
    TBranch *br_PIp_dist = t21->GetBranch("PIp_dist");
    TBranch *br_P_time = t21->GetBranch("P_time");
    TBranch *br_P_dist = t21->GetBranch("P_dist");
    TBranch *br_PIm_time = t21->GetBranch("PIm_time");
    TBranch *br_PIm_dist = t21->GetBranch("PIm_dist");  
//    TBranch *br_pf_x = t21->GetBranch("pf_x");  
//    TBranch *br_pf_y = t21->GetBranch("pf_y"); 
 //   TBranch *br_pf_z = t21->GetBranch("pf_z"); 

  Bool_t adc_cut_switch,tdc_cut_switch;
  
  Qdiff = 0.;
  Qcurr = 0.;
  Qprev = 0.;
  Qtotal = 0.;
  k = 0;
  block = 0;
  last_k = 0;
  nstart = 0;
  nstop = 0;
  n_incl = 0;
  n_elast = 0;
  
  
  for (i=0; i<br_sector->GetEntries(); i++) { 
  
  Qprev = Qcurr;

br_segment->GetEntry(i);
br_pmt_hit->GetEntry(i);
 br_deltaQ->GetEntry(i);
  br_n_incl->GetEntry(i);
  br_n_elast->GetEntry(i);  
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_W->GetEntry(i);
  br_Q2->GetEntry(i); 
  br_npart->GetEntry(i);
  br_n_PIm->GetEntry(i);
  br_n_PIp->GetEntry(i);
  br_n_P->GetEntry(i);
  br_sector->GetEntry(i);
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_ph_cc->GetEntry(i);
  br_theta_cc->GetEntry(i);  
  br_p_el->GetEntry(i);  
  br_p_p->GetEntry(i); 
  br_p_pim->GetEntry(i);
  br_p_pip->GetEntry(i); 
  br_ectot_el->GetEntry(i); 
  br_ecin_el->GetEntry(i); 
  br_ecout_el->GetEntry(i);
  br_nphe_el->GetEntry(i);
  br_nphe_pip->GetEntry(i);
  br_th_EL->GetEntry(i);
  br_ph_EL->GetEntry(i); 
  br_th_P->GetEntry(i);
  br_ph_P->GetEntry(i);
  br_th_PIp->GetEntry(i);
  br_ph_PIp->GetEntry(i);
  br_th_PIm->GetEntry(i);
  br_ph_PIm->GetEntry(i);
  br_beta_PIp->GetEntry(i);
  br_beta_PIm->GetEntry(i);
  br_beta_P->GetEntry(i);
  br_beta_PIp_time->GetEntry(i);
  br_beta_PIm_time->GetEntry(i);
  br_beta_P_time->GetEntry(i);  
  br_PdHit_EL->GetEntry(i);
  br_PdHit_PIp->GetEntry(i);
  br_PdHit_PIm->GetEntry(i);
  br_PdHit_P->GetEntry(i);
  br_sc_x->GetEntry(i); 
  br_sc_y->GetEntry(i);  
  br_sc_z->GetEntry(i);     
  br_z_EL->GetEntry(i);
  br_dc_z_EL->GetEntry(i);
  br_dc_z_PIm->GetEntry(i);  
  br_dc_z_PIp->GetEntry(i); 
  br_dc_z_P->GetEntry(i);      
  br_W->GetEntry(i);
  br_Q2->GetEntry(i);
  br_indtype->GetEntry(i);
  
  
  br_PIp_time->GetEntry(i);
  br_PIm_time->GetEntry(i);
  br_P_time->GetEntry(i); 
  br_PIp_dist->GetEntry(i);
  br_PIm_dist->GetEntry(i);
  br_P_dist->GetEntry(i);
//   br_pf_x->GetEntry(i);
 // br_pf_y->GetEntry(i);
//  br_pf_z->GetEntry(i);
  
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  P_P = br_p_p->GetLeaf("P_P")->GetValue();
  P_PIp = br_p_pip->GetLeaf("P_PIp")->GetValue();
  P_PIm = br_p_pim->GetLeaf("P_PIm")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();
  th_P = br_th_P->GetLeaf("th_P")->GetValue();
  ph_P = br_ph_P->GetLeaf("ph_P")->GetValue();
  th_PIp = br_th_PIp->GetLeaf("th_PIp")->GetValue();
  ph_PIp = br_ph_PIp->GetLeaf("ph_PIp")->GetValue();
  npart = br_npart->GetLeaf("npart")->GetValue();
  n_P = br_n_P->GetLeaf("n_P")->GetValue();
  n_PIp = br_n_PIp->GetLeaf("n_PIp")->GetValue();
  Nphe_pip = br_nphe_pip->GetLeaf("NpheCC_PIp")->GetValue();
  PdHit_EL = br_PdHit_EL->GetLeaf("PdHit_EL")->GetValue();
  PdHit_PIp = br_PdHit_PIp->GetLeaf("PdHit_PIp")->GetValue();
  PdHit_PIm = br_PdHit_PIm->GetLeaf("PdHit_PIm")->GetValue();
  PdHit_P = br_PdHit_P->GetLeaf("PdHit_P")->GetValue();
  sc_x = br_sc_x->GetLeaf("sc_x")->GetValue(); 
  sc_y = br_sc_y->GetLeaf("sc_y")->GetValue();   
  sc_z = br_sc_z->GetLeaf("sc_z")->GetValue();   
  beta_PIp = br_beta_PIp_time->GetLeaf("beta_PIp_time")->GetValue();
  beta_PIm = br_beta_PIm_time->GetLeaf("beta_PIm_time")->GetValue();
  beta_P = br_beta_P_time->GetLeaf("beta_P_time")->GetValue();
  n_PIm = br_n_PIm->GetLeaf("n_PIm")->GetValue();
  ph_PIm = br_ph_PIm->GetLeaf("ph_PIm")->GetValue();
  th_PIm = br_th_PIm->GetLeaf("th_PIm")->GetValue();
  z_EL = br_z_EL->GetLeaf("z_EL")->GetValue();
  dc_z_EL = br_dc_z_EL->GetLeaf("dc_z_EL")->GetValue();
  dc_z_PIm = br_dc_z_PIm->GetLeaf("dc_z_PIm")->GetValue();  
  dc_z_PIp = br_dc_z_PIp->GetLeaf("dc_z_PIp")->GetValue();   
  dc_z_P = br_dc_z_P->GetLeaf("dc_z_P")->GetValue();    
  W = br_W->GetLeaf("W")->GetValue();
  Q2 = br_Q2->GetLeaf("Q2")->GetValue();
  segment = br_segment->GetLeaf("segment")->GetValue();
  pmt_hit = br_pmt_hit->GetLeaf("pmt_hit")->GetValue();
  
 
  PIp_time = br_PIp_time->GetLeaf("PIp_time")->GetValue();
  PIm_time = br_PIm_time->GetLeaf("PIm_time")->GetValue();
  P_time = br_P_time->GetLeaf("P_time")->GetValue();
  PIp_dist = br_PIp_dist->GetLeaf("PIp_dist")->GetValue();
  PIm_dist = br_PIm_dist->GetLeaf("PIm_dist")->GetValue();
  P_dist = br_P_dist->GetLeaf("P_dist")->GetValue();
  indtype = br_indtype->GetLeaf("indtype")->GetValue();

  
  if (br_ectot_el->GetLeaf("ECtot_EL")->GetValue() > (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue())) { 
          ECT  = br_ectot_el->GetLeaf("ECtot_EL")->GetValue();  
         } else {
          ECT  = (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue());    
          };
	  
	  //quatity check   
	     
	     block_curr = br_block->GetLeaf("block")->GetValue(); 

    if ((m_old != 0) && (m != m_old))block_total=block_total+block_prev;
    
    m_old=m;


   LiveTime=br_LiveTime->GetLeaf("LiveTime")->GetValue();
   
   inclusive = br_n_incl->GetLeaf("n_incl")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue();
   elastic =  br_n_elast->GetLeaf("n_elast")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue(); 

   if (block_curr != block_prev) {


 
 Qfull_sim = Qfull_sim + br_deltaQ->GetLeaf("deltaQ")->GetValue();
 
   
  
   
   block_prev = block_curr;
   };   
	  
  
 
sector =  int(br_sector->GetLeaf("sector")->GetValue());
nphe = br_nphe_el->GetLeaf("NpheCC_EL")->GetValue();
theta_cc = br_theta_cc->GetLeaf("theta_cc")->GetValue();
ph_cc = br_ph_cc->GetLeaf("ph_cc")->GetValue();  
  
//fermi_bonn();

P4_P.SetXYZT(0,0,0,m_proton);

P4_ELP_reg.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),P_EL);

P4_PP_reg.SetXYZT(P_P*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*cos(th_P*M_PI/180.),sqrt(m_proton*m_proton+P_P*P_P));

P4_PIp_reg.SetXYZT(P_PIp*cos(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*sin(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*cos(th_PIp*M_PI/180.),sqrt(m_pip*m_pip+P_PIp*P_PIp));


P4_PIm_reg.SetXYZT(P_PIm*cos(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*sin(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*cos(th_PIm*M_PI/180.),sqrt(m_pip*m_pip+P_PIm*P_PIm));






P4_PIm_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;

P4_PIp_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_PIp_miss_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_P_miss = P4_EL + P4_P - P4_ELP_reg - P4_PIp_reg - P4_PIm_reg;
P4_miss_0 = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;
P4_miss_0_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;

//W_old = W;
//W = sqrt((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2());


//if ((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2()>0) W=sqrt((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2());
//if ((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2()<=0) W=0;




th_PIm_miss = (180./M_PI)*acos(P4_PIm_miss[2]/sqrt(P4_PIm_miss[0]*P4_PIm_miss[0]+P4_PIm_miss[1]*P4_PIm_miss[1]+P4_PIm_miss[2]*P4_PIm_miss[2]));
/*if(P4_PIm_miss[0] == 0.) {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;

};


if ((P4_PIm_miss[1]>0)&&(P4_PIm_miss[0]>0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
};


if ((P4_PIm_miss[1]>0)&&(P4_PIm_miss[0]<0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 180. - (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};

if ((P4_PIm_miss[1]<0)&&(P4_PIm_miss[0]<0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 180. + (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};

if ((P4_PIm_miss[1]<0)&&(P4_PIm_miss[0]>0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 360. - (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};*/

if (P4_PIm_miss[0] != 0.) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
}
 else {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;
};
				   

if ((P4_PIm_miss[0] < 0.) && (P4_PIm_miss[1] > 0)) ph_PIm_miss = ph_PIm_miss+180.;
if (( P4_PIm_miss[0]< 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+180.;
if ((P4_PIm_miss[0] > 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+360.;

if (indtype==2) h_delta_w_vs_w_old_gen->Fill(W,W-W_old,1.);


if (sector ==1) {
if (indtype ==2) hist_z_el_1_sim_2->Fill(z_EL,1.);
 };

if (sector ==2)  {
if (indtype ==2) hist_z_el_2_sim_2->Fill(z_EL,1.);
};

if (sector ==3) {
if (indtype ==2)  hist_z_el_3_sim_2->Fill(z_EL,1.);
 };

if (sector ==4)  {
if (indtype ==2) hist_z_el_4_sim_2->Fill(z_EL,1.);
};

 if (sector ==5) {
if (indtype ==2)   hist_z_el_5_sim_2->Fill(z_EL,1.);
 };

if (sector ==6) {
if (indtype ==2)  hist_z_el_6_sim_2->Fill(z_EL,1.);
};



if((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if (indtype==2){
rot_boost_cmsyst();

Var_1[0] = inv_m_pip_p;
Var_1[1] = inv_m_pip_pim;
Var_1[2] = theta_P_cm;
Var_1[3] = phi_P_cm;
Var_1[4] = alpha_PIpPIm_pipf;

Var_2[0] = inv_m_pip_p;
Var_2[1] = inv_m_pip_pim;
Var_2[2] = theta_PIm_cm;
Var_2[3] = phi_PIm_cm;
Var_2[4] = alpha_PPIp_piPIm;

Var_3[0] = inv_m_pim_p;
Var_3[1] = inv_m_pip_pim;
Var_3[2] = theta_PIp_cm;
Var_3[3] = phi_PIp_cm;
Var_3[4] = alpha_PPIm_piPIp;

h_5dim_1_sim_2[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,1.);
h_5dim_2_sim_2[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,1.);
h_5dim_3_sim_2[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,1.); 

/*if((W > 1.4)&&(W < 1.625)&&(Q2 > 0.4)&&(Q2 < 0.9)){

h_5dim_1_sim_excl_2[int((W-1.4)/0.025)]->Fill(Var_1,1.);
h_5dim_2_sim_excl_2[int((W-1.4)/0.025)]->Fill(Var_2,1.);
h_5dim_3_sim_excl_2[int((W-1.4)/0.025)]->Fill(Var_3,1.);
};*/

};
};


bool_el_id_sim=particle_ID_sim.Electron_cuts_sim();
bool_proton_id_sim=particle_ID_sim.Proton_cuts_sim();
bool_pip_id_sim=particle_ID_sim.PIp_cuts_sim();
bool_pim_id_sim=particle_ID_sim.PIm_cuts_sim();

selection = false;
selection_pim_miss_sim = false;
selection_pip_miss_sim = false;
selection_proton_miss_sim = false;
selection_0_miss_sim = false;

if (bool_el_id_sim) {





if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)) {
if ((bool_pip_id_sim)&&(bool_pim_id_sim)&&(!bool_proton_id_sim)){


if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) < W - m_pip + 0.05) {


if ((abs(dc_z_EL - dc_z_PIp) < 4.) && (abs(dc_z_EL - dc_z_PIm) < 4.) && (abs(dc_z_PIp - dc_z_PIm) < 4.)) {
if (P4_P_miss.Mag2() > m_proton*m_proton - 0.04 ){
if (P4_P_miss.Mag2() < m_proton*m_proton + 0.04 ){

if ((W > 1.58) && (W < 1.63) && (Q2 > 0.5) && (Q2 < 0.55)) hist_P_miss_en_sim->Fill(P4_P_miss[3],1.);

if (P4_P_miss[3] > m_proton  - 0.05) {

selection_proton_miss_sim = true;
P4_PP_reg = P4_P_miss;
};
};

};
};
};
};
};
};
};
};

};
};
//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=4)) {
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)) {
if ((bool_pip_id_sim)&&(bool_proton_id_sim)&&(bool_pim_id_sim)){

if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {

//P4_inprot_miss =-(P4_EL - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg);

//h_mm_0_vs_npart_sim->Fill(npart,P4_miss_0.Mag2(),1.);



if ((abs(dc_z_EL - dc_z_PIp) < 4.) && (abs(dc_z_EL - dc_z_PIm) < 4.) && (abs(dc_z_PIp - dc_z_PIm) < 4.) && (abs(dc_z_PIp - dc_z_P) < 4.) && (abs(dc_z_PIm - dc_z_P) < 4.) && (abs(dc_z_EL - dc_z_P) < 4.)) {



h_miss_mass_0_d_sim-> Fill(P4_miss_0_d.Mag2(),1.);


if ((P4_miss_0.Mag2()>-0.04)&&(P4_miss_0.Mag2()<0.04)){

if ((W > 1.58) && (W < 1.63) && (Q2 > 0.5) && (Q2 < 0.55)) h_miss_en_0_sim->Fill(P4_miss_0[3],1.);

if ((P4_miss_0[3] >  - 0.05)) {

selection_0_miss_sim = true;

if ((Q2>0.4)&&(Q2<0.45)) h_pr_mom_sim_test[int((W-1.3)/0.025)]->Fill((P4_PP_reg.Vect()).Mag());
if ((Q2>0.4)&&(Q2<0.45)) h_pip_mom_sim_test[int((W-1.3)/0.025)]->Fill((P4_PIp_reg.Vect()).Mag());
};
};
};
};
};
};
};
};
};
};

};

//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=3)){
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if ((bool_pip_id_sim)&&(bool_proton_id_sim)&&(!bool_pim_id_sim)){
//if ((bool_pip_id_sim)&&(bool_proton_id_sim)){

if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) < W - m_pip + 0.05) {

h_mm_pim_vs_npart_sim->Fill(npart,P4_PIm_miss.Mag2(),1.);





if ((abs(dc_z_EL - dc_z_PIp) < 4.) && (abs(dc_z_EL - dc_z_P) < 4.) && (abs(dc_z_PIp - dc_z_P) < 4.)) {
h_pim_mm_q2_w_sim[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]-> Fill(P4_PIm_miss.Mag2(),1.);

h_PIm_miss_sim-> Fill(P4_PIm_miss.Mag2(),1.);
if ((P4_PIm_miss.Mag2()>-0.04)&&(P4_PIm_miss.Mag2()<0.06)){ 

if ((W > 1.58) && (W < 1.63) && (Q2 > 0.5) && (Q2 < 0.55)) h_PIm_miss_en_sim->Fill(P4_PIm_miss[3],1.);

if (P4_PIm_miss[3] > m_pip  - 0.05) {

selection_pim_miss_sim= true;
P4_PIm_reg = P4_PIm_miss;

if ((Q2>0.4)&&(Q2<0.45)) h_pr_mom_sim_test[int((W-1.3)/0.025)]->Fill((P4_PP_reg.Vect()).Mag());
if ((Q2>0.4)&&(Q2<0.45)) h_pip_mom_sim_test[int((W-1.3)/0.025)]->Fill((P4_PIp_reg.Vect()).Mag());
};
};
};


};
};
};
};
};
};
};

if ((bool_pim_id_sim)&&(bool_proton_id_sim)&&(!bool_pip_id_sim)){

if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {




if ((abs(dc_z_EL - dc_z_PIm) < 4.) && (abs(dc_z_EL - dc_z_P) < 4.) && (abs(dc_z_PIm - dc_z_P) < 4.)) {
if ((P4_PIp_miss.Mag2()>-0.04)&&(P4_PIp_miss.Mag2()<0.06)){ 


if ((W > 1.58) && (W < 1.63) && (Q2 > 0.5) && (Q2 < 0.55)) h_PIp_miss_en_sim->Fill(P4_PIp_miss[3],1.);

if (P4_PIp_miss[3] > m_pip  - 0.05) {

selection_pip_miss_sim= true;
P4_PIp_reg = P4_PIp_miss;
};
};
};


};
};
};
};
};
};
};




};





};//konets ifa electronnix cutov


if ((selection_pim_miss_sim)||(selection_0_miss_sim)||(selection_pip_miss_sim)||(selection_proton_miss_sim)) {
//if ((selection_pim_miss_sim)||(selection_0_miss_sim)||(selection_proton_miss_sim)) {
//if ((selection_pim_miss_sim)) {
//if ((selection_pim_miss_sim)) {

/*
if (sector ==1) hist_z_el_1_sim_1->Fill(dc_z_EL,1.);
if (sector ==2) hist_z_el_2_sim_1->Fill(dc_z_EL,1.);
if (sector ==3) hist_z_el_3_sim_1->Fill(dc_z_EL,1.);
if (sector ==4) hist_z_el_4_sim_1->Fill(dc_z_EL,1.);
if (sector ==5) hist_z_el_5_sim_1->Fill(dc_z_EL,1.);
if (sector ==6) hist_z_el_6_sim_1->Fill(dc_z_EL,1.);
*/


rot_boost_cmsyst();

if((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
Var_1[0] = inv_m_pip_p;
Var_1[1] = inv_m_pip_pim;
Var_1[2] = theta_P_cm;
Var_1[3] = phi_P_cm;
Var_1[4] = alpha_PIpPIm_pipf;

Var_2[0] = inv_m_pip_p;
Var_2[1] = inv_m_pip_pim;
Var_2[2] = theta_PIm_cm;
Var_2[3] = phi_PIm_cm;
Var_2[4] = alpha_PPIp_piPIm;

Var_3[0] = inv_m_pim_p;
Var_3[1] = inv_m_pip_pim;
Var_3[2] = theta_PIp_cm;
Var_3[3] = phi_PIp_cm;
Var_3[4] = alpha_PPIm_piPIp;

if ((indtype==1)||(indtype==4)){
//if ((selection_0_miss_sim)||(selection_pim_miss_sim)||(selection_pip_miss_sim)||(selection_proton_miss_sim)){
//if ((selection_0_miss_sim)||(selection_pim_miss_sim)||(selection_proton_miss_sim)){
//if ((selection_pim_miss_sim)){
//if ((selection_pim_miss_sim)){

if (selection_0_miss_sim) {
if (indtype==1) {
h_5dim_0_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,1.);
h_5dim_0_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,1.);
h_5dim_0_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,1.); 
};

if ((indtype==4)&&(W>1.65)&&(W<1.8)) {
h_5dim_0_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,10.);
h_5dim_0_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,10.);
h_5dim_0_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,10.); 
};
if ((indtype==4)&&(W>1.8)) {
h_5dim_0_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,7.2);
h_5dim_0_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,7.2);
h_5dim_0_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,7.2); 
};

if ((Q2>0.45)&&(Q2<0.5)) {

if ((indtype==1))h_0_mis_radcorr_sim[int((W-1.3)/0.025)]->Fill(P4_miss_0.Mag2(),1.);
if ((indtype==4))h_0_mis_radcorr_sim_bckg[int((W-1.3)/0.025)]->Fill(P4_miss_0.Mag2(),1.);

};
};
if (selection_pim_miss_sim) {
if (indtype==1) {
h_5dim_pim_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,1.);
h_5dim_pim_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,1.);
h_5dim_pim_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,1.); 
};
if ((indtype==4)&&(W>1.65)&&(W<1.8)) {
h_5dim_pim_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,10.);
h_5dim_pim_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,10.);
h_5dim_pim_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,10.); 
};
if ((indtype==4)&&(W>1.8)) {
h_5dim_pim_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,7.2);
h_5dim_pim_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,7.2);
h_5dim_pim_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,7.2); 
};

if ((Q2>0.45)&&(Q2<0.5)) {

if ((n_PIm == 0)&&(indtype==1))h_pim_mis_radcorr_sim[int((W-1.3)/0.025)]->Fill(P4_PIm_miss.Mag2(),1.);
if ((n_PIm == 0)&&(indtype==4))h_pim_mis_radcorr_sim_bckg[int((W-1.3)/0.025)]->Fill(P4_PIm_miss.Mag2(),1.);

};

};
if (selection_pip_miss_sim) {
if (indtype==1) {
h_5dim_pip_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,1.);
h_5dim_pip_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,1.);
h_5dim_pip_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,1.); 
};
if ((indtype==4)&&(W>1.65)&&(W<1.8)) {
h_5dim_pip_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,10.);
h_5dim_pip_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,10.);
h_5dim_pip_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,10.); 
};
if ((indtype==4)&&(W>1.8)) {
h_5dim_pip_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,7.2);
h_5dim_pip_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,7.2);
h_5dim_pip_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,7.2); 
};
if ((Q2>0.45)&&(Q2<0.5)) {

if ((n_PIp == 0)&&(indtype==1))h_pip_mis_radcorr_sim[int((W-1.3)/0.025)]->Fill(P4_PIp_miss.Mag2(),1.);
if ((n_PIp == 0)&&(indtype==4))h_pip_mis_radcorr_sim_bckg[int((W-1.3)/0.025)]->Fill(P4_PIp_miss.Mag2(),1.);

};


};
if (selection_proton_miss_sim) {
if (indtype==1) {
h_5dim_pr_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,1.);
h_5dim_pr_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,1.);
h_5dim_pr_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,1.); 
};
if ((indtype==4)&&(W>1.65)&&(W<1.8)) {
h_5dim_pr_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,10.5);
h_5dim_pr_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,10.5);
h_5dim_pr_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,10.5); 
};
if ((indtype==4)&&(W>1.8)) {
h_5dim_pr_1_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,7.2);
h_5dim_pr_2_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,7.2);
h_5dim_pr_3_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,7.2); 
};

if ((Q2>0.45)&&(Q2<0.5)) {

if ((n_P == 0)&&(indtype==1))h_pr_mis_radcorr_sim[int((W-1.3)/0.025)]->Fill(P4_P_miss.Mag2(),1.);
if ((n_P == 0)&&(indtype==4))h_pr_mis_radcorr_sim_bckg[int((W-1.3)/0.025)]->Fill(P4_P_miss.Mag2(),1.);

};

};



//};
};

};





}; //selection


 //}; //konec ifa electronnih cutov
    };// konets tsicla po sobitiyam (po i)
  
  
  t21->Delete();

 finp->Close();  
 
     };//konets tsicla po failam (po m)  
     
       
     
//ERRORS


    
  
     
     outFile = new TFile(outfile_inp.c_str(),"recreate");

h_cos_th = new TH1D("h_cos_th","h_cos_th",10,0.,180.);
     
Double_t temp;    
for (j=1; j<=10; j++) {
//temp = cos((h_cos_th->GetBinLowEdge(j))*M_PI/180.)-cos(M_PI/180.*(h_cos_th->GetBinLowEdge(j)+h_cos_th->GetBinWidth(j)));
//temp = cos(((h_cos_th->GetBinLowEdge(j))+h_cos_th->GetBinWidth(j)/2.)*M_PI/180.);
temp = sin(((h_cos_th->GetBinLowEdge(j))+h_cos_th->GetBinWidth(j)/2.)*M_PI/180.);
h_cos_th->SetBinContent(j,temp);
//h_cos_th->SetBinError(j,0.);
};     
     
    
       
//VIRTUAL PHOTON FLUX CALCULATION
Double_t W_bin[21];
Double_t Q2_bin[12];
Double_t omega[12][21],en_elp[12][21],th_elp[12][21],epsilon[12][21],flux[12][21],factor,L_0;
Long64_t Entr_sim_rec;
THnSparseD *ha[12][21];

for (j=0; j<12; j++) {
Q2_bin[j] = 0.425+0.05*j;
 for (i=0; i<21;i++) {
 W_bin[i] = 1.3125+0.025*i;
omega[j][i] = (W_bin[i]*W_bin[i] + Q2_bin[j] - m_proton*m_proton)/2./m_proton ;
en_elp[j][i] = 2.039 - omega[j][i];
th_elp[j][i]  = 2*asin(sqrt(Q2_bin[j]/4./2.039/en_elp[j][i]));

epsilon[j][i] = 1/(1. + 2.*(1. + omega[j][i]*omega[j][i]/Q2_bin[j])*(tan(th_elp[j][i]/2.))*(tan(th_elp[j][i]/2.)));
flux[j][i] = (omega[j][i]-Q2_bin[j]/2./m_proton)/137.;

flux[j][i]= flux[j][i]/2./(M_PI)/2.039/Q2_bin[j]/(1-epsilon[j][i]);
flux[j][i] = flux[j][i]*W_bin[i]/2.039/m_proton; 



//ERRORS





//cout << "j = " << j << " i = " << i << " flux = " << flux[j][i] << "\n";
factor = 1./flux[j][i];
L_0 = 0.63E12;
//factor = factor/L_0;





 };
  };


gStyle->SetTitleSize(0.08,"t");
gStyle->SetOptStat("e");
gStyle->SetStatY(0.88); 
gStyle->SetTitleY(0.99);
gStyle->SetTitleX(0.445);
 gStyle->SetErrorX(0);
//gStyle->SetError(0); 










 
 

 
 
 
 
 
 
 nphefile->Close(); 
// twopicutfile->Close(); 
 
  output();
     
 	outFile->Close();

	/*hist_sector1->Delete();
	hist_nphe_sector1->Delete();
	hist_sector2->Delete();
	hist_nphe_sector2->Delete();
	hist_sector3->Delete();
	hist_nphe_sector3->Delete();
	hist_sector4->Delete();
	hist_nphe_sector4->Delete();
	hist_sector5->Delete();
	hist_nphe_sector5->Delete();
	hist_sector6->Delete();
	hist_nphe_sector6->Delete();*/
	
	cout << "Qfull = " << Qfull << "\n";
	cout << "Qfull_empty = " << Qfull_empty << "\n";
	cout << "Qfull_sim = " << Qfull_sim << "\n";
	
    };

    
     

















