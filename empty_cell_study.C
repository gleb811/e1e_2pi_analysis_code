//TH1D *h_cos_th;
# define M_PI           3.14159265358979323846  /* pi */

THnSparseD *h_data_1[21],*h_data_2[21],*h_data_3[21];
THnSparseD *h_data_1_tmp[21],*h_data_2_tmp[21],*h_data_3_tmp[21];
THnSparseD *h_empty_1[21],*h_empty_2[21],*h_empty_3[21];
THnSparseD *h_empty_fill_1[21],*h_empty_fill_2[21],*h_empty_fill_3[21];
THnSparseD *h_rec_1[21],*h_rec_2[21],*h_rec_3[21];

THnSparseD *h_gen_1[21],*h_gen_2[21],*h_gen_3[21];
THnSparseD *h_eff_1[21],*h_eff_2[21],*h_eff_3[21];
THnSparseD *h_cr_sect_1[21],*h_cr_sect_2[21],*h_cr_sect_3[21];


THnSparseD *h_cr_sect_noemptcells_1[21],*h_cr_sect_noemptcells_2[21],*h_cr_sect_noemptcells_3[21];
THnSparseD *h_rec_unit_1[21],*h_rec_unit_2[21],*h_rec_unit_3[21];
THnSparseD *h_gen_unit_1[21],*h_gen_unit_2[21],*h_gen_unit_3[21];
THnSparseD *model_1[21],*model_2[21],*model_3[21];

TH2F *h_eff_err = new TH2F("h_eff_err","h_eff_err",100,0.001,0.2,100,0.,1.);

//Open output file
TFile *out_file = new TFile("eff_err.root","RECREATE");
/////// 

Float_t avrg_eff,avrg_eff_err;
Float_t avrg_data,avrg_data_err;
Float_t n_eff;

Float_t rad_corr[25];
Float_t empty_cell_persentage;
ostringstream qqq;
 Float_t W_bin[21];
 Float_t Q2_bin = 0.475;
Float_t m_proton = 0.938272;
 TLegend *leg;
 TLegend *leg_w_int;
 
void empty_cell_study() {
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
gStyle->SetTitleSize(0.07,"t");
gStyle->SetTitleY(1.01);
gStyle->SetOptStat(0);
gErrorIgnoreLevel = kError;
gStyle->SetStatY(0.88); 
TCanvas *c = new TCanvas("c","c",700,700);

//TH1F *h_eff_err = new TH1F("h_eff_err","h_eff_err",39,0.001-0.0005/2.,0.02+0.0005/2.);
//TH1F *h_eff_err = new TH1F("h_eff_err","h_eff_err",39,0.15-0.01/2.,0.53+0.01/2.);


h_eff_err->SetMarkerStyle(20);

TCanvas *c2 = new TCanvas("c2","c2",500,500);
TCanvas *c3 = new TCanvas("c3","c3",500,500);
c2->cd();
leg_w_int = new TLegend(0.11,0.7,0.89,0.89); 
leg_w_int->SetNColumns(3);
leg_w_int->SetFillStyle(0);

Int_t k,i,j;
Int_t number_entr;
Float_t eff_threshould = 0.00;
Float_t rec_threshould = 0.;
Float_t eff_err_threshould = 1000.15;
Float_t W_selected = 1.5125;
Float_t factor;




Double_t w_err;
Float_t Int[21];






//FILLING OUT -d(cos) HISTOGRAM
/*h_cos_th = new TH1D("h_cos_th","h_cos_th",10,0.,180.);  
Double_t temp;   
for (j=1; j<=10; j++) {
temp = cos((h_cos_th->GetBinLowEdge(j))*M_PI/180.)-cos(M_PI/180.*(h_cos_th->GetBinLowEdge(j)+h_cos_th->GetBinWidth(j)));
h_cos_th->SetBinContent(j,temp);
h_cos_th->SetBinError(j,0.);
};
*/

TH1D *h_w_int[16];

//Reading radiative correction factor

read_rad_corr("avrg_rad_corr.dat");


//Define input files

TFile *file_sim = new TFile("haddmy/combine_small_bin_all_top.root","READ");
TFile *file_data = new TFile("haddmy/combine_data_small_bin_all_top.root","READ");
TFile *file_empty = new TFile("haddmy/combine_data_small_bin_all_top.root","READ");
TFile *file_model = new TFile("empty_cells_new_bin.root","READ");

for (Int_t err_cyc=0; err_cyc<1;err_cyc++) {

eff_err_threshould = 0.15 + 0.01*err_cyc;
//eff_threshould = 0;
//eff_err_threshould = 0.15 + 0.005*err_cyc;
avrg_eff = 0.;
avrg_eff_err = 0.;
avrg_data = 0.;
avrg_data_err = 0.;
n_eff = 0.;

 for (Int_t qq2=2; qq2<3;qq2++) {
 Q2_bin = 0.425 + 0.05*qq2;

qqq.str("");
qqq << "h_w_int_" << qq2;
 h_w_int[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str("");

 for (i=13; i<14;i++) {
 W_bin[i] = 1.3125+0.025*i; 


//cout << i << "\n";

TH1D *h_prj_crs;

//Reading input files
read_sim(file_sim,i);
read_data(file_data,i);
read_empty(file_empty,i);
read_model(file_model,i);

set_eff_errors(i);
//Setting up error bars for model cross sections
model_errors(i);

//Applying cut on efficioncy error
eff_err_cut (eff_err_threshould, eff_threshould, rec_threshould, i);




set_all_errors(i);


//UNIT RECONSTRUCTED

qqq.str("");
qqq << "h_rec_unit_1_" << i;
h_rec_unit_1[i] = (THnSparseD*)h_rec_1[i]->Clone(qqq.str().c_str());
h_rec_unit_1[i]->Divide(h_rec_1[i]);

qqq.str("");
qqq << "h_rec_unit_2_" << i;
h_rec_unit_2[i] = (THnSparseD*)h_rec_2[i]->Clone(qqq.str().c_str());
h_rec_unit_2[i]->Divide(h_rec_2[i]);

qqq.str("");
qqq << "h_rec_unit_3_" << i;
h_rec_unit_3[i] = (THnSparseD*)h_rec_3[i]->Clone(qqq.str().c_str());
h_rec_unit_3[i]->Divide(h_rec_3[i]);

//UNIT GENERATED

qqq.str("");
qqq << "h_gen_unit_1_" << i;
h_gen_unit_1[i] = (THnSparseD*)h_gen_1[i]->Clone(qqq.str().c_str());
h_gen_unit_1[i]->Divide(h_gen_1[i]);

qqq.str("");
qqq << "h_gen_unit_2_" << i;
h_gen_unit_2[i] = (THnSparseD*)h_gen_2[i]->Clone(qqq.str().c_str());
h_gen_unit_2[i]->Divide(h_gen_2[i]);

qqq.str("");
qqq << "h_gen_unit_3_" << i;
h_gen_unit_3[i] = (THnSparseD*)h_gen_3[i]->Clone(qqq.str().c_str());
h_gen_unit_3[i]->Divide(h_gen_3[i]);

//1 IN BINS WHERE Ngen>1, Nrec=0 or efficiency errors > eff_err_threshould ("map" of empty cells)

h_gen_unit_1[i]->Add(h_rec_unit_1[i],-1);
h_gen_unit_2[i]->Add(h_rec_unit_2[i],-1);
h_gen_unit_3[i]->Add(h_rec_unit_3[i],-1);


//DIVIDE DATA AND EMPTY ON CHARGE ON FARADAY CUP

h_data_1_tmp[i] = (THnSparseD*)h_data_1[i]->Clone();
h_data_2_tmp[i] = (THnSparseD*)h_data_2[i]->Clone();
h_data_3_tmp[i] = (THnSparseD*)h_data_3[i]->Clone();


h_data_1[i]->Scale(1/5999.63);
h_data_2[i]->Scale(1/5999.63);
h_data_3[i]->Scale(1/5999.63);

h_empty_1[i]->Scale(1/334.602);
h_empty_2[i]->Scale(1/334.602);
h_empty_3[i]->Scale(1/334.602);


//SUBSTRUCT EMPTY TARGET


h_data_1[i]->Add(h_empty_1[i],-1.);
h_data_2[i]->Add(h_empty_2[i],-1.);
h_data_3[i]->Add(h_empty_3[i],-1.);

/*
h_data_1[i]->Scale(1. -((h_empty_1[i]->Projection(4,""))->Integral())/ ((h_data_1[i]->Projection(4,""))->Integral()));
h_data_2[i]->Scale(1. -((h_empty_2[i]->Projection(4,""))->Integral())/ ((h_data_2[i]->Projection(4,""))->Integral()));
h_data_3[i]->Scale(1. -((h_empty_3[i]->Projection(4,""))->Integral())/ ((h_data_3[i]->Projection(4,""))->Integral()));
*/

//cout << (h_data_1[i]->Projection(4,""))->Integral()  <<(h_empty_1[i]->Projection(4,""))->Integral() << "  " <<  1. -((h_empty_1[i]->Projection(4,""))->Integral())/ ((h_data_1[i]->Projection(4,""))->Integral()) <<  endl;


// Devide model crossection to the W and Q2 bin width

h_empty_fill_1[i]->Scale(0.025*0.05);
h_empty_fill_2[i]->Scale(0.025*0.05);
h_empty_fill_3[i]->Scale(0.025*0.05);


//Kepping histigrams with model cross section

qqq.str("");
qqq << "model_1_" << i;
model_1[i] = (THnSparseD*)h_empty_fill_1[i]->Clone(qqq.str().c_str());


qqq.str("");
qqq << "model_2_" << i;
model_2[i] = (THnSparseD*)h_empty_fill_2[i]->Clone(qqq.str().c_str());


qqq.str("");
qqq << "model_3_" << i;
model_3[i] = (THnSparseD*)h_empty_fill_3[i]->Clone(qqq.str().c_str());
qqq.str("");



//Kepp model cross section only in bins where  Ngen>1, Nrec=0 or efficiency error > eff_err_threshould


h_empty_fill_1[i]->Divide(h_gen_unit_1[i]);
h_empty_fill_2[i]->Divide(h_gen_unit_2[i]);
h_empty_fill_3[i]->Divide(h_gen_unit_3[i]);




//Deviding data to the efficiency

qqq.str("");
qqq << "h_cr_sect_1_" << i;
h_cr_sect_1[i] = (THnSparseD*)h_data_1[i]->Clone(qqq.str().c_str());
h_cr_sect_1[i]->Divide(h_eff_1[i]);

qqq.str("");
qqq << "h_cr_sect_2_" << i;
h_cr_sect_2[i] = (THnSparseD*)h_data_2[i]->Clone(qqq.str().c_str());
h_cr_sect_2[i]->Divide(h_eff_2[i]);

qqq.str("");
qqq << "h_cr_sect_3_" << i;
h_cr_sect_3[i] = (THnSparseD*)h_data_3[i]->Clone(qqq.str().c_str());
h_cr_sect_3[i]->Divide(h_eff_3[i]);


/*
h_empty_1[i]->Divide(h_eff_1[i]);
cout << (h_cr_sect_1[i]->Projection(4,""))->Integral()  <<(h_empty_1[i]->Projection(4,""))->Integral() << "  " <<  1. -((h_empty_1[i]->Projection(4,""))->Integral())/ ((h_cr_sect_1[i]->Projection(4,""))->Integral()) <<  endl;
*/

//Cross sections with no empty cells filling

qqq.str("");
qqq << "h_cr_sect_noemptcells_1_" << i;
h_cr_sect_noemptcells_1[i] = (THnSparseD*)h_cr_sect_1[i]->Clone(qqq.str().c_str());

qqq.str("");
qqq << "h_cr_sect_noemptcells_2_" << i;
h_cr_sect_noemptcells_2[i] = (THnSparseD*)h_cr_sect_2[i]->Clone(qqq.str().c_str());

qqq.str("");
qqq << "h_cr_sect_noemptcells_3_" << i;
h_cr_sect_noemptcells_3[i] = (THnSparseD*)h_cr_sect_3[i]->Clone(qqq.str().c_str());
qqq.str("");


//Scaling with virtual photon flux and luminosity

h_cr_sect_1[i]->Scale(1/flux(i));
h_cr_sect_2[i]->Scale(1/flux(i));
h_cr_sect_3[i]->Scale(1/flux(i));

h_cr_sect_1[i]->Scale(1.893785248e-6);
h_cr_sect_2[i]->Scale(1.893785248e-6);
h_cr_sect_3[i]->Scale(1.893785248e-6);

h_cr_sect_noemptcells_1[i]->Scale(1/flux(i));
h_cr_sect_noemptcells_2[i]->Scale(1/flux(i));
h_cr_sect_noemptcells_3[i]->Scale(1/flux(i));

h_cr_sect_noemptcells_1[i]->Scale(1.893785248e-6);
h_cr_sect_noemptcells_2[i]->Scale(1.893785248e-6);
h_cr_sect_noemptcells_3[i]->Scale(1.893785248e-6);

//Applying RAD CORR factor
h_cr_sect_1[i]->Scale(rad_corr[i+1]);
h_cr_sect_2[i]->Scale(rad_corr[i+1]);
h_cr_sect_3[i]->Scale(rad_corr[i+1]);

h_cr_sect_noemptcells_1[i]->Scale(rad_corr[i+1]);
h_cr_sect_noemptcells_2[i]->Scale(rad_corr[i+1]);
h_cr_sect_noemptcells_3[i]->Scale(rad_corr[i+1]);

//Filling out empty cells

h_cr_sect_1[i]->Add(h_empty_fill_1[i],1.);
h_cr_sect_2[i]->Add(h_empty_fill_2[i],1.);
h_cr_sect_3[i]->Add(h_empty_fill_3[i],1.);


//cr_sec_errors(i);


//h_prj_crs = h_cr_sect_2[i]->Projection(4,"");
h_prj_crs = empty_fill_err(h_cr_sect_2[i]->Projection(4,""),h_cr_sect_noemptcells_2[i]->Projection(4,""));
h_prj_crs->Scale(1/(0.025*0.05));


 Int[i]= h_prj_crs->Integral();
 h_prj_crs->IntegralAndError(1,h_prj_crs->GetNbinsX(),w_err);
 h_w_int[qq2]->Fill(W_bin[i],Int[i]); 
h_w_int[qq2]->SetBinError(h_w_int[qq2]->FindBin(W_bin[i]),w_err);

//empty_cell_persentage = (((h_cr_sect_2[i]->Projection(4,""))->Integral())-((h_cr_sect_noemptcells_2[i]->Projection(4,""))->Integral()))/((h_cr_sect_2[i]->Projection(4,""))->Integral())*100.;

//cout << "eff cut " << eff_threshould << " persentage of empty cells " << empty_cell_persentage << endl;
h_data_2_tmp[i]->Divide(h_gen_unit_2[i]);
//h_eff_err->Fill(eff_threshould,empty_cell_persentage);
//h_eff_err->Fill(eff_threshould,((h_gen_unit_2[i]->Projection(4,""))->Integral()));
//h_eff_err->Fill(eff_threshould,((h_data_2_tmp[i]->Projection(4,""))->Integral())*avrg_eff_err);
//h_eff_err->Fill(eff_err_threshould,avrg_eff);
//cout << eff_threshould << "  " << avrg_data_err << endl;
//h_eff_err->Fill(eff_threshould,avrg_data_err);
//(h_data_2_tmp[i]->Projection(4,""))->IntegralAndError(1,(h_data_2_tmp[i]->Projection(4,""))->GetNbinsX(),w_err);
//h_eff_err->Fill(eff_threshould,w_err);
//if (W_bin[i] == W_selected){

draw_1d_canvas(i,qq2);

//};




}; //end of W loop


c2->cd();

if (qq2<6) {
h_w_int[qq2]->SetMarkerStyle(20);
h_w_int[qq2]->SetMarkerColor(qq2+1);
};
if (qq2>=6) {
h_w_int[qq2]->SetMarkerStyle(29);
h_w_int[qq2]->SetMarkerColor(qq2+1-6);
};
if (qq2>=12) {
h_w_int[qq2]->SetMarkerStyle(29);
h_w_int[qq2]->SetMarkerColor(qq2+1-12);
};
 h_w_int[qq2]->GetXaxis()->SetTitle("W, GeV");
 h_w_int[qq2]->GetXaxis()->SetNdivisions(8);
 h_w_int[qq2]->GetXaxis()->SetLabelSize(0.04);
 h_w_int[qq2]->GetYaxis()->SetLabelSize(0.04);
 Double_t max = h_w_int[qq2]->GetMaximum(); 
 h_w_int[qq2]->SetAxisRange(0.,1.5*max,"Y"); 
 h_w_int[qq2]->SetAxisRange(1.3,1.8,"X");
  h_w_int[qq2]->SetTitle("Intergated cross section");
 h_w_int[qq2] ->GetYaxis()->SetTitle("#sigma, #mubn");
 if(qq2==0) {
 h_w_int[qq2]->Draw("e1pX0");
 } else {
 h_w_int[qq2]->Draw("e1pX0 same"); 
 };
 qqq.str("");
 qqq << "Q^{2} = " << Q2_bin << " GeV^{2}";
 leg_w_int->AddEntry(h_w_int[qq2],qqq.str().c_str(),"p");
 qqq.str("");
 c2->Update();
// ostringstream qqq3; 
/* for  (j=1; j<12; j++) { 
 h_w_int[j]->SetTitle("  ");
 h_w_int[j]->SetMarkerStyle(20);
  h_w_int[j]->SetMarkerColor(j+1); 
 if (j == 9) h_w_int[j]->SetMarkerColor(46); 
 if (j == 8) h_w_int[j]->SetMarkerColor(41);
 h_w_int[j]->Draw("same e1pX0"); 
 };
h_w_int[2]->SetAxisRange(1.3,1.79,"X"); 
 h_w_int[3]->SetAxisRange(1.3,1.79,"X"); 
 h_w_int[4]->SetAxisRange(1.3,1.774,"X"); 
 h_w_int[5]->SetAxisRange(1.3,1.74,"X"); 
 h_w_int[6]->SetAxisRange(1.3,1.724,"X"); 
 h_w_int[7]->SetAxisRange(1.3,1.69,"X"); 
 h_w_int[8]->SetAxisRange(1.3,1.64,"X"); 
 h_w_int[21]->SetAxisRange(1.3,1.624,"X"); 
 h_w_int[10]->SetAxisRange(1.3,1.574,"X"); 
 h_w_int[11]->SetAxisRange(1.3,1.54,"X"); 
TLegend*leg = new TLegend(0.15,0.6,0.45,0.9); 
 for  (j=0; j<12; j++) { 
 qqq3 <<"Q^{2} = " << 0.475+0.05*j << "GeV^{2}" ; 
 leg->AddEntry(h_w_int[j],qqq3.str().c_str(),"p"); qqq3.str(""); 
 }; 
 leg->Draw();*/
/*
for (k=1; k<=h_gen_cor->GetNbins(); k++) {
if (h_gen_cor->GetBinContent(k)>50) cout << h_gen_cor->GetBinContent(k) << "\n";
};
*/
//gDirectory->ls("theta*");
//MyFile->Close();

}; //end of Q2 loop

c3->cd();
//h_eff_err->Draw("P");
h_eff_err->Draw("colz");
c3->Update();

}; //end err_cyc loop

c2->cd();
leg_w_int->Draw();
c2->SaveAs("w_int.eps");
c2->SaveAs("w_int.png");
out_file->cd();
h_eff_err->Write();
out_file->Close();
}; //end of main program

void model_errors (Int_t i) {
Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t y_max = 8;
Int_t t_max = 8;

if ((i==0)||(i==1)) {
o_max = p_max = 8;
r_max = 6;
y_max = 5; 
t_max = 5;
};
if ((i==2)||(i==3)) {
o_max = p_max = 10;
r_max = 8;
y_max = 6;
t_max = 5;
}; 
if ((i>=4)&&(i<=6)){
o_max = p_max =12;
r_max =10;
t_max = 5;
y_max = 8;
};
 
for (Int_t o=1; o<=o_max; o++) {
for (Int_t p=1; p<=p_max; p++) {
for (Int_t r=1; r<=r_max; r++) {
for (Int_t t=1; t<=t_max; t++) {
for (Int_t y=1; y<=y_max; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;


Long64_t tmp_model1_bin = h_empty_fill_1[i]->GetBin(bins);
Long64_t tmp_model2_bin = h_empty_fill_2[i]->GetBin(bins);
Long64_t tmp_model3_bin = h_empty_fill_3[i]->GetBin(bins);



h_empty_fill_1[i]->SetBinError(tmp_model1_bin,(h_empty_fill_1[i]->GetBinContent(tmp_model1_bin)));
h_empty_fill_2[i]->SetBinError(tmp_model2_bin,(h_empty_fill_2[i]->GetBinContent(tmp_model2_bin)));
h_empty_fill_3[i]->SetBinError(tmp_model3_bin,(h_empty_fill_3[i]->GetBinContent(tmp_model3_bin)));


};
};
};
}; 
};

};

void set_all_errors  (Int_t i) {
Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t y_max = 8;
Int_t t_max = 8;

if ((i==0)||(i==1)) {
o_max = p_max = 8;
r_max = 6;
y_max = 5; 
t_max = 5;
};
if ((i==2)||(i==3)) {
o_max = p_max = 10;
r_max = 8;
y_max = 6;
t_max = 5;
}; 
if ((i>=4)&&(i<=6)){
o_max = p_max =12;
r_max =10;
t_max = 5;
y_max = 8;
};



 
for (Int_t o=1; o<=o_max; o++) {
for (Int_t p=1; p<=p_max; p++) {
for (Int_t r=1; r<=r_max; r++) {
for (Int_t t=1; t<=t_max; t++) {
for (Int_t y=1; y<=y_max; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;





if (h_data_1[i]->GetBinContent(bins) > 0.) {
h_data_1[i]->SetBinError(bins,sqrt(h_data_1[i]->GetBinContent(bins)));
} else {
h_data_1[i]->SetBinError(bins,0.);
};

if (h_data_2[i]->GetBinContent(bins) > 0.) {
h_data_2[i]->SetBinError(bins,sqrt(h_data_2[i]->GetBinContent(bins)));
} else {
h_data_2[i]->SetBinError(bins,0.);
};
if (h_data_3[i]->GetBinContent(bins) > 0.) {
h_data_3[i]->SetBinError(bins,sqrt(h_data_3[i]->GetBinContent(bins)));
} else {
h_data_3[i]->SetBinError(bins,0.);
};


h_empty_1[i]->SetBinError(bins,sqrt(h_empty_1[i]->GetBinContent(bins)));
h_empty_2[i]->SetBinError(bins,sqrt(h_empty_2[i]->GetBinContent(bins)));
h_empty_3[i]->SetBinError(bins,sqrt(h_empty_3[i]->GetBinContent(bins)));

/*
h_rec_1[i]->SetBinError(bins,0.);
h_rec_2[i]->SetBinError(bins,0.);
h_rec_3[i]->SetBinError(bins,0.);

h_gen_1[i]->SetBinError(bins,0.);
h_gen_2[i]->SetBinError(bins,0.);
h_gen_3[i]->SetBinError(bins,0.);


h_eff_1[i]->SetBinError(bins,0.);
h_eff_2[i]->SetBinError(bins,0.);
h_eff_3[i]->SetBinError(bins,0.);
*/
h_empty_fill_1[i]->SetBinError(bins,0.);
h_empty_fill_2[i]->SetBinError(bins,0.);
h_empty_fill_3[i]->SetBinError(bins,0.);


};
};
};
}; 
};



};


void set_eff_errors (Int_t i) {
Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t y_max = 8;
Int_t t_max = 8;

if ((i==0)||(i==1)) {
o_max = p_max = 8;
r_max = 6;
y_max = 5; 
t_max = 5;
};
if ((i==2)||(i==3)) {
o_max = p_max = 10;
r_max = 8;
y_max = 6;
t_max = 5;
}; 
if ((i>=4)&&(i<=6)){
o_max = p_max =12;
r_max =10;
t_max = 5;
y_max = 8;
};
 
for (Int_t o=1; o<=o_max; o++) {
for (Int_t p=1; p<=p_max; p++) {
for (Int_t r=1; r<=r_max; r++) {
for (Int_t t=1; t<=t_max; t++) {
for (Int_t y=1; y<=y_max; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;


Long64_t tmp_rec1_bin = h_rec_1[i]->GetBin(bins);
Long64_t tmp_rec2_bin = h_rec_2[i]->GetBin(bins);
Long64_t tmp_rec3_bin = h_rec_3[i]->GetBin(bins);

Long64_t tmp_gen1_bin = h_gen_1[i]->GetBin(bins);
Long64_t tmp_gen2_bin = h_gen_2[i]->GetBin(bins);
Long64_t tmp_gen3_bin = h_gen_3[i]->GetBin(bins);

Long64_t tmp_eff1_bin = h_eff_1[i]->GetBin(bins);
Long64_t tmp_eff2_bin = h_eff_2[i]->GetBin(bins);
Long64_t tmp_eff3_bin = h_eff_3[i]->GetBin(bins);

Double_t err1,err2;

if ((h_eff_1[i]->GetBinContent(tmp_eff1_bin) > 0.) && ((h_gen_1[i]->GetBinContent(tmp_gen1_bin)) > (h_rec_1[i]->GetBinContent(tmp_rec1_bin)))) {

err1 = (h_rec_1[i]->GetBinContent(tmp_rec1_bin))*((h_gen_1[i]->GetBinContent(tmp_gen1_bin))-(h_rec_1[i]->GetBinContent(tmp_rec1_bin)));
err1 = err1/(h_gen_1[i]->GetBinContent(tmp_gen1_bin))/(h_gen_1[i]->GetBinContent(tmp_gen1_bin))/(h_gen_1[i]->GetBinContent(tmp_gen1_bin));
err1 = sqrt(err1);
h_eff_1[i]->SetBinError(tmp_eff1_bin,err1);

};

if ((h_eff_2[i]->GetBinContent(tmp_eff2_bin) > 0.) && ((h_gen_2[i]->GetBinContent(tmp_gen2_bin)) > (h_rec_2[i]->GetBinContent(tmp_rec2_bin)))) {

err1 = (h_rec_2[i]->GetBinContent(tmp_rec2_bin))*((h_gen_2[i]->GetBinContent(tmp_gen2_bin))-(h_rec_2[i]->GetBinContent(tmp_rec2_bin)));
err1 = err1/(h_gen_2[i]->GetBinContent(tmp_gen2_bin))/(h_gen_2[i]->GetBinContent(tmp_gen2_bin))/(h_gen_2[i]->GetBinContent(tmp_gen2_bin));
err1 = sqrt(err1);
h_eff_2[i]->SetBinError(tmp_eff2_bin,err1);

};

if ((h_eff_3[i]->GetBinContent(tmp_eff3_bin) > 0.) && ((h_gen_3[i]->GetBinContent(tmp_gen3_bin)) > (h_rec_3[i]->GetBinContent(tmp_rec3_bin)))) {

err1 = (h_rec_3[i]->GetBinContent(tmp_rec3_bin))*((h_gen_3[i]->GetBinContent(tmp_gen3_bin))-(h_rec_3[i]->GetBinContent(tmp_rec3_bin)));
err1 = err1/(h_gen_3[i]->GetBinContent(tmp_gen3_bin))/(h_gen_3[i]->GetBinContent(tmp_gen3_bin))/(h_gen_3[i]->GetBinContent(tmp_gen3_bin));
err1 = sqrt(err1);
err1 = err1/(h_gen_3[i]->GetBinContent(tmp_gen3_bin));

h_eff_3[i]->SetBinError(tmp_eff3_bin,err1);

};





};
};
};
}; 
};



};



void cr_sec_errors (Int_t i) {
Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t y_max = 8;
Int_t t_max = 8;

if ((i==0)||(i==1)) {
o_max = p_max = 8;
r_max = 6;
y_max = 5; 
t_max = 5;
};
if ((i==2)||(i==3)) {
o_max = p_max = 10;
r_max = 8;
y_max = 6;
t_max = 5;
}; 
if ((i>=4)&&(i<=6)){
o_max = p_max =12;
r_max =10;
t_max = 5;
y_max = 8;
};
 
for (Int_t o=1; o<=o_max; o++) {
for (Int_t p=1; p<=p_max; p++) {
for (Int_t r=1; r<=r_max; r++) {
for (Int_t t=1; t<=t_max; t++) {
for (Int_t y=1; y<=y_max; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;


Long64_t tmp_data1_bin = h_data_1_tmp[i]->GetBin(bins);
Long64_t tmp_data2_bin = h_data_2_tmp[i]->GetBin(bins);
Long64_t tmp_data3_bin = h_data_3_tmp[i]->GetBin(bins);

Long64_t tmp_cr_sect1_bin = h_cr_sect_1[i]->GetBin(bins);
Long64_t tmp_cr_sect2_bin = h_cr_sect_2[i]->GetBin(bins);
Long64_t tmp_cr_sect3_bin = h_cr_sect_3[i]->GetBin(bins);

Long64_t tmp_eff1_bin = h_eff_1[i]->GetBin(bins);
Long64_t tmp_eff2_bin = h_eff_2[i]->GetBin(bins);
Long64_t tmp_eff3_bin = h_eff_3[i]->GetBin(bins);

Double_t err1,err2;

if (((h_data_1[i]->GetBinContent(tmp_data1_bin)) > 0.) &&(h_eff_1[i]->GetBinContent(tmp_eff1_bin) > 0.)) {

err1 = 1./sqrt(h_data_1_tmp[i]->GetBinContent(tmp_data1_bin));
err2 = (h_eff_1[i]->GetBinError(tmp_eff1_bin))/(h_eff_1[i]->GetBinContent(tmp_eff1_bin));

h_cr_sect_1[i]->SetBinError(tmp_cr_sect1_bin,(h_cr_sect_1[i]->GetBinContent(tmp_cr_sect1_bin))*sqrt(err1*err1+err2*err2));

} else {

h_cr_sect_1[i]->SetBinError(tmp_cr_sect1_bin,0.);
};


if (((h_data_2[i]->GetBinContent(tmp_data2_bin)) > 0.) &&(h_eff_2[i]->GetBinContent(tmp_eff2_bin) > 0.)) {

err1 = 1./sqrt(h_data_2_tmp[i]->GetBinContent(tmp_data2_bin));
err2 = (h_eff_2[i]->GetBinError(tmp_eff2_bin))/(h_eff_2[i]->GetBinContent(tmp_eff2_bin));

h_cr_sect_2[i]->SetBinError(tmp_cr_sect2_bin,(h_cr_sect_2[i]->GetBinContent(tmp_cr_sect2_bin))*sqrt(err1*err1+err2*err2));

} else {

h_cr_sect_2[i]->SetBinError(tmp_cr_sect2_bin,0.);
};


if (((h_data_3[i]->GetBinContent(tmp_data3_bin)) > 0.) &&(h_eff_3[i]->GetBinContent(tmp_eff3_bin) > 0.)) {

err1 = 1./sqrt(h_data_3_tmp[i]->GetBinContent(tmp_data3_bin));
err2 = (h_eff_3[i]->GetBinError(tmp_eff3_bin))/(h_eff_3[i]->GetBinContent(tmp_eff3_bin));

h_cr_sect_3[i]->SetBinError(tmp_cr_sect3_bin,(h_cr_sect_3[i]->GetBinContent(tmp_cr_sect3_bin))*sqrt(err1*err1+err2*err2));

} else {

h_cr_sect_3[i]->SetBinError(tmp_cr_sect3_bin,0.);
};


};
};
};
}; 
};



};

void draw_1d_hist (Int_t canvas, TH1D *h, string title, string name, string ytitle, string xtitle, Int_t color, string draw_options, string distr_flag,Int_t i) {

TH1D *h_cos_th;
//FILLING OUT -d(cos) HISTOGRAM
Int_t n_theta_bins;
if ((i==0)||(i==1)) n_theta_bins = 6;
if ((i==2)||(i==3)) n_theta_bins = 8;
if (i>=4) n_theta_bins = 10; 
h_cos_th = new TH1D("h_cos_th","h_cos_th",n_theta_bins,0.,180.);  
Double_t temp;   
for (Int_t j=1; j<=n_theta_bins; j++) {
temp = cos((h_cos_th->GetBinLowEdge(j))*M_PI/180.)-cos(M_PI/180.*(h_cos_th->GetBinLowEdge(j)+h_cos_th->GetBinWidth(j)));
h_cos_th->SetBinContent(j,temp);
h_cos_th->SetBinError(j,0.);
};

c->cd(canvas);
c->cd(canvas)->SetBottomMargin(0.2);
c->cd(canvas)->SetTopMargin(0.1);
c->cd(canvas)->SetLeftMargin(0.17);
gPad->SetFillStyle(0);

h->Sumw2();
h->SetMarkerStyle(20);
h->SetMarkerColor(color);
h->SetLineColor(color);
h->SetOption("pX0");
h->SetTitle(title.c_str());
h->SetTitleSize(0.1);


h->SetName(name.c_str());


 h->GetYaxis()->SetTitle(ytitle.c_str());
 h ->GetXaxis()->SetTitle(xtitle.c_str());
 h ->GetXaxis()->SetTitleSize(0.08);
 h->GetYaxis()->SetTitleSize(0.08);
 h->GetXaxis()->SetLabelSize(0.08);
 h->GetXaxis()->SetNdivisions(6);
 h->GetYaxis()->SetLabelSize(0.07);
 h->GetYaxis()->SetNdivisions(5);
 if (distr_flag == "mass") Float_t factor = 0.025*0.05*(h->GetBinWidth(5));
 if (distr_flag == "theta") {
 h->Divide(h_cos_th);
 Float_t factor = 0.025*0.05;
 }; 
  if (distr_flag == "alpha") Float_t factor = M_PI*0.025*0.05*(h->GetBinWidth(5))/180;
h->Scale(1./factor);
h->SetAxisRange(0, h-> GetMaximum()+(h-> GetMaximum()), "Y");
if ((name == "model_thetapip_"))h->SetAxisRange(0, 1.2*(h-> GetMaximum()), "Y");
if (name == "h1prj_inv_m_pip_p_1_") leg->AddEntry(h,"empty cells from model","p");
if (name == "h1prj_inv_m_pip_p_1_noempty") leg->AddEntry(h,"zero in empty cells","p");
if (name == "model_pipP_1_") leg->AddEntry(h,"model","p");
h->Draw(draw_options.c_str());
};

void eff_err_cut (Float_t  eff_err_threshould, Float_t  eff_threshould, Float_t  rec_threshould, Int_t i) {

Int_t *bins = new Int_t[5];
Int_t o_max =12;
Int_t p_max =12;
Int_t r_max = 10;
Int_t y_max = 8;
Int_t t_max = 8;

if ((i==0)||(i==1)) {
o_max = p_max = 8;
r_max = 6;
y_max = 5; 
t_max = 5;
};
if ((i==2)||(i==3)) {
o_max = p_max = 10;
r_max = 8;
y_max = 6;
t_max = 5;
}; 
if ((i>=4)&&(i<=6)){
o_max = p_max =12;
r_max =10;
t_max = 5;
y_max = 8;
};
 
for (Int_t o=1; o<=o_max; o++) {
for (Int_t p=1; p<=p_max; p++) {
for (Int_t r=1; r<=r_max; r++) {
for (Int_t t=1; t<=t_max; t++) {
for (Int_t y=1; y<=y_max; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;


Long64_t tmp_eff1_bin = h_eff_1[i]->GetBin(bins);
Long64_t tmp_rec1_bin = h_rec_1[i]->GetBin(bins);
Long64_t tmp_gen1_bin = h_gen_1[i]->GetBin(bins);

Long64_t tmp_eff2_bin = h_eff_2[i]->GetBin(bins);
Long64_t tmp_rec2_bin = h_rec_2[i]->GetBin(bins);
Long64_t tmp_gen2_bin = h_gen_2[i]->GetBin(bins);

Long64_t tmp_eff3_bin = h_eff_3[i]->GetBin(bins);
Long64_t tmp_rec3_bin = h_rec_3[i]->GetBin(bins);
Long64_t tmp_gen3_bin = h_gen_3[i]->GetBin(bins);


if ((((h_eff_1[i]->GetBinError(tmp_eff1_bin))/(h_eff_1[i]->GetBinContent(tmp_eff1_bin))) > eff_err_threshould) || ((h_gen_1[i]->GetBinContent(tmp_gen1_bin)) <= (h_rec_1[i]->GetBinContent(tmp_rec1_bin)))) {
h_eff_1[i]->SetBinContent(tmp_eff1_bin,0);
h_eff_1[i]->SetBinError(tmp_eff1_bin,0);
h_rec_1[i]->SetBinContent(tmp_rec1_bin,0);
h_rec_1[i]->SetBinError(tmp_rec1_bin,0);
};

if (h_eff_1[i]->GetBinContent(tmp_eff1_bin) < eff_threshould) {
h_eff_1[i]->SetBinContent(tmp_eff1_bin,0);
h_eff_1[i]->SetBinError(tmp_eff1_bin,0);
h_rec_1[i]->SetBinContent(tmp_rec1_bin,0);
h_rec_1[i]->SetBinError(tmp_rec1_bin,0);
};
if (h_rec_1[i]->GetBinContent(tmp_rec1_bin) < rec_threshould) {
h_eff_1[i]->SetBinContent(tmp_eff1_bin,0);
h_eff_1[i]->SetBinError(tmp_eff1_bin,0);
h_rec_1[i]->SetBinContent(tmp_rec1_bin,0);
h_rec_1[i]->SetBinError(tmp_rec1_bin,0);
};

if (h_eff_2[i]->GetBinContent(tmp_eff2_bin) > 0.) {
//if (h_data_2[i]->GetBinContent(bins) > 0.) {
h_eff_err->Fill(h_eff_2[i]->GetBinContent(tmp_eff2_bin),((h_eff_2[i]->GetBinError(tmp_eff2_bin))/(h_eff_2[i]->GetBinContent(tmp_eff2_bin))));
//};
};


if ((((h_eff_2[i]->GetBinError(tmp_eff2_bin))/(h_eff_2[i]->GetBinContent(tmp_eff2_bin))) > eff_err_threshould) || ((h_gen_2[i]->GetBinContent(tmp_gen2_bin)) <= (h_rec_2[i]->GetBinContent(tmp_rec2_bin)))) {
if (((h_eff_2[i]->GetBinError(tmp_eff2_bin))/(h_eff_2[i]->GetBinContent(tmp_eff2_bin))) < eff_err_threshould + 0.01) {
avrg_eff = avrg_eff + h_eff_2[i]->GetBinContent(tmp_eff2_bin);
n_eff = n_eff + 1.;
};


h_eff_2[i]->SetBinContent(tmp_eff2_bin,0);
h_eff_2[i]->SetBinError(tmp_eff2_bin,0);
h_rec_2[i]->SetBinContent(tmp_rec2_bin,0);
h_rec_2[i]->SetBinError(tmp_rec2_bin,0);
};



if (h_eff_2[i]->GetBinContent(tmp_eff2_bin) < eff_threshould) {
/*
if (h_eff_2[i]->GetBinContent(tmp_eff2_bin) > 0.) {
if (h_eff_2[i]->GetBinContent(tmp_eff2_bin) > eff_threshould - 0.0005) {
avrg_eff = avrg_eff + h_eff_2[i]->GetBinContent(tmp_eff2_bin);
avrg_eff_err = avrg_eff_err + (h_eff_2[i]->GetBinError(tmp_eff2_bin));
};
};
if (h_data_2[i]->GetBinContent(bins) > 0.) {
avrg_data = avrg_data + h_data_2[i]->GetBinContent(bins);
avrg_data_err = avrg_data_err + (h_data_2[i]->GetBinError(bins));
};
*/
h_eff_2[i]->SetBinContent(tmp_eff2_bin,0);
h_eff_2[i]->SetBinError(tmp_eff2_bin,0);
h_rec_2[i]->SetBinContent(tmp_rec2_bin,0);
h_rec_2[i]->SetBinError(tmp_rec2_bin,0);
//n_eff = n_eff +1;
};
if (h_rec_2[i]->GetBinContent(tmp_rec2_bin) < rec_threshould) {
h_eff_2[i]->SetBinContent(tmp_eff2_bin,0);
h_eff_2[i]->SetBinError(tmp_eff2_bin,0);
h_rec_2[i]->SetBinContent(tmp_rec2_bin,0);
h_rec_2[i]->SetBinError(tmp_rec2_bin,0);
};


if ((((h_eff_3[i]->GetBinError(tmp_eff3_bin))/(h_eff_3[i]->GetBinContent(tmp_eff3_bin))) > eff_err_threshould) || ((h_gen_3[i]->GetBinContent(tmp_gen3_bin)) <= (h_rec_3[i]->GetBinContent(tmp_rec3_bin)))) {
h_eff_3[i]->SetBinContent(tmp_eff3_bin,0);
h_eff_3[i]->SetBinError(tmp_eff3_bin,0);
h_rec_3[i]->SetBinContent(tmp_rec3_bin,0);
h_rec_3[i]->SetBinError(tmp_rec3_bin,0);
};

if (h_eff_3[i]->GetBinContent(tmp_eff3_bin) < eff_threshould) {
h_eff_3[i]->SetBinContent(tmp_eff3_bin,0);
h_eff_3[i]->SetBinError(tmp_eff3_bin,0);
h_rec_3[i]->SetBinContent(tmp_rec3_bin,0);
h_rec_3[i]->SetBinError(tmp_rec3_bin,0);
};
if (h_rec_3[i]->GetBinContent(tmp_rec3_bin) < rec_threshould) {
h_eff_3[i]->SetBinContent(tmp_eff3_bin,0);
h_eff_3[i]->SetBinError(tmp_eff3_bin,0);
h_rec_3[i]->SetBinContent(tmp_rec3_bin,0);
h_rec_3[i]->SetBinError(tmp_rec3_bin,0);
};

};
};
};
}; 
};

//avrg_eff_err = avrg_eff_err/avrg_eff;
//avrg_data_err = avrg_data_err/avrg_data;
avrg_eff = avrg_eff/n_eff;
};

Int_t get_max_w (Float_t Q2_bin) {
Int_t get_max_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.5))get_max_w = 21;
if ((Q2_bin> 0.5)&&(Q2_bin< 0.6))get_max_w = 20;
if ((Q2_bin> 0.6)&&(Q2_bin< 0.65))get_max_w = 19;
if ((Q2_bin> 0.65)&&(Q2_bin< 0.7))get_max_w = 18;
if ((Q2_bin> 0.7)&&(Q2_bin< 0.75))get_max_w = 17;
if ((Q2_bin> 0.75)&&(Q2_bin< 0.8))get_max_w = 16;
if ((Q2_bin> 0.8)&&(Q2_bin< 0.85))get_max_w = 14;
if ((Q2_bin> 0.85)&&(Q2_bin< 0.9))get_max_w = 13;
if ((Q2_bin> 0.9)&&(Q2_bin< 0.95))get_max_w = 12;
if ((Q2_bin> 0.95)&&(Q2_bin< 1.0))get_max_w = 10;
if ((Q2_bin> 1.0)&&(Q2_bin< 1.05))get_max_w = 9;
if ((Q2_bin> 1.05)&&(Q2_bin< 1.1))get_max_w = 7;
if ((Q2_bin> 1.1)&&(Q2_bin< 1.15))get_max_w = 6;
if ((Q2_bin> 1.15)&&(Q2_bin< 1.2))get_max_w = 4;
if ((Q2_bin> 1.2)&&(Q2_bin< 1.25))get_max_w = 2;
if ((Q2_bin> 1.25)&&(Q2_bin< 1.3))get_max_w = 1;
return get_max_w;
};

Int_t get_min_w (Float_t Q2_bin) {
Int_t get_min_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 8;
return get_min_w;
};


void read_rad_corr (string file_name) {

ifstream input(file_name.c_str());
if(input.is_open()){
Int_t i=0;
    while(!input.eof()){
          string line1,t_str,e_str;
           getline(input,line1); //read number
	   if (line1.length() != 0){ 
              t_str= line1.substr(0,line1.find("   "));
	    e_str = line1.substr(t_str.length()+1,t_str.length()+7);  
            rad_corr[i] = atof(e_str.c_str());
		   
	   
	    	    	    
//	    cout << i << "   " << rad_corr[i] <<" \n";
	  
	    i=i+1;
	    	    };
	    
    };
};


};


void read_sim(TFile *MyFile,Int_t i) {

MyFile->cd();
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_rec_3[i]);
qqq.str("");


qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_gen_3[i]);
qqq.str("");






qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_1_" << i;
gDirectory->GetObject(qqq.str().c_str(),h_eff_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_2_" << i;
gDirectory->GetObject(qqq.str().c_str(),h_eff_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_3_" << i;
gDirectory->GetObject(qqq.str().c_str(),h_eff_3[i]);
qqq.str("");

};


void read_data(TFile *MyFile,Int_t i) {

MyFile->cd();
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_data_3[i]);
qqq.str("");



};

void read_empty(TFile *MyFile,Int_t i) {

MyFile->cd();
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empty_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empty_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empty_3[i]);
qqq.str("");

};


void read_model(TFile *MyFile,Int_t i) {

MyFile->cd();

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_empty1_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empty_fill_1[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_empty2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empty_fill_2[i]);
qqq.str("");

qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_empty3_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),h_empty_fill_3[i]);
qqq.str("");

};

Float_t flux(Int_t i) {

 //VIRTUAL PHOTON FLUX CALCULATION
 
 
Float_t omega = (W_bin[i]*W_bin[i] + Q2_bin - m_proton*m_proton)/2./m_proton ;
Float_t en_elp = 2.039 - omega;
Float_t th_elp = 2*asin(sqrt(Q2_bin/4./2.039/en_elp));

Float_t epsilon = 1/(1. + 2.*(1. + omega*omega/Q2_bin)*(tan(th_elp/2.))*(tan(th_elp/2.)));
Float_t flux= (omega-Q2_bin/2./m_proton)/137.;

flux= flux/2./(M_PI)/2.039/Q2_bin/(1-epsilon);
flux = flux*W_bin[i]/2.039/m_proton; 

return flux;

};

void draw_1d_canvas( Int_t i, Int_t qq2 ) {

 c->Divide(3,3); 
 
leg = new TLegend(0.3,0.965,0.99,1.0); 
leg->SetNColumns(3);
leg->SetFillStyle(0);

// pi+ proton mass model

qqq << "Q^{2} = " << Q2_bin << " GeV^{2}" <<", W = " << W_bin[i] <<" GeV" ;

draw_1d_hist(1,model_1[i]->Projection(0,""),qqq.str(),"model_pipP_1_","d#sigma/dM (#mubn/GeV)","m_{#pi^{+}p} (GeV)",3,"e1","mass",i);
qqq.str("");

 
// pi+ proton mass


draw_1d_hist(1,empty_fill_err(h_cr_sect_1[i]->Projection(0,""),h_cr_sect_noemptcells_1[i]->Projection(0,"")),"","h1prj_inv_m_pip_p_1_","d#sigma/dM (#mubn/GeV)","m_{#pi^{+}p} (GeV)",1,"e1AP same","mass",i);
qqq.str("");


// pi+ proton mass no empty cells filling

draw_1d_hist(1,h_cr_sect_noemptcells_1[i]->Projection(0,""),"","h1prj_inv_m_pip_p_1_noempty","d#sigma/dM (#mubn/GeV)","m_{#pi^{+}p} (GeV)",2,"e1AP same","mass",i);

// pi+ pi- mass model

     draw_1d_hist(2,model_1[i]->Projection(1,""),"","model_pippim_1_","d#sigma/dM (#mubn/GeV)","m_{#pi^{+}#pi^{-}} (GeV)",3,"e1","mass",i);

// pi+ pi- mass

draw_1d_hist(2,empty_fill_err(h_cr_sect_1[i]->Projection(1,""),h_cr_sect_noemptcells_1[i]->Projection(1,"")),"","h1prj_inv_m_pip_pim_1_","d#sigma/dM (#mubn/GeV)","m_{#pi^{+}#pi^{-}} (GeV)",1,"e1AP same","mass",i);

// pi+ pi- mass  no empty cells filling

draw_1d_hist(2,h_cr_sect_noemptcells_1[i]->Projection(1,""),"","h1prj_inv_m_pip_pim_1_noempty","d#sigma/dM (#mubn/GeV)","m_{#pi^{+}#pi^{-}} (GeV)",2,"e1AP same","mass",i);

// pi- proton mass model

draw_1d_hist(3,model_3[i]->Projection(0,""),"","model_pimP_1_","d#sigma/dM (#mubn/GeV)","m_{#pi^{-}p} (GeV)",3,"e1","mass",i);

// pi- proton mass

draw_1d_hist(3,empty_fill_err(h_cr_sect_3[i]->Projection(0,""),h_cr_sect_noemptcells_3[i]->Projection(0,"")),"","h3prj_inv_m_pim_p_1_","d#sigma/dM (#mubn/GeV)","m_{#pi^{-}p} (GeV)",1,"e1AP same","mass",i);

// pi- proton mass  no empty cells filling

draw_1d_hist(3,h_cr_sect_noemptcells_3[i]->Projection(0,""),"","h3prj_inv_m_pim_p_1_noempty","d#sigma/dM (#mubn/GeV)","m_{#pi^{-}p}, GeV",2,"e1AP same","mass",i);

// theta proton model

draw_1d_hist(4,model_1[i]->Projection(2,""),"","model_thetaP_","d#sigma/d(-cos#theta) (#mubn/rad)","#theta_{p'} in c.m. (deg)",3,"e1","theta",i);

// theta proton

draw_1d_hist(4,empty_fill_err(h_cr_sect_1[i]->Projection(2,""),h_cr_sect_noemptcells_1[i]->Projection(2,"")),"","h1prj_th_P_","d#sigma/d(-cos#theta) (#mubn/rad)","#theta_{p'} in c.m. (deg)",1,"e1AP same","theta",i);

// theta proton  no empty cells filling

draw_1d_hist(4,h_cr_sect_noemptcells_1[i]->Projection(2,""),"","h1prj_th_P_noempt","d#sigma/d(-cos#theta) (#mubn/rad)","#theta_{p'} in c.m. (deg)",2,"e1AP same","theta",i);

// theta pi- model

draw_1d_hist(5,model_2[i]->Projection(2,""),"","model_thetapim_","d#sigma/d(-cos#theta) (#mubn/rad)","#theta_{#pi^{-}} in c.m. (deg)",3,"e1","theta",i);

// theta pi-

draw_1d_hist(5,empty_fill_err(h_cr_sect_2[i]->Projection(2,""),h_cr_sect_noemptcells_2[i]->Projection(2,"")),"","h1prj_th_PIm_","d#sigma/d(-cos#theta) (#mubn/rad)","#theta_{#pi^{-}} in c.m. (deg)",1,"e1AP same","theta",i);


// theta pi-  no empty cells filling

draw_1d_hist(5,h_cr_sect_noemptcells_2[i]->Projection(2,""),"","h1prj_th_PIm_noempt","d#sigma/d(-cos#theta) (#mubn/rad)","#theta_{#pi^{-}} in c.m. (deg)",2,"e1AP same","theta",i);

// theta pi+ model

draw_1d_hist(6,model_3[i]->Projection(2,""),"","model_thetapip_","d#sigma/d(-cos#theta) (#mubn/rad)","#theta_{#pi^{+}} in c.m. (deg)",3,"e1","theta",i);

// theta pi+


draw_1d_hist(6,empty_fill_err(h_cr_sect_3[i]->Projection(2,""),h_cr_sect_noemptcells_3[i]->Projection(2,"")),"","h1prj_th_PIp_","d#sigma/d(-cos#theta) (#mubn/rad)","#theta_{#pi^{+}} in c.m. (deg)",1,"e1AP same","theta",i);

// theta pi+  no empty cells filling

draw_1d_hist(6,h_cr_sect_noemptcells_3[i]->Projection(2,""),"","h1prj_th_PIp_noempt","d#sigma/d(-cos#theta) (#mubn/rad)","#theta_{#pi^{+}} in c.m. (deg)",2,"e1AP same","theta",i);

// alpha proton model

draw_1d_hist(7,model_1[i]->Projection(4,""),"","model_alpha_proton_","d#sigma/d#alpha (#mubn/rad)","#alpha_{p'} (deg)",3,"e1","alpha",i);

// alpha proton

       draw_1d_hist(7,empty_fill_err(h_cr_sect_1[i]->Projection(4,""), h_cr_sect_noemptcells_1[i]->Projection(4,"")),"","h1prj_alpha_PIpPIm_pipf_","d#sigma/d#alpha (#mubn/rad)","#alpha_{p'} (deg)",1,"e1AP same","alpha",i);

// alpha proton  no empty cells filling

draw_1d_hist(7, h_cr_sect_noemptcells_1[i]->Projection(4,""),"","h1prj_alpha_PIpPIm_pipf_noempt","d#sigma/d#alpha (#mubn/rad)","#alpha_{p'} (deg)",2,"e1AP same","alpha",i);

// alpha pi+model

draw_1d_hist(8,model_2[i]->Projection(4,""),"","model_alpha_pip_","d#sigma/d#alpha (#mubn/rad)","alpha_{#pi^{-}} (deg)",3,"e1","alpha",i);

// alpha pi+

draw_1d_hist(8,empty_fill_err( h_cr_sect_2[i]->Projection(4,""),h_cr_sect_noemptcells_2[i]->Projection(4,"")),"","h2prj_alpha_PPIp_piPIm_","d#sigma/d#alpha (#mubn/rad)","alpha_{#pi^{-}} (deg)",1,"e1AP same","alpha",i);

// alpha pi+  no empty cells filling

draw_1d_hist(8,h_cr_sect_noemptcells_2[i]->Projection(4,""),"","h2prj_alpha_PPIp_piPIm_noemptcells","d#sigma/d#alpha (#mubn/rad)","alpha_{#pi^{-}} (deg)",2,"e1AP same","alpha",i);

// alpha pi-

draw_1d_hist(9,model_3[i]->Projection(4,""),"","model_alpha_pim_","d#sigma/d#alpha (#mubn/rad)","#alpha_{#pi^{+}} (deg)",3,"e1","alpha",i);

// alpha pi-

draw_1d_hist(9,empty_fill_err(h_cr_sect_3[i]->Projection(4,""),h_cr_sect_noemptcells_3[i]->Projection(4,"")),"","h3prj_alpha_PPIm_piPIp_","d#sigma/d#alpha (#mubn/rad)","#alpha_{#pi^{+}} (deg)",1,"e1AP same","alpha",i);

// alpha pi-  no empty cells filling

draw_1d_hist(9,h_cr_sect_noemptcells_3[i]->Projection(4,""),"","h3prj_alpha_PPIm_piPIp_noemptcells","d#sigma/d#alpha (#mubn/rad)","#alpha_{#pi^{+}} (deg)",2,"e1AP same","alpha",i);

 c->cd();
 leg->Draw();
 c->Update();

//saving canvas to the file

//qqq << "q2_" << Q2_bin << "_w_" << W_bin[i] << ".eps" ;
qqq.str("");
qqq << "Q2_0" << ((10*(4 + qq2/2) + (qq2%2)*5)/10.)*(1+(qq2%2)*9) << "_0" << ((10*(4 + (qq2+1)/2) + ((qq2+1)%2)*5)/10.)*(1+((qq2+1)%2)*9) << "/" << "w_" << W_bin[i] << ".png" ;
if (qq2 == 11) {
qqq.str("");
qqq << "Q2_095_10/"  << "w_" << W_bin[i] << ".png" ;
};
//qqq.str("");
//gSystem->cd(qqq.str().c_str());
//qqq << "w_" << W_bin[i] << ".png" ;
c->SaveAs(qqq.str().c_str());
//gSystem->cd("/volatile/clas/clase1/gleb/test/2pi_analysis_fermi2/");

qqq.str("");

c->Clear();


};

TH1D *empty_fill_err(TH1D *h1, TH1D *h2) {


Double_t err1,err2;


for (Int_t zzz=1; zzz<=h1->GetNbinsX(); zzz++) {
err1=h1->GetBinError(zzz);
err2= (h1->GetBinContent(zzz) - h2->GetBinContent(zzz))/2.;
h1->SetBinError(zzz,(sqrt(err1*err1+err2*err2)));



};

return h1;


};
