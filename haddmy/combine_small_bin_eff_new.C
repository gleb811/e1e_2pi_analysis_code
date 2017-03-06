void combine_small_bin_eff_new() {

ostringstream qqq;
Float_t Q2_bin = 0.975;
Float_t W_bin[21];
Float_t rec1[12][12][10][8][8];
Float_t gen1[12][12][10][8][8];
Double_t eff_err;
THnSparseD *tmp_rec1[4];
THnSparseD *tmp_rec2[4];
THnSparseD *tmp_rec3[4];
THnSparseD *tmp_gen1[4];
THnSparseD *tmp_gen2[4];
THnSparseD *tmp_gen3[4];
THnSparseD *tmp_eff1[4];
THnSparseD *tmp_eff2[4];
THnSparseD *tmp_eff3[4];
THnSparseD *tmp_data1[4];
THnSparseD *tmp_data2[4];
THnSparseD *tmp_data3[4];
THnSparseD *tmp_empty1[4];
THnSparseD *tmp_empty2[4];
THnSparseD *tmp_empty3[4];
TH1D *h_eff_err  = new TH1D("h_eff_err","h_eff_err",1200, 0.,1.2);
TH1D *h_eff1  = new TH1D("h_eff1","h_eff1",1000, 0.,1.);
TH1D *h_eff2  = new TH1D("h_eff2","h_eff2",200, 0.,0.2);
TH1D *h_rec1  = new TH1D("h_rec1","h_rec",2000, 0.5,2000.5);
TH2D *h_eff_vs_rec  = new TH2D("h_eff_vs_rec","h_eff_vs_rec",200, 0.,0.2.,100,0.5,100.5);
TH2D *h_gen_vs_rec  = new TH2D("h_gen_vs_rec","h_gen_vs_rec",2000, 0.5,2000.5,100,0.5,100.5);
TH2D *h_eff_vs_gen  = new TH2D("h_eff_vs_gen","h_eff_vs_gen",200, 0.,0.2.,4000,0.5,4000.5);
Long64_t tmp_rec1_bin,tmp_rec2_bin,tmp_rec3_bin;
Long64_t tmp_empty1_bin,tmp_empty2_bin,tmp_empty3_bin;
Long64_t tmp_gen1_bin,tmp_gen2_bin,tmp_gen3_bin;
Long64_t tmp_eff1_bin,tmp_eff2_bin,tmp_eff3_bin;
Long64_t tmp_data1_bin,tmp_data2_bin,tmp_data3_bin;
Long64_t rec1_bin,rec2_bin,rec3_bin;
Long64_t gen1_bin,gen2_bin,gen3_bin;
Long64_t eff1_bin,eff2_bin,eff3_bin;
Long64_t data1_bin,data2_bin,data3_bin;
Long64_t empty1_bin,empty2_bin,empty3_bin;

Int_t *bins = new Int_t[5];

Double_t xmin[5] = {1.028, 0.229,0.,0.,0.};
Double_t xmax[5];
xmax[2] = 180.;
xmax[3] = 360.;
xmax[4] = 360.;

bins[0]=12;
bins[1]=12;
bins[2]=10;
bins[3]=8;
bins[4]=8;

for (Int_t i=0; i<21; i++) {
W_bin[i] = 1.3125+0.025*i; 
xmax[0] =  W_bin[i] - 0.13957 + 0.05;
xmax[1] =  W_bin[i] - 0.938272 + 0.05;
};

qqq.str("");
qqq << "all_w_pim.root";
TFile *MyFile_pim = new TFile(qqq.str().c_str(),"READ");
qqq.str("");
qqq << "all_w_pip.root";
TFile *MyFile_pip = new TFile(qqq.str().c_str(),"READ");
qqq.str("");
qqq << "all_w_pr.root";
TFile *MyFile_pr = new TFile(qqq.str().c_str(),"READ");
qqq.str("");
qqq << "all_w_0.root";
TFile *MyFile_0 = new TFile(qqq.str().c_str(),"READ");
qqq.str("");
qqq << "../out_data.root";
TFile *MyFile_data = new TFile(qqq.str().c_str(),"READ");
qqq.str("");

TFile *MyFile_out = new TFile("combine_small_bin_eff_time_corr1.root","RECREATE");
TFile *MyFile_data_out = new TFile("combine_data_small_bin_eff_time_corr1.root","RECREATE");

for (Int_t k=0; k<5;k++) {
Q2_bin = 0.425 + k*0.05; 
for (Int_t i=0; i<21;i++) { 
W_bin[i] = 1.3125+0.025*i; 

//READ pim

MyFile_pim->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec1[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec2[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec3[0]);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen1[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen2[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen3[0]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_1_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff1[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_2_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff2[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_3_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff3[0]);

MyFile_data->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data1[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data2[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data3[0]);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_1_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty1[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_2_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty2[0]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pim_3_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty3[0]);

// READ pip

MyFile_pip->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_1_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec1[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_2_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec2[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_3_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec3[1]);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen1[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen2[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen3[1]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_1_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff1[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_2_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff2[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_3_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff3[1]);


MyFile_data->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data1[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data2[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_3_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data3[1]);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_1_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty1[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_2_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty2[1]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pip_3_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty3[1]);

// READ proton


MyFile_pr->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_1_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec1[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_2_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec2[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_3_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec3[2]);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen1[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen2[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen3[2]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_1_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff1[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_2_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff2[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_3_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff3[2]);


MyFile_data->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data1[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data2[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_3_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data3[2]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_1_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty1[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_2_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty2[2]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_3_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty3[2]);


// READ exclusive


MyFile_0->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_1_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec1[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_2_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec2[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_3_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec3[3]);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen1[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen2[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen3[3]);


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_1_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff1[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_2_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff2[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_eff_3_" <<i;
gDirectory->GetObject(qqq.str().c_str(),tmp_eff3[3]);

MyFile_data->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data1[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data2[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_3_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_data3[3]);



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_1_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty1[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_2_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty2[3]);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_0_3_empty_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_empty3[3]);

// Begin merging



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

for (Int_t top=1; top<4; top++) {

eff1_bin = tmp_eff1[0]->GetBin(bins);
tmp_eff1_bin = tmp_eff1[top]->GetBin(bins);

rec1_bin = tmp_rec1[0]->GetBin(bins);
tmp_rec1_bin = tmp_rec1[top]->GetBin(bins);

gen1_bin = tmp_gen1[0]->GetBin(bins);
tmp_gen1_bin = tmp_gen1[top]->GetBin(bins);

data1_bin = tmp_data1[0]->GetBin(bins);
tmp_data1_bin = tmp_data1[top]->GetBin(bins);

empty1_bin = tmp_empty1[0]->GetBin(bins);
tmp_empty1_bin = tmp_empty1[top]->GetBin(bins);

if ((tmp_rec1[0]->GetBinContent(rec1_bin)) < (tmp_rec1[top]->GetBinContent(tmp_rec1_bin))) {

tmp_rec1[0]->SetBinContent(rec1_bin,tmp_rec1[top]->GetBinContent(tmp_rec1_bin));
tmp_eff1[0]->SetBinContent(eff1_bin,tmp_eff1[top]->GetBinContent(tmp_eff1_bin));
tmp_data1[0]->SetBinContent(data1_bin,tmp_data1[top]->GetBinContent(tmp_data1_bin));
tmp_empty1[0]->SetBinContent(empty1_bin,tmp_empty1[top]->GetBinContent(tmp_empty1_bin));
};

/*
if ( ((tmp_data1[0]->GetBinContent(data1_bin)) <=0.1 ) &&((tmp_eff1[0]->GetBinContent(eff1_bin)) < (tmp_eff1[top]->GetBinContent(tmp_eff1_bin))) ) {

tmp_rec1[0]->SetBinContent(rec1_bin,tmp_rec1[top]->GetBinContent(tmp_rec1_bin));
tmp_eff1[0]->SetBinContent(eff1_bin,tmp_eff1[top]->GetBinContent(tmp_eff1_bin));
tmp_data1[0]->SetBinContent(data1_bin,tmp_data1[top]->GetBinContent(tmp_data1_bin));
tmp_empty1[0]->SetBinContent(empty1_bin,tmp_empty1[top]->GetBinContent(tmp_empty1_bin));
};
*/







eff2_bin = tmp_eff2[0]->GetBin(bins);
tmp_eff2_bin = tmp_eff2[top]->GetBin(bins);

rec2_bin = tmp_rec2[0]->GetBin(bins);
tmp_rec2_bin = tmp_rec2[top]->GetBin(bins);

gen2_bin = tmp_gen2[0]->GetBin(bins);
tmp_gen2_bin = tmp_gen2[top]->GetBin(bins);

data2_bin = tmp_data2[0]->GetBin(bins);
tmp_data2_bin = tmp_data2[top]->GetBin(bins);

empty2_bin = tmp_empty2[0]->GetBin(bins);
tmp_empty2_bin = tmp_empty2[top]->GetBin(bins);

if ((tmp_rec2[0]->GetBinContent(rec2_bin)) < (tmp_rec2[top]->GetBinContent(tmp_rec2_bin))) {

tmp_rec2[0]->SetBinContent(rec2_bin,tmp_rec2[top]->GetBinContent(tmp_rec2_bin));
tmp_eff2[0]->SetBinContent(eff2_bin,tmp_eff2[top]->GetBinContent(tmp_eff2_bin));
tmp_data2[0]->SetBinContent(data2_bin,tmp_data2[top]->GetBinContent(tmp_data2_bin));
tmp_empty2[0]->SetBinContent(empty2_bin,tmp_empty2[top]->GetBinContent(tmp_empty2_bin));
};

/*
if ( ((tmp_data2[0]->GetBinContent(data2_bin)) <=0.1 ) &&((tmp_eff2[0]->GetBinContent(eff2_bin)) < (tmp_eff2[top]->GetBinContent(tmp_eff2_bin))) ) {

tmp_rec2[0]->SetBinContent(rec2_bin,tmp_rec2[top]->GetBinContent(tmp_rec2_bin));
tmp_eff2[0]->SetBinContent(eff2_bin,tmp_eff2[top]->GetBinContent(tmp_eff2_bin));
tmp_data2[0]->SetBinContent(data2_bin,tmp_data2[top]->GetBinContent(tmp_data2_bin));
tmp_empty2[0]->SetBinContent(empty2_bin,tmp_empty2[top]->GetBinContent(tmp_empty2_bin));
};
*/






eff3_bin = tmp_eff3[0]->GetBin(bins);
tmp_eff3_bin = tmp_eff3[top]->GetBin(bins);

rec3_bin = tmp_rec3[0]->GetBin(bins);
tmp_rec3_bin = tmp_rec3[top]->GetBin(bins);

gen3_bin = tmp_gen3[0]->GetBin(bins);
tmp_gen3_bin = tmp_gen3[top]->GetBin(bins);

data3_bin = tmp_data3[0]->GetBin(bins);
tmp_data3_bin = tmp_data3[top]->GetBin(bins);

empty3_bin = tmp_empty3[0]->GetBin(bins);
tmp_empty3_bin = tmp_empty3[top]->GetBin(bins);

if ((tmp_rec3[0]->GetBinContent(rec3_bin)) < (tmp_rec3[top]->GetBinContent(tmp_rec3_bin))) {

tmp_rec3[0]->SetBinContent(rec3_bin,tmp_rec3[top]->GetBinContent(tmp_rec3_bin));
tmp_eff3[0]->SetBinContent(eff3_bin,tmp_eff3[top]->GetBinContent(tmp_eff3_bin));
tmp_data3[0]->SetBinContent(data3_bin,tmp_data3[top]->GetBinContent(tmp_data3_bin));
tmp_empty3[0]->SetBinContent(empty3_bin,tmp_empty3[top]->GetBinContent(tmp_empty3_bin));
};

/*
if ( ((tmp_data3[0]->GetBinContent(data3_bin)) <=0.1 ) &&((tmp_eff3[0]->GetBinContent(eff3_bin)) < (tmp_eff3[top]->GetBinContent(tmp_eff3_bin))) ) {

tmp_rec3[0]->SetBinContent(rec3_bin,tmp_rec3[top]->GetBinContent(tmp_rec3_bin));
tmp_eff3[0]->SetBinContent(eff3_bin,tmp_eff3[top]->GetBinContent(tmp_eff3_bin));
tmp_data3[0]->SetBinContent(data3_bin,tmp_data3[top]->GetBinContent(tmp_data3_bin));
tmp_empty3[0]->SetBinContent(empty3_bin,tmp_empty3[top]->GetBinContent(tmp_empty3_bin));
};

*/


}; // loop over topologies

};
};
};
};
};


// Output to the file
MyFile_out->cd();
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile_out->mkdir(qqq.str().c_str());
MyFile_out->cd(qqq.str().c_str());

tmp_rec1[0]->Write();
tmp_rec2[0]->Write();
tmp_rec3[0]->Write();

tmp_gen1[0]->Write();
tmp_gen2[0]->Write();
tmp_gen3[0]->Write();

tmp_eff1[0]->Write();
tmp_eff2[0]->Write();
tmp_eff3[0]->Write();

tmp_rec1[0]->Delete();
tmp_rec2[0]->Delete();
tmp_rec3[0]->Delete();

tmp_gen1[0]->Delete();
tmp_gen2[0]->Delete();
tmp_gen3[0]->Delete();

tmp_eff1[0]->Delete();
tmp_eff2[0]->Delete();
tmp_eff3[0]->Delete();

// Output to the file
MyFile_data_out->cd();
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile_data_out->mkdir(qqq.str().c_str());
MyFile_data_out->cd(qqq.str().c_str());

tmp_data1[0]->Write();
tmp_data2[0]->Write();
tmp_data3[0]->Write();

tmp_empty1[0]->Write();
tmp_empty2[0]->Write();
tmp_empty3[0]->Write();

tmp_data1[0]->Delete();
tmp_data2[0]->Delete();
tmp_data3[0]->Delete();

tmp_empty1[0]->Delete();
tmp_empty2[0]->Delete();
tmp_empty3[0]->Delete();


//

};
};

MyFile_pim->Close();
MyFile_pip->Close();
MyFile_pr->Close();
MyFile_0->Close();
MyFile_out->Close();
MyFile_data_out->Close();
};

