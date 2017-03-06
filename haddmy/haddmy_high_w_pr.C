void haddmy_high_w_pr() {

ostringstream qqq;
Float_t Q2_bin = 0.975;
Float_t W_bin[9];
Float_t rec1[12][12][10][8][8];
Float_t gen1[12][12][10][8][8];
Double_t eff_err;
THnSparseD *h_rec_1[9][12];
THnSparseD *h_rec_2[9][12];
THnSparseD *h_rec_3[9][12];
THnSparseD *h_gen_1[9][12];
THnSparseD *h_gen_2[9][12];
THnSparseD *h_gen_3[9][12];
THnSparseD *h_eff_1[9][12];
THnSparseD *h_eff_2[9][12];
THnSparseD *h_eff_3[9][12];
THnSparseD *tmp_rec1;
THnSparseD *tmp_rec2;
THnSparseD *tmp_rec3;
THnSparseD *tmp_gen1;
THnSparseD *tmp_gen2;
THnSparseD *tmp_gen3;
TH1D *h_eff_err  = new TH1D("h_eff_err","h_eff_err",1200, 0.,1.2);
TH1D *h_eff1  = new TH1D("h_eff1","h_eff1",1000, 0.,1.);
TH1D *h_eff2  = new TH1D("h_eff2","h_eff2",200, 0.,0.2);
TH1D *h_rec1  = new TH1D("h_rec1","h_rec",2000, 0.5,2000.5);
TH2D *h_eff_vs_rec  = new TH2D("h_eff_vs_rec","h_eff_vs_rec",200, 0.,0.2.,100,0.5,100.5);
TH2D *h_gen_vs_rec  = new TH2D("h_gen_vs_rec","h_gen_vs_rec",2000, 0.5,2000.5,100,0.5,100.5);
TH2D *h_eff_vs_gen  = new TH2D("h_eff_vs_gen","h_eff_vs_gen",200, 0.,0.2.,4000,0.5,4000.5);
Long64_t tmp_rec1_bin,tmp_rec2_bin,tmp_rec3_bin;
Long64_t tmp_gen1_bin,tmp_gen2_bin,tmp_gen3_bin;
Long64_t rec1_bin,rec2_bin,rec3_bin;
Long64_t gen1_bin,gen2_bin,gen3_bin;
Long64_t eff1_bin,eff2_bin,eff3_bin;

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

for (Int_t i=0; i<6; i++) {
W_bin[i] = 1.6875+0.025*i; 
xmax[0] =  W_bin[i] - 0.13957 + 0.05;
xmax[1] =  W_bin[i] - 0.938272 + 0.05;
/*qqq <<  "h_5dim_1_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
h_rec_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),5,bins,xmin,xmax);
qqq.str("");
qqq <<"h_5dim_2_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
h_rec_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),5,bins,xmin,xmax);
qqq.str("");
qqq << "h_5dim_3_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
h_rec_3[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),5,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_1_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
h_gen_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),5,bins,xmin,xmax);
qqq.str("");
qqq <<"h_5dim_2_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
h_gen_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),5,bins,xmin,xmax);
qqq.str("");
qqq << "h_5dim_3_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
h_gen_3[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),5,bins,xmin,xmax);
qqq.str("");*/
};


for (Int_t j=0; j<=99; j++) {
qqq.str("");
qqq << "../w_high_" << j << ".root";
TFile *MyFile = new TFile(qqq.str().c_str(),"READ");
MyFile->cd();
qqq.str("");
for (Int_t k=0; k<12;k++) {
Q2_bin = 0.425 + k*0.05; 
for (Int_t i=0; i<6;i++) { 
W_bin[i] = 1.6875+0.025*i; 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_1_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec1);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_2_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec2);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_pr_3_sim_1_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_rec3);

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_1_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen1);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_2_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen2);
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h_5dim_3_sim_2_q2_" << Q2_bin*1000 << "_w_" << 1000*W_bin[i];
gDirectory->GetObject(qqq.str().c_str(),tmp_gen3);

if (j == 0) {
h_rec_1[i][k] = (THnSparseD*)tmp_rec1->Clone();
h_rec_2[i][k] = (THnSparseD*)tmp_rec2->Clone();
h_rec_3[i][k] = (THnSparseD*)tmp_rec3->Clone();

h_gen_1[i][k] = (THnSparseD*)tmp_gen1->Clone();
h_gen_2[i][k] = (THnSparseD*)tmp_gen2->Clone();
h_gen_3[i][k] = (THnSparseD*)tmp_gen3->Clone();


/*for (Int_t o=1; o<=12; o++) {
for (Int_t p=1; p<=12; p++) {
for (Int_t r=1; r<=10; r++) {
for (Int_t t=1; t<=8; t++) {
for (Int_t y=1; y<=8; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;
tmp_rec1_bin = tmp_rec1->GetBin(bins);
tmp_gen1_bin = tmp_gen1->GetBin(bins);
rec1_bin = h_rec_1[i]->GetBin(bins);
gen1_bin = h_gen_1[i]->GetBin(bins);
tmp_rec2_bin = tmp_rec2->GetBin(bins);
tmp_gen2_bin = tmp_gen2->GetBin(bins);
rec2_bin = h_rec_2[i]->GetBin(bins);
gen2_bin = h_gen_2[i]->GetBin(bins);
tmp_rec3_bin = tmp_rec3->GetBin(bins);
tmp_gen3_bin = tmp_gen3->GetBin(bins);
rec3_bin = h_rec_3[i]->GetBin(bins);
gen3_bin = h_gen_3[i]->GetBin(bins);

if ((tmp_rec1->GetBinContent(tmp_rec1_bin) > 0.0001)) {
h_rec_1[i]->SetBinContent(rec1_bin,tmp_rec1->GetBinContent(tmp_rec1_bin));
};
if ((tmp_gen1->GetBinContent(tmp_gen1_bin) > 0.0001)) {
h_gen_1[i]->SetBinContent(gen1_bin,tmp_gen1->GetBinContent(tmp_gen1_bin));
};

if ((tmp_rec2->GetBinContent(tmp_rec2_bin) > 0.0001)) {
h_rec_2[i]->SetBinContent(rec2_bin,tmp_rec2->GetBinContent(tmp_rec2_bin));
};
if ((tmp_gen2->GetBinContent(tmp_gen2_bin) > 0.0001)) {
h_gen_2[i]->SetBinContent(gen2_bin,tmp_gen2->GetBinContent(tmp_gen2_bin));
};

if ((tmp_rec3->GetBinContent(tmp_rec3_bin) > 0.0001)) {
h_rec_3[i]->SetBinContent(rec3_bin,tmp_rec3->GetBinContent(tmp_rec3_bin));
};
if ((tmp_gen3->GetBinContent(tmp_gen3_bin) > 0.0001)) {
h_gen_3[i]->SetBinContent(gen3_bin,tmp_gen3->GetBinContent(tmp_gen3_bin));
};


if (i ==7){
if (tmp_rec1->GetBinContent(tmp_rec1_bin) > 0) {
rec1[o-1][p-1][r-1][t-1][y-1] = rec1[o-1][p-1][r-1][t-1][y-1] + tmp_rec1->GetBinContent(tmp_rec1_bin);
};
if (tmp_gen1->GetBinContent(tmp_gen1_bin) > 0) {
gen1[o-1][p-1][r-1][t-1][y-1] = gen1[o-1][p-1][r-1][t-1][y-1] + tmp_gen1->GetBinContent(tmp_gen1_bin);
};
//cout << h_gen_1[i]->GetBinContent(gen1_bin) << "          " << gen1[o-1][p-1][r-1][t-1][y-1] << "\n";
};
/*if(tmp_rec1->GetBinContent(bins)>0.1)h_rec_1[i]->SetBinContent(bins,h_rec_1[i]->GetBinContent(bins)+tmp_rec1->GetBinContent(bins));
h_rec_2[i]->SetBinContent(bins,h_rec_2[i]->GetBinContent(bins)+tmp_rec2->GetBinContent(bins));
h_rec_3[i]->SetBinContent(bins,h_rec_3[i]->GetBinContent(bins)+tmp_rec3->GetBinContent(bins));

if(tmp_gen1->GetBinContent(bins)>0.1)h_gen_1[i]->SetBinContent(bins,h_gen_1[i]->GetBinContent(bins)+tmp_gen1->GetBinContent(bins));
h_gen_2[i]->SetBinContent(bins,h_gen_2[i]->GetBinContent(bins)+tmp_gen2->GetBinContent(bins));
h_gen_3[i]->SetBinContent(bins,h_gen_3[i]->GetBinContent(bins)+tmp_gen3->GetBinContent(bins));


};
};
};
};
};
*/

} else {


h_rec_1[i][k]->Add(tmp_rec1);
h_rec_2[i][k]->Add(tmp_rec2);
h_rec_3[i][k]->Add(tmp_rec3);

h_gen_1[i][k]->Add(tmp_gen1);
h_gen_2[i][k]->Add(tmp_gen2);
h_gen_3[i][k]->Add(tmp_gen3);



/*for (Int_t o=1; o<=12; o++) {
for (Int_t p=1; p<=12; p++) {
for (Int_t r=1; r<=10; r++) {
for (Int_t t=1; t<=8; t++) {
for (Int_t y=1; y<=8; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;
tmp_rec1_bin = tmp_rec1->GetBin(bins);
tmp_gen1_bin = tmp_gen1->GetBin(bins);
rec1_bin = h_rec_1[i]->GetBin(bins);
gen1_bin = h_gen_1[i]->GetBin(bins);
tmp_rec2_bin = tmp_rec2->GetBin(bins);
tmp_gen2_bin = tmp_gen2->GetBin(bins);
rec2_bin = h_rec_2[i]->GetBin(bins);
gen2_bin = h_gen_2[i]->GetBin(bins);
tmp_rec3_bin = tmp_rec3->GetBin(bins);
tmp_gen3_bin = tmp_gen3->GetBin(bins);
rec3_bin = h_rec_3[i]->GetBin(bins);
gen3_bin = h_gen_3[i]->GetBin(bins);

if ((tmp_rec1->GetBinContent(tmp_rec1_bin) > 0.0001)) {
h_rec_1[i]->SetBinContent(rec1_bin,h_rec_1[i]->GetBinContent(rec1_bin)+tmp_rec1->GetBinContent(tmp_rec1_bin));
};
if ((tmp_gen1->GetBinContent(tmp_gen1_bin) > 0.0001)) {
h_gen_1[i]->SetBinContent(gen1_bin,h_gen_1[i]->GetBinContent(gen1_bin)+tmp_gen1->GetBinContent(tmp_gen1_bin));
};

if ((tmp_rec2->GetBinContent(tmp_rec2_bin) > 0.0001)) {
h_rec_2[i]->AddBinContent(rec2_bin,tmp_rec2->GetBinContent(tmp_rec2_bin));
};
if ((tmp_gen2->GetBinContent(tmp_gen2_bin) > 0.0001)) {
h_gen_2[i]->AddBinContent(gen2_bin,tmp_gen2->GetBinContent(tmp_gen2_bin));
};

if ((tmp_rec3->GetBinContent(tmp_rec3_bin) > 0.0001))) {
h_rec_3[i]->AddBinContent(rec3_bin,tmp_rec3->GetBinContent(tmp_rec3_bin));
};
if ((tmp_gen3->GetBinContent(tmp_gen3_bin) > 0.0001)) {
h_gen_3[i]->AddBinContent(gen3_bin,tmp_gen3->GetBinContent(tmp_gen3_bin));
};
//rec1_bin = h_rec_1[i]->GetBin(bins);
//cout << tmp_rec1->GetBinContent(tmp_rec1_bin) << "\n";
//gen1_bin = h_gen_1[i]->GetBin(bins);

//h_rec_1[i]->AddBinContent(rec1_bin,tmp_rec1->GetBinContent(tmp_rec1_bin));
//h_gen_1[i]->AddBinContent(gen1_bin,tmp_gen1->GetBinContent(tmp_gen1_bin));

if (i == 7) {
if (tmp_rec1->GetBinContent(tmp_rec1_bin) > 0) {
rec1[o-1][p-1][r-1][t-1][y-1] = rec1[o-1][p-1][r-1][t-1][y-1] + tmp_rec1->GetBinContent(tmp_rec1_bin);
};
if (tmp_gen1->GetBinContent(tmp_gen1_bin) > 0){
gen1[o-1][p-1][r-1][t-1][y-1] = gen1[o-1][p-1][r-1][t-1][y-1] + tmp_gen1->GetBinContent(tmp_gen1_bin);
};
};
/*if(tmp_rec1->GetBinContent(bins)>0.1)h_rec_1[i]->SetBinContent(bins,h_rec_1[i]->GetBinContent(bins)+tmp_rec1->GetBinContent(bins));
h_rec_2[i]->SetBinContent(bins,h_rec_2[i]->GetBinContent(bins)+tmp_rec2->GetBinContent(bins));
h_rec_3[i]->SetBinContent(bins,h_rec_3[i]->GetBinContent(bins)+tmp_rec3->GetBinContent(bins));

if(tmp_gen1->GetBinContent(bins)>0.1)h_gen_1[i]->SetBinContent(bins,h_gen_1[i]->GetBinContent(bins)+tmp_gen1->GetBinContent(bins));
h_gen_2[i]->SetBinContent(bins,h_gen_2[i]->GetBinContent(bins)+tmp_gen2->GetBinContent(bins));
h_gen_3[i]->SetBinContent(bins,h_gen_3[i]->GetBinContent(bins)+tmp_gen3->GetBinContent(bins));


};
};
};
};
};
*/

};
}; 
};
MyFile->Close();
};

for (k=0; k<12;k++) { 
for (i=0; i<6;i++) { 
qqq.str("");
qqq << "h_eff_1_" << i+15;
h_eff_1[i][k]= (THnSparseD*)h_rec_1[i][k]->Clone(qqq.str().c_str());
h_eff_1[i][k]->Divide(h_gen_1[i][k]);
qqq.str("");
qqq << "h_eff_2_" << i+15;
h_eff_2[i][k]= (THnSparseD*)h_rec_2[i][k]->Clone(qqq.str().c_str());
h_eff_2[i][k]->Divide(h_gen_2[i][k]);
qqq.str("");
qqq << "h_eff_3_" << i+15;
h_eff_3[i][k]= (THnSparseD*)h_rec_3[i][k]->Clone(qqq.str().c_str());
h_eff_3[i][k]->Divide(h_gen_3[i][k]);
qqq.str("");

for (Int_t o=1; o<=12; o++) {
for (Int_t p=1; p<=12; p++) {
for (Int_t r=1; r<=10; r++) {
for (Int_t t=1; t<=8; t++) {
for (Int_t y=1; y<=8; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;

rec1_bin = h_rec_1[i][k]->GetBin(bins);
gen1_bin = h_gen_1[i][k]->GetBin(bins);
eff1_bin = h_eff_1[i][k]->GetBin(bins);

rec2_bin = h_rec_2[i][k]->GetBin(bins);
gen2_bin = h_gen_2[i][k]->GetBin(bins);
eff2_bin = h_eff_2[i][k]->GetBin(bins);

rec3_bin = h_rec_3[i][k]->GetBin(bins);
gen3_bin = h_gen_3[i][k]->GetBin(bins);
eff3_bin = h_eff_3[i][k]->GetBin(bins);


if ((h_eff_1[i][k]->GetBinContent(eff1_bin)>0)){
eff_err = abs(h_gen_1[i][k]->GetBinContent(gen1_bin) - h_rec_1[i][k]->GetBinContent(rec1_bin));
eff_err = eff_err/(h_gen_1[i][k]->GetBinContent(gen1_bin));
eff_err = eff_err/(h_gen_1[i][k]->GetBinContent(gen1_bin));
eff_err = eff_err/(h_gen_1[i][k]->GetBinContent(gen1_bin));

eff_err = eff_err*(h_rec_1[i][k]->GetBinContent(rec1_bin));
eff_err = sqrt(eff_err);
h_eff_1[i][k]->SetBinError(eff1_bin,eff_err);
} else {
h_eff_1[i][k]->SetBinError(eff1_bin,0);
};


if ((h_eff_2[i][k]->GetBinContent(eff2_bin)>0)){
eff_err = abs(h_gen_2[i][k]->GetBinContent(gen2_bin) - h_rec_2[i][k]->GetBinContent(rec2_bin));
eff_err = eff_err/(h_gen_2[i][k]->GetBinContent(gen2_bin));
eff_err = eff_err/(h_gen_2[i][k]->GetBinContent(gen2_bin));
eff_err = eff_err/(h_gen_2[i][k]->GetBinContent(gen2_bin));

eff_err = eff_err*(h_rec_2[i][k]->GetBinContent(rec2_bin));
eff_err = sqrt(eff_err);
h_eff_2[i][k]->SetBinError(eff2_bin,eff_err);
} else {
h_eff_2[i][k]->SetBinError(eff2_bin,0);
};


if ((h_eff_3[i][k]->GetBinContent(eff3_bin)>0)){
eff_err = abs(h_gen_3[i][k]->GetBinContent(gen3_bin) - h_rec_3[i][k]->GetBinContent(rec3_bin));
eff_err = eff_err/(h_gen_3[i][k]->GetBinContent(gen3_bin));
eff_err = eff_err/(h_gen_3[i][k]->GetBinContent(gen3_bin));
eff_err = eff_err/(h_gen_3[i][k]->GetBinContent(gen3_bin));

eff_err = eff_err*(h_rec_3[i][k]->GetBinContent(rec3_bin));
eff_err = sqrt(eff_err);
h_eff_3[i][k]->SetBinError(eff3_bin,eff_err);
} else {
h_eff_3[i][k]->SetBinError(eff3_bin,0);
};




};
};
};
};
};

};
};

TFile *MyFile = new TFile("high_w_pr.root","RECREATE");
MyFile->cd();

for (k=0; k<12;k++) { 
Q2_bin = 0.425 + 0.05*k;
for (i=0; i<6;i++) { 
qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i];
MyFile->mkdir(qqq.str().c_str());
MyFile->cd(qqq.str().c_str());
h_rec_1[i][k]->Write();
h_rec_2[i][k]->Write();
h_rec_3[i][k]->Write();

h_gen_1[i][k]->Write();
h_gen_2[i][k]->Write();
h_gen_3[i][k]->Write();
/*qqq.str("");
qqq << "h_eff_1_" << i;
h_eff_1[i]= (THnSparseD*)h_rec_1[i]->Clone(qqq.str().c_str());
h_eff_1[i]->Divide(h_gen_1[i]);
qqq.str("");
qqq << "h_eff_2_" << i;
h_eff_2[i]= (THnSparseD*)h_rec_2[i]->Clone(qqq.str().c_str());
h_eff_2[i]->Divide(h_gen_2[i]);
qqq.str("");
qqq << "h_eff_3_" << i;
h_eff_3[i]= (THnSparseD*)h_rec_3[i]->Clone(qqq.str().c_str());
h_eff_3[i]->Divide(h_gen_3[i]);
qqq.str("");
*/
h_eff_1[i][k]->Write();
h_eff_2[i][k]->Write();
h_eff_3[i][k]->Write();
};
};
MyFile->Close();


//h_eff_1[0]= (THnSparseD*)h_rec_1[0]->Clone("h_eff_1");
//h_eff_1[0]->Divide(h_gen_1[0]);

/*for (Int_t o=1; o<=12; o++) {
for (Int_t p=1; p<=12; p++) {
for (Int_t r=1; r<=10; r++) {
for (Int_t t=1; t<=8; t++) {
for (Int_t y=1; y<=8; y++) {
bins[0] = o;
bins[1] = p;
bins[2] = r;
bins[3] = t;
bins[4] = y;

rec1_bin = h_rec_1[7]->GetBin(bins);
gen1_bin = h_gen_1[7]->GetBin(bins);
eff1_bin = h_eff_1[7]->GetBin(bins);

if ((h_eff_1[7]->GetBinContent(eff1_bin)>0)&&(h_eff_1[7]->GetBinContent(eff1_bin)<1.)){
eff_err = (h_eff_1[7]->GetBinError(eff1_bin))/(h_eff_1[7]->GetBinContent(eff1_bin));
h_eff_err-> Fill(eff_err);
//if (eff_err<0.2) {
h_eff1->Fill(h_eff_1[7]->GetBinContent(eff1_bin));
h_rec1->Fill(h_rec_1[7]->GetBinContent(rec1_bin));
//};
};

};
};
};
};
};
*/

/*TCanvas *c = new TCanvas("c","c",600,600);
c->Divide(2,2);
//TCanvas *c1;
c->cd(1);
h_eff_err->Draw();
c->cd(2);
h_eff1->Draw();
c->cd(3);
h_rec1->Draw();
*/
//TCanvas *c2;
//c2->cd();
//h_eff_vs_rec->Draw("colz");
//h_gen_vs_rec->Draw("colz");
//h_eff_vs_gen->Draw("colz");
}
