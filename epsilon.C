void epsilon(Float_t Q2){

Float_t W;

TF1 *eps = new TF1("eps_func",eps_func,1.22,2.5,2);
eps->SetParameter(0,Q2);
eps->SetParameter(1,2.);

TCanvas *c = new TCanvas();
c->cd();
//TPad *pad1 = new TPad("pad1","",1.2,1.,1.9,2.);
//pad1->Draw();
//pad1->cd();
//c->cd()->SetPad(1.2,1.,1.9,2.);

eps->Eval(2000.);
eps->SetMaximum(15.);
eps->SetMinimum(0.3);
eps->SetTitle("Q^{2} = 0.325 GeV^{2}");
eps->Draw();
eps->GetXaxis()->SetTitle("W (GeV)");
eps->GetXaxis()->SetTitleSize(0.06);
eps->GetXaxis()->SetTitleOffset(0.7);
eps->GetYaxis()->SetTitle("#varepsilon_{L}");
eps->GetYaxis()->SetTitleSize(0.07);
eps->GetYaxis()->SetTitleOffset(0.5);

TF1 *eps_1 = new TF1("eps_func1",eps_func,1.22,2.5,2);
eps_1->SetParameter(0,Q2);
eps_1->SetParameter(1,2.);
eps_1->SetLineColor(1);
eps_1->Eval(2000.);
eps_1->Draw("same");

TF1 *eps_2 = new TF1("eps_func2",eps_func,1.22,2.5,2);
eps_2->SetParameter(0,Q2);
eps_2->SetParameter(1,1.5);
eps_2->SetLineColor(3);
eps_2->Eval(2000.);
eps_2->Draw("same");

 leg = new TLegend(0.15,0.3,0.55,0.5);
 leg->SetFillColor(kWhite);
    leg->AddEntry(eps_2,"E_{beam} = 1.5 GeV","l");
    leg->AddEntry(eps,"E_{beam} = 6.0 GeV","l");
    leg->AddEntry(eps_1,"E_{beam} = 2.0 GeV","l");
   
    leg->SetTextSize(0.05);
    
    leg->Draw();

c->Update();

TCanvas *c1 = new TCanvas();
c1->cd();
c1->cd()->SetBottomMargin(0.15);
c1->cd()->SetLeftMargin(0.15);
TF1 *eps_ratio_func = new TF1("eps_ratio_func",eps_ratio,1.22,1.85,4);
eps_ratio_func->SetParameter(0,Q2);
eps_ratio_func->SetParameter(1,2.);
eps_ratio_func->SetParameter(2,Q2);
eps_ratio_func->SetParameter(3,1.5);
eps_ratio_func->SetLineColor(3);
eps_ratio_func->Eval(2000.);
eps_ratio_func->SetTitle("Q^{2} = 0.475 GeV^{2}");
eps_ratio_func->GetXaxis()->SetTitle("W (GeV)");
eps_ratio_func->GetXaxis()->SetTitleSize(0.06);
//eps_ratio_func->GetXaxis()->SetTitleOffset(0.7);
eps_ratio_func->GetYaxis()->SetTitle("#varepsilon_{2GeV}/#varepsilon_{1.5GeV}");
eps_ratio_func->GetYaxis()->SetTitleSize(0.07);
//eps_ratio_func->GetYaxis()->SetTitleOffset(0.5);
eps_ratio_func->Draw();

c1->Update();

for (Int_t i=1; i<=20; i++) {
cout << "W = " << W << "Q2 = " << Q2 << "\n";
};



};

Double_t eps_func(Double_t *x, Double_t *par) {
/*
x = W
par[0] = Q2
par[1] = beam energy


*/
Double_t epsilon;

Double_t nu = (x[0]*x[0]+par[0]-0.938*0.938)/2./0.938;
Double_t eprime = par[1] - nu;
Double_t sin = par[0]/4./par[1]/eprime;
Double_t tan = sin/(1-sin);
epsilon = 1 + 2.*(1 + nu*nu/par[0])*tan;
epsilon = 1./epsilon;
epsilon = epsilon*par[0]/nu/nu;

return epsilon;


};



Double_t eps_ratio(Double_t *x, Double_t *par) {
/*
x = W
par[0] = Q2
par[1] = beam energy


*/
Double_t epsilon_1,epsilon_2;

Double_t nu_1 = (x[0]*x[0]+par[0]-0.938*0.938)/2./0.938;
Double_t nu_2 = (x[0]*x[0]+par[2]-0.938*0.938)/2./0.938;
Double_t eprime_1 = par[1] - nu_1;
Double_t eprime_2 = par[3] - nu_2;
Double_t sin_1 = par[0]/4./par[1]/eprime_1;
Double_t sin_2 = par[2]/4./par[3]/eprime_2;
Double_t tan_1 = sin_1/(1-sin_1);
Double_t tan_2 = sin_2/(1-sin_2);
epsilon_1 = 1 + 2.*(1 + nu_1*nu_1/par[0])*tan_1;
epsilon_2 = 1 + 2.*(1 + nu_2*nu_2/par[2])*tan_2;
epsilon_1 = 1./epsilon_1;
epsilon_1 = epsilon_1*par[0]/nu_1/nu_1; 
epsilon_2 = 1./epsilon_2;
epsilon_2 = epsilon_2*par[2]/nu_2/nu_2; 

return epsilon_1/epsilon_2;


};

