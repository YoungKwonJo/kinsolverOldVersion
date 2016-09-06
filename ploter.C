#include <iostream>
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TFile.h"
#include <iostream>
#include "TROOT.h"

void ploter()
{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0); //remove statistics box
    gStyle->SetOptTitle(0); //remove title
    gStyle->SetPalette(1);
    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.06, "XYZ");
    gStyle->SetTitleXOffset(1.2);
    gStyle->SetTitleYOffset(1.2);
    gStyle->SetTitleX(0.2);
    gStyle->SetTitleY(0.985);


  gROOT->ProcessLine(".L tdrstyle.C");
  defaultStyle();

  ploter("MuMu");
  ploter("ElEl");
  ploter("MuEl");
  ploter(11);

}

void ploter(TString decay)
{

    TFile* f1 = new TFile(Form("result_%s_TTbarFullLepMGDecays.root",decay.Data()  ) );
    TFile* f2 = new TFile(Form("result_%s_Run2012%s.root",decay.Data(),decay.Data()) );

    TString hh[16] = {"htM","htbarM","httM","hweight","htM2","htbarM2","httM2","hweight2","hlnuM","hlnuM2","hlnubarM","hlnubarM2","hCSVD11","hCSVD12","hCSVD21","hCSVD22" };
    TString yx[16] = {"M_{t} (GeV/c^{2})", "M_{#bar{t}} (GeV/c^{2})","M_{t#bar{t}} (GeV/c^{2})","weight",
                     "M_{t} (GeV/c^{2})", "M_{#bar{t}} (GeV/c^{2})","M_{t#bar{t}} (GeV/c^{2})","weight",
                     "M_{l#nu} (GeV/c^{2})","M_{l#nu} (GeV/c^{2})", "M_{l#bar{#nu}} (GeV/c^{2})","M_{l#bar{#nu}} (GeV/c^{2})",
                     "b-Discriminator (CSV)","b-Discriminator (CSV)","b-Discriminator (CSV)","b-Discriminator (CSV)"  };

    for(int i=0;i<16;i++)
    {
        TH1F* h1 = (TH1*)f1->Get(Form("%s",hh[i].Data()));
        TH1F* h2 = (TH1*)f2->Get(Form("%s",hh[i].Data()));

        h2->Sumw2();
        int nBin = h1->GetNbinsX();
        h1->AddBinContent(nBin,h1->GetBinContent(nBin+1));
        h2->AddBinContent(nBin,h2->GetBinContent(nBin+1));

        h1->Scale(1/h1->Integral() );
        h2->Scale(1/h2->Integral() );
        h1->GetYaxis()->SetTitle("Normalized Entries");
        h1->GetXaxis()->SetTitle(yx[i].Data());

        TCanvas * c_zmass = new TCanvas(Form("c_zmass%d%s",i,decay.Data()),"c_zmass",500,500);
        
        double ymax = h1->GetMaximum(); 
        if(ymax< h2->GetMaximum()) ymax=h2->GetMaximum();
 
        ymax=ymax*1.5;
        h1->SetMaximum(ymax);   
        
        h1->SetLineColor(2);
        h2->SetLineColor(4);
        h1->Draw();
        h2->SetMarkerStyle(20);
        h2->SetMarkerSize(1);
        h2->Draw("esame");
        
        TLegend *l = new TLegend(0.68,0.76,0.89,0.87);
        l->AddEntry(h1,"TTbar","L");
        l->AddEntry(h2,"DATA","LP");
        SetLegend(l);
        c_zmass->Print(Form("plots/%s%s.eps",decay.Data(),hh[i].Data() ));
    }

}
void ploter(int ii)
{

    TFile* f1 = new TFile(Form("result_MuMu_TTbarFullLepMGDecays.root") );
    TFile* f2 = new TFile(Form("result_MuMu_Run2012MuMu.root") );
    TFile* f12 = new TFile(Form("result_ElEl_TTbarFullLepMGDecays.root") );
    TFile* f22 = new TFile(Form("result_ElEl_Run2012ElEl.root") );
    TFile* f13 = new TFile(Form("result_MuEl_TTbarFullLepMGDecays.root") );
    TFile* f23 = new TFile(Form("result_MuEl_Run2012MuEl.root") );

    TString hh[16] = {"htM","htbarM","httM","hweight","htM2","htbarM2","httM2","hweight2","hlnuM","hlnuM2","hlnubarM","hlnubarM2","hCSVD11","hCSVD12","hCSVD21","hCSVD22" };
    TString yx[16] = {"M_{t} (GeV/c^{2})", "M_{#bar{t}} (GeV/c^{2})","M_{t#bar{t}} (GeV/c^{2})","weight",
                     "M_{t} (GeV/c^{2})", "M_{#bar{t}} (GeV/c^{2})","M_{t#bar{t}} (GeV/c^{2})","weight",
                     "M_{l#nu} (GeV/c^{2})","M_{l#nu} (GeV/c^{2})", "M_{l#bar{#nu}} (GeV/c^{2})","M_{l#bar{#nu}} (GeV/c^{2})",
                     "b-Discriminator (CSV)","b-Discriminator (CSV)","b-Discriminator (CSV)","b-Discriminator (CSV)"  };

    for(int i=0;i<16;i++)
    {
        TH1F* h1 = (TH1*)f1->Get(Form("%s",hh[i].Data()));
        TH1F* h2 = (TH1*)f2->Get(Form("%s",hh[i].Data()));
        TH1F* h12 = (TH1*)f12->Get(Form("%s",hh[i].Data()));
        TH1F* h22 = (TH1*)f22->Get(Form("%s",hh[i].Data()));
        TH1F* h13 = (TH1*)f13->Get(Form("%s",hh[i].Data()));
        TH1F* h23 = (TH1*)f23->Get(Form("%s",hh[i].Data()));

        h2->Sumw2(); h22->Sumw2(); h23->Sumw2();
        h1->Add(h12);  h1->Add(h13);
        h2->Add(h22);  h2->Add(h23);

        int nBin = h1->GetNbinsX();
        h1->AddBinContent(nBin,h1->GetBinContent(nBin+1));
        h2->AddBinContent(nBin,h2->GetBinContent(nBin+1));

        h1->Scale(1/h1->Integral() );
        h2->Scale(1/h2->Integral() );
        h1->GetYaxis()->SetTitle("Normalized Entries");
        h1->GetXaxis()->SetTitle(yx[i].Data());

        TCanvas * c_zmass = new TCanvas(Form("c_zmass%dall",i),"c_zmass",500,500);
        
        double ymax = h1->GetMaximum(); 
        if(ymax< h2->GetMaximum()) ymax=h2->GetMaximum();
 
        ymax=ymax*1.5;
        h1->SetMaximum(ymax);   
        
        h1->SetLineColor(2);
        h2->SetLineColor(4);
        h1->Draw();
        h2->SetMarkerStyle(20);
        h2->SetMarkerSize(1);
        h2->Draw("esame");
        
        TLegend *l = new TLegend(0.68,0.76,0.89,0.87);
        l->AddEntry(h1,"TTbar","L");
        l->AddEntry(h2,"DATA","LP");
        SetLegend(l);
        c_zmass->Print(Form("plots/all%s.eps",hh[i].Data() ));
    }

}


void SetLegend(TLegend* l){
  l->SetTextSize(0.04);
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->Draw();
}
