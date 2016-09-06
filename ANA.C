#define ANA_cxx
#include "ANA.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
//#include "TMath.h"
#include "TtFullLepKinSolver.h"

void ANA::Loop()
{

    double tmassbegin_= 100;
    double tmassend_  = 300;
    double tmassstep_ = 1;

    double myints[] = {30.7137,56.2880,23.0744,59.1015,24.9145};
    std::vector<double> nupars_ (myints, myints + sizeof(myints) / sizeof(int) );

    TtFullLepKinSolver* solver = new TtFullLepKinSolver(tmassbegin_, tmassend_, tmassstep_, nupars_);

    TH1F *htM, *htbarM, *httM, *hlnuM, *hlnubarM, *hweight, *hCSVD11, *hCSVD12;         

    htM       = new TH1F(Form("htM"),      "t_{M}  ",      70, 0, 350    );
    htbarM    = new TH1F(Form("htbarM"),   "#bar{t}_{M}  ", 70, 0, 350    );
    httM      = new TH1F(Form("httM"),     "t#bar{t}_{M}  ",     80, 0, 800    );
    hlnuM     = new TH1F(Form("hlnuM"),    "lnu_{M}  ",    30, 0, 150    );
    hlnubarM  = new TH1F(Form("hlnubarM"), "lnubar_{M}  ", 30, 0, 150    );
    hweight   = new TH1F(Form("hweight"),  "weight  ",     52, 0, 1.2    );

    hCSVD11    =  new TH1F(Form("hCSVD11"),  "CSVD  ",     10, 0, 1    );
    hCSVD12    =  new TH1F(Form("hCSVD12"),  "CSVD  ",     10, 0, 1    );

    TH1F *htM2, *htbarM2, *httM2, *hlnuM2, *hlnubarM2, *hweight2, *hCSVD21, *hCSVD22;         
    htM2       = new TH1F(Form("htM2"),      "t_{M}  ",      70, 0, 350    );
    htbarM2    = new TH1F(Form("htbarM2"),   "#bar{t}_{M}  ", 70, 0, 350    );
    httM2      = new TH1F(Form("httM2"),     "t#bar{t}_{M}  ",     80, 0, 800    );
    hlnuM2     = new TH1F(Form("hlnuM2"),    "lnu_{M}  ",    30, 0, 150    );
    hlnubarM2  = new TH1F(Form("hlnubarM2"), "lnubar_{M}  ", 30, 0, 150    );
    hweight2   = new TH1F(Form("hweight2"),  "weight  ",     52, 0, 1.2    );

    hCSVD21    =  new TH1F(Form("hCSVD21"),  "CSVD  ",     10, 0, 1    );
    hCSVD22    =  new TH1F(Form("hCSVD22"),  "CSVD  ",     10, 0, 1    );

   if (fChain == 0) return;

   TLorentzVector lep1, lep2, met, jet1, jet2, nu1, nu2;
   TLorentzVector nu11[30][30], nu12[30][30], nu21[30][30], nu22[30][30];

   TLorentzVector jets[30];
   double weight1[30][30], weight2[30][30];
   
   //int iiii=0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
///////////

     if(ZMass > 12 && lep1_relIso03 < 0.15 && lep2_relIso03 < 0.15 && PairSign < 0 
             && fabs(ZMass - 91.2) > 15 && MET > 30 && nJet30 >= 2 ) ;
     else  continue;

     lep1.SetPtEtaPhiM(lep1_pt,lep1_eta,lep1_phi,0);
     lep2.SetPtEtaPhiM(lep2_pt,lep2_eta,lep2_phi,0);
     int jid1 = csvd_jetid->at(0);
     int jid2 = csvd_jetid->at(1);
     //cout << "jet id : "<< jid1 << ", " << jid2 << " : " << nJet30 << endl;

     met.SetPtEtaPhiM(MET,0,metphi,0);

     //jet1.SetPtEtaPhiM(jets_pt->at(jid1),jets_eta->at(jid1),jets_phi->at(jid1),jets_secvtxmass->at(jid1) );
     //jet2.SetPtEtaPhiM(jets_pt->at(jid2),jets_eta->at(jid2),jets_phi->at(jid2),jets_secvtxmass->at(jid2) );
     for(int i=0;i < (int)nJet30;i++) jets[i].SetPtEtaPhiM(jets_pt->at(i),jets_eta->at(i),jets_phi->at(i),jets_secvtxmass->at(i) );

     int j1id=-1, j2id=-1;
     double weight=-1;
     double csvd1=-1, csvd2=-1; 

     for(int i=0;i<(int)nJet30-1;i++) for(int j=i;j<(int)nJet30;j++)
     {
          //selecttion jet of CSVM 
          if(jets_bDiscriminatorCSV->at(i) >0.679 && jets_bDiscriminatorCSV->at(j) >0.679 ) ;
          else continue; 

          double xconstraint = lep1.Px()+lep2.Px()+ jets[i].Px() + jets[i].Px() +met.Px();
          double yconstraint = lep1.Py()+lep2.Py()+ jets[j].Py() + jets[j].Py() +met.Py();
        
          solver->SetConstraints(xconstraint, yconstraint);
          TtFullLepKinSolver::NeutrinoSolution nuSol= solver->getNuSolution( lep1, lep2 , jets[i], jets[j]);
          weight1[i][j] = nuSol.weight;
          nu11[i][j] = nuSol.neutrino;
          nu12[i][j] = nuSol.neutrinoBar;
        
          TtFullLepKinSolver::NeutrinoSolution nuSol2= solver->getNuSolution( lep1, lep2 , jets[j], jets[i]);
          weight2[j][i] = nuSol2.weight;
          nu21[j][i] = nuSol2.neutrino;
          nu22[j][i] = nuSol2.neutrinoBar;     

          if(weight1[i][j] > weight2[j][i])
          {
              if(weight<weight1[i][j]) 
              {
                 weight=weight1[i][j];
                 j1id=i; j2id=j;
                 jet1=jets[i]; jet2=jets[j];
                 nu1=nu11[i][j]; nu2=nu12[i][j];
                 csvd1=jets_bDiscriminatorCSV->at(i);
                 csvd2=jets_bDiscriminatorCSV->at(j);
              }
          }
          else
          {          
              if(weight<weight2[j][i])
              {
                 weight=weight2[j][i]; 
                 j1id=j; j2id=i;
                 jet1=jets[j]; jet2=jets[i];
                 nu1=nu21[j][i]; nu2=nu22[j][i];
                 csvd1=jets_bDiscriminatorCSV->at(j);
                 csvd2=jets_bDiscriminatorCSV->at(i);
              }
          }
     }

     if(weight < 0) continue;
     //if(iiii>1000) break;
     //iiii++;

      htM      ->Fill( (lep1+jet1+nu1).M() ,puweight);
      htbarM   ->Fill( (lep2+jet2+nu2).M() ,puweight);
      httM     ->Fill( (lep1+jet1+nu1+lep2+jet2+nu2).M() ,puweight);
      hlnuM    ->Fill( (lep1+nu1).M() ,puweight);
      hlnubarM ->Fill( (lep2+nu2).M() ,puweight);
      hweight  ->Fill( weight         ,puweight);
     
 
      hCSVD11  ->Fill( csvd1         ,puweight);
      hCSVD12  ->Fill( csvd2         ,puweight);

    
      if(nbjets30_CSVM>=2)
      {
           htM2      ->Fill( (lep1+jet1+nu1).M() ,puweight);
           htbarM2   ->Fill( (lep2+jet2+nu2).M() ,puweight);
           httM2     ->Fill( (lep1+jet1+nu1+lep2+jet2+nu2).M() ,puweight);
           hlnuM2    ->Fill( (lep1+nu1).M() ,puweight);
           hlnubarM2 ->Fill( (lep2+nu2).M() ,puweight);
           hweight2  ->Fill( weight         ,puweight);

          hCSVD21  ->Fill( csvd1         ,puweight);
          hCSVD22  ->Fill( csvd2         ,puweight);


      }  

////////////
   }
}
