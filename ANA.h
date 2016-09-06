//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 25 10:40:39 2013 by ROOT version 5.32/00
// from TTree tree/Tree for Top quark study
// found on file: /afs/cern.ch/work/y/youngjo/public/For8Tev/v20130321_V00-00-07/vallot_TTbarFullLepMGDecays.root
//////////////////////////////////////////////////////////

#ifndef ANA_h
#define ANA_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class ANA {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          EVENT;
   UInt_t          RUN;
   UInt_t          LUMI;
   UInt_t          npileup;
   UInt_t          nvertex;
   Double_t        puweight;
   Double_t        puweightplus;
   Double_t        puweightminus;
   Double_t        bweight30CSVT;
   Double_t        bweight30CSVTup;
   Double_t        bweight30CSVTdw;
   Double_t        bweight30CSVTuplight;
   Double_t        bweight30CSVTdwlight;
   Double_t        ZMass;
   Double_t        genZMass;
   UInt_t          ZtauDecay;
   Double_t        isIso;
   Double_t        PairSign;
   Double_t        lep1_relIso04;
   Double_t        lep2_relIso04;
   Double_t        lep1_chIso03;
   Double_t        lep2_chIso03;
   Double_t        lep1_nhIso03;
   Double_t        lep2_nhIso03;
   Double_t        lep1_phIso03;
   Double_t        lep2_phIso03;
   Double_t        lep1_relIso03;
   Double_t        lep2_relIso03;
   Double_t        lep1_relIso03db;
   Double_t        lep2_relIso03db;
   Double_t        lep1_relIso03rho;
   Double_t        lep2_relIso03rho;
   Double_t        lep1_pt;
   Double_t        lep2_pt;
   Double_t        lep1_eta;
   Double_t        lep2_eta;
   Double_t        lep1_phi;
   Double_t        lep2_phi;
   Double_t        lep1_charge;
   Double_t        lep2_charge;
   Double_t        lepweight;
   vector<double>  *jets_secvtxmass;
   vector<double>  *jets_pt;
   vector<double>  *jets_eta;
   vector<double>  *jets_phi;
   vector<int>     *jets_flavor;
   vector<int>     *jets_fromtop;
   vector<double>  *jets_bDiscriminatorJP;
   vector<double>  *jets_bDiscriminatorCSV;
   vector<double>  *jets_bDisCSVweight;
   vector<int>     *csvd_jetid;
   UInt_t          nJet30;
   UInt_t          nJet30Up;
   UInt_t          nJet30Dw;
   UInt_t          nGenJet20;
   UInt_t          nGenbJet20;
   UInt_t          nGencJet20;
   Double_t        genLep1_pt;
   Double_t        genLep2_pt;
   Double_t        genLep1_eta;
   Double_t        genLep2_eta;
   UInt_t          ttbarGen_dileptonic;
   UInt_t          visible;
   UInt_t          ttbb;
   UInt_t          ttcc;
   UInt_t          ttLF;
   UInt_t          nbjets30_CSVL;
   UInt_t          nbjets30_CSVM;
   UInt_t          nbjets30_CSVT;
   Double_t        MET;
   Double_t        metphi;
   Double_t        METUp;
   Double_t        METDw;

   // List of branches
   TBranch        *b_EVENT;   //!
   TBranch        *b_RUN;   //!
   TBranch        *b_LUMI;   //!
   TBranch        *b_npileup;   //!
   TBranch        *b_nvertex;   //!
   TBranch        *b_puweight;   //!
   TBranch        *b_puweightplus;   //!
   TBranch        *b_puweightminus;   //!
   TBranch        *b_bweight30CSVT;   //!
   TBranch        *b_bweight30CSVTup;   //!
   TBranch        *b_bweight30CSVTdw;   //!
   TBranch        *b_bweight30CSVTuplight;   //!
   TBranch        *b_bweight30CSVTdwlight;   //!
   TBranch        *b_ZMass;   //!
   TBranch        *b_ZMass2;   //!
   TBranch        *b_ZtauDecay;   //!
   TBranch        *b_isIso;   //!
   TBranch        *b_PairSign;   //!
   TBranch        *b_lep1_relIso04;   //!
   TBranch        *b_lep2_relIso04;   //!
   TBranch        *b_lep1_chIso03;   //!
   TBranch        *b_lep2_chIso03;   //!
   TBranch        *b_lep1_nhIso03;   //!
   TBranch        *b_lep2_nhIso03;   //!
   TBranch        *b_lep1_phIso03;   //!
   TBranch        *b_lep2_phIso03;   //!
   TBranch        *b_lep1_relIso03;   //!
   TBranch        *b_lep2_relIso03;   //!
   TBranch        *b_lep1_relIso03db;   //!
   TBranch        *b_lep2_relIso03db;   //!
   TBranch        *b_lep1_relIso03rho;   //!
   TBranch        *b_lep2_relIso03rho;   //!
   TBranch        *b_lep1_pt;   //!
   TBranch        *b_lep2_pt;   //!
   TBranch        *b_lep1_eta;   //!
   TBranch        *b_lep2_eta;   //!
   TBranch        *b_lep1_phi;   //!
   TBranch        *b_lep2_phi;   //!
   TBranch        *b_lep1_charge;   //!
   TBranch        *b_lep2_charge;   //!
   TBranch        *b_lepweight;   //!
   TBranch        *b_jets_secvtxmass;   //!
   TBranch        *b_jets_pt;   //!
   TBranch        *b_jets_eta;   //!
   TBranch        *b_jets_phi;   //!
   TBranch        *b_jets_flavor;   //!
   TBranch        *b_jets_fromtop;   //!
   TBranch        *b_jets_bDiscriminatorJP;   //!
   TBranch        *b_jets_bDiscriminatorCSV;   //!
   TBranch        *b_jets_bDisCSVweight;   //!
   TBranch        *b_csvd_jetid;   //!
   TBranch        *b_nJet30;   //!
   TBranch        *b_nJet30Up;   //!
   TBranch        *b_nJet30Dw;   //!
   TBranch        *b_nGenJet20;   //!
   TBranch        *b_nGenbJet20;   //!
   TBranch        *b_nGencJet20;   //!
   TBranch        *b_genLep1_pt;   //!
   TBranch        *b_genLep2_pt;   //!
   TBranch        *b_genLep1_eta;   //!
   TBranch        *b_genLep2_eta;   //!
   TBranch        *b_ttbarGen_dileptonic;   //!
   TBranch        *b_visible;   //!
   TBranch        *b_ttbb;   //!
   TBranch        *b_ttcc;   //!
   TBranch        *b_ttLF;   //!
   TBranch        *b_nbjets30_CSVL;   //!
   TBranch        *b_nbjets30_CSVM;   //!
   TBranch        *b_nbjets30_CSVT;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_METUp;   //!
   TBranch        *b_METDw;   //!

   ANA(TTree *tree=0);
   virtual ~ANA();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ANA_cxx
ANA::ANA(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/work/y/youngjo/public/For8Tev/v20130321_V00-00-07/vallot_TTbarFullLepMGDecays.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/work/y/youngjo/public/For8Tev/v20130321_V00-00-07/vallot_TTbarFullLepMGDecays.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/afs/cern.ch/work/y/youngjo/public/For8Tev/v20130321_V00-00-07/vallot_TTbarFullLepMGDecays.root:/MuMu");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

ANA::~ANA()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ANA::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ANA::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ANA::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jets_secvtxmass = 0;
   jets_pt = 0;
   jets_eta = 0;
   jets_phi = 0;
   jets_flavor = 0;
   jets_fromtop = 0;
   jets_bDiscriminatorJP = 0;
   jets_bDiscriminatorCSV = 0;
   jets_bDisCSVweight = 0;
   csvd_jetid = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EVENT", &EVENT, &b_EVENT);
   fChain->SetBranchAddress("RUN", &RUN, &b_RUN);
   fChain->SetBranchAddress("LUMI", &LUMI, &b_LUMI);
   fChain->SetBranchAddress("npileup", &npileup, &b_npileup);
   fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   fChain->SetBranchAddress("puweightplus", &puweightplus, &b_puweightplus);
   fChain->SetBranchAddress("puweightminus", &puweightminus, &b_puweightminus);
   fChain->SetBranchAddress("bweight30CSVT", &bweight30CSVT, &b_bweight30CSVT);
   fChain->SetBranchAddress("bweight30CSVTup", &bweight30CSVTup, &b_bweight30CSVTup);
   fChain->SetBranchAddress("bweight30CSVTdw", &bweight30CSVTdw, &b_bweight30CSVTdw);
   fChain->SetBranchAddress("bweight30CSVTuplight", &bweight30CSVTuplight, &b_bweight30CSVTuplight);
   fChain->SetBranchAddress("bweight30CSVTdwlight", &bweight30CSVTdwlight, &b_bweight30CSVTdwlight);
   fChain->SetBranchAddress("ZMass", &ZMass, &b_ZMass);
   fChain->SetBranchAddress("genZMass", &genZMass, &b_ZMass2);
   fChain->SetBranchAddress("ZtauDecay", &ZtauDecay, &b_ZtauDecay);
   fChain->SetBranchAddress("isIso", &isIso, &b_isIso);
   fChain->SetBranchAddress("PairSign", &PairSign, &b_PairSign);
   fChain->SetBranchAddress("lep1_relIso04", &lep1_relIso04, &b_lep1_relIso04);
   fChain->SetBranchAddress("lep2_relIso04", &lep2_relIso04, &b_lep2_relIso04);
   fChain->SetBranchAddress("lep1_chIso03", &lep1_chIso03, &b_lep1_chIso03);
   fChain->SetBranchAddress("lep2_chIso03", &lep2_chIso03, &b_lep2_chIso03);
   fChain->SetBranchAddress("lep1_nhIso03", &lep1_nhIso03, &b_lep1_nhIso03);
   fChain->SetBranchAddress("lep2_nhIso03", &lep2_nhIso03, &b_lep2_nhIso03);
   fChain->SetBranchAddress("lep1_phIso03", &lep1_phIso03, &b_lep1_phIso03);
   fChain->SetBranchAddress("lep2_phIso03", &lep2_phIso03, &b_lep2_phIso03);
   fChain->SetBranchAddress("lep1_relIso03", &lep1_relIso03, &b_lep1_relIso03);
   fChain->SetBranchAddress("lep2_relIso03", &lep2_relIso03, &b_lep2_relIso03);
   fChain->SetBranchAddress("lep1_relIso03db", &lep1_relIso03db, &b_lep1_relIso03db);
   fChain->SetBranchAddress("lep2_relIso03db", &lep2_relIso03db, &b_lep2_relIso03db);
   fChain->SetBranchAddress("lep1_relIso03rho", &lep1_relIso03rho, &b_lep1_relIso03rho);
   fChain->SetBranchAddress("lep2_relIso03rho", &lep2_relIso03rho, &b_lep2_relIso03rho);
   fChain->SetBranchAddress("lep1_pt", &lep1_pt, &b_lep1_pt);
   fChain->SetBranchAddress("lep2_pt", &lep2_pt, &b_lep2_pt);
   fChain->SetBranchAddress("lep1_eta", &lep1_eta, &b_lep1_eta);
   fChain->SetBranchAddress("lep2_eta", &lep2_eta, &b_lep2_eta);
   fChain->SetBranchAddress("lep1_phi", &lep1_phi, &b_lep1_phi);
   fChain->SetBranchAddress("lep2_phi", &lep2_phi, &b_lep2_phi);
   fChain->SetBranchAddress("lep1_charge", &lep1_charge, &b_lep1_charge);
   fChain->SetBranchAddress("lep2_charge", &lep2_charge, &b_lep2_charge);
   fChain->SetBranchAddress("lepweight", &lepweight, &b_lepweight);
   fChain->SetBranchAddress("jets_secvtxmass", &jets_secvtxmass, &b_jets_secvtxmass);
   fChain->SetBranchAddress("jets_pt", &jets_pt, &b_jets_pt);
   fChain->SetBranchAddress("jets_eta", &jets_eta, &b_jets_eta);
   fChain->SetBranchAddress("jets_phi", &jets_phi, &b_jets_phi);
   fChain->SetBranchAddress("jets_flavor", &jets_flavor, &b_jets_flavor);
   fChain->SetBranchAddress("jets_fromtop", &jets_fromtop, &b_jets_fromtop);
   fChain->SetBranchAddress("jets_bDiscriminatorJP", &jets_bDiscriminatorJP, &b_jets_bDiscriminatorJP);
   fChain->SetBranchAddress("jets_bDiscriminatorCSV", &jets_bDiscriminatorCSV, &b_jets_bDiscriminatorCSV);
   fChain->SetBranchAddress("jets_bDisCSVweight", &jets_bDisCSVweight, &b_jets_bDisCSVweight);
   fChain->SetBranchAddress("csvd_jetid", &csvd_jetid, &b_csvd_jetid);
   fChain->SetBranchAddress("nJet30", &nJet30, &b_nJet30);
   fChain->SetBranchAddress("nJet30Up", &nJet30Up, &b_nJet30Up);
   fChain->SetBranchAddress("nJet30Dw", &nJet30Dw, &b_nJet30Dw);
   fChain->SetBranchAddress("nGenJet20", &nGenJet20, &b_nGenJet20);
   fChain->SetBranchAddress("nGenbJet20", &nGenbJet20, &b_nGenbJet20);
   fChain->SetBranchAddress("nGencJet20", &nGencJet20, &b_nGencJet20);
   fChain->SetBranchAddress("genLep1_pt", &genLep1_pt, &b_genLep1_pt);
   fChain->SetBranchAddress("genLep2_pt", &genLep2_pt, &b_genLep2_pt);
   fChain->SetBranchAddress("genLep1_eta", &genLep1_eta, &b_genLep1_eta);
   fChain->SetBranchAddress("genLep2_eta", &genLep2_eta, &b_genLep2_eta);
   fChain->SetBranchAddress("ttbarGen_dileptonic", &ttbarGen_dileptonic, &b_ttbarGen_dileptonic);
   fChain->SetBranchAddress("visible", &visible, &b_visible);
   fChain->SetBranchAddress("ttbb", &ttbb, &b_ttbb);
   fChain->SetBranchAddress("ttcc", &ttcc, &b_ttcc);
   fChain->SetBranchAddress("ttLF", &ttLF, &b_ttLF);
   fChain->SetBranchAddress("nbjets30_CSVL", &nbjets30_CSVL, &b_nbjets30_CSVL);
   fChain->SetBranchAddress("nbjets30_CSVM", &nbjets30_CSVM, &b_nbjets30_CSVM);
   fChain->SetBranchAddress("nbjets30_CSVT", &nbjets30_CSVT, &b_nbjets30_CSVT);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("METUp", &METUp, &b_METUp);
   fChain->SetBranchAddress("METDw", &METDw, &b_METDw);
   Notify();
}

Bool_t ANA::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ANA::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ANA::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ANA_cxx
