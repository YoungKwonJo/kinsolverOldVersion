void run()
{

// if(i==0)
 {
  run("MuMu",0);
  run("MuMu",1);
 }
// if(i==1)
 {
  run("ElEl",0);
  run("ElEl",1);
 }
// if(i==2)
 {
  run("MuEl",0);
  run("MuEl",1);
 }
}
void run(TString decay,int i){

   gROOT->ProcessLine(".L TtFullLepKinSolver.C+g");
   gROOT->ProcessLine(".L ANA.C+g");

   TString location = "/data/CMSDATA/v20130321_V00-00-07/";
   TString channel[] = {"TTbarFullLepMGDecays",Form("Run2012%s",decay.Data())};


   TFile *f1 =  new TFile(Form("%svallot_%s.root",location.Data(),channel[i].Data() ));

   TTree *atree1= dynamic_cast<TTree *>(f1->Get( Form("%s/tree",decay.Data())  ));
   TFile fout(Form("result_%s_%s.root",decay.Data(),channel[i].Data() ), "RECREATE");

   ANA t(atree1);
   t.Loop();

   fout.Write();
   fout.Close();


}
