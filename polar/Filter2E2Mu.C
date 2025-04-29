#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TSystem.h>
#include "ExRootAnalysis/ExRootAnalysis/ExRootClasses.h"

using namespace std;

void Filter2E2Mu(string rootin = "for_ran_jhugen_0p3m.root", string rootout = "")
{
  gSystem->Load("ExRootAnalysis/libExRootAnalysis.so");
  if(rootout == "") rootout = rootin.substr(0, rootin.size() - 5) + "_2e2mu.root";

  TFile *filein = new TFile(rootin.c_str());
  TTree *LHEF = (TTree *)filein->Get("LHEF");
  TFile *fileout = new TFile(rootout.c_str(), "RECREATE");
  TTree *LHEF_new = LHEF->CloneTree(0);

  TClonesArray *Particle = NULL;
  LHEF->SetBranchAddress("Particle", &Particle);
  for(Long64_t i = 0; LHEF->GetEntry(i); ++i) {
    Int_t np = Particle->GetEntries();
    Int_t cE = 0, cMu = 0, cTau = 0;
    for(int ip = 0; ip < np; ++ip) {
      TRootLHEFParticle *particle = (TRootLHEFParticle *)Particle->UncheckedAt(ip);
      if(abs(particle->PID) == 11) ++cE;
      if(abs(particle->PID) == 13) ++cMu;
      if(abs(particle->PID) == 15) ++cTau;
    }
    if(cE != 2 || cMu != 2 || cTau != 0) continue;
    LHEF_new->Fill();
  }
  fileout->cd();
  LHEF_new->Write(NULL, TObject::kOverwrite);

  fileout->Close();
  filein->Close();
}
