//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  7 22:34:16 2016 by ROOT version 6.02/13
// from TTree HltTree/
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root
//////////////////////////////////////////////////////////

#ifndef Skim_h
#define Skim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Skim {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ana_step;
   Int_t           pHBHENoiseFilterResultProducer;
   Int_t           HBHENoiseFilterResult;
   Int_t           HBHENoiseFilterResultRun1;
   Int_t           HBHENoiseFilterResultRun2Loose;
   Int_t           HBHENoiseFilterResultRun2Tight;
   Int_t           HBHEIsoNoiseFilterResult;
   Int_t           pPAprimaryVertexFilter;
   Int_t           pBeamScrapingFilter;
   Int_t           pVertexFilterCutG;
   Int_t           pVertexFilterCutGloose;
   Int_t           pVertexFilterCutGtight;
   Int_t           pVertexFilterCutGplus;
   Int_t           pVertexFilterCutE;
   Int_t           pVertexFilterCutEandG;

   Int_t           pcollisionEventSelection;
   Int_t           pprimaryVertexFilter;
   Int_t           phfCoincFilter1;
   Int_t           phfCoincFilter2;
   Int_t           phfCoincFilter3;
   Int_t           phfCoincFilter4;
   Int_t           phfCoincFilter5;
   Int_t           pclusterCompatibilityFilter;

   // List of branches
   TBranch        *b_ana_step;   //!
   TBranch        *b_pHBHENoiseFilterResultProducer;   //!
   TBranch        *b_HBHENoiseFilterResult;   //!
   TBranch        *b_HBHENoiseFilterResultRun1;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Loose;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Tight;   //!
   TBranch        *b_HBHEIsoNoiseFilterResult;   //!
   TBranch        *b_pPAprimaryVertexFilter;   //!
   TBranch        *b_pBeamScrapingFilter;   //!
   TBranch        *b_pVertexFilterCutG;   //!
   TBranch        *b_pVertexFilterCutGloose;   //!
   TBranch        *b_pVertexFilterCutGtight;   //!
   TBranch        *b_pVertexFilterCutGplus;   //!
   TBranch        *b_pVertexFilterCutE;   //!
   TBranch        *b_pVertexFilterCutEandG;   //!

   TBranch        *b_pcollisionEventSelection;   //!
   TBranch        *b_pprimaryVertexFilter;   //!
   TBranch        *b_phfCoincFilter1;   //!
   TBranch        *b_phfCoincFilter2;   //!
   TBranch        *b_phfCoincFilter3;   //!
   TBranch        *b_phfCoincFilter4;   //!
   TBranch        *b_phfCoincFilter5;   //!
   TBranch        *b_pclusterCompatibilityFilter;   //!


   Skim(TTree *tree=0, bool ispp);
   virtual ~Skim();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, bool ispp);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

Skim::Skim(TTree *tree, bool ispp) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root:/skimanalysis");
      dir->GetObject("HltTree",tree);

   }
   Init(tree);
}

Skim::~Skim()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Skim::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Skim::LoadTree(Long64_t entry)
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

void Skim::Init(TTree *tree, bool ispp)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
   fChain->SetBranchAddress("pHBHENoiseFilterResultProducer", &pHBHENoiseFilterResultProducer, &b_pHBHENoiseFilterResultProducer);
   fChain->SetBranchAddress("HBHENoiseFilterResult", &HBHENoiseFilterResult, &b_HBHENoiseFilterResult);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun1", &HBHENoiseFilterResultRun1, &b_HBHENoiseFilterResultRun1);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose, &b_HBHENoiseFilterResultRun2Loose);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &HBHENoiseFilterResultRun2Tight, &b_HBHENoiseFilterResultRun2Tight);
   fChain->SetBranchAddress("HBHEIsoNoiseFilterResult", &HBHEIsoNoiseFilterResult, &b_HBHEIsoNoiseFilterResult);
   if(ispp){
      fChain->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter, &b_pPAprimaryVertexFilter);
      fChain->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter, &b_pBeamScrapingFilter);
      fChain->SetBranchAddress("pVertexFilterCutG", &pVertexFilterCutG, &b_pVertexFilterCutG);
      fChain->SetBranchAddress("pVertexFilterCutGloose", &pVertexFilterCutGloose, &b_pVertexFilterCutGloose);
      fChain->SetBranchAddress("pVertexFilterCutGtight", &pVertexFilterCutGtight, &b_pVertexFilterCutGtight);
      fChain->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus, &b_pVertexFilterCutGplus);
      fChain->SetBranchAddress("pVertexFilterCutE", &pVertexFilterCutE, &b_pVertexFilterCutE);
      fChain->SetBranchAddress("pVertexFilterCutEandG", &pVertexFilterCutEandG, &b_pVertexFilterCutEandG);
   }
   else{
      fChain->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection, &b_pcollisionEventSelection);
      fChain->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter, &b_pprimaryVertexFilter);
      fChain->SetBranchAddress("phfCoincFilter1", &phfCoincFilter1, &b_phfCoincFilter1);
      fChain->SetBranchAddress("phfCoincFilter2", &phfCoincFilter2, &b_phfCoincFilter2);
      fChain->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3, &b_phfCoincFilter3);
      fChain->SetBranchAddress("phfCoincFilter4", &phfCoincFilter4, &b_phfCoincFilter4);
      fChain->SetBranchAddress("phfCoincFilter5", &phfCoincFilter5, &b_phfCoincFilter5);
      fChain->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter, &b_pclusterCompatibilityFilter);
   }
   Notify();
}

Bool_t Skim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Skim::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Skim::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Skim_cxx
