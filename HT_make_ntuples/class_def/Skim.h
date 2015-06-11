//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 14 12:58:23 2014 by ROOT version 5.14/00e
// from TTree HltTree/
// found on file: /mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root
//////////////////////////////////////////////////////////

#ifndef Skim_h
#define Skim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Skim {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           pcollisionEventSelection;
   Int_t           pHBHENoiseFilter;
   Int_t           phfCoincFilter;
   Int_t           phfCoincFilter3;
   Int_t           pprimaryVertexFilter;
   Int_t           phltPixelClusterShapeFilter;
   Int_t           phiEcalRecHitSpikeFilter;
   Int_t           ana_step;
   Int_t           pPAcollisionEventSelectionPA;

   // List of branches
   TBranch        *b_pcollisionEventSelection;   //!
   TBranch        *b_pPAcollisionEventSelectionPA;
   TBranch        *b_pHBHENoiseFilter;   //!
   TBranch        *b_phfCoincFilter;   //!
   TBranch        *b_phfCoincFilter3;   //!
   TBranch        *b_pprimaryVertexFilter;   //!
   TBranch        *b_phltPixelClusterShapeFilter;   //!
   TBranch        *b_phiEcalRecHitSpikeFilter;   //!
   TBranch        *b_ana_step;   //!

   Skim(TTree *tree=0);
   virtual ~Skim();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Skim_cxx
Skim::Skim(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
  /*    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root");
      if (!f) {
         f = new TFile("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root");
      }
      tree = (TTree*)gDirectory->Get("HltTree");
*/
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
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Skim::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection, &b_pcollisionEventSelection);
   fChain->SetBranchAddress("pPAcollisionEventSelectionPA", & pPAcollisionEventSelectionPA, &b_pPAcollisionEventSelectionPA);
   fChain->SetBranchAddress("pHBHENoiseFilter", &pHBHENoiseFilter, &b_pHBHENoiseFilter);
   fChain->SetBranchAddress("phfCoincFilter", &phfCoincFilter, &b_phfCoincFilter);
   fChain->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3, &b_phfCoincFilter3);
   fChain->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter, &b_pprimaryVertexFilter);
   fChain->SetBranchAddress("phltPixelClusterShapeFilter", &phltPixelClusterShapeFilter, &b_phltPixelClusterShapeFilter);
   fChain->SetBranchAddress("phiEcalRecHitSpikeFilter", &phiEcalRecHitSpikeFilter, &b_phiEcalRecHitSpikeFilter);
   fChain->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
   Notify();
}

Bool_t Skim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
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
