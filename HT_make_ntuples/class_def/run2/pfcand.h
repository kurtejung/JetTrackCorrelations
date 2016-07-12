//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  7 22:33:41 2016 by ROOT version 6.02/13
// from TTree pfTree/dijet tree
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root
//////////////////////////////////////////////////////////

#ifndef pfcand_h
#define pfcand_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class pfcand {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nPFpart;
   vector<int>     *pfId;
   vector<float>   *pfPt;
   vector<float>   *pfEnergy;
   vector<float>   *pfVsPtInitial;
   vector<float>   *pfEta;
   vector<float>   *pfPhi;
   Float_t         vn[5][15];
   Float_t         psin[5][15];
   Float_t         sumpt[15];

   // List of branches
   TBranch        *b_nPFpart;   //!
   TBranch        *b_pfId;   //!
   TBranch        *b_pfPt;   //!
   TBranch        *b_pfEnergy;   //!
   TBranch        *b_pfVsPtInitial;   //!
   TBranch        *b_pfEta;   //!
   TBranch        *b_pfPhi;   //!
   TBranch        *b_vn;   //!
   TBranch        *b_vpsi;   //!
   TBranch        *b_sumpt;   //!

   pfcand(TTree *tree=0, bool ispp=0);
   virtual ~pfcand();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, bool ispp);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

pfcand::pfcand(TTree *tree, bool ispp) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root:/pfcandAnalyzer");
      dir->GetObject("pfTree",tree);

   }
   Init(tree, ispp);
}

pfcand::~pfcand()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pfcand::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pfcand::LoadTree(Long64_t entry)
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

void pfcand::Init(TTree *tree, bool ispp)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pfId = 0;
   pfPt = 0;
   pfEnergy = 0;
   pfVsPtInitial = 0;
   pfEta = 0;
   pfPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nPFpart", &nPFpart, &b_nPFpart);
   fChain->SetBranchAddress("pfId", &pfId, &b_pfId);
   fChain->SetBranchAddress("pfPt", &pfPt, &b_pfPt);
   fChain->SetBranchAddress("pfEnergy", &pfEnergy, &b_pfEnergy);
   fChain->SetBranchAddress("pfVsPtInitial", &pfVsPtInitial, &b_pfVsPtInitial);
   fChain->SetBranchAddress("pfEta", &pfEta, &b_pfEta);
   fChain->SetBranchAddress("pfPhi", &pfPhi, &b_pfPhi);
   fChain->SetBranchAddress("vn", vn, &b_vn);
   fChain->SetBranchAddress("psin", psin, &b_vpsi);
   fChain->SetBranchAddress("sumpt", sumpt, &b_sumpt);

   ispp = 0; //dummy variable for now - might be used in the future
   Notify();
}

Bool_t pfcand::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pfcand::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pfcand::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pfcand_cxx
