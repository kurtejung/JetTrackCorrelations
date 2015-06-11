//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 14 12:57:43 2014 by ROOT version 5.14/00e
// from TTree t/akVs3CalopatJets Jet Analysis Tree
// found on file: /mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root
//////////////////////////////////////////////////////////

#ifndef JetAna_h
#define JetAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class JetAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           evt;
   Float_t         b;
   Int_t           nref;
   Int_t           ngen;
   Float_t         pthat;  
   Float_t         rawpt[206];   //[nref]
 
   Float_t         jtpt[206];   //[nref]
   Float_t         jteta[206];   //[nref]
   Float_t         jty[206];   //[nref]
   Float_t         jtphi[206];   //[nref]
   Float_t         genpt[206];   //[nref]
   Float_t         geneta[206];   //[nref]
   Float_t         genphi[206];   //[nref]
   Float_t         jtpu[206];   //[nref]
   Float_t         jtm[206];   //[nref]
   Float_t         discr_fr01[206];   //[nref]
   Float_t         trackMax[206];   //[nref]
   Float_t         trackSum[206];   //[nref]
   Int_t           trackN[206];   //[nref]
   Float_t         trackHardSum[206];   //[nref]
   Int_t           trackHardN[206];   //[nref]
   Float_t         chargedMax[206];   //[nref]
   Float_t         chargedSum[206];   //[nref]
   Int_t           chargedN[206];   //[nref]
   Float_t         chargedHardSum[206];   //[nref]
   Int_t           chargedHardN[206];   //[nref]
   Float_t         photonMax[206];   //[nref]
   Float_t         photonSum[206];   //[nref]
   Int_t           photonN[206];   //[nref]
   Float_t         photonHardSum[206];   //[nref]
   Int_t           photonHardN[206];   //[nref]
   Float_t         neutralMax[206];   //[nref]
   Float_t         neutralSum[206];   //[nref]
   Int_t           neutralN[206];   //[nref]
   Float_t         hcalSum[206];   //[nref]
   Float_t         ecalSum[206];   //[nref]
   Float_t         eMax[206];   //[nref]
   Float_t         eSum[206];   //[nref]
   Int_t           eN[206];   //[nref]
   Float_t         muMax[206];   //[nref]
   Float_t         muSum[206];   //[nref]
   Int_t           muN[206];   //[nref]

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_b;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_jtpu;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_discr_fr01;   //!
   TBranch        *b_trackMax;   //!
   TBranch        *b_trackSum;   //!
   TBranch        *b_trackN;   //!
   TBranch        *b_trackHardSum;   //!
   TBranch        *b_trackHardN;   //!
   TBranch        *b_chargedMax;   //!
   TBranch        *b_chargedSum;   //!
   TBranch        *b_chargedN;   //!
   TBranch        *b_chargedHardSum;   //!
   TBranch        *b_chargedHardN;   //!
   TBranch        *b_photonMax;   //!
   TBranch        *b_photonSum;   //!
   TBranch        *b_photonN;   //!
   TBranch        *b_photonHardSum;   //!
   TBranch        *b_photonHardN;   //!
   TBranch        *b_neutralMax;   //!
   TBranch        *b_neutralSum;   //!
   TBranch        *b_neutralN;   //!
   TBranch        *b_hcalSum;   //!
   TBranch        *b_ecalSum;   //!
   TBranch        *b_eMax;   //!
   TBranch        *b_eSum;   //!
   TBranch        *b_eN;   //!
   TBranch        *b_muMax;   //!
   TBranch        *b_muSum;   //!
   TBranch        *b_muN;   //!

   JetAna(TTree *tree=0);
   virtual ~JetAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JetAna_cxx
JetAna::JetAna(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
  /*    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root");
      if (!f) {
         f = new TFile("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root");
      }
      tree = (TTree*)gDirectory->Get("t");
*/
   }
   Init(tree);
}

JetAna::~JetAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetAna::LoadTree(Long64_t entry)
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

void JetAna::Init(TTree *tree)
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

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("nref", &nref, &b_nref);
   fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
   fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
   fChain->SetBranchAddress("pthat",&pthat, &b_pthat);
   fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", jteta, &b_jteta);
   fChain->SetBranchAddress("jty", jty, &b_jty);
   fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
   fChain->SetBranchAddress("genpt", genpt, &b_genpt);
   fChain->SetBranchAddress("geneta", geneta, &b_geneta);
   fChain->SetBranchAddress("genphi", genphi, &b_genphi);
   fChain->SetBranchAddress("jtpu", jtpu, &b_jtpu);
   fChain->SetBranchAddress("jtm", jtm, &b_jtm);
   fChain->SetBranchAddress("discr_fr01", discr_fr01, &b_discr_fr01);
   fChain->SetBranchAddress("trackMax", trackMax, &b_trackMax);
   fChain->SetBranchAddress("trackSum", trackSum, &b_trackSum);
   fChain->SetBranchAddress("trackN", trackN, &b_trackN);
   fChain->SetBranchAddress("trackHardSum", trackHardSum, &b_trackHardSum);
   fChain->SetBranchAddress("trackHardN", trackHardN, &b_trackHardN);
   fChain->SetBranchAddress("chargedMax", chargedMax, &b_chargedMax);
   fChain->SetBranchAddress("chargedSum", chargedSum, &b_chargedSum);
   fChain->SetBranchAddress("chargedN", chargedN, &b_chargedN);
   fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
   fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
   fChain->SetBranchAddress("photonMax", photonMax, &b_photonMax);
   fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
   fChain->SetBranchAddress("photonN", photonN, &b_photonN);
   fChain->SetBranchAddress("photonHardSum", photonHardSum, &b_photonHardSum);
   fChain->SetBranchAddress("photonHardN", photonHardN, &b_photonHardN);
   fChain->SetBranchAddress("neutralMax", neutralMax, &b_neutralMax);
   fChain->SetBranchAddress("neutralSum", neutralSum, &b_neutralSum);
   fChain->SetBranchAddress("neutralN", neutralN, &b_neutralN);
   fChain->SetBranchAddress("hcalSum", hcalSum, &b_hcalSum);
   fChain->SetBranchAddress("ecalSum", ecalSum, &b_ecalSum);
   fChain->SetBranchAddress("eMax", eMax, &b_eMax);
   fChain->SetBranchAddress("eSum", eSum, &b_eSum);
   fChain->SetBranchAddress("eN", eN, &b_eN);
   fChain->SetBranchAddress("muMax", muMax, &b_muMax);
   fChain->SetBranchAddress("muSum", muSum, &b_muSum);
   fChain->SetBranchAddress("muN", muN, &b_muN);
   Notify();
}

Bool_t JetAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JetAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JetAna_cxx
