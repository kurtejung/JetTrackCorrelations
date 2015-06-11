//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb 15 09:55:16 2014 by ROOT version 5.14/00e
// from TTree hi/Tree of Hi gen Event
// found on file: /mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet7_STARTHI53_LV1/merged/HiForest_Pythia_Hydjet_Jet80_Track8_Jet7_STARTHI53_LV1_merged_forest_4.root
//////////////////////////////////////////////////////////

#ifndef GenParticles_h
#define GenParticles_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class GenParticles {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           event;
   Float_t         b;
   Float_t         npart;
   Float_t         ncoll;
   Float_t         nhard;
   Float_t         phi0;
   Float_t         scale;
   Int_t           n[3];
   Float_t         ptav[3];
   Int_t           mult;
   Float_t         pt[49303];   //[mult]
   Float_t         eta[49303];   //[mult]
   Float_t         phi[49303];   //[mult]
   Int_t           pdg[49303];   //[mult]
   Int_t           chg[49303];   //[mult]
   Int_t           sube[49303];   //[mult]
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         vr;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_b;   //!
   TBranch        *b_npart;   //!
   TBranch        *b_ncoll;   //!
   TBranch        *b_nhard;   //!
   TBranch        *b_phi0;   //!
   TBranch        *b_scale;   //!
   TBranch        *b_n;   //!
   TBranch        *b_ptav;   //!
   TBranch        *b_mult;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_chg;   //!
   TBranch        *b_sube;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vr;   //!

   GenParticles(TTree *tree=0);
   virtual ~GenParticles();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GenParticles_cxx
GenParticles::GenParticles(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet7_STARTHI53_LV1/merged/HiForest_Pythia_Hydjet_Jet80_Track8_Jet7_STARTHI53_LV1_merged_forest_4.root");
      if (!f) {
         f = new TFile("/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet7_STARTHI53_LV1/merged/HiForest_Pythia_Hydjet_Jet80_Track8_Jet7_STARTHI53_LV1_merged_forest_4.root");
      }
      tree = (TTree*)gDirectory->Get("hi");

   }
   Init(tree);
}

GenParticles::~GenParticles()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GenParticles::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GenParticles::LoadTree(Long64_t entry)
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

void GenParticles::Init(TTree *tree)
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

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("ncoll", &ncoll, &b_ncoll);
   fChain->SetBranchAddress("nhard", &nhard, &b_nhard);
   fChain->SetBranchAddress("phi0", &phi0, &b_phi0);
   fChain->SetBranchAddress("scale", &scale, &b_scale);
   fChain->SetBranchAddress("n", n, &b_n);
   fChain->SetBranchAddress("ptav", ptav, &b_ptav);
   fChain->SetBranchAddress("mult", &mult, &b_mult);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("pdg", pdg, &b_pdg);
   fChain->SetBranchAddress("chg", chg, &b_chg);
   fChain->SetBranchAddress("sube", sube, &b_sube);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("vr", &vr, &b_vr);
   Notify();
}

Bool_t GenParticles::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GenParticles::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GenParticles::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GenParticles_cxx
