//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  7 22:34:32 2016 by ROOT version 6.02/13
// from TTree trackTree/v1
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root
//////////////////////////////////////////////////////////

#ifndef Tracks_h
#define Tracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Tracks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nEv;
   Int_t           nLumi;
   Int_t           nBX;
   Int_t           nRun;
   Int_t           N;
   Int_t           nVtx;
   Int_t           nTrk;
   Int_t           maxPtVtx;
   Int_t           maxMultVtx;
   Int_t           nTrkVtx[9];   //[nVtx]
   Float_t         normChi2Vtx[9];   //[nVtx]
   Float_t         sumPtVtx[9];   //[nVtx]
   Float_t         xVtx[9];   //[nVtx]
   Float_t         yVtx[9];   //[nVtx]
   Float_t         zVtx[9];   //[nVtx]
   Float_t         xVtxErr[9];   //[nVtx]
   Float_t         yVtxErr[9];   //[nVtx]
   Float_t         zVtxErr[9];   //[nVtx]
   Float_t         vtxDist2D[9];   //[nVtx]
   Float_t         vtxDist2DErr[9];   //[nVtx]
   Float_t         vtxDist2DSig[9];   //[nVtx]
   Float_t         vtxDist3D[9];   //[nVtx]
   Float_t         vtxDist3DErr[9];   //[nVtx]
   Float_t         vtxDist3DSig[9];   //[nVtx]
   Int_t           nVtxSim;
   Float_t         xVtxSim[1];   //[nVtxSim]
   Float_t         yVtxSim[1];   //[nVtxSim]
   Float_t         zVtxSim[1];   //[nVtxSim]
   Float_t         trkPt[20000];   //[nTrk]
   Float_t         trkPtError[20000];   //[nTrk]
   UChar_t         trkNHit[20000];   //[nTrk]
   UChar_t         trkNlayer[20000];   //[nTrk]
   Float_t         trkEta[20000];   //[nTrk]
   Float_t         trkPhi[20000];   //[nTrk]
   Int_t           trkCharge[20000];   //[nTrk]
   UChar_t         trkNVtx[20000];   //[nTrk]
   Int_t           nTrkTimesnVtx;
   Bool_t          trkAssocVtx[1470];   //[nTrkTimesnVtx]
   Float_t         trkDxyOverDxyError[1470];   //[nTrkTimesnVtx]
   Float_t         trkDzOverDzError[1470];   //[nTrkTimesnVtx]
   Bool_t          highPurity[20000];   //[nTrk]
   Bool_t          tight[20000];   //[nTrk]
   Bool_t          loose[20000];   //[nTrk]
   Float_t         trkChi2[20000];   //[nTrk]
   UChar_t         trkNdof[20000];   //[nTrk]
   Float_t         trkDxy1[20000];   //[nTrk]
   Float_t         trkDxyError1[20000];   //[nTrk]
   Float_t         trkDz1[20000];   //[nTrk]
   Float_t         trkDzError1[20000];   //[nTrk]
   Bool_t          trkFake[20000];   //[nTrk]
   UChar_t         trkAlgo[20000];   //[nTrk]
   UChar_t         trkOriginalAlgo[20000];   //[nTrk]
   Float_t         trkMVA[20000];   //[nTrk]
   Bool_t          trkMVALoose[20000];   //[nTrk]
   Bool_t          trkMVATight[20000];   //[nTrk]
   Int_t           pfType[20000];   //[nTrk]
   Float_t         pfCandPt[20000];   //[nTrk]
   Float_t         pfEcal[20000];   //[nTrk]
   Float_t         pfHcal[20000];   //[nTrk]

   // List of branches
   TBranch        *b_nEv;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_N;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_maxPtVtx;   //!
   TBranch        *b_maxMultVtx;   //!
   TBranch        *b_nTrkVtx;   //!
   TBranch        *b_normChi2Vtx;   //!
   TBranch        *b_sumPtVtx;   //!
   TBranch        *b_xVtx;   //!
   TBranch        *b_yVtx;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_xVtxErr;   //!
   TBranch        *b_yVtxErr;   //!
   TBranch        *b_zVtxErr;   //!
   TBranch        *b_vtxDist2D;   //!
   TBranch        *b_vtxDist2DErr;   //!
   TBranch        *b_vtxDist2DSig;   //!
   TBranch        *b_vtxDist3D;   //!
   TBranch        *b_vtxDist3DErr;   //!
   TBranch        *b_vtxDist3DSig;   //!
   TBranch        *b_nVtxSim;   //!
   TBranch        *b_xVtxSim;   //!
   TBranch        *b_yVtxSim;   //!
   TBranch        *b_zVtxSim;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkNHit;   //!
   TBranch        *b_trkNlayer;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkNVtx;   //!
   TBranch        *b_nTrkTimesnVtx;   //!
   TBranch        *b_trkAssocVtx;   //!
   TBranch        *b_trkDxyOverDxyError;   //!
   TBranch        *b_trkDzOverDzError;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_tight;   //!
   TBranch        *b_loose;   //!
   TBranch        *b_trkChi2;   //!
   TBranch        *b_trkNdof;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
   TBranch        *b_trkFake;   //!
   TBranch        *b_trkAlgo;   //!
   TBranch        *b_trkOriginalAlgo;   //!
   TBranch        *b_trkMVA;   //!
   TBranch        *b_trkMVALoose;   //!
   TBranch        *b_trkMVATight;   //!
   TBranch        *b_pfType;   //!
   TBranch        *b_pfCandPt;   //!
   TBranch        *b_pfEcal;   //!
   TBranch        *b_pfHcal;   //!

   Tracks(TTree *tree=0, bool ispp=0);
   virtual ~Tracks();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, bool ispp);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

Tracks::Tracks(TTree *tree, bool ispp) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root:/ppTrack");
      dir->GetObject("trackTree",tree);

   }
   Init(tree, ispp);
}

Tracks::~Tracks()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Tracks::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Tracks::LoadTree(Long64_t entry)
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

void Tracks::Init(TTree *tree, bool ispp)
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

   fChain->SetBranchAddress("nEv", &nEv, &b_nEv);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("N", &N, &b_N);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("maxPtVtx", &maxPtVtx, &b_maxPtVtx);
   fChain->SetBranchAddress("maxMultVtx", &maxMultVtx, &b_maxMultVtx);
   fChain->SetBranchAddress("nTrkVtx", nTrkVtx, &b_nTrkVtx);
   fChain->SetBranchAddress("normChi2Vtx", normChi2Vtx, &b_normChi2Vtx);
   fChain->SetBranchAddress("sumPtVtx", sumPtVtx, &b_sumPtVtx);
   fChain->SetBranchAddress("xVtx", xVtx, &b_xVtx);
   fChain->SetBranchAddress("yVtx", yVtx, &b_yVtx);
   fChain->SetBranchAddress("zVtx", zVtx, &b_zVtx);
   fChain->SetBranchAddress("xVtxErr", xVtxErr, &b_xVtxErr);
   fChain->SetBranchAddress("yVtxErr", yVtxErr, &b_yVtxErr);
   fChain->SetBranchAddress("zVtxErr", zVtxErr, &b_zVtxErr);
   fChain->SetBranchAddress("vtxDist2D", vtxDist2D, &b_vtxDist2D);
   fChain->SetBranchAddress("vtxDist2DErr", vtxDist2DErr, &b_vtxDist2DErr);
   fChain->SetBranchAddress("vtxDist2DSig", vtxDist2DSig, &b_vtxDist2DSig);
   fChain->SetBranchAddress("vtxDist3D", vtxDist3D, &b_vtxDist3D);
   fChain->SetBranchAddress("vtxDist3DErr", vtxDist3DErr, &b_vtxDist3DErr);
   fChain->SetBranchAddress("vtxDist3DSig", vtxDist3DSig, &b_vtxDist3DSig);
   fChain->SetBranchAddress("nVtxSim", &nVtxSim, &b_nVtxSim);
   fChain->SetBranchAddress("xVtxSim", &xVtxSim, &b_xVtxSim);
   fChain->SetBranchAddress("yVtxSim", &yVtxSim, &b_yVtxSim);
   fChain->SetBranchAddress("zVtxSim", &zVtxSim, &b_zVtxSim);
   fChain->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
   fChain->SetBranchAddress("trkNHit", trkNHit, &b_trkNHit);
   fChain->SetBranchAddress("trkNlayer", trkNlayer, &b_trkNlayer);
   fChain->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("trkCharge", trkCharge, &b_trkCharge);
   fChain->SetBranchAddress("trkNVtx", trkNVtx, &b_trkNVtx);
   fChain->SetBranchAddress("nTrkTimesnVtx", &nTrkTimesnVtx, &b_nTrkTimesnVtx);
   fChain->SetBranchAddress("trkAssocVtx", trkAssocVtx, &b_trkAssocVtx);
   fChain->SetBranchAddress("trkDxyOverDxyError", trkDxyOverDxyError, &b_trkDxyOverDxyError);
   fChain->SetBranchAddress("trkDzOverDzError", trkDzOverDzError, &b_trkDzOverDzError);
   fChain->SetBranchAddress("highPurity", highPurity, &b_highPurity);
   fChain->SetBranchAddress("tight", tight, &b_tight);
   fChain->SetBranchAddress("loose", loose, &b_loose);
   fChain->SetBranchAddress("trkChi2", trkChi2, &b_trkChi2);
   fChain->SetBranchAddress("trkNdof", trkNdof, &b_trkNdof);
   fChain->SetBranchAddress("trkDxy1", trkDxy1, &b_trkDxy1);
   fChain->SetBranchAddress("trkDxyError1", trkDxyError1, &b_trkDxyError1);
   fChain->SetBranchAddress("trkDz1", trkDz1, &b_trkDz1);
   fChain->SetBranchAddress("trkDzError1", trkDzError1, &b_trkDzError1);
   fChain->SetBranchAddress("trkFake", trkFake, &b_trkFake);
   fChain->SetBranchAddress("trkAlgo", trkAlgo, &b_trkAlgo);
   fChain->SetBranchAddress("trkOriginalAlgo", trkOriginalAlgo, &b_trkOriginalAlgo);
   fChain->SetBranchAddress("trkMVA", trkMVA, &b_trkMVA);
   if(ispp) fChain->SetBranchAddress("trkMVALoose", trkMVALoose, &b_trkMVALoose);
   fChain->SetBranchAddress("trkMVATight", trkMVATight, &b_trkMVATight);
   fChain->SetBranchAddress("pfType", pfType, &b_pfType);
   fChain->SetBranchAddress("pfCandPt", pfCandPt, &b_pfCandPt);
   fChain->SetBranchAddress("pfEcal", pfEcal, &b_pfEcal);
   fChain->SetBranchAddress("pfHcal", pfHcal, &b_pfHcal);
   Notify();
}

Bool_t Tracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Tracks::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Tracks::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Tracks_cxx
