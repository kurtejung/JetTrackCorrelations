//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 14 12:58:37 2014 by ROOT version 5.14/00e
// from TTree trackTree/v1
// found on file: /mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root
//////////////////////////////////////////////////////////

#ifndef Tracks_h
#define Tracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Tracks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           nEv;
   Int_t           nLumi;
   Int_t           nBX;
   Int_t           nRun;
   Int_t           N;
   Int_t           nv;
   Float_t         vx[2];   //[nv]
   Float_t         vy[2];   //[nv]
   Float_t         vz[2];   //[nv]
   Float_t         vxErr[2];   //[nv]
   Float_t         vyErr[2];   //[nv]
   Float_t         vzErr[2];   //[nv]
   Int_t           nDaugher[2];   //[nv]
   Int_t           nVtx;
   Int_t           maxVtx;
   Int_t           nTrkVtx[1];   //[nVtx]
   Float_t         normChi2Vtx[1];   //[nVtx]
   Float_t         sumPtVtx[1];   //[nVtx]
   Float_t         xVtx[1];   //[nVtx]
   Float_t         yVtx[1];   //[nVtx]
   Float_t         zVtx[1];   //[nVtx]
   Float_t         xVtxErr[1];   //[nVtx]
   Float_t         yVtxErr[1];   //[nVtx]
   Float_t         zVtxErr[1];   //[nVtx]
   Float_t         vtxDist2D[1];   //[nVtx]
   Float_t         vtxDist2DErr[1];   //[nVtx]
   Float_t         vtxDist2DSig[1];   //[nVtx]
   //Float_t         vtxDist2D[1];   //[nVtx]
   Float_t         vtxDist3DErr[1];   //[nVtx]
   Float_t         vtxDist3DSig[1];   //[nVtx]
   Int_t           nTrk;
   Float_t         trkPt[12739];   //[nTrk]
   Float_t         trkPtError[12739];   //[nTrk]
   Int_t           trkNHit[12739];   //[nTrk]
   Int_t           trkNlayer[12739];   //[nTrk]
   Float_t         trkEta[12739];   //[nTrk]
   Float_t         trkPhi[12739];   //[nTrk]
   Float_t         pEta[12739];   //[nTrk]
   Float_t         pPhi[12739];   //[nTrk]
   Float_t         pPt[12739];   //[nTrk]
   Int_t           trkCharge[12739];   //[nTrk]
   Int_t           trkVtxIndex[12739];   //[nTrk]
   Bool_t          highPurity[12739];   //[nTrk]
   Bool_t          highPuritySetWithPV[12739];   //[nTrk]
   Float_t         trkChi2[12739];   //[nTrk]
   Float_t         trkNdof[12739];   //[nTrk]
   Float_t         trkDxy1[12739];   //[nTrk]
   Float_t         trkDxyError1[12739];   //[nTrk]
   Float_t         trkDz1[12739];   //[nTrk]
   Float_t         trkDzError1[12739];   //[nTrk]
   Bool_t          trkFake[12739];   //[nTrk]
   Float_t         trkAlgo[12739];   //[nTrk]

   // List of branches
   TBranch        *b_nEv;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_N;   //!
   TBranch        *b_nv;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vxErr;   //!
   TBranch        *b_vyErr;   //!
   TBranch        *b_vzErr;   //!
   TBranch        *b_nDaugher;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_maxVtx;   //!
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
   //TBranch        *b_vtxDist2D;   //!
   TBranch        *b_vtxDist3DErr;   //!
   TBranch        *b_vtxDist3DSig;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkNHit;   //!
   TBranch        *b_trkNlayer;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_pEta;   //!
   TBranch        *b_pPhi;   //!
   TBranch        *b_pPt;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkVtxIndex;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_highPuritySetWithPV;   //!
   TBranch        *b_trkChi2;   //!
   TBranch        *b_trkNdof;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
   TBranch        *b_trkFake;   //!
   TBranch        *b_trkAlgo;   //!

   Tracks(TTree *tree=0);
   virtual ~Tracks();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Tracks_cxx
Tracks::Tracks(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
 /*     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root");
      if (!f) {
         f = new TFile("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root");
      }
      tree = (TTree*)gDirectory->Get("trackTree");
*/
   }
   Init(tree);
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
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Tracks::Init(TTree *tree)
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

   fChain->SetBranchAddress("nEv", &nEv, &b_nEv);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("N", &N, &b_N);
   fChain->SetBranchAddress("nv", &nv, &b_nv);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("vxErr", vxErr, &b_vxErr);
   fChain->SetBranchAddress("vyErr", vyErr, &b_vyErr);
   fChain->SetBranchAddress("vzErr", vzErr, &b_vzErr);
   fChain->SetBranchAddress("nDaugher", nDaugher, &b_nDaugher);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("maxVtx", &maxVtx, &b_maxVtx);
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
   fChain->SetBranchAddress("vtxDist2D", vtxDist2D, &b_vtxDist2D);
   fChain->SetBranchAddress("vtxDist3DErr", vtxDist3DErr, &b_vtxDist3DErr);
   fChain->SetBranchAddress("vtxDist3DSig", vtxDist3DSig, &b_vtxDist3DSig);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
   fChain->SetBranchAddress("trkNHit", trkNHit, &b_trkNHit);
   fChain->SetBranchAddress("trkNlayer", trkNlayer, &b_trkNlayer);
   fChain->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("pEta", pEta, &b_pEta);
   fChain->SetBranchAddress("pPhi", pPhi, &b_pPhi);
   fChain->SetBranchAddress("pPt", pPt, &b_pPt);
   fChain->SetBranchAddress("trkCharge", trkCharge, &b_trkCharge);
   fChain->SetBranchAddress("trkVtxIndex", trkVtxIndex, &b_trkVtxIndex);
   fChain->SetBranchAddress("highPurity", highPurity, &b_highPurity);
   fChain->SetBranchAddress("highPuritySetWithPV", highPuritySetWithPV, &b_highPuritySetWithPV);
   fChain->SetBranchAddress("trkChi2", trkChi2, &b_trkChi2);
   fChain->SetBranchAddress("trkNdof", trkNdof, &b_trkNdof);
   fChain->SetBranchAddress("trkDxy1", trkDxy1, &b_trkDxy1);
   fChain->SetBranchAddress("trkDxyError1", trkDxyError1, &b_trkDxyError1);
   fChain->SetBranchAddress("trkDz1", trkDz1, &b_trkDz1);
   fChain->SetBranchAddress("trkDzError1", trkDzError1, &b_trkDzError1);
   fChain->SetBranchAddress("trkFake", trkFake, &b_trkFake);
   fChain->SetBranchAddress("trkAlgo", trkAlgo, &b_trkAlgo);
   Notify();
}

Bool_t Tracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
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
