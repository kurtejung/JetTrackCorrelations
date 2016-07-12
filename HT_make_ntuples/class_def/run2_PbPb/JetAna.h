//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 11 14:24:57 2016 by ROOT version 5.34/20
// from TTree t/akPu3CalopatJetsWithBtagging Jet Analysis Tree
// found on file: HiForestAOD_1128.root
//////////////////////////////////////////////////////////

#ifndef JetAna_h
#define JetAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class JetAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           evt;
   Float_t         b;
   Int_t           nref;
   Float_t         rawpt[16];   //[nref]
   Float_t         jtpt[16];   //[nref]
   Float_t         jteta[16];   //[nref]
   Float_t         jty[16];   //[nref]
   Float_t         jtphi[16];   //[nref]
   Float_t         jtpu[16];   //[nref]
   Float_t         jtm[16];   //[nref]
   Float_t         jtarea[16];   //[nref]
   Float_t         discr_fr01[16];   //[nref]
   Float_t         trackMax[16];   //[nref]
   Float_t         trackSum[16];   //[nref]
   Int_t           trackN[16];   //[nref]
   Float_t         trackHardSum[16];   //[nref]
   Int_t           trackHardN[16];   //[nref]
   Float_t         chargedMax[16];   //[nref]
   Float_t         chargedSum[16];   //[nref]
   Int_t           chargedN[16];   //[nref]
   Float_t         chargedHardSum[16];   //[nref]
   Int_t           chargedHardN[16];   //[nref]
   Float_t         photonMax[16];   //[nref]
   Float_t         photonSum[16];   //[nref]
   Int_t           photonN[16];   //[nref]
   Float_t         photonHardSum[16];   //[nref]
   Int_t           photonHardN[16];   //[nref]
   Float_t         neutralMax[16];   //[nref]
   Float_t         neutralSum[16];   //[nref]
   Int_t           neutralN[16];   //[nref]
   Float_t         hcalSum[16];   //[nref]
   Float_t         ecalSum[16];   //[nref]
   Float_t         eMax[16];   //[nref]
   Float_t         eSum[16];   //[nref]
   Int_t           eN[16];   //[nref]
   Float_t         muMax[16];   //[nref]
   Float_t         muSum[16];   //[nref]
   Int_t           muN[16];   //[nref]
   Float_t         discr_ssvHighEff[16];   //[nref]
   Float_t         discr_ssvHighPur[16];   //[nref]
   Float_t         discr_csvMva[16];   //[nref]
   Float_t         discr_csvSimple[16];   //[nref]
   Float_t         discr_muByIp3[16];   //[nref]
   Float_t         discr_muByPt[16];   //[nref]
   Float_t         discr_prob[16];   //[nref]
   Float_t         discr_probb[16];   //[nref]
   Float_t         discr_tcHighEff[16];   //[nref]
   Float_t         discr_tcHighPur[16];   //[nref]
   Float_t         ndiscr_ssvHighEff[16];   //[nref]
   Float_t         ndiscr_ssvHighPur[16];   //[nref]
   Float_t         ndiscr_csvSimple[16];   //[nref]
   Float_t         ndiscr_muByPt[16];   //[nref]
   Float_t         ndiscr_prob[16];   //[nref]
   Float_t         ndiscr_probb[16];   //[nref]
   Float_t         ndiscr_tcHighEff[16];   //[nref]
   Float_t         ndiscr_tcHighPur[16];   //[nref]
   Float_t         pdiscr_csvSimple[16];   //[nref]
   Float_t         pdiscr_prob[16];   //[nref]
   Float_t         pdiscr_probb[16];   //[nref]
   Int_t           nsvtx[16];   //[nref]
   Int_t           svtxntrk[16];   //[nref]
   Float_t         svtxdl[16];   //[nref]
   Float_t         svtxdls[16];   //[nref]
   Float_t         svtxdl2d[16];   //[nref]
   Float_t         svtxdls2d[16];   //[nref]
   Float_t         svtxm[16];   //[nref]
   Float_t         svtxpt[16];   //[nref]
   Float_t         svtxmcorr[16];   //[nref]
   Int_t           nIPtrk[16];   //[nref]
   Int_t           nselIPtrk[16];   //[nref]
   Float_t         mue[16];   //[nref]
   Float_t         mupt[16];   //[nref]
   Float_t         mueta[16];   //[nref]
   Float_t         muphi[16];   //[nref]
   Float_t         mudr[16];   //[nref]
   Float_t         muptrel[16];   //[nref]
   Int_t           muchg[16];   //[nref]

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_b;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtpu;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_jtarea;   //!
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
   TBranch        *b_discr_ssvHighEff;   //!
   TBranch        *b_discr_ssvHighPur;   //!
   TBranch        *b_discr_csvMva;   //!
   TBranch        *b_discr_csvSimple;   //!
   TBranch        *b_discr_muByIp3;   //!
   TBranch        *b_discr_muByPt;   //!
   TBranch        *b_discr_prob;   //!
   TBranch        *b_discr_probb;   //!
   TBranch        *b_discr_tcHighEff;   //!
   TBranch        *b_discr_tcHighPur;   //!
   TBranch        *b_ndiscr_ssvHighEff;   //!
   TBranch        *b_ndiscr_ssvHighPur;   //!
   TBranch        *b_ndiscr_csvSimple;   //!
   TBranch        *b_ndiscr_muByPt;   //!
   TBranch        *b_ndiscr_prob;   //!
   TBranch        *b_ndiscr_probb;   //!
   TBranch        *b_ndiscr_tcHighEff;   //!
   TBranch        *b_ndiscr_tcHighPur;   //!
   TBranch        *b_pdiscr_csvSimple;   //!
   TBranch        *b_pdiscr_prob;   //!
   TBranch        *b_pdiscr_probb;   //!
   TBranch        *b_nsvtx;   //!
   TBranch        *b_svtxntrk;   //!
   TBranch        *b_svtxdl;   //!
   TBranch        *b_svtxdls;   //!
   TBranch        *b_svtxdl2d;   //!
   TBranch        *b_svtxdls2d;   //!
   TBranch        *b_svtxm;   //!
   TBranch        *b_svtxpt;   //!
   TBranch        *b_svtxmcorr;   //!
   TBranch        *b_nIPtrk;   //!
   TBranch        *b_nselIPtrk;   //!
   TBranch        *b_mue;   //!
   TBranch        *b_mupt;   //!
   TBranch        *b_mueta;   //!
   TBranch        *b_muphi;   //!
   TBranch        *b_mudr;   //!
   TBranch        *b_muptrel;   //!
   TBranch        *b_muchg;   //!

   JetAna(TTree *tree=0);
   virtual ~JetAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

JetAna::JetAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("HiForestAOD_1128.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("HiForestAOD_1128.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("HiForestAOD_1128.root:/akPu3CaloJetAnalyzer");
      dir->GetObject("t",tree);

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
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JetAna::Init(TTree *tree)
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

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("nref", &nref, &b_nref);
   fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
   fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", jteta, &b_jteta);
   fChain->SetBranchAddress("jty", jty, &b_jty);
   fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
   fChain->SetBranchAddress("jtpu", jtpu, &b_jtpu);
   fChain->SetBranchAddress("jtm", jtm, &b_jtm);
   fChain->SetBranchAddress("jtarea", jtarea, &b_jtarea);
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
   fChain->SetBranchAddress("discr_ssvHighEff", discr_ssvHighEff, &b_discr_ssvHighEff);
   fChain->SetBranchAddress("discr_ssvHighPur", discr_ssvHighPur, &b_discr_ssvHighPur);
   fChain->SetBranchAddress("discr_csvMva", discr_csvMva, &b_discr_csvMva);
   fChain->SetBranchAddress("discr_csvSimple", discr_csvSimple, &b_discr_csvSimple);
   fChain->SetBranchAddress("discr_muByIp3", discr_muByIp3, &b_discr_muByIp3);
   fChain->SetBranchAddress("discr_muByPt", discr_muByPt, &b_discr_muByPt);
   fChain->SetBranchAddress("discr_prob", discr_prob, &b_discr_prob);
   fChain->SetBranchAddress("discr_probb", discr_probb, &b_discr_probb);
   fChain->SetBranchAddress("discr_tcHighEff", discr_tcHighEff, &b_discr_tcHighEff);
   fChain->SetBranchAddress("discr_tcHighPur", discr_tcHighPur, &b_discr_tcHighPur);
   fChain->SetBranchAddress("ndiscr_ssvHighEff", ndiscr_ssvHighEff, &b_ndiscr_ssvHighEff);
   fChain->SetBranchAddress("ndiscr_ssvHighPur", ndiscr_ssvHighPur, &b_ndiscr_ssvHighPur);
   fChain->SetBranchAddress("ndiscr_csvSimple", ndiscr_csvSimple, &b_ndiscr_csvSimple);
   fChain->SetBranchAddress("ndiscr_muByPt", ndiscr_muByPt, &b_ndiscr_muByPt);
   fChain->SetBranchAddress("ndiscr_prob", ndiscr_prob, &b_ndiscr_prob);
   fChain->SetBranchAddress("ndiscr_probb", ndiscr_probb, &b_ndiscr_probb);
   fChain->SetBranchAddress("ndiscr_tcHighEff", ndiscr_tcHighEff, &b_ndiscr_tcHighEff);
   fChain->SetBranchAddress("ndiscr_tcHighPur", ndiscr_tcHighPur, &b_ndiscr_tcHighPur);
   fChain->SetBranchAddress("pdiscr_csvSimple", pdiscr_csvSimple, &b_pdiscr_csvSimple);
   fChain->SetBranchAddress("pdiscr_prob", pdiscr_prob, &b_pdiscr_prob);
   fChain->SetBranchAddress("pdiscr_probb", pdiscr_probb, &b_pdiscr_probb);
   fChain->SetBranchAddress("nsvtx", nsvtx, &b_nsvtx);
   fChain->SetBranchAddress("svtxntrk", svtxntrk, &b_svtxntrk);
   fChain->SetBranchAddress("svtxdl", svtxdl, &b_svtxdl);
   fChain->SetBranchAddress("svtxdls", svtxdls, &b_svtxdls);
   fChain->SetBranchAddress("svtxdl2d", svtxdl2d, &b_svtxdl2d);
   fChain->SetBranchAddress("svtxdls2d", svtxdls2d, &b_svtxdls2d);
   fChain->SetBranchAddress("svtxm", svtxm, &b_svtxm);
   fChain->SetBranchAddress("svtxpt", svtxpt, &b_svtxpt);
   fChain->SetBranchAddress("svtxmcorr", svtxmcorr, &b_svtxmcorr);
   fChain->SetBranchAddress("nIPtrk", nIPtrk, &b_nIPtrk);
   fChain->SetBranchAddress("nselIPtrk", nselIPtrk, &b_nselIPtrk);
   fChain->SetBranchAddress("mue", mue, &b_mue);
   fChain->SetBranchAddress("mupt", mupt, &b_mupt);
   fChain->SetBranchAddress("mueta", mueta, &b_mueta);
   fChain->SetBranchAddress("muphi", muphi, &b_muphi);
   fChain->SetBranchAddress("mudr", mudr, &b_mudr);
   fChain->SetBranchAddress("muptrel", muptrel, &b_muptrel);
   fChain->SetBranchAddress("muchg", muchg, &b_muchg);
   Notify();
}

Bool_t JetAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
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
