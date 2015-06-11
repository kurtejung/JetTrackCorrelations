//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 14 12:58:51 2014 by ROOT version 5.14/00e
// from TTree HltTree/
// found on file: /mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_08Mar2014_0300CET_Track8_Jet21/0.root
//////////////////////////////////////////////////////////

#ifndef HLT_h
#define HLT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class HLT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           NL1IsolEm;
   Float_t         L1IsolEmEt[4];   //[NL1IsolEm]
   Float_t         L1IsolEmE[4];   //[NL1IsolEm]
   Float_t         L1IsolEmEta[4];   //[NL1IsolEm]
   Float_t         L1IsolEmPhi[4];   //[NL1IsolEm]
   Int_t           NL1NIsolEm;
   Float_t         L1NIsolEmEt[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmE[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmEta[4];   //[NL1NIsolEm]
   Float_t         L1NIsolEmPhi[4];   //[NL1NIsolEm]
   Int_t           NL1Mu;
   Float_t         L1MuPt[4];   //[NL1Mu]
   Float_t         L1MuE[4];   //[NL1Mu]
   Float_t         L1MuEta[4];   //[NL1Mu]
   Float_t         L1MuPhi[4];   //[NL1Mu]
   Int_t           L1MuIsol[4];   //[NL1Mu]
   Int_t           L1MuMip[4];   //[NL1Mu]
   Int_t           L1MuFor[4];   //[NL1Mu]
   Int_t           L1MuRPC[4];   //[NL1Mu]
   Int_t           L1MuQal[4];   //[NL1Mu]
   Int_t           L1MuChg[4];   //[NL1Mu]
   Int_t           NL1CenJet;
   Float_t         L1CenJetEt[4];   //[NL1CenJet]
   Float_t         L1CenJetE[4];   //[NL1CenJet]
   Float_t         L1CenJetEta[4];   //[NL1CenJet]
   Float_t         L1CenJetPhi[4];   //[NL1CenJet]
   Int_t           NL1ForJet;
   Float_t         L1ForJetEt[4];   //[NL1ForJet]
   Float_t         L1ForJetE[4];   //[NL1ForJet]
   Float_t         L1ForJetEta[4];   //[NL1ForJet]
   Float_t         L1ForJetPhi[4];   //[NL1ForJet]
   Int_t           NL1Tau;
   Float_t         L1TauEt[4];   //[NL1Tau]
   Float_t         L1TauE[4];   //[NL1Tau]
   Float_t         L1TauEta[4];   //[NL1Tau]
   Float_t         L1TauPhi[4];   //[NL1Tau]
   Float_t         L1Met;
   Float_t         L1MetPhi;
   Float_t         L1EtTot;
   Float_t         L1Mht;
   Float_t         L1MhtPhi;
   Float_t         L1EtHad;
   Int_t           L1HfRing1EtSumPositiveEta;
   Int_t           L1HfRing2EtSumPositiveEta;
   Int_t           L1HfRing1EtSumNegativeEta;
   Int_t           L1HfRing2EtSumNegativeEta;
   Int_t           L1HfTowerCountPositiveEtaRing1;
   Int_t           L1HfTowerCountNegativeEtaRing1;
   Int_t           L1HfTowerCountPositiveEtaRing2;
   Int_t           L1HfTowerCountNegativeEtaRing2;
   Int_t           NMCpart;
   Int_t           MCpid[1];   //[NMCpart]
   Int_t           MCstatus[1];   //[NMCpart]
   Float_t         MCvtxX[1];   //[NMCpart]
   Float_t         MCvtxY[1];   //[NMCpart]
   Float_t         MCvtxZ[1];   //[NMCpart]
   Float_t         MCpt[1];   //[NMCpart]
   Float_t         MCeta[1];   //[NMCpart]
   Float_t         MCphi[1];   //[NMCpart]
   Float_t         MCPtHat;
   Int_t           MCmu3;
   Int_t           MCel3;
   Int_t           MCbb;
   Int_t           MCab;
   Int_t           MCWenu;
   Int_t           MCWmunu;
   Int_t           MCZee;
   Int_t           MCZmumu;
   Float_t         MCptEleMax;
   Float_t         MCptMuMax;
   Int_t           NPUTrueBX0;
   Int_t           NPUgenBX0;
   Int_t           recoNVrtOffline0;
   Float_t         recoVrtXOffline0[1];   //[NVrtx]
   Float_t         recoVrtYOffline0[1];   //[NVrtx]
   Float_t         recoVrtZOffline0[1];   //[NVrtx]
   Int_t           recoVrtNtrkOffline0[1];   //[NVrtx]
   Float_t         recoVrtChi2Offline0[1];   //[NVrtx]
   Float_t         recoVrtNdofOffline0[1];   //[NVrtx]
   Int_t           Run;
   Int_t           Event;
   Int_t           LumiBlock;
   Int_t           Bx;
   Int_t           Orbit;
   Double_t        AvgInstDelLumi;
   Int_t           HLT_HIMET220_v1;
   Int_t           HLT_HIMET220_v1_Prescl;
   Int_t           HLT_HIMET200_v1;
   Int_t           HLT_HIMET200_v1_Prescl;
   Int_t           HLT_HIMET120_v1;
   Int_t           HLT_HIMET120_v1_Prescl;
   Int_t           HLT_HIPhysics_v1;
   Int_t           HLT_HIPhysics_v1_Prescl;
   Int_t           HLT_HIDTCalibration_v1;
   Int_t           HLT_HIDTCalibration_v1_Prescl;
   Int_t           HLT_HIEcalCalibration_v1;
   Int_t           HLT_HIEcalCalibration_v1_Prescl;
   Int_t           HLT_HIHcalCalibration_v1;
   Int_t           HLT_HIHcalCalibration_v1_Prescl;
   Int_t           HLT_HIZeroBias_v1;
   Int_t           HLT_HIZeroBias_v1_Prescl;
   Int_t           HLT_HIZeroBiasXOR_v1;
   Int_t           HLT_HIZeroBiasXOR_v1_Prescl;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v1;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIMinBiasBSC_v1;
   Int_t           HLT_HIMinBiasBSC_v1_Prescl;
   Int_t           HLT_HIMinBiasBSC_OR_v1;
   Int_t           HLT_HIMinBiasBSC_OR_v1_Prescl;
   Int_t           HLT_HIMinBiasHF_v1;
   Int_t           HLT_HIMinBiasHF_v1_Prescl;
   Int_t           HLT_HIMinBiasHf_OR_v1;
   Int_t           HLT_HIMinBiasHf_OR_v1_Prescl;
   Int_t           HLT_HIMinBiasHfOrBSC_v1;
   Int_t           HLT_HIMinBiasHfOrBSC_v1_Prescl;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v1;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_v1;
   Int_t           HLT_HIMinBiasZDC_Calo_v1_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v1;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIBptxXOR_v1;
   Int_t           HLT_HIBptxXOR_v1_Prescl;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v1;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl;
   Int_t           HLT_HIL1DoubleMuOpen_v2;
   Int_t           HLT_HIL1DoubleMuOpen_v2_Prescl;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v2;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v2_Prescl;
   Int_t           HLT_HIL2Mu3_v2;
   Int_t           HLT_HIL2Mu3_v2_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ_v2;
   Int_t           HLT_HIL2Mu3_NHitQ_v2_Prescl;
   Int_t           HLT_HIL2Mu7_v2;
   Int_t           HLT_HIL2Mu7_v2_Prescl;
   Int_t           HLT_HIL2Mu15_v2;
   Int_t           HLT_HIL2Mu15_v2_Prescl;
   Int_t           HLT_HIL2DoubleMu0_v2;
   Int_t           HLT_HIL2DoubleMu0_v2_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v2;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v2_Prescl;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl;
   Int_t           HLT_HIL2DoubleMu3_v2;
   Int_t           HLT_HIL2DoubleMu3_v2_Prescl;
   Int_t           HLT_HIL3Mu3_v2;
   Int_t           HLT_HIL3Mu3_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_v2;
   Int_t           HLT_HIL3DoubleMuOpen_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v2;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v2;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v2;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl;
   Int_t           HLT_HISinglePhoton15_v1;
   Int_t           HLT_HISinglePhoton15_v1_Prescl;
   Int_t           HLT_HISinglePhoton20_v2;
   Int_t           HLT_HISinglePhoton20_v2_Prescl;
   Int_t           HLT_HISinglePhoton30_v2;
   Int_t           HLT_HISinglePhoton30_v2_Prescl;
   Int_t           HLT_HISinglePhoton40_v2;
   Int_t           HLT_HISinglePhoton40_v2_Prescl;
   Int_t           HLT_HIPhoton10_Photon15_v1;
   Int_t           HLT_HIPhoton10_Photon15_v1_Prescl;
   Int_t           HLT_HIPhoton15_Photon20_v1;
   Int_t           HLT_HIPhoton15_Photon20_v1_Prescl;
   Int_t           HLT_HIDoublePhoton10_v1;
   Int_t           HLT_HIDoublePhoton10_v1_Prescl;
   Int_t           HLT_HIDoublePhoton15_v1;
   Int_t           HLT_HIDoublePhoton15_v1_Prescl;
   Int_t           HLT_HIDoublePhoton20_v1;
   Int_t           HLT_HIDoublePhoton20_v1_Prescl;
   Int_t           HLT_HIJet55_v1;
   Int_t           HLT_HIJet55_v1_Prescl;
   Int_t           HLT_HIJet65_v1;
   Int_t           HLT_HIJet65_v1_Prescl;
   Int_t           HLT_HIJet80_v1;
   Int_t           HLT_HIJet80_v1_Prescl;
   Int_t           HLT_HIJet80_v7;
   Int_t           HLT_HIJet80_v7_Prescl;
   Int_t           HLT_HIJet95_v1;
   Int_t           HLT_HIJet95_v1_Prescl;
   Int_t           HLT_HIDiJet55_v1;
   Int_t           HLT_HIDiJet55_v1_Prescl;
   Int_t           HLT_HIJet65_Jet55_v1;
   Int_t           HLT_HIJet65_Jet55_v1_Prescl;
   Int_t           HLT_HIJetE30_NoBPTX_v1;
   Int_t           HLT_HIJetE30_NoBPTX_v1_Prescl;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v1;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl;
   Int_t           HLT_HIActivityHF_Coincidence3_v1;
   Int_t           HLT_HIActivityHF_Coincidence3_v1_Prescl;
   Int_t           HLT_HIActivityHF_Single3_v1;
   Int_t           HLT_HIActivityHF_Single3_v1_Prescl;
   Int_t           HLT_HIClusterVertexCompatibility_v1;
   Int_t           HLT_HIClusterVertexCompatibility_v1_Prescl;
   Int_t           HLT_HICentralityVeto_v1;
   Int_t           HLT_HICentralityVeto_v1_Prescl;
   Int_t           HLT_HIFullTrack12_L1Central_v1;
   Int_t           HLT_HIFullTrack12_L1Central_v1_Prescl;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v1;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v1_Prescl;
   Int_t           HLT_HIFullTrack14_L1Central_v1;
   Int_t           HLT_HIFullTrack14_L1Central_v1_Prescl;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v1;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v1_Prescl;
   Int_t           HLT_HIFullTrack20_L1Central_v1;
   Int_t           HLT_HIFullTrack20_L1Central_v1_Prescl;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v1;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v1_Prescl;
   Int_t           HLT_HIFullTrack25_L1Central_v1;
   Int_t           HLT_HIFullTrack25_L1Central_v1_Prescl;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v1;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v1_Prescl;
   Int_t           HLT_HIRandom_v1;
   Int_t           HLT_HIRandom_v1_Prescl;
   Int_t           HLT_HIUCC010_v2;
   Int_t           HLT_HIUCC010_v2_Prescl;
   Int_t           HLT_HIUCC015_v2;
   Int_t           HLT_HIUCC015_v2_Prescl;
   Int_t           HLT_HICentral10_v2;
   Int_t           HLT_HICentral10_v2_Prescl;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl;
   Int_t           HLTriggerFinalPath;
   Int_t           HLTriggerFinalPath_Prescl;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v1;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v1_Prescl;
   Int_t           HLT_HIL1DoubleMuOpen_v1;
   Int_t           HLT_HIL1DoubleMuOpen_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v1;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu0_v1;
   Int_t           HLT_HIL2DoubleMu0_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu3_v1;
   Int_t           HLT_HIL2DoubleMu3_v1_Prescl;
   Int_t           HLT_HIL2Mu15_v1;
   Int_t           HLT_HIL2Mu15_v1_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ_v1;
   Int_t           HLT_HIL2Mu3_NHitQ_v1_Prescl;
   Int_t           HLT_HIL2Mu3_v1;
   Int_t           HLT_HIL2Mu3_v1_Prescl;
   Int_t           HLT_HIL2Mu7_v1;
   Int_t           HLT_HIL2Mu7_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v1;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v1;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v1;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_v1;
   Int_t           HLT_HIL3DoubleMuOpen_v1_Prescl;
   Int_t           HLT_HIL3Mu3_v1;
   Int_t           HLT_HIL3Mu3_v1_Prescl;
   Int_t           HLT_HISinglePhoton20_v1;
   Int_t           HLT_HISinglePhoton20_v1_Prescl;
   Int_t           HLT_HISinglePhoton30_v1;
   Int_t           HLT_HISinglePhoton30_v1_Prescl;
   Int_t           HLT_HISinglePhoton40_v1;
   Int_t           HLT_HISinglePhoton40_v1_Prescl;
   Int_t           HLT_HICentral10_v1;
   Int_t           HLT_HICentral10_v1_Prescl;
   Int_t           HLT_HIUCC010_v1;
   Int_t           HLT_HIUCC010_v1_Prescl;
   Int_t           HLT_HIUCC015_v1;
   Int_t           HLT_HIUCC015_v1_Prescl;
   Int_t           HLT_HIDiJet55_v2;
   Int_t           HLT_HIDiJet55_v2_Prescl;
   Int_t           HLT_HIDoublePhoton10_v2;
   Int_t           HLT_HIDoublePhoton10_v2_Prescl;
   Int_t           HLT_HIDoublePhoton15_v2;
   Int_t           HLT_HIDoublePhoton15_v2_Prescl;
   Int_t           HLT_HIDoublePhoton20_v2;
   Int_t           HLT_HIDoublePhoton20_v2_Prescl;
   Int_t           HLT_HIFullTrack12_L1Central_v2;
   Int_t           HLT_HIFullTrack12_L1Central_v2_Prescl;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v2;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v2_Prescl;
   Int_t           HLT_HIFullTrack14_L1Central_v2;
   Int_t           HLT_HIFullTrack14_L1Central_v2_Prescl;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v2;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v2_Prescl;
   Int_t           HLT_HIFullTrack20_L1Central_v2;
   Int_t           HLT_HIFullTrack20_L1Central_v2_Prescl;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v2;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v2_Prescl;
   Int_t           HLT_HIFullTrack25_L1Central_v2;
   Int_t           HLT_HIFullTrack25_L1Central_v2_Prescl;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v2;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v2_Prescl;
   Int_t           HLT_HIJet55_v2;
   Int_t           HLT_HIJet55_v2_Prescl;
   Int_t           HLT_HIJet65_Jet55_v2;
   Int_t           HLT_HIJet65_Jet55_v2_Prescl;
   Int_t           HLT_HIJet65_v2;
   Int_t           HLT_HIJet65_v2_Prescl;
   Int_t           HLT_HIJet80_v2;
   Int_t           HLT_HIJet80_v2_Prescl;
   Int_t           HLT_HIJet95_v2;
   Int_t           HLT_HIJet95_v2_Prescl;
   Int_t           HLT_HIJetE30_NoBPTX_v2;
   Int_t           HLT_HIJetE30_NoBPTX_v2_Prescl;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v2;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl;
   Int_t           HLT_HIMET120_v2;
   Int_t           HLT_HIMET120_v2_Prescl;
   Int_t           HLT_HIMET200_v2;
   Int_t           HLT_HIMET200_v2_Prescl;
   Int_t           HLT_HIMET220_v2;
   Int_t           HLT_HIMET220_v2_Prescl;
   Int_t           HLT_HIPhoton10_Photon15_v2;
   Int_t           HLT_HIPhoton10_Photon15_v2_Prescl;
   Int_t           HLT_HIPhoton15_Photon20_v2;
   Int_t           HLT_HIPhoton15_Photon20_v2_Prescl;
   Int_t           HLT_HISinglePhoton15_v2;
   Int_t           HLT_HISinglePhoton15_v2_Prescl;
   Int_t           HLT_HIBptxXOR_v2;
   Int_t           HLT_HIBptxXOR_v2_Prescl;
   Int_t           HLT_HICentralityVeto_v2;
   Int_t           HLT_HICentralityVeto_v2_Prescl;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v2;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl;
   Int_t           HLT_HIMinBiasBSC_OR_v2;
   Int_t           HLT_HIMinBiasBSC_OR_v2_Prescl;
   Int_t           HLT_HIMinBiasBSC_v2;
   Int_t           HLT_HIMinBiasBSC_v2_Prescl;
   Int_t           HLT_HIMinBiasHF_v2;
   Int_t           HLT_HIMinBiasHF_v2_Prescl;
   Int_t           HLT_HIMinBiasHfOrBSC_v2;
   Int_t           HLT_HIMinBiasHfOrBSC_v2_Prescl;
   Int_t           HLT_HIMinBiasHf_OR_v2;
   Int_t           HLT_HIMinBiasHf_OR_v2_Prescl;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v2;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v2;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_v2;
   Int_t           HLT_HIMinBiasZDC_Calo_v2_Prescl;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIPhysics_v2;
   Int_t           HLT_HIPhysics_v2_Prescl;
   Int_t           HLT_HIRandom_v2;
   Int_t           HLT_HIRandom_v2_Prescl;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v2;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v2;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl;
   Int_t           HLT_HIZeroBiasXOR_v2;
   Int_t           HLT_HIZeroBiasXOR_v2_Prescl;
   Int_t           HLT_HIZeroBias_v2;
   Int_t           HLT_HIZeroBias_v2_Prescl;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v3;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v3_Prescl;
   Int_t           HLT_HIL1DoubleMuOpen_v3;
   Int_t           HLT_HIL1DoubleMuOpen_v3_Prescl;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v3;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v3_Prescl;
   Int_t           HLT_HIL2DoubleMu0_v3;
   Int_t           HLT_HIL2DoubleMu0_v3_Prescl;
   Int_t           HLT_HIL2DoubleMu3_v3;
   Int_t           HLT_HIL2DoubleMu3_v3_Prescl;
   Int_t           HLT_HIL2Mu15_v3;
   Int_t           HLT_HIL2Mu15_v3_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ_v3;
   Int_t           HLT_HIL2Mu3_NHitQ_v3_Prescl;
   Int_t           HLT_HIL2Mu3_v3;
   Int_t           HLT_HIL2Mu3_v3_Prescl;
   Int_t           HLT_HIL2Mu7_v3;
   Int_t           HLT_HIL2Mu7_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v3;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v3;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v3;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_v3;
   Int_t           HLT_HIL3DoubleMuOpen_v3_Prescl;
   Int_t           HLT_HIL3Mu3_v3;
   Int_t           HLT_HIL3Mu3_v3_Prescl;
   Int_t           HLT_HIDiJet55_v3;
   Int_t           HLT_HIDiJet55_v3_Prescl;
   Int_t           HLT_HIDoublePhoton10_v3;
   Int_t           HLT_HIDoublePhoton10_v3_Prescl;
   Int_t           HLT_HIDoublePhoton15_v3;
   Int_t           HLT_HIDoublePhoton15_v3_Prescl;
   Int_t           HLT_HIDoublePhoton20_v3;
   Int_t           HLT_HIDoublePhoton20_v3_Prescl;
   Int_t           HLT_HIFullTrack12_L1Central_v3;
   Int_t           HLT_HIFullTrack12_L1Central_v3_Prescl;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v3;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v3_Prescl;
   Int_t           HLT_HIFullTrack14_L1Central_v3;
   Int_t           HLT_HIFullTrack14_L1Central_v3_Prescl;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v3;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v3_Prescl;
   Int_t           HLT_HIFullTrack20_L1Central_v3;
   Int_t           HLT_HIFullTrack20_L1Central_v3_Prescl;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v3;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v3_Prescl;
   Int_t           HLT_HIFullTrack25_L1Central_v3;
   Int_t           HLT_HIFullTrack25_L1Central_v3_Prescl;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v3;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v3_Prescl;
   Int_t           HLT_HIJet55_v3;
   Int_t           HLT_HIJet55_v3_Prescl;
   Int_t           HLT_HIJet65_Jet55_v3;
   Int_t           HLT_HIJet65_Jet55_v3_Prescl;
   Int_t           HLT_HIJet65_v3;
   Int_t           HLT_HIJet65_v3_Prescl;
   Int_t           HLT_HIJet80_v3;
   Int_t           HLT_HIJet80_v3_Prescl;
   Int_t           HLT_HIJet95_v3;
   Int_t           HLT_HIJet95_v3_Prescl;
   Int_t           HLT_HIJetE30_NoBPTX_v3;
   Int_t           HLT_HIJetE30_NoBPTX_v3_Prescl;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v3;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl;
   Int_t           HLT_HIMET120_v3;
   Int_t           HLT_HIMET120_v3_Prescl;
   Int_t           HLT_HIMET200_v3;
   Int_t           HLT_HIMET200_v3_Prescl;
   Int_t           HLT_HIMET220_v3;
   Int_t           HLT_HIMET220_v3_Prescl;
   Int_t           HLT_HIPhoton10_Photon15_v3;
   Int_t           HLT_HIPhoton10_Photon15_v3_Prescl;
   Int_t           HLT_HIPhoton15_Photon20_v3;
   Int_t           HLT_HIPhoton15_Photon20_v3_Prescl;
   Int_t           HLT_HISinglePhoton15_v3;
   Int_t           HLT_HISinglePhoton15_v3_Prescl;
   Int_t           HLT_HISinglePhoton20_v3;
   Int_t           HLT_HISinglePhoton20_v3_Prescl;
   Int_t           HLT_HISinglePhoton30_v3;
   Int_t           HLT_HISinglePhoton30_v3_Prescl;
   Int_t           HLT_HISinglePhoton40_v3;
   Int_t           HLT_HISinglePhoton40_v3_Prescl;
   Int_t           HLT_HIBptxXOR_v3;
   Int_t           HLT_HIBptxXOR_v3_Prescl;
   Int_t           HLT_HICentral10_v3;
   Int_t           HLT_HICentral10_v3_Prescl;
   Int_t           HLT_HICentralityVeto_v3;
   Int_t           HLT_HICentralityVeto_v3_Prescl;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v3;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl;
   Int_t           HLT_HIMinBiasBSC_OR_v3;
   Int_t           HLT_HIMinBiasBSC_OR_v3_Prescl;
   Int_t           HLT_HIMinBiasBSC_v3;
   Int_t           HLT_HIMinBiasBSC_v3_Prescl;
   Int_t           HLT_HIMinBiasHF_v3;
   Int_t           HLT_HIMinBiasHF_v3_Prescl;
   Int_t           HLT_HIMinBiasHfOrBSC_v3;
   Int_t           HLT_HIMinBiasHfOrBSC_v3_Prescl;
   Int_t           HLT_HIMinBiasHf_OR_v3;
   Int_t           HLT_HIMinBiasHf_OR_v3_Prescl;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v3;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v3;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_v3;
   Int_t           HLT_HIMinBiasZDC_Calo_v3_Prescl;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIPhysics_v3;
   Int_t           HLT_HIPhysics_v3_Prescl;
   Int_t           HLT_HIRandom_v3;
   Int_t           HLT_HIRandom_v3_Prescl;
   Int_t           HLT_HIUCC010_v3;
   Int_t           HLT_HIUCC010_v3_Prescl;
   Int_t           HLT_HIUCC015_v3;
   Int_t           HLT_HIUCC015_v3_Prescl;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v3;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v3;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl;
   Int_t           HLT_HIZeroBiasXOR_v3;
   Int_t           HLT_HIZeroBiasXOR_v3_Prescl;
   Int_t           HLT_HIZeroBias_v3;
   Int_t           HLT_HIZeroBias_v3_Prescl;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v4;
   Int_t           HLT_HIL1DoubleMu0_HighQ_v4_Prescl;
   Int_t           HLT_HIL1DoubleMuOpen_v4;
   Int_t           HLT_HIL1DoubleMuOpen_v4_Prescl;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4;
   Int_t           HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v4;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v4_Prescl;
   Int_t           HLT_HIL2DoubleMu0_v4;
   Int_t           HLT_HIL2DoubleMu0_v4_Prescl;
   Int_t           HLT_HIL2DoubleMu3_v4;
   Int_t           HLT_HIL2DoubleMu3_v4_Prescl;
   Int_t           HLT_HIL2Mu15_v4;
   Int_t           HLT_HIL2Mu15_v4_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ_v4;
   Int_t           HLT_HIL2Mu3_NHitQ_v4_Prescl;
   Int_t           HLT_HIL2Mu3_v4;
   Int_t           HLT_HIL2Mu3_v4_Prescl;
   Int_t           HLT_HIL2Mu7_v4;
   Int_t           HLT_HIL2Mu7_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v4;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v4;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v4;
   Int_t           HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl;
   Int_t           HLT_HIL3DoubleMuOpen_v4;
   Int_t           HLT_HIL3DoubleMuOpen_v4_Prescl;
   Int_t           HLT_HIL3Mu3_v4;
   Int_t           HLT_HIL3Mu3_v4_Prescl;
   Int_t           HLT_HIDiJet55_v4;
   Int_t           HLT_HIDiJet55_v4_Prescl;
   Int_t           HLT_HIDoublePhoton10_v4;
   Int_t           HLT_HIDoublePhoton10_v4_Prescl;
   Int_t           HLT_HIDoublePhoton15_v4;
   Int_t           HLT_HIDoublePhoton15_v4_Prescl;
   Int_t           HLT_HIDoublePhoton20_v4;
   Int_t           HLT_HIDoublePhoton20_v4_Prescl;
   Int_t           HLT_HIFullTrack12_L1Central_v4;
   Int_t           HLT_HIFullTrack12_L1Central_v4_Prescl;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v4;
   Int_t           HLT_HIFullTrack12_L1Peripheral_v4_Prescl;
   Int_t           HLT_HIFullTrack14_L1Central_v4;
   Int_t           HLT_HIFullTrack14_L1Central_v4_Prescl;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v4;
   Int_t           HLT_HIFullTrack14_L1Peripheral_v4_Prescl;
   Int_t           HLT_HIFullTrack20_L1Central_v4;
   Int_t           HLT_HIFullTrack20_L1Central_v4_Prescl;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v4;
   Int_t           HLT_HIFullTrack20_L1Peripheral_v4_Prescl;
   Int_t           HLT_HIFullTrack25_L1Central_v4;
   Int_t           HLT_HIFullTrack25_L1Central_v4_Prescl;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v4;
   Int_t           HLT_HIFullTrack25_L1Peripheral_v4_Prescl;
   Int_t           HLT_HIJet55_v4;
   Int_t           HLT_HIJet55_v4_Prescl;
   Int_t           HLT_HIJet65_Jet55_v4;
   Int_t           HLT_HIJet65_Jet55_v4_Prescl;
   Int_t           HLT_HIJet65_v4;
   Int_t           HLT_HIJet65_v4_Prescl;
   Int_t           HLT_HIJet80_v4;
   Int_t           HLT_HIJet80_v4_Prescl;
   Int_t           HLT_HIJet95_v4;
   Int_t           HLT_HIJet95_v4_Prescl;
   Int_t           HLT_HIJetE30_NoBPTX_v4;
   Int_t           HLT_HIJetE30_NoBPTX_v4_Prescl;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v4;
   Int_t           HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl;
   Int_t           HLT_HIMET120_v4;
   Int_t           HLT_HIMET120_v4_Prescl;
   Int_t           HLT_HIMET200_v4;
   Int_t           HLT_HIMET200_v4_Prescl;
   Int_t           HLT_HIMET220_v4;
   Int_t           HLT_HIMET220_v4_Prescl;
   Int_t           HLT_HIPhoton10_Photon15_v4;
   Int_t           HLT_HIPhoton10_Photon15_v4_Prescl;
   Int_t           HLT_HIPhoton15_Photon20_v4;
   Int_t           HLT_HIPhoton15_Photon20_v4_Prescl;
   Int_t           HLT_HISinglePhoton15_v4;
   Int_t           HLT_HISinglePhoton15_v4_Prescl;
   Int_t           HLT_HISinglePhoton20_v4;
   Int_t           HLT_HISinglePhoton20_v4_Prescl;
   Int_t           HLT_HISinglePhoton30_v4;
   Int_t           HLT_HISinglePhoton30_v4_Prescl;
   Int_t           HLT_HISinglePhoton40_v4;
   Int_t           HLT_HISinglePhoton40_v4_Prescl;
   Int_t           HLT_HIBptxXOR_v4;
   Int_t           HLT_HIBptxXOR_v4_Prescl;
   Int_t           HLT_HICentral10_v4;
   Int_t           HLT_HICentral10_v4_Prescl;
   Int_t           HLT_HICentralityVeto_v4;
   Int_t           HLT_HICentralityVeto_v4_Prescl;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v4;
   Int_t           HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl;
   Int_t           HLT_HIMinBiasBSC_OR_v4;
   Int_t           HLT_HIMinBiasBSC_OR_v4_Prescl;
   Int_t           HLT_HIMinBiasBSC_v4;
   Int_t           HLT_HIMinBiasBSC_v4_Prescl;
   Int_t           HLT_HIMinBiasHF_v4;
   Int_t           HLT_HIMinBiasHF_v4_Prescl;
   Int_t           HLT_HIMinBiasHfOrBSC_v4;
   Int_t           HLT_HIMinBiasHfOrBSC_v4_Prescl;
   Int_t           HLT_HIMinBiasHf_OR_v4;
   Int_t           HLT_HIMinBiasHf_OR_v4_Prescl;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v4;
   Int_t           HLT_HIMinBiasPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v4;
   Int_t           HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4;
   Int_t           HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl;
   Int_t           HLT_HIMinBiasZDC_Calo_v4;
   Int_t           HLT_HIMinBiasZDC_Calo_v4_Prescl;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4;
   Int_t           HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIPhysics_v4;
   Int_t           HLT_HIPhysics_v4_Prescl;
   Int_t           HLT_HIRandom_v4;
   Int_t           HLT_HIRandom_v4_Prescl;
   Int_t           HLT_HIUCC010_v4;
   Int_t           HLT_HIUCC010_v4_Prescl;
   Int_t           HLT_HIUCC015_v4;
   Int_t           HLT_HIUCC015_v4_Prescl;
   Int_t           HLT_PAJet80_NoJetID_v1;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v4;
   Int_t           HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v4;
   Int_t           HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl;
   Int_t           HLT_HIZeroBiasXOR_v4;
   Int_t           HLT_HIZeroBiasXOR_v4_Prescl;
   Int_t           HLT_HIZeroBias_v4;
   Int_t           HLT_HIZeroBias_v4_Prescl;
   Int_t           L1_BptxMinus;
   Int_t           L1_BptxMinus_Prescl;
   Int_t           L1_BptxMinus_NotBptxPlus;
   Int_t           L1_BptxMinus_NotBptxPlus_Prescl;
   Int_t           L1_BptxPlus;
   Int_t           L1_BptxPlus_Prescl;
   Int_t           L1_BptxPlusANDMinus;
   Int_t           L1_BptxPlusANDMinus_Prescl;
   Int_t           L1_BptxPlusORMinus;
   Int_t           L1_BptxPlusORMinus_Prescl;
   Int_t           L1_BptxPlusORMinus_instance1;
   Int_t           L1_BptxPlusORMinus_instance1_Prescl;
   Int_t           L1_BptxPlus_NotBptxMinus;
   Int_t           L1_BptxPlus_NotBptxMinus_Prescl;
   Int_t           L1_BptxXOR;
   Int_t           L1_BptxXOR_Prescl;
   Int_t           L1_BptxXOR_BscMinBiasOR;
   Int_t           L1_BptxXOR_BscMinBiasOR_Prescl;
   Int_t           L1_Bsc2Minus_BptxMinus;
   Int_t           L1_Bsc2Minus_BptxMinus_Prescl;
   Int_t           L1_Bsc2Plus_BptxPlus;
   Int_t           L1_Bsc2Plus_BptxPlus_Prescl;
   Int_t           L1_BscMinBiasInnerThreshold1;
   Int_t           L1_BscMinBiasInnerThreshold1_Prescl;
   Int_t           L1_BscMinBiasInnerThreshold1_BptxAND;
   Int_t           L1_BscMinBiasInnerThreshold1_BptxAND_Prescl;
   Int_t           L1_BscMinBiasInnerThreshold2_BptxAND;
   Int_t           L1_BscMinBiasInnerThreshold2_BptxAND_Prescl;
   Int_t           L1_BscMinBiasOR_BptxAND;
   Int_t           L1_BscMinBiasOR_BptxAND_Prescl;
   Int_t           L1_BscMinBiasOR_BptxPlusORMinus;
   Int_t           L1_BscMinBiasOR_BptxPlusORMinus_Prescl;
   Int_t           L1_BscMinBiasThreshold1;
   Int_t           L1_BscMinBiasThreshold1_Prescl;
   Int_t           L1_BscMinBiasThreshold1_BptxAND;
   Int_t           L1_BscMinBiasThreshold1_BptxAND_Prescl;
   Int_t           L1_BscMinBiasThreshold2_BptxAND;
   Int_t           L1_BscMinBiasThreshold2_BptxAND_Prescl;
   Int_t           L1_DoubleEG10;
   Int_t           L1_DoubleEG10_Prescl;
   Int_t           L1_DoubleEG3;
   Int_t           L1_DoubleEG3_Prescl;
   Int_t           L1_DoubleEG5_BptxAND;
   Int_t           L1_DoubleEG5_BptxAND_Prescl;
   Int_t           L1_DoubleEG8_BptxAND;
   Int_t           L1_DoubleEG8_BptxAND_Prescl;
   Int_t           L1_DoubleForJet32_EtaOpp;
   Int_t           L1_DoubleForJet32_EtaOpp_Prescl;
   Int_t           L1_DoubleForJet44_EtaOpp;
   Int_t           L1_DoubleForJet44_EtaOpp_Prescl;
   Int_t           L1_DoubleHfBitCountsRing1_P1N1;
   Int_t           L1_DoubleHfBitCountsRing1_P1N1_Prescl;
   Int_t           L1_DoubleHfBitCountsRing2_P1N1;
   Int_t           L1_DoubleHfBitCountsRing2_P1N1_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing1_P200N200;
   Int_t           L1_DoubleHfRingEtSumsRing1_P200N200_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing1_P4N4;
   Int_t           L1_DoubleHfRingEtSumsRing1_P4N4_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing2_P200N200;
   Int_t           L1_DoubleHfRingEtSumsRing2_P200N200_Prescl;
   Int_t           L1_DoubleHfRingEtSumsRing2_P4N4;
   Int_t           L1_DoubleHfRingEtSumsRing2_P4N4_Prescl;
   Int_t           L1_DoubleMu3_BptxAND;
   Int_t           L1_DoubleMu3_BptxAND_Prescl;
   Int_t           L1_DoubleMuOpen_BptxAND;
   Int_t           L1_DoubleMuOpen_BptxAND_Prescl;
   Int_t           L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND;
   Int_t           L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;
   Int_t           L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;
   Int_t           L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;
   Int_t           L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND;
   Int_t           L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;
   Int_t           L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;
   Int_t           L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;
   Int_t           L1_ETM30_BptxAND;
   Int_t           L1_ETM30_BptxAND_Prescl;
   Int_t           L1_ETM50_BptxAND;
   Int_t           L1_ETM50_BptxAND_Prescl;
   Int_t           L1_ETT100_BptxAND;
   Int_t           L1_ETT100_BptxAND_Prescl;
   Int_t           L1_ETT140_BptxAND;
   Int_t           L1_ETT140_BptxAND_Prescl;
   Int_t           L1_ETT2000;
   Int_t           L1_ETT2000_Prescl;
   Int_t           L1_ETT220_BptxAND;
   Int_t           L1_ETT220_BptxAND_Prescl;
   Int_t           L1_ETT70_BptxAND;
   Int_t           L1_ETT70_BptxAND_Prescl;
   Int_t           L1_ETT800;
   Int_t           L1_ETT800_Prescl;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND;
   Int_t           L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl;
   Int_t           L1_HcalHfCoincidencePm;
   Int_t           L1_HcalHfCoincidencePm_Prescl;
   Int_t           L1_HcalHfCoincidencePm_BptxAND;
   Int_t           L1_HcalHfCoincidencePm_BptxAND_Prescl;
   Int_t           L1_HcalHfMmOrPpOrPm;
   Int_t           L1_HcalHfMmOrPpOrPm_Prescl;
   Int_t           L1_HcalHfMmOrPpOrPm_BptxAND;
   Int_t           L1_HcalHfMmOrPpOrPm_BptxAND_Prescl;
   Int_t           L1_HcalHfMmpOrMpp;
   Int_t           L1_HcalHfMmpOrMpp_Prescl;
   Int_t           L1_HcalHfMmpOrMpp_BptxAND;
   Int_t           L1_HcalHfMmpOrMpp_BptxAND_Prescl;
   Int_t           L1_HcalHoTotalOR;
   Int_t           L1_HcalHoTotalOR_Prescl;
   Int_t           L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND;
   Int_t           L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;
   Int_t           L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;
   Int_t           L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasOR;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasThreshold1;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasThreshold2;
   Int_t           L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl;
   Int_t           L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1;
   Int_t           L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl;
   Int_t           L1_NotBsc2_BptxAND_HcalHfCoincidencePm;
   Int_t           L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl;
   Int_t           L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm;
   Int_t           L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl;
   Int_t           L1_NotBsc2_BptxAND_HcalHfMmpOrMpp;
   Int_t           L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasInnerThreshold1;
   Int_t           L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasInnerThreshold2;
   Int_t           L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasOR;
   Int_t           L1_NotBsc2_BscMinBiasOR_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasThreshold1;
   Int_t           L1_NotBsc2_BscMinBiasThreshold1_Prescl;
   Int_t           L1_NotBsc2_BscMinBiasThreshold2;
   Int_t           L1_NotBsc2_BscMinBiasThreshold2_Prescl;
   Int_t           L1_NotBsc2_HcalHfCoincidencePm;
   Int_t           L1_NotBsc2_HcalHfCoincidencePm_Prescl;
   Int_t           L1_NotBsc2_HcalHfMmOrPpOrPm;
   Int_t           L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl;
   Int_t           L1_NotBsc2_HcalHfMmpOrMpp;
   Int_t           L1_NotBsc2_HcalHfMmpOrMpp_Prescl;
   Int_t           L1_SingleEG12;
   Int_t           L1_SingleEG12_Prescl;
   Int_t           L1_SingleEG15;
   Int_t           L1_SingleEG15_Prescl;
   Int_t           L1_SingleEG18;
   Int_t           L1_SingleEG18_Prescl;
   Int_t           L1_SingleEG20;
   Int_t           L1_SingleEG20_Prescl;
   Int_t           L1_SingleEG22;
   Int_t           L1_SingleEG22_Prescl;
   Int_t           L1_SingleEG2_BptxAND;
   Int_t           L1_SingleEG2_BptxAND_Prescl;
   Int_t           L1_SingleEG30;
   Int_t           L1_SingleEG30_Prescl;
   Int_t           L1_SingleEG5_BptxAND;
   Int_t           L1_SingleEG5_BptxAND_Prescl;
   Int_t           L1_SingleEG8_BptxAND;
   Int_t           L1_SingleEG8_BptxAND_Prescl;
   Int_t           L1_SingleJet128;
   Int_t           L1_SingleJet128_Prescl;
   Int_t           L1_SingleJet128_BptxAND;
   Int_t           L1_SingleJet128_BptxAND_Prescl;
   Int_t           L1_SingleJet128_Central_BptxAND;
   Int_t           L1_SingleJet128_Central_BptxAND_Prescl;
   Int_t           L1_SingleJet16;
   Int_t           L1_SingleJet16_Prescl;
   Int_t           L1_SingleJet16_Central_NotETT140_BptxAND;
   Int_t           L1_SingleJet16_Central_NotETT140_BptxAND_Prescl;
   Int_t           L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo;
   Int_t           L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl;
   Int_t           L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo;
   Int_t           L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl;
   Int_t           L1_SingleJet36_BptxAND;
   Int_t           L1_SingleJet36_BptxAND_Prescl;
   Int_t           L1_SingleJet52_BptxAND;
   Int_t           L1_SingleJet52_BptxAND_Prescl;
   Int_t           L1_SingleJet52_Central_BptxAND;
   Int_t           L1_SingleJet52_Central_BptxAND_Prescl;
   Int_t           L1_SingleJet68_BptxAND;
   Int_t           L1_SingleJet68_BptxAND_Prescl;
   Int_t           L1_SingleJet80_BptxAND;
   Int_t           L1_SingleJet80_BptxAND_Prescl;
   Int_t           L1_SingleJet80_Central;
   Int_t           L1_SingleJet80_Central_Prescl;
   Int_t           L1_SingleJet80_Central_BptxAND;
   Int_t           L1_SingleJet80_Central_BptxAND_Prescl;
   Int_t           L1_SingleJet92;
   Int_t           L1_SingleJet92_Prescl;
   Int_t           L1_SingleJet92_BptxAND;
   Int_t           L1_SingleJet92_BptxAND_Prescl;
   Int_t           L1_SingleJet92_Central;
   Int_t           L1_SingleJet92_Central_Prescl;
   Int_t           L1_SingleJet92_Central_BptxAND;
   Int_t           L1_SingleJet92_Central_BptxAND_Prescl;
   Int_t           L1_SingleMu0_BptxAND;
   Int_t           L1_SingleMu0_BptxAND_Prescl;
   Int_t           L1_SingleMu10;
   Int_t           L1_SingleMu10_Prescl;
   Int_t           L1_SingleMu14;
   Int_t           L1_SingleMu14_Prescl;
   Int_t           L1_SingleMu20;
   Int_t           L1_SingleMu20_Prescl;
   Int_t           L1_SingleMu3_BptxAND;
   Int_t           L1_SingleMu3_BptxAND_Prescl;
   Int_t           L1_SingleMu5_BptxAND;
   Int_t           L1_SingleMu5_BptxAND_Prescl;
   Int_t           L1_SingleMu7;
   Int_t           L1_SingleMu7_Prescl;
   Int_t           L1_SingleMuBeamHalo;
   Int_t           L1_SingleMuBeamHalo_Prescl;
   Int_t           L1_SingleMuOpen_BptxAND;
   Int_t           L1_SingleMuOpen_BptxAND_Prescl;
   Int_t           L1_ZdcCaloMinus;
   Int_t           L1_ZdcCaloMinus_Prescl;
   Int_t           L1_ZdcCaloMinus_BptxAND;
   Int_t           L1_ZdcCaloMinus_BptxAND_Prescl;
   Int_t           L1_ZdcCaloPlus;
   Int_t           L1_ZdcCaloPlus_Prescl;
   Int_t           L1_ZdcCaloPlus_BptxAND;
   Int_t           L1_ZdcCaloPlus_BptxAND_Prescl;
   Int_t           L1_ZdcCaloPlus_ZdcCaloMinus;
   Int_t           L1_ZdcCaloPlus_ZdcCaloMinus_Prescl;
   Int_t           L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND;
   Int_t           L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl;
   Int_t           L1_ZeroBias;
   Int_t           L1_ZeroBias_Prescl;
   Int_t           L1Tech_BPTX_PreBPTX_v0;
   Int_t           L1Tech_BPTX_PreBPTX_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_v0;
   Int_t           L1Tech_BPTX_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0;
   Int_t           L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_v0;
   Int_t           L1Tech_BPTX_plus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0;
   Int_t           L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0;
   Int_t           L1Tech_BPTX_plus_OR_minus_v0_Prescl;
   Int_t           L1Tech_BPTX_quiet_v0;
   Int_t           L1Tech_BPTX_quiet_v0_Prescl;
   Int_t           L1Tech_BSC_HighMultiplicity_v0;
   Int_t           L1Tech_BSC_HighMultiplicity_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam1_inner_v0;
   Int_t           L1Tech_BSC_halo_beam1_inner_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam1_outer_v0;
   Int_t           L1Tech_BSC_halo_beam1_outer_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam2_inner_v0;
   Int_t           L1Tech_BSC_halo_beam2_inner_v0_Prescl;
   Int_t           L1Tech_BSC_halo_beam2_outer_v0;
   Int_t           L1Tech_BSC_halo_beam2_outer_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_OR_v0;
   Int_t           L1Tech_BSC_minBias_OR_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_inner_threshold1_v0;
   Int_t           L1Tech_BSC_minBias_inner_threshold1_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_inner_threshold2_v0;
   Int_t           L1Tech_BSC_minBias_inner_threshold2_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_threshold1_v0;
   Int_t           L1Tech_BSC_minBias_threshold1_v0_Prescl;
   Int_t           L1Tech_BSC_minBias_threshold2_v0;
   Int_t           L1Tech_BSC_minBias_threshold2_v0_Prescl;
   Int_t           L1Tech_BSC_splash_beam1_v0;
   Int_t           L1Tech_BSC_splash_beam1_v0_Prescl;
   Int_t           L1Tech_BSC_splash_beam2_v0;
   Int_t           L1Tech_BSC_splash_beam2_v0_Prescl;
   Int_t           L1Tech_CASTOR_HaloMuon_v0;
   Int_t           L1Tech_CASTOR_HaloMuon_v0_Prescl;
   Int_t           L1Tech_DT_GlobalOR_v0;
   Int_t           L1Tech_DT_GlobalOR_v0_Prescl;
   Int_t           L1Tech_FSC_St1Sect45_down_v0;
   Int_t           L1Tech_FSC_St1Sect45_down_v0_Prescl;
   Int_t           L1Tech_FSC_St1Sect45_upp_v0;
   Int_t           L1Tech_FSC_St1Sect45_upp_v0_Prescl;
   Int_t           L1Tech_FSC_St1Sect56_down_v0;
   Int_t           L1Tech_FSC_St1Sect56_down_v0_Prescl;
   Int_t           L1Tech_FSC_St1Sect56_upp_v0;
   Int_t           L1Tech_FSC_St1Sect56_upp_v0_Prescl;
   Int_t           L1Tech_FSC_St2Sect45_down_v0;
   Int_t           L1Tech_FSC_St2Sect45_down_v0_Prescl;
   Int_t           L1Tech_FSC_St2Sect45_upp_v0;
   Int_t           L1Tech_FSC_St2Sect45_upp_v0_Prescl;
   Int_t           L1Tech_FSC_St2Sect56_down_v0;
   Int_t           L1Tech_FSC_St2Sect56_down_v0_Prescl;
   Int_t           L1Tech_FSC_St2Sect56_upp_v0;
   Int_t           L1Tech_FSC_St2Sect56_upp_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_downLeft_v0;
   Int_t           L1Tech_FSC_St3Sect45_downLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_downRight_v0;
   Int_t           L1Tech_FSC_St3Sect45_downRight_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_uppLeft_v0;
   Int_t           L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect45_uppRight_v0;
   Int_t           L1Tech_FSC_St3Sect45_uppRight_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect56_downLeft_v0;
   Int_t           L1Tech_FSC_St3Sect56_downLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect56_downRight_v0;
   Int_t           L1Tech_FSC_St3Sect56_downRight_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect56_uppLeft_v0;
   Int_t           L1Tech_FSC_St3Sect56_uppLeft_v0_Prescl;
   Int_t           L1Tech_FSC_St3Sect56_uppRight_v0;
   Int_t           L1Tech_FSC_St3Sect56_uppRight_v0_Prescl;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0;
   Int_t           L1Tech_HCAL_HBHE_totalOR_v0_Prescl;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1;
   Int_t           L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2;
   Int_t           L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;
   Int_t           L1Tech_HCAL_HF_single_channel_v0;
   Int_t           L1Tech_HCAL_HF_single_channel_v0_Prescl;
   Int_t           L1Tech_HCAL_HO_totalOR_v0;
   Int_t           L1Tech_HCAL_HO_totalOR_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RB0_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBminus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBminus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_RBst1_collisions_v0;
   Int_t           L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0;
   Int_t           L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;

   // List of branches
   TBranch        *b_NL1IsolEm;   //!
   TBranch        *b_L1IsolEmEt;   //!
   TBranch        *b_L1IsolEmE;   //!
   TBranch        *b_L1IsolEmEta;   //!
   TBranch        *b_L1IsolEmPhi;   //!
   TBranch        *b_NL1NIsolEm;   //!
   TBranch        *b_L1NIsolEmEt;   //!
   TBranch        *b_L1NIsolEmE;   //!
   TBranch        *b_L1NIsolEmEta;   //!
   TBranch        *b_L1NIsolEmPhi;   //!
   TBranch        *b_NL1Mu;   //!
   TBranch        *b_L1MuPt;   //!
   TBranch        *b_L1MuE;   //!
   TBranch        *b_L1MuEta;   //!
   TBranch        *b_L1MuPhi;   //!
   TBranch        *b_L1MuIsol;   //!
   TBranch        *b_L1MuMip;   //!
   TBranch        *b_L1MuFor;   //!
   TBranch        *b_L1MuRPC;   //!
   TBranch        *b_L1MuQal;   //!
   TBranch        *b_L1MuChg;   //!
   TBranch        *b_NL1CenJet;   //!
   TBranch        *b_L1CenJetEt;   //!
   TBranch        *b_L1CenJetE;   //!
   TBranch        *b_L1CenJetEta;   //!
   TBranch        *b_L1CenJetPhi;   //!
   TBranch        *b_NL1ForJet;   //!
   TBranch        *b_L1ForJetEt;   //!
   TBranch        *b_L1ForJetE;   //!
   TBranch        *b_L1ForJetEta;   //!
   TBranch        *b_L1ForJetPhi;   //!
   TBranch        *b_NL1Tau;   //!
   TBranch        *b_L1TauEt;   //!
   TBranch        *b_L1TauE;   //!
   TBranch        *b_L1TauEta;   //!
   TBranch        *b_L1TauPhi;   //!
   TBranch        *b_L1Met;   //!
   TBranch        *b_L1MetPhi;   //!
   TBranch        *b_L1EtTot;   //!
   TBranch        *b_L1Mht;   //!
   TBranch        *b_L1MhtPhi;   //!
   TBranch        *b_L1EtHad;   //!
   TBranch        *b_L1HfRing1EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing2EtSumPositiveEta;   //!
   TBranch        *b_L1HfRing1EtSumNegativeEta;   //!
   TBranch        *b_L1HfRing2EtSumNegativeEta;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing1;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing1;   //!
   TBranch        *b_L1HfTowerCountPositiveEtaRing2;   //!
   TBranch        *b_L1HfTowerCountNegativeEtaRing2;   //!
   TBranch        *b_NMCpart;   //!
   TBranch        *b_MCpid;   //!
   TBranch        *b_MCstatus;   //!
   TBranch        *b_MCvtxX;   //!
   TBranch        *b_MCvtxY;   //!
   TBranch        *b_MCvtxZ;   //!
   TBranch        *b_MCpt;   //!
   TBranch        *b_MCeta;   //!
   TBranch        *b_MCphi;   //!
   TBranch        *b_MCPtHat;   //!
   TBranch        *b_MCmu3;   //!
   TBranch        *b_MCel3;   //!
   TBranch        *b_MCbb;   //!
   TBranch        *b_MCab;   //!
   TBranch        *b_MCWenu;   //!
   TBranch        *b_MCmunu;   //!
   TBranch        *b_MCZee;   //!
   TBranch        *b_MCZmumu;   //!
   TBranch        *b_MCptEleMax;   //!
   TBranch        *b_MCptMuMax;   //!
   TBranch        *b_NPUTrueBX0;   //!
   TBranch        *b_NPUgenBX0;   //!
   TBranch        *b_NVrtx;   //!
   TBranch        *b_recoVrtXOffline0;   //!
   TBranch        *b_recoVrtYOffline0;   //!
   TBranch        *b_recoVrtZOffline0;   //!
   TBranch        *b_recoVrtNtrkOffline0;   //!
   TBranch        *b_recoVrtChi2Offline0;   //!
   TBranch        *b_recoVrtNdofOffline0;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Bx;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_AvgInstDelLumi;   //!
   TBranch        *b_HLT_HIMET220_v1;   //!
   TBranch        *b_HLT_HIMET220_v1_Prescl;   //!
   TBranch        *b_HLT_HIMET200_v1;   //!
   TBranch        *b_HLT_HIMET200_v1_Prescl;   //!
   TBranch        *b_HLT_HIMET120_v1;   //!
   TBranch        *b_HLT_HIMET120_v1_Prescl;   //!
   TBranch        *b_HLT_HIPhysics_v1;   //!
   TBranch        *b_HLT_HIPhysics_v1_Prescl;   //!
   TBranch        *b_HLT_HIDTCalibration_v1;   //!
   TBranch        *b_HLT_HIDTCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_HIEcalCalibration_v1;   //!
   TBranch        *b_HLT_HIEcalCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_HIHcalCalibration_v1;   //!
   TBranch        *b_HLT_HIHcalCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_HIZeroBias_v1;   //!
   TBranch        *b_HLT_HIZeroBias_v1_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v1;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v1_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v1;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v1;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHF_v1;   //!
   TBranch        *b_HLT_HIMinBiasHF_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v1;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v1;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v1;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIBptxXOR_v1;   //!
   TBranch        *b_HLT_HIBptxXOR_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v1;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v2;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v2_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v2;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_v2;   //!
   TBranch        *b_HLT_HIL2Mu3_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v2;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_v2;   //!
   TBranch        *b_HLT_HIL2Mu7_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v2;   //!
   TBranch        *b_HLT_HIL2Mu15_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v2;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v2;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v2;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_v2;   //!
   TBranch        *b_HLT_HIL3Mu3_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_v1;   //!
   TBranch        *b_HLT_HISinglePhoton15_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_v2;   //!
   TBranch        *b_HLT_HISinglePhoton20_v2_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_v2;   //!
   TBranch        *b_HLT_HISinglePhoton30_v2_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_v2;   //!
   TBranch        *b_HLT_HISinglePhoton40_v2_Prescl;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v1;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v1_Prescl;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v1;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v1_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v1;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v1_Prescl;   //!
   TBranch        *b_HLT_HIJet55_v1;   //!
   TBranch        *b_HLT_HIJet55_v1_Prescl;   //!
   TBranch        *b_HLT_HIJet65_v1;   //!
   TBranch        *b_HLT_HIJet65_v1_Prescl;   //!
   TBranch        *b_HLT_HIJet80_v1;   //!
   TBranch        *b_HLT_HIJet80_v1_Prescl;   //!
   TBranch        *b_HLT_HIJet80_v7;   //!
   TBranch        *b_HLT_HIJet80_v7_Prescl;   //!
   TBranch        *b_HLT_HIJet95_v1;   //!
   TBranch        *b_HLT_HIJet95_v1_Prescl;   //!
   TBranch        *b_HLT_HIDiJet55_v1;   //!
   TBranch        *b_HLT_HIDiJet55_v1_Prescl;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v1;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v1_Prescl;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v1;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v1_Prescl;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v1;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl;   //!
   TBranch        *b_HLT_HIActivityHF_Coincidence3_v1;   //!
   TBranch        *b_HLT_HIActivityHF_Coincidence3_v1_Prescl;   //!
   TBranch        *b_HLT_HIActivityHF_Single3_v1;   //!
   TBranch        *b_HLT_HIActivityHF_Single3_v1_Prescl;   //!
   TBranch        *b_HLT_HIClusterVertexCompatibility_v1;   //!
   TBranch        *b_HLT_HIClusterVertexCompatibility_v1_Prescl;   //!
   TBranch        *b_HLT_HICentralityVeto_v1;   //!
   TBranch        *b_HLT_HICentralityVeto_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v1;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v1;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v1;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v1;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v1;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v1;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v1;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v1_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v1;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v1_Prescl;   //!
   TBranch        *b_HLT_HIRandom_v1;   //!
   TBranch        *b_HLT_HIRandom_v1_Prescl;   //!
   TBranch        *b_HLT_HIUCC010_v2;   //!
   TBranch        *b_HLT_HIUCC010_v2_Prescl;   //!
   TBranch        *b_HLT_HIUCC015_v2;   //!
   TBranch        *b_HLT_HIUCC015_v2_Prescl;   //!
   TBranch        *b_HLT_HICentral10_v2;   //!
   TBranch        *b_HLT_HICentral10_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFinalPath_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v1;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v1;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v1;   //!
   TBranch        *b_HLT_HIL2Mu15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v1;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_v1;   //!
   TBranch        *b_HLT_HIL2Mu3_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_v1;   //!
   TBranch        *b_HLT_HIL2Mu7_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_v1;   //!
   TBranch        *b_HLT_HIL3Mu3_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_v1;   //!
   TBranch        *b_HLT_HISinglePhoton20_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_v1;   //!
   TBranch        *b_HLT_HISinglePhoton30_v1_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_v1;   //!
   TBranch        *b_HLT_HISinglePhoton40_v1_Prescl;   //!
   TBranch        *b_HLT_HICentral10_v1;   //!
   TBranch        *b_HLT_HICentral10_v1_Prescl;   //!
   TBranch        *b_HLT_HIUCC010_v1;   //!
   TBranch        *b_HLT_HIUCC010_v1_Prescl;   //!
   TBranch        *b_HLT_HIUCC015_v1;   //!
   TBranch        *b_HLT_HIUCC015_v1_Prescl;   //!
   TBranch        *b_HLT_HIDiJet55_v2;   //!
   TBranch        *b_HLT_HIDiJet55_v2_Prescl;   //!
   TBranch        *b_HLT_PAJet80_NoJetID_v1;
   TBranch        *b_HLT_HIDoublePhoton10_v2;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v2_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v2;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v2_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v2;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v2;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v2;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v2;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v2;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v2;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v2;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v2;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v2_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v2;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet55_v2;   //!
   TBranch        *b_HLT_HIJet55_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v2;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet65_v2;   //!
   TBranch        *b_HLT_HIJet65_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet80_v2;   //!
   TBranch        *b_HLT_HIJet80_v2_Prescl;   //!
   TBranch        *b_HLT_HIJet95_v2;   //!
   TBranch        *b_HLT_HIJet95_v2_Prescl;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v2;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v2_Prescl;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v2;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl;   //!
   TBranch        *b_HLT_HIMET120_v2;   //!
   TBranch        *b_HLT_HIMET120_v2_Prescl;   //!
   TBranch        *b_HLT_HIMET200_v2;   //!
   TBranch        *b_HLT_HIMET200_v2_Prescl;   //!
   TBranch        *b_HLT_HIMET220_v2;   //!
   TBranch        *b_HLT_HIMET220_v2_Prescl;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v2;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v2_Prescl;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v2;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v2_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_v2;   //!
   TBranch        *b_HLT_HISinglePhoton15_v2_Prescl;   //!
   TBranch        *b_HLT_HIBptxXOR_v2;   //!
   TBranch        *b_HLT_HIBptxXOR_v2_Prescl;   //!
   TBranch        *b_HLT_HICentralityVeto_v2;   //!
   TBranch        *b_HLT_HICentralityVeto_v2_Prescl;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v2;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v2;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v2;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHF_v2;   //!
   TBranch        *b_HLT_HIMinBiasHF_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v2;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v2;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v2;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v2_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIPhysics_v2;   //!
   TBranch        *b_HLT_HIPhysics_v2_Prescl;   //!
   TBranch        *b_HLT_HIRandom_v2;   //!
   TBranch        *b_HLT_HIRandom_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v2;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v2;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v2_Prescl;   //!
   TBranch        *b_HLT_HIZeroBias_v2;   //!
   TBranch        *b_HLT_HIZeroBias_v2_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v3;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v3_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v3;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v3;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v3;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v3;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v3;   //!
   TBranch        *b_HLT_HIL2Mu15_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v3;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_v3;   //!
   TBranch        *b_HLT_HIL2Mu3_v3_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_v3;   //!
   TBranch        *b_HLT_HIL2Mu7_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v3;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v3_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_v3;   //!
   TBranch        *b_HLT_HIL3Mu3_v3_Prescl;   //!
   TBranch        *b_HLT_HIDiJet55_v3;   //!
   TBranch        *b_HLT_HIDiJet55_v3_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v3;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v3_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v3;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v3_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v3;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v3;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v3;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v3;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v3;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v3;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v3;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v3;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v3_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v3;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet55_v3;   //!
   TBranch        *b_HLT_HIJet55_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v3;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet65_v3;   //!
   TBranch        *b_HLT_HIJet65_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet80_v3;   //!
   TBranch        *b_HLT_HIJet80_v3_Prescl;   //!
   TBranch        *b_HLT_HIJet95_v3;   //!
   TBranch        *b_HLT_HIJet95_v3_Prescl;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v3;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v3_Prescl;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v3;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl;   //!
   TBranch        *b_HLT_HIMET120_v3;   //!
   TBranch        *b_HLT_HIMET120_v3_Prescl;   //!
   TBranch        *b_HLT_HIMET200_v3;   //!
   TBranch        *b_HLT_HIMET200_v3_Prescl;   //!
   TBranch        *b_HLT_HIMET220_v3;   //!
   TBranch        *b_HLT_HIMET220_v3_Prescl;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v3;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v3_Prescl;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v3;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v3_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_v3;   //!
   TBranch        *b_HLT_HISinglePhoton15_v3_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_v3;   //!
   TBranch        *b_HLT_HISinglePhoton20_v3_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_v3;   //!
   TBranch        *b_HLT_HISinglePhoton30_v3_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_v3;   //!
   TBranch        *b_HLT_HISinglePhoton40_v3_Prescl;   //!
   TBranch        *b_HLT_HIBptxXOR_v3;   //!
   TBranch        *b_HLT_HIBptxXOR_v3_Prescl;   //!
   TBranch        *b_HLT_HICentral10_v3;   //!
   TBranch        *b_HLT_HICentral10_v3_Prescl;   //!
   TBranch        *b_HLT_HICentralityVeto_v3;   //!
   TBranch        *b_HLT_HICentralityVeto_v3_Prescl;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v3;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v3;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v3;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHF_v3;   //!
   TBranch        *b_HLT_HIMinBiasHF_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v3;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v3;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v3;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v3_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIPhysics_v3;   //!
   TBranch        *b_HLT_HIPhysics_v3_Prescl;   //!
   TBranch        *b_HLT_HIRandom_v3;   //!
   TBranch        *b_HLT_HIRandom_v3_Prescl;   //!
   TBranch        *b_HLT_HIUCC010_v3;   //!
   TBranch        *b_HLT_HIUCC010_v3_Prescl;   //!
   TBranch        *b_HLT_HIUCC015_v3;   //!
   TBranch        *b_HLT_HIUCC015_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v3;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v3;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v3_Prescl;   //!
   TBranch        *b_HLT_HIZeroBias_v3;   //!
   TBranch        *b_HLT_HIZeroBias_v3_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v4;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_HighQ_v4_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v4;   //!
   TBranch        *b_HLT_HIL1DoubleMuOpen_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v4;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v4;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v4;   //!
   TBranch        *b_HLT_HIL2DoubleMu3_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v4;   //!
   TBranch        *b_HLT_HIL2Mu15_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v4;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_v4;   //!
   TBranch        *b_HLT_HIL2Mu3_v4_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_v4;   //!
   TBranch        *b_HLT_HIL2Mu7_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v4;   //!
   TBranch        *b_HLT_HIL3DoubleMuOpen_v4_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_v4;   //!
   TBranch        *b_HLT_HIL3Mu3_v4_Prescl;   //!
   TBranch        *b_HLT_HIDiJet55_v4;   //!
   TBranch        *b_HLT_HIDiJet55_v4_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v4;   //!
   TBranch        *b_HLT_HIDoublePhoton10_v4_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v4;   //!
   TBranch        *b_HLT_HIDoublePhoton15_v4_Prescl;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v4;   //!
   TBranch        *b_HLT_HIDoublePhoton20_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v4;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Central_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v4;   //!
   TBranch        *b_HLT_HIFullTrack12_L1Peripheral_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v4;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Central_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v4;   //!
   TBranch        *b_HLT_HIFullTrack14_L1Peripheral_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v4;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Central_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v4;   //!
   TBranch        *b_HLT_HIFullTrack20_L1Peripheral_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v4;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Central_v4_Prescl;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v4;   //!
   TBranch        *b_HLT_HIFullTrack25_L1Peripheral_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet55_v4;   //!
   TBranch        *b_HLT_HIJet55_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v4;   //!
   TBranch        *b_HLT_HIJet65_Jet55_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet65_v4;   //!
   TBranch        *b_HLT_HIJet65_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet80_v4;   //!
   TBranch        *b_HLT_HIJet80_v4_Prescl;   //!
   TBranch        *b_HLT_HIJet95_v4;   //!
   TBranch        *b_HLT_HIJet95_v4_Prescl;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v4;   //!
   TBranch        *b_HLT_HIJetE30_NoBPTX_v4_Prescl;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v4;   //!
   TBranch        *b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl;   //!
   TBranch        *b_HLT_HIMET120_v4;   //!
   TBranch        *b_HLT_HIMET120_v4_Prescl;   //!
   TBranch        *b_HLT_HIMET200_v4;   //!
   TBranch        *b_HLT_HIMET200_v4_Prescl;   //!
   TBranch        *b_HLT_HIMET220_v4;   //!
   TBranch        *b_HLT_HIMET220_v4_Prescl;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v4;   //!
   TBranch        *b_HLT_HIPhoton10_Photon15_v4_Prescl;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v4;   //!
   TBranch        *b_HLT_HIPhoton15_Photon20_v4_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton15_v4;   //!
   TBranch        *b_HLT_HISinglePhoton15_v4_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton20_v4;   //!
   TBranch        *b_HLT_HISinglePhoton20_v4_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton30_v4;   //!
   TBranch        *b_HLT_HISinglePhoton30_v4_Prescl;   //!
   TBranch        *b_HLT_HISinglePhoton40_v4;   //!
   TBranch        *b_HLT_HISinglePhoton40_v4_Prescl;   //!
   TBranch        *b_HLT_HIBptxXOR_v4;   //!
   TBranch        *b_HLT_HIBptxXOR_v4_Prescl;   //!
   TBranch        *b_HLT_HICentral10_v4;   //!
   TBranch        *b_HLT_HICentral10_v4_Prescl;   //!
   TBranch        *b_HLT_HICentralityVeto_v4;   //!
   TBranch        *b_HLT_HICentralityVeto_v4_Prescl;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v4;   //!
   TBranch        *b_HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v4;   //!
   TBranch        *b_HLT_HIMinBiasBSC_OR_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v4;   //!
   TBranch        *b_HLT_HIMinBiasBSC_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHF_v4;   //!
   TBranch        *b_HLT_HIMinBiasHF_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v4;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v4;   //!
   TBranch        *b_HLT_HIMinBiasHf_OR_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIMinBiasPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v4;   //!
   TBranch        *b_HLT_HIMinBiasZDC_Calo_v4_Prescl;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIPhysics_v4;   //!
   TBranch        *b_HLT_HIPhysics_v4_Prescl;   //!
   TBranch        *b_HLT_HIRandom_v4;   //!
   TBranch        *b_HLT_HIRandom_v4_Prescl;   //!
   TBranch        *b_HLT_HIUCC010_v4;   //!
   TBranch        *b_HLT_HIUCC010_v4_Prescl;   //!
   TBranch        *b_HLT_HIUCC015_v4;   //!
   TBranch        *b_HLT_HIUCC015_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v4;   //!
   TBranch        *b_HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v4;   //!
   TBranch        *b_HLT_HIZeroBiasXOR_v4_Prescl;   //!
   TBranch        *b_HLT_HIZeroBias_v4;   //!
   TBranch        *b_HLT_HIZeroBias_v4_Prescl;   //!
   TBranch        *b_L1_BptxMinus;   //!
   TBranch        *b_L1_BptxMinus_Prescl;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus;   //!
   TBranch        *b_L1_BptxMinus_NotBptxPlus_Prescl;   //!
   TBranch        *b_L1_BptxPlus;   //!
   TBranch        *b_L1_BptxPlus_Prescl;   //!
   TBranch        *b_L1_BptxPlusANDMinus;   //!
   TBranch        *b_L1_BptxPlusANDMinus_Prescl;   //!
   TBranch        *b_L1_BptxPlusORMinus;   //!
   TBranch        *b_L1_BptxPlusORMinus_Prescl;   //!
   TBranch        *b_L1_BptxPlusORMinus_instance1;   //!
   TBranch        *b_L1_BptxPlusORMinus_instance1_Prescl;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus;   //!
   TBranch        *b_L1_BptxPlus_NotBptxMinus_Prescl;   //!
   TBranch        *b_L1_BptxXOR;   //!
   TBranch        *b_L1_BptxXOR_Prescl;   //!
   TBranch        *b_L1_BptxXOR_BscMinBiasOR;   //!
   TBranch        *b_L1_BptxXOR_BscMinBiasOR_Prescl;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus;   //!
   TBranch        *b_L1_Bsc2Minus_BptxMinus_Prescl;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus;   //!
   TBranch        *b_L1_Bsc2Plus_BptxPlus_Prescl;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1_Prescl;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold1_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold2_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasInnerThreshold2_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusORMinus;   //!
   TBranch        *b_L1_BscMinBiasOR_BptxPlusORMinus_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold1;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasThreshold1_BptxAND_Prescl;   //!
   TBranch        *b_L1_BscMinBiasThreshold2_BptxAND;   //!
   TBranch        *b_L1_BscMinBiasThreshold2_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleEG10;   //!
   TBranch        *b_L1_DoubleEG10_Prescl;   //!
   TBranch        *b_L1_DoubleEG3;   //!
   TBranch        *b_L1_DoubleEG3_Prescl;   //!
   TBranch        *b_L1_DoubleEG5_BptxAND;   //!
   TBranch        *b_L1_DoubleEG5_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleEG8_BptxAND;   //!
   TBranch        *b_L1_DoubleEG8_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleForJet32_EtaOpp;   //!
   TBranch        *b_L1_DoubleForJet32_EtaOpp_Prescl;   //!
   TBranch        *b_L1_DoubleForJet44_EtaOpp;   //!
   TBranch        *b_L1_DoubleForJet44_EtaOpp_Prescl;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing1_P1N1;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing1_P1N1_Prescl;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing2_P1N1;   //!
   TBranch        *b_L1_DoubleHfBitCountsRing2_P1N1_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P200N200;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P200N200_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P4N4;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing1_P4N4_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P200N200;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P200N200_Prescl;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P4N4;   //!
   TBranch        *b_L1_DoubleHfRingEtSumsRing2_P4N4_Prescl;   //!
   TBranch        *b_L1_DoubleMu3_BptxAND;   //!
   TBranch        *b_L1_DoubleMu3_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND_Prescl;   //!
   TBranch        *b_L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND;   //!
   TBranch        *b_L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;   //!
   TBranch        *b_L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;   //!
   TBranch        *b_L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;   //!
   TBranch        *b_L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND;   //!
   TBranch        *b_L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;   //!
   TBranch        *b_L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;   //!
   TBranch        *b_L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETM30_BptxAND;   //!
   TBranch        *b_L1_ETM30_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETM50_BptxAND;   //!
   TBranch        *b_L1_ETM50_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT100_BptxAND;   //!
   TBranch        *b_L1_ETT100_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT140_BptxAND;   //!
   TBranch        *b_L1_ETT140_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT2000;   //!
   TBranch        *b_L1_ETT2000_Prescl;   //!
   TBranch        *b_L1_ETT220_BptxAND;   //!
   TBranch        *b_L1_ETT220_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT70_BptxAND;   //!
   TBranch        *b_L1_ETT70_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT800;   //!
   TBranch        *b_L1_ETT800_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND;   //!
   TBranch        *b_L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincidencePm;   //!
   TBranch        *b_L1_HcalHfCoincidencePm_Prescl;   //!
   TBranch        *b_L1_HcalHfCoincidencePm_BptxAND;   //!
   TBranch        *b_L1_HcalHfCoincidencePm_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHfMmOrPpOrPm;   //!
   TBranch        *b_L1_HcalHfMmOrPpOrPm_Prescl;   //!
   TBranch        *b_L1_HcalHfMmOrPpOrPm_BptxAND;   //!
   TBranch        *b_L1_HcalHfMmOrPpOrPm_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHfMmpOrMpp;   //!
   TBranch        *b_L1_HcalHfMmpOrMpp_Prescl;   //!
   TBranch        *b_L1_HcalHfMmpOrMpp_BptxAND;   //!
   TBranch        *b_L1_HcalHfMmpOrMpp_BptxAND_Prescl;   //!
   TBranch        *b_L1_HcalHoTotalOR;   //!
   TBranch        *b_L1_HcalHoTotalOR_Prescl;   //!
   TBranch        *b_L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND;   //!
   TBranch        *b_L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl;   //!
   TBranch        *b_L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND;   //!
   TBranch        *b_L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasOR;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasThreshold1;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasThreshold2;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfCoincidencePm;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfMmpOrMpp;   //!
   TBranch        *b_L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasInnerThreshold1;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasInnerThreshold2;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasOR;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasOR_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasThreshold1;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasThreshold1_Prescl;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasThreshold2;   //!
   TBranch        *b_L1_NotBsc2_BscMinBiasThreshold2_Prescl;   //!
   TBranch        *b_L1_NotBsc2_HcalHfCoincidencePm;   //!
   TBranch        *b_L1_NotBsc2_HcalHfCoincidencePm_Prescl;   //!
   TBranch        *b_L1_NotBsc2_HcalHfMmOrPpOrPm;   //!
   TBranch        *b_L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl;   //!
   TBranch        *b_L1_NotBsc2_HcalHfMmpOrMpp;   //!
   TBranch        *b_L1_NotBsc2_HcalHfMmpOrMpp_Prescl;   //!
   TBranch        *b_L1_SingleEG12;   //!
   TBranch        *b_L1_SingleEG12_Prescl;   //!
   TBranch        *b_L1_SingleEG15;   //!
   TBranch        *b_L1_SingleEG15_Prescl;   //!
   TBranch        *b_L1_SingleEG18;   //!
   TBranch        *b_L1_SingleEG18_Prescl;   //!
   TBranch        *b_L1_SingleEG20;   //!
   TBranch        *b_L1_SingleEG20_Prescl;   //!
   TBranch        *b_L1_SingleEG22;   //!
   TBranch        *b_L1_SingleEG22_Prescl;   //!
   TBranch        *b_L1_SingleEG2_BptxAND;   //!
   TBranch        *b_L1_SingleEG2_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG30;   //!
   TBranch        *b_L1_SingleEG30_Prescl;   //!
   TBranch        *b_L1_SingleEG5_BptxAND;   //!
   TBranch        *b_L1_SingleEG5_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG8_BptxAND;   //!
   TBranch        *b_L1_SingleEG8_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet128;   //!
   TBranch        *b_L1_SingleJet128_Prescl;   //!
   TBranch        *b_L1_SingleJet128_BptxAND;   //!
   TBranch        *b_L1_SingleJet128_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet128_Central_BptxAND;   //!
   TBranch        *b_L1_SingleJet128_Central_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet16;   //!
   TBranch        *b_L1_SingleJet16_Prescl;   //!
   TBranch        *b_L1_SingleJet16_Central_NotETT140_BptxAND;   //!
   TBranch        *b_L1_SingleJet16_Central_NotETT140_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo;   //!
   TBranch        *b_L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo;   //!
   TBranch        *b_L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleJet36_BptxAND;   //!
   TBranch        *b_L1_SingleJet36_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet52_BptxAND;   //!
   TBranch        *b_L1_SingleJet52_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet52_Central_BptxAND;   //!
   TBranch        *b_L1_SingleJet52_Central_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet68_BptxAND;   //!
   TBranch        *b_L1_SingleJet68_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet80_BptxAND;   //!
   TBranch        *b_L1_SingleJet80_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet80_Central;   //!
   TBranch        *b_L1_SingleJet80_Central_Prescl;   //!
   TBranch        *b_L1_SingleJet80_Central_BptxAND;   //!
   TBranch        *b_L1_SingleJet80_Central_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet92;   //!
   TBranch        *b_L1_SingleJet92_Prescl;   //!
   TBranch        *b_L1_SingleJet92_BptxAND;   //!
   TBranch        *b_L1_SingleJet92_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet92_Central;   //!
   TBranch        *b_L1_SingleJet92_Central_Prescl;   //!
   TBranch        *b_L1_SingleJet92_Central_BptxAND;   //!
   TBranch        *b_L1_SingleJet92_Central_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu0_BptxAND;   //!
   TBranch        *b_L1_SingleMu0_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu10;   //!
   TBranch        *b_L1_SingleMu10_Prescl;   //!
   TBranch        *b_L1_SingleMu14;   //!
   TBranch        *b_L1_SingleMu14_Prescl;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu20_Prescl;   //!
   TBranch        *b_L1_SingleMu3_BptxAND;   //!
   TBranch        *b_L1_SingleMu3_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu5_BptxAND;   //!
   TBranch        *b_L1_SingleMu5_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMu7_Prescl;   //!
   TBranch        *b_L1_SingleMuBeamHalo;   //!
   TBranch        *b_L1_SingleMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen_BptxAND;   //!
   TBranch        *b_L1_SingleMuOpen_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZdcCaloMinus;   //!
   TBranch        *b_L1_ZdcCaloMinus_Prescl;   //!
   TBranch        *b_L1_ZdcCaloMinus_BptxAND;   //!
   TBranch        *b_L1_ZdcCaloMinus_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZdcCaloPlus;   //!
   TBranch        *b_L1_ZdcCaloPlus_Prescl;   //!
   TBranch        *b_L1_ZdcCaloPlus_BptxAND;   //!
   TBranch        *b_L1_ZdcCaloPlus_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZdcCaloPlus_ZdcCaloMinus;   //!
   TBranch        *b_L1_ZdcCaloPlus_ZdcCaloMinus_Prescl;   //!
   TBranch        *b_L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND;   //!
   TBranch        *b_L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0;   //!
   TBranch        *b_L1Tech_BPTX_PreBPTX_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0;   //!
   TBranch        *b_L1Tech_BPTX_plus_OR_minus_v0_Prescl;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0;   //!
   TBranch        *b_L1Tech_BPTX_quiet_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_HighMultiplicity_v0;   //!
   TBranch        *b_L1Tech_BSC_HighMultiplicity_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_inner_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_inner_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_outer_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam1_outer_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_inner_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_inner_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_outer_v0;   //!
   TBranch        *b_L1Tech_BSC_halo_beam2_outer_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_OR_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_OR_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold1_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold2_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_inner_threshold2_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold1_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold2_v0;   //!
   TBranch        *b_L1Tech_BSC_minBias_threshold2_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_splash_beam1_v0;   //!
   TBranch        *b_L1Tech_BSC_splash_beam1_v0_Prescl;   //!
   TBranch        *b_L1Tech_BSC_splash_beam2_v0;   //!
   TBranch        *b_L1Tech_BSC_splash_beam2_v0_Prescl;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0;   //!
   TBranch        *b_L1Tech_CASTOR_HaloMuon_v0_Prescl;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0;   //!
   TBranch        *b_L1Tech_DT_GlobalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St1Sect45_down_v0;   //!
   TBranch        *b_L1Tech_FSC_St1Sect45_down_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St1Sect45_upp_v0;   //!
   TBranch        *b_L1Tech_FSC_St1Sect45_upp_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St1Sect56_down_v0;   //!
   TBranch        *b_L1Tech_FSC_St1Sect56_down_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St1Sect56_upp_v0;   //!
   TBranch        *b_L1Tech_FSC_St1Sect56_upp_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St2Sect45_down_v0;   //!
   TBranch        *b_L1Tech_FSC_St2Sect45_down_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St2Sect45_upp_v0;   //!
   TBranch        *b_L1Tech_FSC_St2Sect45_upp_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St2Sect56_down_v0;   //!
   TBranch        *b_L1Tech_FSC_St2Sect56_down_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St2Sect56_upp_v0;   //!
   TBranch        *b_L1Tech_FSC_St2Sect56_upp_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_downRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect45_uppRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_downLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_downLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_downRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_downRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_uppLeft_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_uppLeft_v0_Prescl;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_uppRight_v0;   //!
   TBranch        *b_L1Tech_FSC_St3Sect56_uppRight_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1;   //!
   TBranch        *b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2;   //!
   TBranch        *b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0;   //!
   TBranch        *b_L1Tech_HCAL_HF_single_channel_v0_Prescl;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0;   //!
   TBranch        *b_L1Tech_HCAL_HO_totalOR_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RB0_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBst1_collisions_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0;   //!
   TBranch        *b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl;   //!

   HLT(TTree *tree=0);
   virtual ~HLT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HLT_cxx
HLT::HLT(TTree *tree)
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

HLT::~HLT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HLT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HLT::LoadTree(Long64_t entry)
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

void HLT::Init(TTree *tree)
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

   fChain->SetBranchAddress("NL1IsolEm", &NL1IsolEm, &b_NL1IsolEm);
   fChain->SetBranchAddress("L1IsolEmEt", L1IsolEmEt, &b_L1IsolEmEt);
   fChain->SetBranchAddress("L1IsolEmE", L1IsolEmE, &b_L1IsolEmE);
   fChain->SetBranchAddress("L1IsolEmEta", L1IsolEmEta, &b_L1IsolEmEta);
   fChain->SetBranchAddress("L1IsolEmPhi", L1IsolEmPhi, &b_L1IsolEmPhi);
   fChain->SetBranchAddress("NL1NIsolEm", &NL1NIsolEm, &b_NL1NIsolEm);
   fChain->SetBranchAddress("L1NIsolEmEt", L1NIsolEmEt, &b_L1NIsolEmEt);
   fChain->SetBranchAddress("L1NIsolEmE", L1NIsolEmE, &b_L1NIsolEmE);
   fChain->SetBranchAddress("L1NIsolEmEta", L1NIsolEmEta, &b_L1NIsolEmEta);
   fChain->SetBranchAddress("L1NIsolEmPhi", L1NIsolEmPhi, &b_L1NIsolEmPhi);
   fChain->SetBranchAddress("NL1Mu", &NL1Mu, &b_NL1Mu);
   fChain->SetBranchAddress("L1MuPt", L1MuPt, &b_L1MuPt);
   fChain->SetBranchAddress("L1MuE", L1MuE, &b_L1MuE);
   fChain->SetBranchAddress("L1MuEta", L1MuEta, &b_L1MuEta);
   fChain->SetBranchAddress("L1MuPhi", L1MuPhi, &b_L1MuPhi);
   fChain->SetBranchAddress("L1MuIsol", L1MuIsol, &b_L1MuIsol);
   fChain->SetBranchAddress("L1MuMip", L1MuMip, &b_L1MuMip);
   fChain->SetBranchAddress("L1MuFor", L1MuFor, &b_L1MuFor);
   fChain->SetBranchAddress("L1MuRPC", L1MuRPC, &b_L1MuRPC);
   fChain->SetBranchAddress("L1MuQal", L1MuQal, &b_L1MuQal);
   fChain->SetBranchAddress("L1MuChg", L1MuChg, &b_L1MuChg);
   fChain->SetBranchAddress("NL1CenJet", &NL1CenJet, &b_NL1CenJet);
   fChain->SetBranchAddress("L1CenJetEt", L1CenJetEt, &b_L1CenJetEt);
   fChain->SetBranchAddress("L1CenJetE", L1CenJetE, &b_L1CenJetE);
   fChain->SetBranchAddress("L1CenJetEta", L1CenJetEta, &b_L1CenJetEta);
   fChain->SetBranchAddress("L1CenJetPhi", L1CenJetPhi, &b_L1CenJetPhi);
   fChain->SetBranchAddress("NL1ForJet", &NL1ForJet, &b_NL1ForJet);
   fChain->SetBranchAddress("L1ForJetEt", L1ForJetEt, &b_L1ForJetEt);
   fChain->SetBranchAddress("L1ForJetE", L1ForJetE, &b_L1ForJetE);
   fChain->SetBranchAddress("L1ForJetEta", L1ForJetEta, &b_L1ForJetEta);
   fChain->SetBranchAddress("L1ForJetPhi", L1ForJetPhi, &b_L1ForJetPhi);
   fChain->SetBranchAddress("NL1Tau", &NL1Tau, &b_NL1Tau);
   fChain->SetBranchAddress("L1TauEt", L1TauEt, &b_L1TauEt);
   fChain->SetBranchAddress("L1TauE", L1TauE, &b_L1TauE);
   fChain->SetBranchAddress("L1TauEta", L1TauEta, &b_L1TauEta);
   fChain->SetBranchAddress("L1TauPhi", L1TauPhi, &b_L1TauPhi);
   fChain->SetBranchAddress("L1Met", &L1Met, &b_L1Met);
   fChain->SetBranchAddress("L1MetPhi", &L1MetPhi, &b_L1MetPhi);
   fChain->SetBranchAddress("L1EtTot", &L1EtTot, &b_L1EtTot);
   fChain->SetBranchAddress("L1Mht", &L1Mht, &b_L1Mht);
   fChain->SetBranchAddress("L1MhtPhi", &L1MhtPhi, &b_L1MhtPhi);
   fChain->SetBranchAddress("L1EtHad", &L1EtHad, &b_L1EtHad);
   fChain->SetBranchAddress("L1HfRing1EtSumPositiveEta", &L1HfRing1EtSumPositiveEta, &b_L1HfRing1EtSumPositiveEta);
   fChain->SetBranchAddress("L1HfRing2EtSumPositiveEta", &L1HfRing2EtSumPositiveEta, &b_L1HfRing2EtSumPositiveEta);
   fChain->SetBranchAddress("L1HfRing1EtSumNegativeEta", &L1HfRing1EtSumNegativeEta, &b_L1HfRing1EtSumNegativeEta);
   fChain->SetBranchAddress("L1HfRing2EtSumNegativeEta", &L1HfRing2EtSumNegativeEta, &b_L1HfRing2EtSumNegativeEta);
   fChain->SetBranchAddress("L1HfTowerCountPositiveEtaRing1", &L1HfTowerCountPositiveEtaRing1, &b_L1HfTowerCountPositiveEtaRing1);
   fChain->SetBranchAddress("L1HfTowerCountNegativeEtaRing1", &L1HfTowerCountNegativeEtaRing1, &b_L1HfTowerCountNegativeEtaRing1);
   fChain->SetBranchAddress("L1HfTowerCountPositiveEtaRing2", &L1HfTowerCountPositiveEtaRing2, &b_L1HfTowerCountPositiveEtaRing2);
   fChain->SetBranchAddress("L1HfTowerCountNegativeEtaRing2", &L1HfTowerCountNegativeEtaRing2, &b_L1HfTowerCountNegativeEtaRing2);
   fChain->SetBranchAddress("NMCpart", &NMCpart, &b_NMCpart);
   fChain->SetBranchAddress("MCpid", &MCpid, &b_MCpid);
   fChain->SetBranchAddress("MCstatus", &MCstatus, &b_MCstatus);
   fChain->SetBranchAddress("MCvtxX", &MCvtxX, &b_MCvtxX);
   fChain->SetBranchAddress("MCvtxY", &MCvtxY, &b_MCvtxY);
   fChain->SetBranchAddress("MCvtxZ", &MCvtxZ, &b_MCvtxZ);
   fChain->SetBranchAddress("MCpt", &MCpt, &b_MCpt);
   fChain->SetBranchAddress("MCeta", &MCeta, &b_MCeta);
   fChain->SetBranchAddress("MCphi", &MCphi, &b_MCphi);
   fChain->SetBranchAddress("MCPtHat", &MCPtHat, &b_MCPtHat);
   fChain->SetBranchAddress("MCmu3", &MCmu3, &b_MCmu3);
   fChain->SetBranchAddress("MCel3", &MCel3, &b_MCel3);
   fChain->SetBranchAddress("MCbb", &MCbb, &b_MCbb);
   fChain->SetBranchAddress("MCab", &MCab, &b_MCab);
   fChain->SetBranchAddress("MCWenu", &MCWenu, &b_MCWenu);
   fChain->SetBranchAddress("MCWmunu", &MCWmunu, &b_MCmunu);
   fChain->SetBranchAddress("MCZee", &MCZee, &b_MCZee);
   fChain->SetBranchAddress("MCZmumu", &MCZmumu, &b_MCZmumu);
   fChain->SetBranchAddress("MCptEleMax", &MCptEleMax, &b_MCptEleMax);
   fChain->SetBranchAddress("MCptMuMax", &MCptMuMax, &b_MCptMuMax);
   fChain->SetBranchAddress("NPUTrueBX0", &NPUTrueBX0, &b_NPUTrueBX0);
   fChain->SetBranchAddress("NPUgenBX0", &NPUgenBX0, &b_NPUgenBX0);
   fChain->SetBranchAddress("recoNVrtOffline0", &recoNVrtOffline0, &b_NVrtx);
   fChain->SetBranchAddress("recoVrtXOffline0", &recoVrtXOffline0, &b_recoVrtXOffline0);
   fChain->SetBranchAddress("recoVrtYOffline0", &recoVrtYOffline0, &b_recoVrtYOffline0);
   fChain->SetBranchAddress("recoVrtZOffline0", &recoVrtZOffline0, &b_recoVrtZOffline0);
   fChain->SetBranchAddress("recoVrtNtrkOffline0", &recoVrtNtrkOffline0, &b_recoVrtNtrkOffline0);
   fChain->SetBranchAddress("recoVrtChi2Offline0", &recoVrtChi2Offline0, &b_recoVrtChi2Offline0);
   fChain->SetBranchAddress("recoVrtNdofOffline0", &recoVrtNdofOffline0, &b_recoVrtNdofOffline0);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("Bx", &Bx, &b_Bx);
   fChain->SetBranchAddress("Orbit", &Orbit, &b_Orbit);
   fChain->SetBranchAddress("AvgInstDelLumi", &AvgInstDelLumi, &b_AvgInstDelLumi);
   fChain->SetBranchAddress("HLT_HIMET220_v1", &HLT_HIMET220_v1, &b_HLT_HIMET220_v1);
   fChain->SetBranchAddress("HLT_HIMET220_v1_Prescl", &HLT_HIMET220_v1_Prescl, &b_HLT_HIMET220_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMET200_v1", &HLT_HIMET200_v1, &b_HLT_HIMET200_v1);
   fChain->SetBranchAddress("HLT_HIMET200_v1_Prescl", &HLT_HIMET200_v1_Prescl, &b_HLT_HIMET200_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMET120_v1", &HLT_HIMET120_v1, &b_HLT_HIMET120_v1);
   fChain->SetBranchAddress("HLT_HIMET120_v1_Prescl", &HLT_HIMET120_v1_Prescl, &b_HLT_HIMET120_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPhysics_v1", &HLT_HIPhysics_v1, &b_HLT_HIPhysics_v1);
   fChain->SetBranchAddress("HLT_HIPhysics_v1_Prescl", &HLT_HIPhysics_v1_Prescl, &b_HLT_HIPhysics_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDTCalibration_v1", &HLT_HIDTCalibration_v1, &b_HLT_HIDTCalibration_v1);
   fChain->SetBranchAddress("HLT_HIDTCalibration_v1_Prescl", &HLT_HIDTCalibration_v1_Prescl, &b_HLT_HIDTCalibration_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIEcalCalibration_v1", &HLT_HIEcalCalibration_v1, &b_HLT_HIEcalCalibration_v1);
   fChain->SetBranchAddress("HLT_HIEcalCalibration_v1_Prescl", &HLT_HIEcalCalibration_v1_Prescl, &b_HLT_HIEcalCalibration_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIHcalCalibration_v1", &HLT_HIHcalCalibration_v1, &b_HLT_HIHcalCalibration_v1);
   fChain->SetBranchAddress("HLT_HIHcalCalibration_v1_Prescl", &HLT_HIHcalCalibration_v1_Prescl, &b_HLT_HIHcalCalibration_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBias_v1", &HLT_HIZeroBias_v1, &b_HLT_HIZeroBias_v1);
   fChain->SetBranchAddress("HLT_HIZeroBias_v1_Prescl", &HLT_HIZeroBias_v1_Prescl, &b_HLT_HIZeroBias_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBiasXOR_v1", &HLT_HIZeroBiasXOR_v1, &b_HLT_HIZeroBiasXOR_v1);
   fChain->SetBranchAddress("HLT_HIZeroBiasXOR_v1_Prescl", &HLT_HIZeroBiasXOR_v1_Prescl, &b_HLT_HIZeroBiasXOR_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v1", &HLT_HIZeroBiasPixel_SingleTrack_v1, &b_HLT_HIZeroBiasPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl", &HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl, &b_HLT_HIZeroBiasPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_v1", &HLT_HIMinBiasBSC_v1, &b_HLT_HIMinBiasBSC_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_v1_Prescl", &HLT_HIMinBiasBSC_v1_Prescl, &b_HLT_HIMinBiasBSC_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_OR_v1", &HLT_HIMinBiasBSC_OR_v1, &b_HLT_HIMinBiasBSC_OR_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_OR_v1_Prescl", &HLT_HIMinBiasBSC_OR_v1_Prescl, &b_HLT_HIMinBiasBSC_OR_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHF_v1", &HLT_HIMinBiasHF_v1, &b_HLT_HIMinBiasHF_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasHF_v1_Prescl", &HLT_HIMinBiasHF_v1_Prescl, &b_HLT_HIMinBiasHF_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHf_OR_v1", &HLT_HIMinBiasHf_OR_v1, &b_HLT_HIMinBiasHf_OR_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasHf_OR_v1_Prescl", &HLT_HIMinBiasHf_OR_v1_Prescl, &b_HLT_HIMinBiasHf_OR_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1", &HLT_HIMinBiasHfOrBSC_v1, &b_HLT_HIMinBiasHfOrBSC_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1_Prescl", &HLT_HIMinBiasHfOrBSC_v1_Prescl, &b_HLT_HIMinBiasHfOrBSC_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v1", &HLT_HIMinBiasPixel_SingleTrack_v1, &b_HLT_HIMinBiasPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v1_Prescl", &HLT_HIMinBiasPixel_SingleTrack_v1_Prescl, &b_HLT_HIMinBiasPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v1", &HLT_HIMinBiasZDC_Calo_v1, &b_HLT_HIMinBiasZDC_Calo_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v1_Prescl", &HLT_HIMinBiasZDC_Calo_v1_Prescl, &b_HLT_HIMinBiasZDC_Calo_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1", &HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1, &b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl", &HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl, &b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v1", &HLT_HIMinBiasZDCPixel_SingleTrack_v1, &b_HLT_HIMinBiasZDCPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl", &HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl, &b_HLT_HIMinBiasZDCPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1", &HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1, &b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl", &HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl, &b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIBptxXOR_v1", &HLT_HIBptxXOR_v1, &b_HLT_HIBptxXOR_v1);
   fChain->SetBranchAddress("HLT_HIBptxXOR_v1_Prescl", &HLT_HIBptxXOR_v1_Prescl, &b_HLT_HIBptxXOR_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v1", &HLT_HIL1Algo_BptxXOR_BSC_OR_v1, &b_HLT_HIL1Algo_BptxXOR_BSC_OR_v1);
   fChain->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl", &HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl, &b_HLT_HIL1Algo_BptxXOR_BSC_OR_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMuOpen_v2", &HLT_HIL1DoubleMuOpen_v2, &b_HLT_HIL1DoubleMuOpen_v2);
   fChain->SetBranchAddress("HLT_HIL1DoubleMuOpen_v2_Prescl", &HLT_HIL1DoubleMuOpen_v2_Prescl, &b_HLT_HIL1DoubleMuOpen_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v2", &HLT_HIL1DoubleMu0_HighQ_v2, &b_HLT_HIL1DoubleMu0_HighQ_v2);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v2_Prescl", &HLT_HIL1DoubleMu0_HighQ_v2_Prescl, &b_HLT_HIL1DoubleMu0_HighQ_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_v2", &HLT_HIL2Mu3_v2, &b_HLT_HIL2Mu3_v2);
   fChain->SetBranchAddress("HLT_HIL2Mu3_v2_Prescl", &HLT_HIL2Mu3_v2_Prescl, &b_HLT_HIL2Mu3_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v2", &HLT_HIL2Mu3_NHitQ_v2, &b_HLT_HIL2Mu3_NHitQ_v2);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v2_Prescl", &HLT_HIL2Mu3_NHitQ_v2_Prescl, &b_HLT_HIL2Mu3_NHitQ_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu7_v2", &HLT_HIL2Mu7_v2, &b_HLT_HIL2Mu7_v2);
   fChain->SetBranchAddress("HLT_HIL2Mu7_v2_Prescl", &HLT_HIL2Mu7_v2_Prescl, &b_HLT_HIL2Mu7_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v2", &HLT_HIL2Mu15_v2, &b_HLT_HIL2Mu15_v2);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v2_Prescl", &HLT_HIL2Mu15_v2_Prescl, &b_HLT_HIL2Mu15_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_v2", &HLT_HIL2DoubleMu0_v2, &b_HLT_HIL2DoubleMu0_v2);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_v2_Prescl", &HLT_HIL2DoubleMu0_v2_Prescl, &b_HLT_HIL2DoubleMu0_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v2", &HLT_HIL2DoubleMu0_NHitQ_v2, &b_HLT_HIL2DoubleMu0_NHitQ_v2);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v2_Prescl", &HLT_HIL2DoubleMu0_NHitQ_v2_Prescl, &b_HLT_HIL2DoubleMu0_NHitQ_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2", &HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2, &b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl", &HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl, &b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu3_v2", &HLT_HIL2DoubleMu3_v2, &b_HLT_HIL2DoubleMu3_v2);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu3_v2_Prescl", &HLT_HIL2DoubleMu3_v2_Prescl, &b_HLT_HIL2DoubleMu3_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu3_v2", &HLT_HIL3Mu3_v2, &b_HLT_HIL3Mu3_v2);
   fChain->SetBranchAddress("HLT_HIL3Mu3_v2_Prescl", &HLT_HIL3Mu3_v2_Prescl, &b_HLT_HIL3Mu3_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_v2", &HLT_HIL3DoubleMuOpen_v2, &b_HLT_HIL3DoubleMuOpen_v2);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_v2_Prescl", &HLT_HIL3DoubleMuOpen_v2_Prescl, &b_HLT_HIL3DoubleMuOpen_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v2", &HLT_HIL3DoubleMuOpen_Mgt2_v2, &b_HLT_HIL3DoubleMuOpen_Mgt2_v2);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v2", &HLT_HIL3DoubleMuOpen_Mgt2_SS_v2, &b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v2);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v2", &HLT_HIL3DoubleMuOpen_Mgt2_OS_v2, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v2);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2", &HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v2_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_v1", &HLT_HISinglePhoton15_v1, &b_HLT_HISinglePhoton15_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_v1_Prescl", &HLT_HISinglePhoton15_v1_Prescl, &b_HLT_HISinglePhoton15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_v2", &HLT_HISinglePhoton20_v2, &b_HLT_HISinglePhoton20_v2);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_v2_Prescl", &HLT_HISinglePhoton20_v2_Prescl, &b_HLT_HISinglePhoton20_v2_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_v2", &HLT_HISinglePhoton30_v2, &b_HLT_HISinglePhoton30_v2);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_v2_Prescl", &HLT_HISinglePhoton30_v2_Prescl, &b_HLT_HISinglePhoton30_v2_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_v2", &HLT_HISinglePhoton40_v2, &b_HLT_HISinglePhoton40_v2);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_v2_Prescl", &HLT_HISinglePhoton40_v2_Prescl, &b_HLT_HISinglePhoton40_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIPhoton10_Photon15_v1", &HLT_HIPhoton10_Photon15_v1, &b_HLT_HIPhoton10_Photon15_v1);
   fChain->SetBranchAddress("HLT_HIPhoton10_Photon15_v1_Prescl", &HLT_HIPhoton10_Photon15_v1_Prescl, &b_HLT_HIPhoton10_Photon15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPhoton15_Photon20_v1", &HLT_HIPhoton15_Photon20_v1, &b_HLT_HIPhoton15_Photon20_v1);
   fChain->SetBranchAddress("HLT_HIPhoton15_Photon20_v1_Prescl", &HLT_HIPhoton15_Photon20_v1_Prescl, &b_HLT_HIPhoton15_Photon20_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton10_v1", &HLT_HIDoublePhoton10_v1, &b_HLT_HIDoublePhoton10_v1);
   fChain->SetBranchAddress("HLT_HIDoublePhoton10_v1_Prescl", &HLT_HIDoublePhoton10_v1_Prescl, &b_HLT_HIDoublePhoton10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_v1", &HLT_HIDoublePhoton15_v1, &b_HLT_HIDoublePhoton15_v1);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_v1_Prescl", &HLT_HIDoublePhoton15_v1_Prescl, &b_HLT_HIDoublePhoton15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton20_v1", &HLT_HIDoublePhoton20_v1, &b_HLT_HIDoublePhoton20_v1);
   fChain->SetBranchAddress("HLT_HIDoublePhoton20_v1_Prescl", &HLT_HIDoublePhoton20_v1_Prescl, &b_HLT_HIDoublePhoton20_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIJet55_v1", &HLT_HIJet55_v1, &b_HLT_HIJet55_v1);
   fChain->SetBranchAddress("HLT_HIJet55_v1_Prescl", &HLT_HIJet55_v1_Prescl, &b_HLT_HIJet55_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIJet65_v1", &HLT_HIJet65_v1, &b_HLT_HIJet65_v1);
   fChain->SetBranchAddress("HLT_HIJet65_v1_Prescl", &HLT_HIJet65_v1_Prescl, &b_HLT_HIJet65_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIJet80_v1", &HLT_HIJet80_v1, &b_HLT_HIJet80_v1);
   fChain->SetBranchAddress("HLT_HIJet80_v_Prescl", &HLT_HIJet80_v1_Prescl, &b_HLT_HIJet80_v1_Prescl);
  fChain->SetBranchAddress("HLT_HIJet80_v7", &HLT_HIJet80_v7, &b_HLT_HIJet80_v7);
   fChain->SetBranchAddress("HLT_HIJet80_v7_Prescl", &HLT_HIJet80_v7_Prescl, &b_HLT_HIJet80_v7_Prescl);
   fChain->SetBranchAddress("HLT_HIJet95_v1", &HLT_HIJet95_v1, &b_HLT_HIJet95_v1);
   fChain->SetBranchAddress("HLT_HIJet95_v1_Prescl", &HLT_HIJet95_v1_Prescl, &b_HLT_HIJet95_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDiJet55_v1", &HLT_HIDiJet55_v1, &b_HLT_HIDiJet55_v1);
   fChain->SetBranchAddress("HLT_HIDiJet55_v1_Prescl", &HLT_HIDiJet55_v1_Prescl, &b_HLT_HIDiJet55_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIJet65_Jet55_v1", &HLT_HIJet65_Jet55_v1, &b_HLT_HIJet65_Jet55_v1);
   fChain->SetBranchAddress("HLT_HIJet65_Jet55_v1_Prescl", &HLT_HIJet65_Jet55_v1_Prescl, &b_HLT_HIJet65_Jet55_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&HLT_PAJet80_NoJetID_v1,&b_HLT_PAJet80_NoJetID_v1);
  fChain->SetBranchAddress("HLT_HIJetE30_NoBPTX_v1", &HLT_HIJetE30_NoBPTX_v1, &b_HLT_HIJetE30_NoBPTX_v1);
   fChain->SetBranchAddress("HLT_HIJetE30_NoBPTX_v1_Prescl", &HLT_HIJetE30_NoBPTX_v1_Prescl, &b_HLT_HIJetE30_NoBPTX_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v1", &HLT_HIJetE50_NoBPTX3BX_NoHalo_v1, &b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v1);
   fChain->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl", &HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl, &b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIActivityHF_Coincidence3_v1", &HLT_HIActivityHF_Coincidence3_v1, &b_HLT_HIActivityHF_Coincidence3_v1);
   fChain->SetBranchAddress("HLT_HIActivityHF_Coincidence3_v1_Prescl", &HLT_HIActivityHF_Coincidence3_v1_Prescl, &b_HLT_HIActivityHF_Coincidence3_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIActivityHF_Single3_v1", &HLT_HIActivityHF_Single3_v1, &b_HLT_HIActivityHF_Single3_v1);
   fChain->SetBranchAddress("HLT_HIActivityHF_Single3_v1_Prescl", &HLT_HIActivityHF_Single3_v1_Prescl, &b_HLT_HIActivityHF_Single3_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIClusterVertexCompatibility_v1", &HLT_HIClusterVertexCompatibility_v1, &b_HLT_HIClusterVertexCompatibility_v1);
   fChain->SetBranchAddress("HLT_HIClusterVertexCompatibility_v1_Prescl", &HLT_HIClusterVertexCompatibility_v1_Prescl, &b_HLT_HIClusterVertexCompatibility_v1_Prescl);
   fChain->SetBranchAddress("HLT_HICentralityVeto_v1", &HLT_HICentralityVeto_v1, &b_HLT_HICentralityVeto_v1);
   fChain->SetBranchAddress("HLT_HICentralityVeto_v1_Prescl", &HLT_HICentralityVeto_v1_Prescl, &b_HLT_HICentralityVeto_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Central_v1", &HLT_HIFullTrack12_L1Central_v1, &b_HLT_HIFullTrack12_L1Central_v1);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Central_v1_Prescl", &HLT_HIFullTrack12_L1Central_v1_Prescl, &b_HLT_HIFullTrack12_L1Central_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v1", &HLT_HIFullTrack12_L1Peripheral_v1, &b_HLT_HIFullTrack12_L1Peripheral_v1);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v1_Prescl", &HLT_HIFullTrack12_L1Peripheral_v1_Prescl, &b_HLT_HIFullTrack12_L1Peripheral_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Central_v1", &HLT_HIFullTrack14_L1Central_v1, &b_HLT_HIFullTrack14_L1Central_v1);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Central_v1_Prescl", &HLT_HIFullTrack14_L1Central_v1_Prescl, &b_HLT_HIFullTrack14_L1Central_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v1", &HLT_HIFullTrack14_L1Peripheral_v1, &b_HLT_HIFullTrack14_L1Peripheral_v1);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v1_Prescl", &HLT_HIFullTrack14_L1Peripheral_v1_Prescl, &b_HLT_HIFullTrack14_L1Peripheral_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Central_v1", &HLT_HIFullTrack20_L1Central_v1, &b_HLT_HIFullTrack20_L1Central_v1);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Central_v1_Prescl", &HLT_HIFullTrack20_L1Central_v1_Prescl, &b_HLT_HIFullTrack20_L1Central_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v1", &HLT_HIFullTrack20_L1Peripheral_v1, &b_HLT_HIFullTrack20_L1Peripheral_v1);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v1_Prescl", &HLT_HIFullTrack20_L1Peripheral_v1_Prescl, &b_HLT_HIFullTrack20_L1Peripheral_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Central_v1", &HLT_HIFullTrack25_L1Central_v1, &b_HLT_HIFullTrack25_L1Central_v1);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Central_v1_Prescl", &HLT_HIFullTrack25_L1Central_v1_Prescl, &b_HLT_HIFullTrack25_L1Central_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v1", &HLT_HIFullTrack25_L1Peripheral_v1, &b_HLT_HIFullTrack25_L1Peripheral_v1);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v1_Prescl", &HLT_HIFullTrack25_L1Peripheral_v1_Prescl, &b_HLT_HIFullTrack25_L1Peripheral_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIRandom_v1", &HLT_HIRandom_v1, &b_HLT_HIRandom_v1);
   fChain->SetBranchAddress("HLT_HIRandom_v1_Prescl", &HLT_HIRandom_v1_Prescl, &b_HLT_HIRandom_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUCC010_v2", &HLT_HIUCC010_v2, &b_HLT_HIUCC010_v2);
   fChain->SetBranchAddress("HLT_HIUCC010_v2_Prescl", &HLT_HIUCC010_v2_Prescl, &b_HLT_HIUCC010_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIUCC015_v2", &HLT_HIUCC015_v2, &b_HLT_HIUCC015_v2);
   fChain->SetBranchAddress("HLT_HIUCC015_v2_Prescl", &HLT_HIUCC015_v2_Prescl, &b_HLT_HIUCC015_v2_Prescl);
   fChain->SetBranchAddress("HLT_HICentral10_v2", &HLT_HICentral10_v2, &b_HLT_HICentral10_v2);
   fChain->SetBranchAddress("HLT_HICentral10_v2_Prescl", &HLT_HICentral10_v2_Prescl, &b_HLT_HICentral10_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v1", &HLT_HIUPCNeuMuPixel_SingleTrack_v1, &b_HLT_HIUPCNeuMuPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl", &HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCNeuMuPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v1", &HLT_HIUPCNeuEG2Pixel_SingleTrack_v1, &b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl", &HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v1", &HLT_HIUPCNeuEG5Pixel_SingleTrack_v1, &b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl", &HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1", &HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1, &b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl", &HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1", &HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1, &b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl", &HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1", &HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1, &b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl", &HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl, &b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("HLTriggerFinalPath_Prescl", &HLTriggerFinalPath_Prescl, &b_HLTriggerFinalPath_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v1", &HLT_HIL1DoubleMu0_HighQ_v1, &b_HLT_HIL1DoubleMu0_HighQ_v1);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v1_Prescl", &HLT_HIL1DoubleMu0_HighQ_v1_Prescl, &b_HLT_HIL1DoubleMu0_HighQ_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMuOpen_v1", &HLT_HIL1DoubleMuOpen_v1, &b_HLT_HIL1DoubleMuOpen_v1);
   fChain->SetBranchAddress("HLT_HIL1DoubleMuOpen_v1_Prescl", &HLT_HIL1DoubleMuOpen_v1_Prescl, &b_HLT_HIL1DoubleMuOpen_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1", &HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1, &b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl", &HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl, &b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v1", &HLT_HIL2DoubleMu0_NHitQ_v1, &b_HLT_HIL2DoubleMu0_NHitQ_v1);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v1_Prescl", &HLT_HIL2DoubleMu0_NHitQ_v1_Prescl, &b_HLT_HIL2DoubleMu0_NHitQ_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_v1", &HLT_HIL2DoubleMu0_v1, &b_HLT_HIL2DoubleMu0_v1);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_v1_Prescl", &HLT_HIL2DoubleMu0_v1_Prescl, &b_HLT_HIL2DoubleMu0_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu3_v1", &HLT_HIL2DoubleMu3_v1, &b_HLT_HIL2DoubleMu3_v1);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu3_v1_Prescl", &HLT_HIL2DoubleMu3_v1_Prescl, &b_HLT_HIL2DoubleMu3_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v1", &HLT_HIL2Mu15_v1, &b_HLT_HIL2Mu15_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v1_Prescl", &HLT_HIL2Mu15_v1_Prescl, &b_HLT_HIL2Mu15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v1", &HLT_HIL2Mu3_NHitQ_v1, &b_HLT_HIL2Mu3_NHitQ_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v1_Prescl", &HLT_HIL2Mu3_NHitQ_v1_Prescl, &b_HLT_HIL2Mu3_NHitQ_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_v1", &HLT_HIL2Mu3_v1, &b_HLT_HIL2Mu3_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3_v1_Prescl", &HLT_HIL2Mu3_v1_Prescl, &b_HLT_HIL2Mu3_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu7_v1", &HLT_HIL2Mu7_v1, &b_HLT_HIL2Mu7_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu7_v1_Prescl", &HLT_HIL2Mu7_v1_Prescl, &b_HLT_HIL2Mu7_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1", &HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v1", &HLT_HIL3DoubleMuOpen_Mgt2_OS_v1, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v1);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v1", &HLT_HIL3DoubleMuOpen_Mgt2_SS_v1, &b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v1);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v1", &HLT_HIL3DoubleMuOpen_Mgt2_v1, &b_HLT_HIL3DoubleMuOpen_Mgt2_v1);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_v1", &HLT_HIL3DoubleMuOpen_v1, &b_HLT_HIL3DoubleMuOpen_v1);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_v1_Prescl", &HLT_HIL3DoubleMuOpen_v1_Prescl, &b_HLT_HIL3DoubleMuOpen_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu3_v1", &HLT_HIL3Mu3_v1, &b_HLT_HIL3Mu3_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu3_v1_Prescl", &HLT_HIL3Mu3_v1_Prescl, &b_HLT_HIL3Mu3_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_v1", &HLT_HISinglePhoton20_v1, &b_HLT_HISinglePhoton20_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_v1_Prescl", &HLT_HISinglePhoton20_v1_Prescl, &b_HLT_HISinglePhoton20_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_v1", &HLT_HISinglePhoton30_v1, &b_HLT_HISinglePhoton30_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_v1_Prescl", &HLT_HISinglePhoton30_v1_Prescl, &b_HLT_HISinglePhoton30_v1_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_v1", &HLT_HISinglePhoton40_v1, &b_HLT_HISinglePhoton40_v1);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_v1_Prescl", &HLT_HISinglePhoton40_v1_Prescl, &b_HLT_HISinglePhoton40_v1_Prescl);
   fChain->SetBranchAddress("HLT_HICentral10_v1", &HLT_HICentral10_v1, &b_HLT_HICentral10_v1);
   fChain->SetBranchAddress("HLT_HICentral10_v1_Prescl", &HLT_HICentral10_v1_Prescl, &b_HLT_HICentral10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUCC010_v1", &HLT_HIUCC010_v1, &b_HLT_HIUCC010_v1);
   fChain->SetBranchAddress("HLT_HIUCC010_v1_Prescl", &HLT_HIUCC010_v1_Prescl, &b_HLT_HIUCC010_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIUCC015_v1", &HLT_HIUCC015_v1, &b_HLT_HIUCC015_v1);
   fChain->SetBranchAddress("HLT_HIUCC015_v1_Prescl", &HLT_HIUCC015_v1_Prescl, &b_HLT_HIUCC015_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIDiJet55_v2", &HLT_HIDiJet55_v2, &b_HLT_HIDiJet55_v2);
   fChain->SetBranchAddress("HLT_HIDiJet55_v2_Prescl", &HLT_HIDiJet55_v2_Prescl, &b_HLT_HIDiJet55_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton10_v2", &HLT_HIDoublePhoton10_v2, &b_HLT_HIDoublePhoton10_v2);
   fChain->SetBranchAddress("HLT_HIDoublePhoton10_v2_Prescl", &HLT_HIDoublePhoton10_v2_Prescl, &b_HLT_HIDoublePhoton10_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_v2", &HLT_HIDoublePhoton15_v2, &b_HLT_HIDoublePhoton15_v2);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_v2_Prescl", &HLT_HIDoublePhoton15_v2_Prescl, &b_HLT_HIDoublePhoton15_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton20_v2", &HLT_HIDoublePhoton20_v2, &b_HLT_HIDoublePhoton20_v2);
   fChain->SetBranchAddress("HLT_HIDoublePhoton20_v2_Prescl", &HLT_HIDoublePhoton20_v2_Prescl, &b_HLT_HIDoublePhoton20_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Central_v2", &HLT_HIFullTrack12_L1Central_v2, &b_HLT_HIFullTrack12_L1Central_v2);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Central_v2_Prescl", &HLT_HIFullTrack12_L1Central_v2_Prescl, &b_HLT_HIFullTrack12_L1Central_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v2", &HLT_HIFullTrack12_L1Peripheral_v2, &b_HLT_HIFullTrack12_L1Peripheral_v2);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v2_Prescl", &HLT_HIFullTrack12_L1Peripheral_v2_Prescl, &b_HLT_HIFullTrack12_L1Peripheral_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Central_v2", &HLT_HIFullTrack14_L1Central_v2, &b_HLT_HIFullTrack14_L1Central_v2);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Central_v2_Prescl", &HLT_HIFullTrack14_L1Central_v2_Prescl, &b_HLT_HIFullTrack14_L1Central_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v2", &HLT_HIFullTrack14_L1Peripheral_v2, &b_HLT_HIFullTrack14_L1Peripheral_v2);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v2_Prescl", &HLT_HIFullTrack14_L1Peripheral_v2_Prescl, &b_HLT_HIFullTrack14_L1Peripheral_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Central_v2", &HLT_HIFullTrack20_L1Central_v2, &b_HLT_HIFullTrack20_L1Central_v2);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Central_v2_Prescl", &HLT_HIFullTrack20_L1Central_v2_Prescl, &b_HLT_HIFullTrack20_L1Central_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v2", &HLT_HIFullTrack20_L1Peripheral_v2, &b_HLT_HIFullTrack20_L1Peripheral_v2);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v2_Prescl", &HLT_HIFullTrack20_L1Peripheral_v2_Prescl, &b_HLT_HIFullTrack20_L1Peripheral_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Central_v2", &HLT_HIFullTrack25_L1Central_v2, &b_HLT_HIFullTrack25_L1Central_v2);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Central_v2_Prescl", &HLT_HIFullTrack25_L1Central_v2_Prescl, &b_HLT_HIFullTrack25_L1Central_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v2", &HLT_HIFullTrack25_L1Peripheral_v2, &b_HLT_HIFullTrack25_L1Peripheral_v2);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v2_Prescl", &HLT_HIFullTrack25_L1Peripheral_v2_Prescl, &b_HLT_HIFullTrack25_L1Peripheral_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIJet55_v2", &HLT_HIJet55_v2, &b_HLT_HIJet55_v2);
   fChain->SetBranchAddress("HLT_HIJet55_v2_Prescl", &HLT_HIJet55_v2_Prescl, &b_HLT_HIJet55_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIJet65_Jet55_v2", &HLT_HIJet65_Jet55_v2, &b_HLT_HIJet65_Jet55_v2);
   fChain->SetBranchAddress("HLT_HIJet65_Jet55_v2_Prescl", &HLT_HIJet65_Jet55_v2_Prescl, &b_HLT_HIJet65_Jet55_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIJet65_v2", &HLT_HIJet65_v2, &b_HLT_HIJet65_v2);
   fChain->SetBranchAddress("HLT_HIJet65_v2_Prescl", &HLT_HIJet65_v2_Prescl, &b_HLT_HIJet65_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIJet80_v2", &HLT_HIJet80_v2, &b_HLT_HIJet80_v2);
   fChain->SetBranchAddress("HLT_HIJet80_v2_Prescl", &HLT_HIJet80_v2_Prescl, &b_HLT_HIJet80_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIJet95_v2", &HLT_HIJet95_v2, &b_HLT_HIJet95_v2);
   fChain->SetBranchAddress("HLT_HIJet95_v2_Prescl", &HLT_HIJet95_v2_Prescl, &b_HLT_HIJet95_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIJetE30_NoBPTX_v2", &HLT_HIJetE30_NoBPTX_v2, &b_HLT_HIJetE30_NoBPTX_v2);
   fChain->SetBranchAddress("HLT_HIJetE30_NoBPTX_v2_Prescl", &HLT_HIJetE30_NoBPTX_v2_Prescl, &b_HLT_HIJetE30_NoBPTX_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v2", &HLT_HIJetE50_NoBPTX3BX_NoHalo_v2, &b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v2);
   fChain->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl", &HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl, &b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMET120_v2", &HLT_HIMET120_v2, &b_HLT_HIMET120_v2);
   fChain->SetBranchAddress("HLT_HIMET120_v2_Prescl", &HLT_HIMET120_v2_Prescl, &b_HLT_HIMET120_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMET200_v2", &HLT_HIMET200_v2, &b_HLT_HIMET200_v2);
   fChain->SetBranchAddress("HLT_HIMET200_v2_Prescl", &HLT_HIMET200_v2_Prescl, &b_HLT_HIMET200_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMET220_v2", &HLT_HIMET220_v2, &b_HLT_HIMET220_v2);
   fChain->SetBranchAddress("HLT_HIMET220_v2_Prescl", &HLT_HIMET220_v2_Prescl, &b_HLT_HIMET220_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIPhoton10_Photon15_v2", &HLT_HIPhoton10_Photon15_v2, &b_HLT_HIPhoton10_Photon15_v2);
   fChain->SetBranchAddress("HLT_HIPhoton10_Photon15_v2_Prescl", &HLT_HIPhoton10_Photon15_v2_Prescl, &b_HLT_HIPhoton10_Photon15_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIPhoton15_Photon20_v2", &HLT_HIPhoton15_Photon20_v2, &b_HLT_HIPhoton15_Photon20_v2);
   fChain->SetBranchAddress("HLT_HIPhoton15_Photon20_v2_Prescl", &HLT_HIPhoton15_Photon20_v2_Prescl, &b_HLT_HIPhoton15_Photon20_v2_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_v2", &HLT_HISinglePhoton15_v2, &b_HLT_HISinglePhoton15_v2);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_v2_Prescl", &HLT_HISinglePhoton15_v2_Prescl, &b_HLT_HISinglePhoton15_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIBptxXOR_v2", &HLT_HIBptxXOR_v2, &b_HLT_HIBptxXOR_v2);
   fChain->SetBranchAddress("HLT_HIBptxXOR_v2_Prescl", &HLT_HIBptxXOR_v2_Prescl, &b_HLT_HIBptxXOR_v2_Prescl);
   fChain->SetBranchAddress("HLT_HICentralityVeto_v2", &HLT_HICentralityVeto_v2, &b_HLT_HICentralityVeto_v2);
   fChain->SetBranchAddress("HLT_HICentralityVeto_v2_Prescl", &HLT_HICentralityVeto_v2_Prescl, &b_HLT_HICentralityVeto_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v2", &HLT_HIL1Algo_BptxXOR_BSC_OR_v2, &b_HLT_HIL1Algo_BptxXOR_BSC_OR_v2);
   fChain->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl", &HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl, &b_HLT_HIL1Algo_BptxXOR_BSC_OR_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_OR_v2", &HLT_HIMinBiasBSC_OR_v2, &b_HLT_HIMinBiasBSC_OR_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_OR_v2_Prescl", &HLT_HIMinBiasBSC_OR_v2_Prescl, &b_HLT_HIMinBiasBSC_OR_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_v2", &HLT_HIMinBiasBSC_v2, &b_HLT_HIMinBiasBSC_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_v2_Prescl", &HLT_HIMinBiasBSC_v2_Prescl, &b_HLT_HIMinBiasBSC_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHF_v2", &HLT_HIMinBiasHF_v2, &b_HLT_HIMinBiasHF_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasHF_v2_Prescl", &HLT_HIMinBiasHF_v2_Prescl, &b_HLT_HIMinBiasHF_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v2", &HLT_HIMinBiasHfOrBSC_v2, &b_HLT_HIMinBiasHfOrBSC_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v2_Prescl", &HLT_HIMinBiasHfOrBSC_v2_Prescl, &b_HLT_HIMinBiasHfOrBSC_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHf_OR_v2", &HLT_HIMinBiasHf_OR_v2, &b_HLT_HIMinBiasHf_OR_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasHf_OR_v2_Prescl", &HLT_HIMinBiasHf_OR_v2_Prescl, &b_HLT_HIMinBiasHf_OR_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v2", &HLT_HIMinBiasPixel_SingleTrack_v2, &b_HLT_HIMinBiasPixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v2_Prescl", &HLT_HIMinBiasPixel_SingleTrack_v2_Prescl, &b_HLT_HIMinBiasPixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v2", &HLT_HIMinBiasZDCPixel_SingleTrack_v2, &b_HLT_HIMinBiasZDCPixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl", &HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl, &b_HLT_HIMinBiasZDCPixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2", &HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2, &b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl", &HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl, &b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v2", &HLT_HIMinBiasZDC_Calo_v2, &b_HLT_HIMinBiasZDC_Calo_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v2_Prescl", &HLT_HIMinBiasZDC_Calo_v2_Prescl, &b_HLT_HIMinBiasZDC_Calo_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2", &HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2, &b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl", &HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl, &b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIPhysics_v2", &HLT_HIPhysics_v2, &b_HLT_HIPhysics_v2);
   fChain->SetBranchAddress("HLT_HIPhysics_v2_Prescl", &HLT_HIPhysics_v2_Prescl, &b_HLT_HIPhysics_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIRandom_v2", &HLT_HIRandom_v2, &b_HLT_HIRandom_v2);
   fChain->SetBranchAddress("HLT_HIRandom_v2_Prescl", &HLT_HIRandom_v2_Prescl, &b_HLT_HIRandom_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v2", &HLT_HIUPCNeuEG2Pixel_SingleTrack_v2, &b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl", &HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl, &b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v2", &HLT_HIUPCNeuEG5Pixel_SingleTrack_v2, &b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl", &HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl, &b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2", &HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2, &b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl", &HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl, &b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2", &HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2, &b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl", &HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl, &b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2", &HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2, &b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl", &HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl, &b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v2", &HLT_HIUPCNeuMuPixel_SingleTrack_v2, &b_HLT_HIUPCNeuMuPixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl", &HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl, &b_HLT_HIUPCNeuMuPixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v2", &HLT_HIZeroBiasPixel_SingleTrack_v2, &b_HLT_HIZeroBiasPixel_SingleTrack_v2);
   fChain->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl", &HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl, &b_HLT_HIZeroBiasPixel_SingleTrack_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBiasXOR_v2", &HLT_HIZeroBiasXOR_v2, &b_HLT_HIZeroBiasXOR_v2);
   fChain->SetBranchAddress("HLT_HIZeroBiasXOR_v2_Prescl", &HLT_HIZeroBiasXOR_v2_Prescl, &b_HLT_HIZeroBiasXOR_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBias_v2", &HLT_HIZeroBias_v2, &b_HLT_HIZeroBias_v2);
   fChain->SetBranchAddress("HLT_HIZeroBias_v2_Prescl", &HLT_HIZeroBias_v2_Prescl, &b_HLT_HIZeroBias_v2_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v3", &HLT_HIL1DoubleMu0_HighQ_v3, &b_HLT_HIL1DoubleMu0_HighQ_v3);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v3_Prescl", &HLT_HIL1DoubleMu0_HighQ_v3_Prescl, &b_HLT_HIL1DoubleMu0_HighQ_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMuOpen_v3", &HLT_HIL1DoubleMuOpen_v3, &b_HLT_HIL1DoubleMuOpen_v3);
   fChain->SetBranchAddress("HLT_HIL1DoubleMuOpen_v3_Prescl", &HLT_HIL1DoubleMuOpen_v3_Prescl, &b_HLT_HIL1DoubleMuOpen_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3", &HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3, &b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl", &HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl, &b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v3", &HLT_HIL2DoubleMu0_NHitQ_v3, &b_HLT_HIL2DoubleMu0_NHitQ_v3);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v3_Prescl", &HLT_HIL2DoubleMu0_NHitQ_v3_Prescl, &b_HLT_HIL2DoubleMu0_NHitQ_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_v3", &HLT_HIL2DoubleMu0_v3, &b_HLT_HIL2DoubleMu0_v3);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_v3_Prescl", &HLT_HIL2DoubleMu0_v3_Prescl, &b_HLT_HIL2DoubleMu0_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu3_v3", &HLT_HIL2DoubleMu3_v3, &b_HLT_HIL2DoubleMu3_v3);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu3_v3_Prescl", &HLT_HIL2DoubleMu3_v3_Prescl, &b_HLT_HIL2DoubleMu3_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v3", &HLT_HIL2Mu15_v3, &b_HLT_HIL2Mu15_v3);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v3_Prescl", &HLT_HIL2Mu15_v3_Prescl, &b_HLT_HIL2Mu15_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v3", &HLT_HIL2Mu3_NHitQ_v3, &b_HLT_HIL2Mu3_NHitQ_v3);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v3_Prescl", &HLT_HIL2Mu3_NHitQ_v3_Prescl, &b_HLT_HIL2Mu3_NHitQ_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_v3", &HLT_HIL2Mu3_v3, &b_HLT_HIL2Mu3_v3);
   fChain->SetBranchAddress("HLT_HIL2Mu3_v3_Prescl", &HLT_HIL2Mu3_v3_Prescl, &b_HLT_HIL2Mu3_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu7_v3", &HLT_HIL2Mu7_v3, &b_HLT_HIL2Mu7_v3);
   fChain->SetBranchAddress("HLT_HIL2Mu7_v3_Prescl", &HLT_HIL2Mu7_v3_Prescl, &b_HLT_HIL2Mu7_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3", &HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v3", &HLT_HIL3DoubleMuOpen_Mgt2_OS_v3, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v3);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v3", &HLT_HIL3DoubleMuOpen_Mgt2_SS_v3, &b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v3);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v3", &HLT_HIL3DoubleMuOpen_Mgt2_v3, &b_HLT_HIL3DoubleMuOpen_Mgt2_v3);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_v3", &HLT_HIL3DoubleMuOpen_v3, &b_HLT_HIL3DoubleMuOpen_v3);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_v3_Prescl", &HLT_HIL3DoubleMuOpen_v3_Prescl, &b_HLT_HIL3DoubleMuOpen_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu3_v3", &HLT_HIL3Mu3_v3, &b_HLT_HIL3Mu3_v3);
   fChain->SetBranchAddress("HLT_HIL3Mu3_v3_Prescl", &HLT_HIL3Mu3_v3_Prescl, &b_HLT_HIL3Mu3_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIDiJet55_v3", &HLT_HIDiJet55_v3, &b_HLT_HIDiJet55_v3);
   fChain->SetBranchAddress("HLT_HIDiJet55_v3_Prescl", &HLT_HIDiJet55_v3_Prescl, &b_HLT_HIDiJet55_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton10_v3", &HLT_HIDoublePhoton10_v3, &b_HLT_HIDoublePhoton10_v3);
   fChain->SetBranchAddress("HLT_HIDoublePhoton10_v3_Prescl", &HLT_HIDoublePhoton10_v3_Prescl, &b_HLT_HIDoublePhoton10_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_v3", &HLT_HIDoublePhoton15_v3, &b_HLT_HIDoublePhoton15_v3);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_v3_Prescl", &HLT_HIDoublePhoton15_v3_Prescl, &b_HLT_HIDoublePhoton15_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton20_v3", &HLT_HIDoublePhoton20_v3, &b_HLT_HIDoublePhoton20_v3);
   fChain->SetBranchAddress("HLT_HIDoublePhoton20_v3_Prescl", &HLT_HIDoublePhoton20_v3_Prescl, &b_HLT_HIDoublePhoton20_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Central_v3", &HLT_HIFullTrack12_L1Central_v3, &b_HLT_HIFullTrack12_L1Central_v3);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Central_v3_Prescl", &HLT_HIFullTrack12_L1Central_v3_Prescl, &b_HLT_HIFullTrack12_L1Central_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v3", &HLT_HIFullTrack12_L1Peripheral_v3, &b_HLT_HIFullTrack12_L1Peripheral_v3);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v3_Prescl", &HLT_HIFullTrack12_L1Peripheral_v3_Prescl, &b_HLT_HIFullTrack12_L1Peripheral_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Central_v3", &HLT_HIFullTrack14_L1Central_v3, &b_HLT_HIFullTrack14_L1Central_v3);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Central_v3_Prescl", &HLT_HIFullTrack14_L1Central_v3_Prescl, &b_HLT_HIFullTrack14_L1Central_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v3", &HLT_HIFullTrack14_L1Peripheral_v3, &b_HLT_HIFullTrack14_L1Peripheral_v3);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v3_Prescl", &HLT_HIFullTrack14_L1Peripheral_v3_Prescl, &b_HLT_HIFullTrack14_L1Peripheral_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Central_v3", &HLT_HIFullTrack20_L1Central_v3, &b_HLT_HIFullTrack20_L1Central_v3);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Central_v3_Prescl", &HLT_HIFullTrack20_L1Central_v3_Prescl, &b_HLT_HIFullTrack20_L1Central_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v3", &HLT_HIFullTrack20_L1Peripheral_v3, &b_HLT_HIFullTrack20_L1Peripheral_v3);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v3_Prescl", &HLT_HIFullTrack20_L1Peripheral_v3_Prescl, &b_HLT_HIFullTrack20_L1Peripheral_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Central_v3", &HLT_HIFullTrack25_L1Central_v3, &b_HLT_HIFullTrack25_L1Central_v3);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Central_v3_Prescl", &HLT_HIFullTrack25_L1Central_v3_Prescl, &b_HLT_HIFullTrack25_L1Central_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v3", &HLT_HIFullTrack25_L1Peripheral_v3, &b_HLT_HIFullTrack25_L1Peripheral_v3);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v3_Prescl", &HLT_HIFullTrack25_L1Peripheral_v3_Prescl, &b_HLT_HIFullTrack25_L1Peripheral_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIJet55_v3", &HLT_HIJet55_v3, &b_HLT_HIJet55_v3);
   fChain->SetBranchAddress("HLT_HIJet55_v3_Prescl", &HLT_HIJet55_v3_Prescl, &b_HLT_HIJet55_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIJet65_Jet55_v3", &HLT_HIJet65_Jet55_v3, &b_HLT_HIJet65_Jet55_v3);
   fChain->SetBranchAddress("HLT_HIJet65_Jet55_v3_Prescl", &HLT_HIJet65_Jet55_v3_Prescl, &b_HLT_HIJet65_Jet55_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIJet65_v3", &HLT_HIJet65_v3, &b_HLT_HIJet65_v3);
   fChain->SetBranchAddress("HLT_HIJet65_v3_Prescl", &HLT_HIJet65_v3_Prescl, &b_HLT_HIJet65_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIJet80_v3", &HLT_HIJet80_v3, &b_HLT_HIJet80_v3);
   fChain->SetBranchAddress("HLT_HIJet80_v3_Prescl", &HLT_HIJet80_v3_Prescl, &b_HLT_HIJet80_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIJet95_v3", &HLT_HIJet95_v3, &b_HLT_HIJet95_v3);
   fChain->SetBranchAddress("HLT_HIJet95_v3_Prescl", &HLT_HIJet95_v3_Prescl, &b_HLT_HIJet95_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIJetE30_NoBPTX_v3", &HLT_HIJetE30_NoBPTX_v3, &b_HLT_HIJetE30_NoBPTX_v3);
   fChain->SetBranchAddress("HLT_HIJetE30_NoBPTX_v3_Prescl", &HLT_HIJetE30_NoBPTX_v3_Prescl, &b_HLT_HIJetE30_NoBPTX_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v3", &HLT_HIJetE50_NoBPTX3BX_NoHalo_v3, &b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v3);
   fChain->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl", &HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl, &b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMET120_v3", &HLT_HIMET120_v3, &b_HLT_HIMET120_v3);
   fChain->SetBranchAddress("HLT_HIMET120_v3_Prescl", &HLT_HIMET120_v3_Prescl, &b_HLT_HIMET120_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMET200_v3", &HLT_HIMET200_v3, &b_HLT_HIMET200_v3);
   fChain->SetBranchAddress("HLT_HIMET200_v3_Prescl", &HLT_HIMET200_v3_Prescl, &b_HLT_HIMET200_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMET220_v3", &HLT_HIMET220_v3, &b_HLT_HIMET220_v3);
   fChain->SetBranchAddress("HLT_HIMET220_v3_Prescl", &HLT_HIMET220_v3_Prescl, &b_HLT_HIMET220_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIPhoton10_Photon15_v3", &HLT_HIPhoton10_Photon15_v3, &b_HLT_HIPhoton10_Photon15_v3);
   fChain->SetBranchAddress("HLT_HIPhoton10_Photon15_v3_Prescl", &HLT_HIPhoton10_Photon15_v3_Prescl, &b_HLT_HIPhoton10_Photon15_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIPhoton15_Photon20_v3", &HLT_HIPhoton15_Photon20_v3, &b_HLT_HIPhoton15_Photon20_v3);
   fChain->SetBranchAddress("HLT_HIPhoton15_Photon20_v3_Prescl", &HLT_HIPhoton15_Photon20_v3_Prescl, &b_HLT_HIPhoton15_Photon20_v3_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_v3", &HLT_HISinglePhoton15_v3, &b_HLT_HISinglePhoton15_v3);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_v3_Prescl", &HLT_HISinglePhoton15_v3_Prescl, &b_HLT_HISinglePhoton15_v3_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_v3", &HLT_HISinglePhoton20_v3, &b_HLT_HISinglePhoton20_v3);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_v3_Prescl", &HLT_HISinglePhoton20_v3_Prescl, &b_HLT_HISinglePhoton20_v3_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_v3", &HLT_HISinglePhoton30_v3, &b_HLT_HISinglePhoton30_v3);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_v3_Prescl", &HLT_HISinglePhoton30_v3_Prescl, &b_HLT_HISinglePhoton30_v3_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_v3", &HLT_HISinglePhoton40_v3, &b_HLT_HISinglePhoton40_v3);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_v3_Prescl", &HLT_HISinglePhoton40_v3_Prescl, &b_HLT_HISinglePhoton40_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIBptxXOR_v3", &HLT_HIBptxXOR_v3, &b_HLT_HIBptxXOR_v3);
   fChain->SetBranchAddress("HLT_HIBptxXOR_v3_Prescl", &HLT_HIBptxXOR_v3_Prescl, &b_HLT_HIBptxXOR_v3_Prescl);
   fChain->SetBranchAddress("HLT_HICentral10_v3", &HLT_HICentral10_v3, &b_HLT_HICentral10_v3);
   fChain->SetBranchAddress("HLT_HICentral10_v3_Prescl", &HLT_HICentral10_v3_Prescl, &b_HLT_HICentral10_v3_Prescl);
   fChain->SetBranchAddress("HLT_HICentralityVeto_v3", &HLT_HICentralityVeto_v3, &b_HLT_HICentralityVeto_v3);
   fChain->SetBranchAddress("HLT_HICentralityVeto_v3_Prescl", &HLT_HICentralityVeto_v3_Prescl, &b_HLT_HICentralityVeto_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v3", &HLT_HIL1Algo_BptxXOR_BSC_OR_v3, &b_HLT_HIL1Algo_BptxXOR_BSC_OR_v3);
   fChain->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl", &HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl, &b_HLT_HIL1Algo_BptxXOR_BSC_OR_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_OR_v3", &HLT_HIMinBiasBSC_OR_v3, &b_HLT_HIMinBiasBSC_OR_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_OR_v3_Prescl", &HLT_HIMinBiasBSC_OR_v3_Prescl, &b_HLT_HIMinBiasBSC_OR_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_v3", &HLT_HIMinBiasBSC_v3, &b_HLT_HIMinBiasBSC_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_v3_Prescl", &HLT_HIMinBiasBSC_v3_Prescl, &b_HLT_HIMinBiasBSC_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHF_v3", &HLT_HIMinBiasHF_v3, &b_HLT_HIMinBiasHF_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasHF_v3_Prescl", &HLT_HIMinBiasHF_v3_Prescl, &b_HLT_HIMinBiasHF_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v3", &HLT_HIMinBiasHfOrBSC_v3, &b_HLT_HIMinBiasHfOrBSC_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v3_Prescl", &HLT_HIMinBiasHfOrBSC_v3_Prescl, &b_HLT_HIMinBiasHfOrBSC_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHf_OR_v3", &HLT_HIMinBiasHf_OR_v3, &b_HLT_HIMinBiasHf_OR_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasHf_OR_v3_Prescl", &HLT_HIMinBiasHf_OR_v3_Prescl, &b_HLT_HIMinBiasHf_OR_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v3", &HLT_HIMinBiasPixel_SingleTrack_v3, &b_HLT_HIMinBiasPixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v3_Prescl", &HLT_HIMinBiasPixel_SingleTrack_v3_Prescl, &b_HLT_HIMinBiasPixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v3", &HLT_HIMinBiasZDCPixel_SingleTrack_v3, &b_HLT_HIMinBiasZDCPixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl", &HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl, &b_HLT_HIMinBiasZDCPixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3", &HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3, &b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl", &HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl, &b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v3", &HLT_HIMinBiasZDC_Calo_v3, &b_HLT_HIMinBiasZDC_Calo_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v3_Prescl", &HLT_HIMinBiasZDC_Calo_v3_Prescl, &b_HLT_HIMinBiasZDC_Calo_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3", &HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3, &b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl", &HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl, &b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIPhysics_v3", &HLT_HIPhysics_v3, &b_HLT_HIPhysics_v3);
   fChain->SetBranchAddress("HLT_HIPhysics_v3_Prescl", &HLT_HIPhysics_v3_Prescl, &b_HLT_HIPhysics_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIRandom_v3", &HLT_HIRandom_v3, &b_HLT_HIRandom_v3);
   fChain->SetBranchAddress("HLT_HIRandom_v3_Prescl", &HLT_HIRandom_v3_Prescl, &b_HLT_HIRandom_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIUCC010_v3", &HLT_HIUCC010_v3, &b_HLT_HIUCC010_v3);
   fChain->SetBranchAddress("HLT_HIUCC010_v3_Prescl", &HLT_HIUCC010_v3_Prescl, &b_HLT_HIUCC010_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIUCC015_v3", &HLT_HIUCC015_v3, &b_HLT_HIUCC015_v3);
   fChain->SetBranchAddress("HLT_HIUCC015_v3_Prescl", &HLT_HIUCC015_v3_Prescl, &b_HLT_HIUCC015_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v3", &HLT_HIUPCNeuEG2Pixel_SingleTrack_v3, &b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl", &HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl, &b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v3", &HLT_HIUPCNeuEG5Pixel_SingleTrack_v3, &b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl", &HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl, &b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3", &HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3, &b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl", &HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl, &b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3", &HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3, &b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl", &HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl, &b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3", &HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3, &b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl", &HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl, &b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v3", &HLT_HIUPCNeuMuPixel_SingleTrack_v3, &b_HLT_HIUPCNeuMuPixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl", &HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl, &b_HLT_HIUPCNeuMuPixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v3", &HLT_HIZeroBiasPixel_SingleTrack_v3, &b_HLT_HIZeroBiasPixel_SingleTrack_v3);
   fChain->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl", &HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl, &b_HLT_HIZeroBiasPixel_SingleTrack_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBiasXOR_v3", &HLT_HIZeroBiasXOR_v3, &b_HLT_HIZeroBiasXOR_v3);
   fChain->SetBranchAddress("HLT_HIZeroBiasXOR_v3_Prescl", &HLT_HIZeroBiasXOR_v3_Prescl, &b_HLT_HIZeroBiasXOR_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBias_v3", &HLT_HIZeroBias_v3, &b_HLT_HIZeroBias_v3);
   fChain->SetBranchAddress("HLT_HIZeroBias_v3_Prescl", &HLT_HIZeroBias_v3_Prescl, &b_HLT_HIZeroBias_v3_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v4", &HLT_HIL1DoubleMu0_HighQ_v4, &b_HLT_HIL1DoubleMu0_HighQ_v4);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_HighQ_v4_Prescl", &HLT_HIL1DoubleMu0_HighQ_v4_Prescl, &b_HLT_HIL1DoubleMu0_HighQ_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMuOpen_v4", &HLT_HIL1DoubleMuOpen_v4, &b_HLT_HIL1DoubleMuOpen_v4);
   fChain->SetBranchAddress("HLT_HIL1DoubleMuOpen_v4_Prescl", &HLT_HIL1DoubleMuOpen_v4_Prescl, &b_HLT_HIL1DoubleMuOpen_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4", &HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4, &b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl", &HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl, &b_HLT_HIL2DoubleMu0_L1HighQL2NHitQ_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v4", &HLT_HIL2DoubleMu0_NHitQ_v4, &b_HLT_HIL2DoubleMu0_NHitQ_v4);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v4_Prescl", &HLT_HIL2DoubleMu0_NHitQ_v4_Prescl, &b_HLT_HIL2DoubleMu0_NHitQ_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_v4", &HLT_HIL2DoubleMu0_v4, &b_HLT_HIL2DoubleMu0_v4);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_v4_Prescl", &HLT_HIL2DoubleMu0_v4_Prescl, &b_HLT_HIL2DoubleMu0_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu3_v4", &HLT_HIL2DoubleMu3_v4, &b_HLT_HIL2DoubleMu3_v4);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu3_v4_Prescl", &HLT_HIL2DoubleMu3_v4_Prescl, &b_HLT_HIL2DoubleMu3_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v4", &HLT_HIL2Mu15_v4, &b_HLT_HIL2Mu15_v4);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v4_Prescl", &HLT_HIL2Mu15_v4_Prescl, &b_HLT_HIL2Mu15_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v4", &HLT_HIL2Mu3_NHitQ_v4, &b_HLT_HIL2Mu3_NHitQ_v4);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ_v4_Prescl", &HLT_HIL2Mu3_NHitQ_v4_Prescl, &b_HLT_HIL2Mu3_NHitQ_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_v4", &HLT_HIL2Mu3_v4, &b_HLT_HIL2Mu3_v4);
   fChain->SetBranchAddress("HLT_HIL2Mu3_v4_Prescl", &HLT_HIL2Mu3_v4_Prescl, &b_HLT_HIL2Mu3_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu7_v4", &HLT_HIL2Mu7_v4, &b_HLT_HIL2Mu7_v4);
   fChain->SetBranchAddress("HLT_HIL2Mu7_v4_Prescl", &HLT_HIL2Mu7_v4_Prescl, &b_HLT_HIL2Mu7_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4", &HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v4", &HLT_HIL3DoubleMuOpen_Mgt2_OS_v4, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v4);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_OS_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v4", &HLT_HIL3DoubleMuOpen_Mgt2_SS_v4, &b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v4);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_SS_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v4", &HLT_HIL3DoubleMuOpen_Mgt2_v4, &b_HLT_HIL3DoubleMuOpen_Mgt2_v4);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl", &HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl, &b_HLT_HIL3DoubleMuOpen_Mgt2_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_v4", &HLT_HIL3DoubleMuOpen_v4, &b_HLT_HIL3DoubleMuOpen_v4);
   fChain->SetBranchAddress("HLT_HIL3DoubleMuOpen_v4_Prescl", &HLT_HIL3DoubleMuOpen_v4_Prescl, &b_HLT_HIL3DoubleMuOpen_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu3_v4", &HLT_HIL3Mu3_v4, &b_HLT_HIL3Mu3_v4);
   fChain->SetBranchAddress("HLT_HIL3Mu3_v4_Prescl", &HLT_HIL3Mu3_v4_Prescl, &b_HLT_HIL3Mu3_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIDiJet55_v4", &HLT_HIDiJet55_v4, &b_HLT_HIDiJet55_v4);
   fChain->SetBranchAddress("HLT_HIDiJet55_v4_Prescl", &HLT_HIDiJet55_v4_Prescl, &b_HLT_HIDiJet55_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton10_v4", &HLT_HIDoublePhoton10_v4, &b_HLT_HIDoublePhoton10_v4);
   fChain->SetBranchAddress("HLT_HIDoublePhoton10_v4_Prescl", &HLT_HIDoublePhoton10_v4_Prescl, &b_HLT_HIDoublePhoton10_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_v4", &HLT_HIDoublePhoton15_v4, &b_HLT_HIDoublePhoton15_v4);
   fChain->SetBranchAddress("HLT_HIDoublePhoton15_v4_Prescl", &HLT_HIDoublePhoton15_v4_Prescl, &b_HLT_HIDoublePhoton15_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIDoublePhoton20_v4", &HLT_HIDoublePhoton20_v4, &b_HLT_HIDoublePhoton20_v4);
   fChain->SetBranchAddress("HLT_HIDoublePhoton20_v4_Prescl", &HLT_HIDoublePhoton20_v4_Prescl, &b_HLT_HIDoublePhoton20_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Central_v4", &HLT_HIFullTrack12_L1Central_v4, &b_HLT_HIFullTrack12_L1Central_v4);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Central_v4_Prescl", &HLT_HIFullTrack12_L1Central_v4_Prescl, &b_HLT_HIFullTrack12_L1Central_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v4", &HLT_HIFullTrack12_L1Peripheral_v4, &b_HLT_HIFullTrack12_L1Peripheral_v4);
   fChain->SetBranchAddress("HLT_HIFullTrack12_L1Peripheral_v4_Prescl", &HLT_HIFullTrack12_L1Peripheral_v4_Prescl, &b_HLT_HIFullTrack12_L1Peripheral_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Central_v4", &HLT_HIFullTrack14_L1Central_v4, &b_HLT_HIFullTrack14_L1Central_v4);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Central_v4_Prescl", &HLT_HIFullTrack14_L1Central_v4_Prescl, &b_HLT_HIFullTrack14_L1Central_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v4", &HLT_HIFullTrack14_L1Peripheral_v4, &b_HLT_HIFullTrack14_L1Peripheral_v4);
   fChain->SetBranchAddress("HLT_HIFullTrack14_L1Peripheral_v4_Prescl", &HLT_HIFullTrack14_L1Peripheral_v4_Prescl, &b_HLT_HIFullTrack14_L1Peripheral_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Central_v4", &HLT_HIFullTrack20_L1Central_v4, &b_HLT_HIFullTrack20_L1Central_v4);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Central_v4_Prescl", &HLT_HIFullTrack20_L1Central_v4_Prescl, &b_HLT_HIFullTrack20_L1Central_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v4", &HLT_HIFullTrack20_L1Peripheral_v4, &b_HLT_HIFullTrack20_L1Peripheral_v4);
   fChain->SetBranchAddress("HLT_HIFullTrack20_L1Peripheral_v4_Prescl", &HLT_HIFullTrack20_L1Peripheral_v4_Prescl, &b_HLT_HIFullTrack20_L1Peripheral_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Central_v4", &HLT_HIFullTrack25_L1Central_v4, &b_HLT_HIFullTrack25_L1Central_v4);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Central_v4_Prescl", &HLT_HIFullTrack25_L1Central_v4_Prescl, &b_HLT_HIFullTrack25_L1Central_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v4", &HLT_HIFullTrack25_L1Peripheral_v4, &b_HLT_HIFullTrack25_L1Peripheral_v4);
   fChain->SetBranchAddress("HLT_HIFullTrack25_L1Peripheral_v4_Prescl", &HLT_HIFullTrack25_L1Peripheral_v4_Prescl, &b_HLT_HIFullTrack25_L1Peripheral_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIJet55_v4", &HLT_HIJet55_v4, &b_HLT_HIJet55_v4);
   fChain->SetBranchAddress("HLT_HIJet55_v4_Prescl", &HLT_HIJet55_v4_Prescl, &b_HLT_HIJet55_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIJet65_Jet55_v4", &HLT_HIJet65_Jet55_v4, &b_HLT_HIJet65_Jet55_v4);
   fChain->SetBranchAddress("HLT_HIJet65_Jet55_v4_Prescl", &HLT_HIJet65_Jet55_v4_Prescl, &b_HLT_HIJet65_Jet55_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIJet65_v4", &HLT_HIJet65_v4, &b_HLT_HIJet65_v4);
   fChain->SetBranchAddress("HLT_HIJet65_v4_Prescl", &HLT_HIJet65_v4_Prescl, &b_HLT_HIJet65_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIJet80_v4", &HLT_HIJet80_v4, &b_HLT_HIJet80_v4);
   fChain->SetBranchAddress("HLT_HIJet80_v4_Prescl", &HLT_HIJet80_v4_Prescl, &b_HLT_HIJet80_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIJet95_v4", &HLT_HIJet95_v4, &b_HLT_HIJet95_v4);
   fChain->SetBranchAddress("HLT_HIJet95_v4_Prescl", &HLT_HIJet95_v4_Prescl, &b_HLT_HIJet95_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIJetE30_NoBPTX_v4", &HLT_HIJetE30_NoBPTX_v4, &b_HLT_HIJetE30_NoBPTX_v4);
   fChain->SetBranchAddress("HLT_HIJetE30_NoBPTX_v4_Prescl", &HLT_HIJetE30_NoBPTX_v4_Prescl, &b_HLT_HIJetE30_NoBPTX_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v4", &HLT_HIJetE50_NoBPTX3BX_NoHalo_v4, &b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v4);
   fChain->SetBranchAddress("HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl", &HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl, &b_HLT_HIJetE50_NoBPTX3BX_NoHalo_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMET120_v4", &HLT_HIMET120_v4, &b_HLT_HIMET120_v4);
   fChain->SetBranchAddress("HLT_HIMET120_v4_Prescl", &HLT_HIMET120_v4_Prescl, &b_HLT_HIMET120_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMET200_v4", &HLT_HIMET200_v4, &b_HLT_HIMET200_v4);
   fChain->SetBranchAddress("HLT_HIMET200_v4_Prescl", &HLT_HIMET200_v4_Prescl, &b_HLT_HIMET200_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMET220_v4", &HLT_HIMET220_v4, &b_HLT_HIMET220_v4);
   fChain->SetBranchAddress("HLT_HIMET220_v4_Prescl", &HLT_HIMET220_v4_Prescl, &b_HLT_HIMET220_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIPhoton10_Photon15_v4", &HLT_HIPhoton10_Photon15_v4, &b_HLT_HIPhoton10_Photon15_v4);
   fChain->SetBranchAddress("HLT_HIPhoton10_Photon15_v4_Prescl", &HLT_HIPhoton10_Photon15_v4_Prescl, &b_HLT_HIPhoton10_Photon15_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIPhoton15_Photon20_v4", &HLT_HIPhoton15_Photon20_v4, &b_HLT_HIPhoton15_Photon20_v4);
   fChain->SetBranchAddress("HLT_HIPhoton15_Photon20_v4_Prescl", &HLT_HIPhoton15_Photon20_v4_Prescl, &b_HLT_HIPhoton15_Photon20_v4_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_v4", &HLT_HISinglePhoton15_v4, &b_HLT_HISinglePhoton15_v4);
   fChain->SetBranchAddress("HLT_HISinglePhoton15_v4_Prescl", &HLT_HISinglePhoton15_v4_Prescl, &b_HLT_HISinglePhoton15_v4_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_v4", &HLT_HISinglePhoton20_v4, &b_HLT_HISinglePhoton20_v4);
   fChain->SetBranchAddress("HLT_HISinglePhoton20_v4_Prescl", &HLT_HISinglePhoton20_v4_Prescl, &b_HLT_HISinglePhoton20_v4_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_v4", &HLT_HISinglePhoton30_v4, &b_HLT_HISinglePhoton30_v4);
   fChain->SetBranchAddress("HLT_HISinglePhoton30_v4_Prescl", &HLT_HISinglePhoton30_v4_Prescl, &b_HLT_HISinglePhoton30_v4_Prescl);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_v4", &HLT_HISinglePhoton40_v4, &b_HLT_HISinglePhoton40_v4);
   fChain->SetBranchAddress("HLT_HISinglePhoton40_v4_Prescl", &HLT_HISinglePhoton40_v4_Prescl, &b_HLT_HISinglePhoton40_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIBptxXOR_v4", &HLT_HIBptxXOR_v4, &b_HLT_HIBptxXOR_v4);
   fChain->SetBranchAddress("HLT_HIBptxXOR_v4_Prescl", &HLT_HIBptxXOR_v4_Prescl, &b_HLT_HIBptxXOR_v4_Prescl);
   fChain->SetBranchAddress("HLT_HICentral10_v4", &HLT_HICentral10_v4, &b_HLT_HICentral10_v4);
   fChain->SetBranchAddress("HLT_HICentral10_v4_Prescl", &HLT_HICentral10_v4_Prescl, &b_HLT_HICentral10_v4_Prescl);
   fChain->SetBranchAddress("HLT_HICentralityVeto_v4", &HLT_HICentralityVeto_v4, &b_HLT_HICentralityVeto_v4);
   fChain->SetBranchAddress("HLT_HICentralityVeto_v4_Prescl", &HLT_HICentralityVeto_v4_Prescl, &b_HLT_HICentralityVeto_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v4", &HLT_HIL1Algo_BptxXOR_BSC_OR_v4, &b_HLT_HIL1Algo_BptxXOR_BSC_OR_v4);
   fChain->SetBranchAddress("HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl", &HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl, &b_HLT_HIL1Algo_BptxXOR_BSC_OR_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_OR_v4", &HLT_HIMinBiasBSC_OR_v4, &b_HLT_HIMinBiasBSC_OR_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_OR_v4_Prescl", &HLT_HIMinBiasBSC_OR_v4_Prescl, &b_HLT_HIMinBiasBSC_OR_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_v4", &HLT_HIMinBiasBSC_v4, &b_HLT_HIMinBiasBSC_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasBSC_v4_Prescl", &HLT_HIMinBiasBSC_v4_Prescl, &b_HLT_HIMinBiasBSC_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHF_v4", &HLT_HIMinBiasHF_v4, &b_HLT_HIMinBiasHF_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasHF_v4_Prescl", &HLT_HIMinBiasHF_v4_Prescl, &b_HLT_HIMinBiasHF_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v4", &HLT_HIMinBiasHfOrBSC_v4, &b_HLT_HIMinBiasHfOrBSC_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v4_Prescl", &HLT_HIMinBiasHfOrBSC_v4_Prescl, &b_HLT_HIMinBiasHfOrBSC_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasHf_OR_v4", &HLT_HIMinBiasHf_OR_v4, &b_HLT_HIMinBiasHf_OR_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasHf_OR_v4_Prescl", &HLT_HIMinBiasHf_OR_v4_Prescl, &b_HLT_HIMinBiasHf_OR_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v4", &HLT_HIMinBiasPixel_SingleTrack_v4, &b_HLT_HIMinBiasPixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasPixel_SingleTrack_v4_Prescl", &HLT_HIMinBiasPixel_SingleTrack_v4_Prescl, &b_HLT_HIMinBiasPixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v4", &HLT_HIMinBiasZDCPixel_SingleTrack_v4, &b_HLT_HIMinBiasZDCPixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl", &HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl, &b_HLT_HIMinBiasZDCPixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4", &HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4, &b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl", &HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl, &b_HLT_HIMinBiasZDC_Calo_PlusOrMinus_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v4", &HLT_HIMinBiasZDC_Calo_v4, &b_HLT_HIMinBiasZDC_Calo_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_Calo_v4_Prescl", &HLT_HIMinBiasZDC_Calo_v4_Prescl, &b_HLT_HIMinBiasZDC_Calo_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4", &HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4, &b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl", &HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl, &b_HLT_HIMinBiasZDC_PlusOrMinusPixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIPhysics_v4", &HLT_HIPhysics_v4, &b_HLT_HIPhysics_v4);
   fChain->SetBranchAddress("HLT_HIPhysics_v4_Prescl", &HLT_HIPhysics_v4_Prescl, &b_HLT_HIPhysics_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIRandom_v4", &HLT_HIRandom_v4, &b_HLT_HIRandom_v4);
   fChain->SetBranchAddress("HLT_HIRandom_v4_Prescl", &HLT_HIRandom_v4_Prescl, &b_HLT_HIRandom_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIUCC010_v4", &HLT_HIUCC010_v4, &b_HLT_HIUCC010_v4);
   fChain->SetBranchAddress("HLT_HIUCC010_v4_Prescl", &HLT_HIUCC010_v4_Prescl, &b_HLT_HIUCC010_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIUCC015_v4", &HLT_HIUCC015_v4, &b_HLT_HIUCC015_v4);
   fChain->SetBranchAddress("HLT_HIUCC015_v4_Prescl", &HLT_HIUCC015_v4_Prescl, &b_HLT_HIUCC015_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v4", &HLT_HIUPCNeuEG2Pixel_SingleTrack_v4, &b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl", &HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl, &b_HLT_HIUPCNeuEG2Pixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v4", &HLT_HIUPCNeuEG5Pixel_SingleTrack_v4, &b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl", &HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl, &b_HLT_HIUPCNeuEG5Pixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4", &HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4, &b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl", &HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl, &b_HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4", &HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4, &b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl", &HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl, &b_HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4", &HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4, &b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl", &HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl, &b_HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v4", &HLT_HIUPCNeuMuPixel_SingleTrack_v4, &b_HLT_HIUPCNeuMuPixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl", &HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl, &b_HLT_HIUPCNeuMuPixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v4", &HLT_HIZeroBiasPixel_SingleTrack_v4, &b_HLT_HIZeroBiasPixel_SingleTrack_v4);
   fChain->SetBranchAddress("HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl", &HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl, &b_HLT_HIZeroBiasPixel_SingleTrack_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBiasXOR_v4", &HLT_HIZeroBiasXOR_v4, &b_HLT_HIZeroBiasXOR_v4);
   fChain->SetBranchAddress("HLT_HIZeroBiasXOR_v4_Prescl", &HLT_HIZeroBiasXOR_v4_Prescl, &b_HLT_HIZeroBiasXOR_v4_Prescl);
   fChain->SetBranchAddress("HLT_HIZeroBias_v4", &HLT_HIZeroBias_v4, &b_HLT_HIZeroBias_v4);
   fChain->SetBranchAddress("HLT_HIZeroBias_v4_Prescl", &HLT_HIZeroBias_v4_Prescl, &b_HLT_HIZeroBias_v4_Prescl);
   fChain->SetBranchAddress("L1_BptxMinus", &L1_BptxMinus, &b_L1_BptxMinus);
   fChain->SetBranchAddress("L1_BptxMinus_Prescl", &L1_BptxMinus_Prescl, &b_L1_BptxMinus_Prescl);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus", &L1_BptxMinus_NotBptxPlus, &b_L1_BptxMinus_NotBptxPlus);
   fChain->SetBranchAddress("L1_BptxMinus_NotBptxPlus_Prescl", &L1_BptxMinus_NotBptxPlus_Prescl, &b_L1_BptxMinus_NotBptxPlus_Prescl);
   fChain->SetBranchAddress("L1_BptxPlus", &L1_BptxPlus, &b_L1_BptxPlus);
   fChain->SetBranchAddress("L1_BptxPlus_Prescl", &L1_BptxPlus_Prescl, &b_L1_BptxPlus_Prescl);
   fChain->SetBranchAddress("L1_BptxPlusANDMinus", &L1_BptxPlusANDMinus, &b_L1_BptxPlusANDMinus);
   fChain->SetBranchAddress("L1_BptxPlusANDMinus_Prescl", &L1_BptxPlusANDMinus_Prescl, &b_L1_BptxPlusANDMinus_Prescl);
   fChain->SetBranchAddress("L1_BptxPlusORMinus", &L1_BptxPlusORMinus, &b_L1_BptxPlusORMinus);
   fChain->SetBranchAddress("L1_BptxPlusORMinus_Prescl", &L1_BptxPlusORMinus_Prescl, &b_L1_BptxPlusORMinus_Prescl);
   fChain->SetBranchAddress("L1_BptxPlusORMinus_instance1", &L1_BptxPlusORMinus_instance1, &b_L1_BptxPlusORMinus_instance1);
   fChain->SetBranchAddress("L1_BptxPlusORMinus_instance1_Prescl", &L1_BptxPlusORMinus_instance1_Prescl, &b_L1_BptxPlusORMinus_instance1_Prescl);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus", &L1_BptxPlus_NotBptxMinus, &b_L1_BptxPlus_NotBptxMinus);
   fChain->SetBranchAddress("L1_BptxPlus_NotBptxMinus_Prescl", &L1_BptxPlus_NotBptxMinus_Prescl, &b_L1_BptxPlus_NotBptxMinus_Prescl);
   fChain->SetBranchAddress("L1_BptxXOR", &L1_BptxXOR, &b_L1_BptxXOR);
   fChain->SetBranchAddress("L1_BptxXOR_Prescl", &L1_BptxXOR_Prescl, &b_L1_BptxXOR_Prescl);
   fChain->SetBranchAddress("L1_BptxXOR_BscMinBiasOR", &L1_BptxXOR_BscMinBiasOR, &b_L1_BptxXOR_BscMinBiasOR);
   fChain->SetBranchAddress("L1_BptxXOR_BscMinBiasOR_Prescl", &L1_BptxXOR_BscMinBiasOR_Prescl, &b_L1_BptxXOR_BscMinBiasOR_Prescl);
   fChain->SetBranchAddress("L1_Bsc2Minus_BptxMinus", &L1_Bsc2Minus_BptxMinus, &b_L1_Bsc2Minus_BptxMinus);
   fChain->SetBranchAddress("L1_Bsc2Minus_BptxMinus_Prescl", &L1_Bsc2Minus_BptxMinus_Prescl, &b_L1_Bsc2Minus_BptxMinus_Prescl);
   fChain->SetBranchAddress("L1_Bsc2Plus_BptxPlus", &L1_Bsc2Plus_BptxPlus, &b_L1_Bsc2Plus_BptxPlus);
   fChain->SetBranchAddress("L1_Bsc2Plus_BptxPlus_Prescl", &L1_Bsc2Plus_BptxPlus_Prescl, &b_L1_Bsc2Plus_BptxPlus_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold1", &L1_BscMinBiasInnerThreshold1, &b_L1_BscMinBiasInnerThreshold1);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold1_Prescl", &L1_BscMinBiasInnerThreshold1_Prescl, &b_L1_BscMinBiasInnerThreshold1_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold1_BptxAND", &L1_BscMinBiasInnerThreshold1_BptxAND, &b_L1_BscMinBiasInnerThreshold1_BptxAND);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold1_BptxAND_Prescl", &L1_BscMinBiasInnerThreshold1_BptxAND_Prescl, &b_L1_BscMinBiasInnerThreshold1_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold2_BptxAND", &L1_BscMinBiasInnerThreshold2_BptxAND, &b_L1_BscMinBiasInnerThreshold2_BptxAND);
   fChain->SetBranchAddress("L1_BscMinBiasInnerThreshold2_BptxAND_Prescl", &L1_BscMinBiasInnerThreshold2_BptxAND_Prescl, &b_L1_BscMinBiasInnerThreshold2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxAND", &L1_BscMinBiasOR_BptxAND, &b_L1_BscMinBiasOR_BptxAND);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxAND_Prescl", &L1_BscMinBiasOR_BptxAND_Prescl, &b_L1_BscMinBiasOR_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusORMinus", &L1_BscMinBiasOR_BptxPlusORMinus, &b_L1_BscMinBiasOR_BptxPlusORMinus);
   fChain->SetBranchAddress("L1_BscMinBiasOR_BptxPlusORMinus_Prescl", &L1_BscMinBiasOR_BptxPlusORMinus_Prescl, &b_L1_BscMinBiasOR_BptxPlusORMinus_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold1", &L1_BscMinBiasThreshold1, &b_L1_BscMinBiasThreshold1);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold1_Prescl", &L1_BscMinBiasThreshold1_Prescl, &b_L1_BscMinBiasThreshold1_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold1_BptxAND", &L1_BscMinBiasThreshold1_BptxAND, &b_L1_BscMinBiasThreshold1_BptxAND);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold1_BptxAND_Prescl", &L1_BscMinBiasThreshold1_BptxAND_Prescl, &b_L1_BscMinBiasThreshold1_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold2_BptxAND", &L1_BscMinBiasThreshold2_BptxAND, &b_L1_BscMinBiasThreshold2_BptxAND);
   fChain->SetBranchAddress("L1_BscMinBiasThreshold2_BptxAND_Prescl", &L1_BscMinBiasThreshold2_BptxAND_Prescl, &b_L1_BscMinBiasThreshold2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG10", &L1_DoubleEG10, &b_L1_DoubleEG10);
   fChain->SetBranchAddress("L1_DoubleEG10_Prescl", &L1_DoubleEG10_Prescl, &b_L1_DoubleEG10_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG3", &L1_DoubleEG3, &b_L1_DoubleEG3);
   fChain->SetBranchAddress("L1_DoubleEG3_Prescl", &L1_DoubleEG3_Prescl, &b_L1_DoubleEG3_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG5_BptxAND", &L1_DoubleEG5_BptxAND, &b_L1_DoubleEG5_BptxAND);
   fChain->SetBranchAddress("L1_DoubleEG5_BptxAND_Prescl", &L1_DoubleEG5_BptxAND_Prescl, &b_L1_DoubleEG5_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleEG8_BptxAND", &L1_DoubleEG8_BptxAND, &b_L1_DoubleEG8_BptxAND);
   fChain->SetBranchAddress("L1_DoubleEG8_BptxAND_Prescl", &L1_DoubleEG8_BptxAND_Prescl, &b_L1_DoubleEG8_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleForJet32_EtaOpp", &L1_DoubleForJet32_EtaOpp, &b_L1_DoubleForJet32_EtaOpp);
   fChain->SetBranchAddress("L1_DoubleForJet32_EtaOpp_Prescl", &L1_DoubleForJet32_EtaOpp_Prescl, &b_L1_DoubleForJet32_EtaOpp_Prescl);
   fChain->SetBranchAddress("L1_DoubleForJet44_EtaOpp", &L1_DoubleForJet44_EtaOpp, &b_L1_DoubleForJet44_EtaOpp);
   fChain->SetBranchAddress("L1_DoubleForJet44_EtaOpp_Prescl", &L1_DoubleForJet44_EtaOpp_Prescl, &b_L1_DoubleForJet44_EtaOpp_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing1_P1N1", &L1_DoubleHfBitCountsRing1_P1N1, &b_L1_DoubleHfBitCountsRing1_P1N1);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing1_P1N1_Prescl", &L1_DoubleHfBitCountsRing1_P1N1_Prescl, &b_L1_DoubleHfBitCountsRing1_P1N1_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing2_P1N1", &L1_DoubleHfBitCountsRing2_P1N1, &b_L1_DoubleHfBitCountsRing2_P1N1);
   fChain->SetBranchAddress("L1_DoubleHfBitCountsRing2_P1N1_Prescl", &L1_DoubleHfBitCountsRing2_P1N1_Prescl, &b_L1_DoubleHfBitCountsRing2_P1N1_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P200N200", &L1_DoubleHfRingEtSumsRing1_P200N200, &b_L1_DoubleHfRingEtSumsRing1_P200N200);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P200N200_Prescl", &L1_DoubleHfRingEtSumsRing1_P200N200_Prescl, &b_L1_DoubleHfRingEtSumsRing1_P200N200_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P4N4", &L1_DoubleHfRingEtSumsRing1_P4N4, &b_L1_DoubleHfRingEtSumsRing1_P4N4);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing1_P4N4_Prescl", &L1_DoubleHfRingEtSumsRing1_P4N4_Prescl, &b_L1_DoubleHfRingEtSumsRing1_P4N4_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P200N200", &L1_DoubleHfRingEtSumsRing2_P200N200, &b_L1_DoubleHfRingEtSumsRing2_P200N200);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P200N200_Prescl", &L1_DoubleHfRingEtSumsRing2_P200N200_Prescl, &b_L1_DoubleHfRingEtSumsRing2_P200N200_Prescl);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P4N4", &L1_DoubleHfRingEtSumsRing2_P4N4, &b_L1_DoubleHfRingEtSumsRing2_P4N4);
   fChain->SetBranchAddress("L1_DoubleHfRingEtSumsRing2_P4N4_Prescl", &L1_DoubleHfRingEtSumsRing2_P4N4_Prescl, &b_L1_DoubleHfRingEtSumsRing2_P4N4_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu3_BptxAND", &L1_DoubleMu3_BptxAND, &b_L1_DoubleMu3_BptxAND);
   fChain->SetBranchAddress("L1_DoubleMu3_BptxAND_Prescl", &L1_DoubleMu3_BptxAND_Prescl, &b_L1_DoubleMu3_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMuOpen_BptxAND", &L1_DoubleMuOpen_BptxAND, &b_L1_DoubleMuOpen_BptxAND);
   fChain->SetBranchAddress("L1_DoubleMuOpen_BptxAND_Prescl", &L1_DoubleMuOpen_BptxAND_Prescl, &b_L1_DoubleMuOpen_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND", &L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND, &b_L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND);
   fChain->SetBranchAddress("L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl", &L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl, &b_L1_EG2_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND", &L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND, &b_L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND);
   fChain->SetBranchAddress("L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl", &L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl, &b_L1_EG2_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND", &L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND, &b_L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND);
   fChain->SetBranchAddress("L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl", &L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl, &b_L1_EG5_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND", &L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND, &b_L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND);
   fChain->SetBranchAddress("L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl", &L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl, &b_L1_EG5_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETM30_BptxAND", &L1_ETM30_BptxAND, &b_L1_ETM30_BptxAND);
   fChain->SetBranchAddress("L1_ETM30_BptxAND_Prescl", &L1_ETM30_BptxAND_Prescl, &b_L1_ETM30_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETM50_BptxAND", &L1_ETM50_BptxAND, &b_L1_ETM50_BptxAND);
   fChain->SetBranchAddress("L1_ETM50_BptxAND_Prescl", &L1_ETM50_BptxAND_Prescl, &b_L1_ETM50_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT100_BptxAND", &L1_ETT100_BptxAND, &b_L1_ETT100_BptxAND);
   fChain->SetBranchAddress("L1_ETT100_BptxAND_Prescl", &L1_ETT100_BptxAND_Prescl, &b_L1_ETT100_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT140_BptxAND", &L1_ETT140_BptxAND, &b_L1_ETT140_BptxAND);
   fChain->SetBranchAddress("L1_ETT140_BptxAND_Prescl", &L1_ETT140_BptxAND_Prescl, &b_L1_ETT140_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT2000", &L1_ETT2000, &b_L1_ETT2000);
   fChain->SetBranchAddress("L1_ETT2000_Prescl", &L1_ETT2000_Prescl, &b_L1_ETT2000_Prescl);
   fChain->SetBranchAddress("L1_ETT220_BptxAND", &L1_ETT220_BptxAND, &b_L1_ETT220_BptxAND);
   fChain->SetBranchAddress("L1_ETT220_BptxAND_Prescl", &L1_ETT220_BptxAND_Prescl, &b_L1_ETT220_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT70_BptxAND", &L1_ETT70_BptxAND, &b_L1_ETT70_BptxAND);
   fChain->SetBranchAddress("L1_ETT70_BptxAND_Prescl", &L1_ETT70_BptxAND_Prescl, &b_L1_ETT70_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT800", &L1_ETT800, &b_L1_ETT800);
   fChain->SetBranchAddress("L1_ETT800_Prescl", &L1_ETT800_Prescl, &b_L1_ETT800_Prescl);
   fChain->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND", &L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND, &b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND);
   fChain->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl", &L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl, &b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1", &L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1, &b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1);
   fChain->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl", &L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl, &b_L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND_instance1_Prescl);
   fChain->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND", &L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND, &b_L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND);
   fChain->SetBranchAddress("L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl", &L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl, &b_L1_HcalHfCoincPmORBscMinBiasThresh2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_HcalHfCoincidencePm", &L1_HcalHfCoincidencePm, &b_L1_HcalHfCoincidencePm);
   fChain->SetBranchAddress("L1_HcalHfCoincidencePm_Prescl", &L1_HcalHfCoincidencePm_Prescl, &b_L1_HcalHfCoincidencePm_Prescl);
   fChain->SetBranchAddress("L1_HcalHfCoincidencePm_BptxAND", &L1_HcalHfCoincidencePm_BptxAND, &b_L1_HcalHfCoincidencePm_BptxAND);
   fChain->SetBranchAddress("L1_HcalHfCoincidencePm_BptxAND_Prescl", &L1_HcalHfCoincidencePm_BptxAND_Prescl, &b_L1_HcalHfCoincidencePm_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_HcalHfMmOrPpOrPm", &L1_HcalHfMmOrPpOrPm, &b_L1_HcalHfMmOrPpOrPm);
   fChain->SetBranchAddress("L1_HcalHfMmOrPpOrPm_Prescl", &L1_HcalHfMmOrPpOrPm_Prescl, &b_L1_HcalHfMmOrPpOrPm_Prescl);
   fChain->SetBranchAddress("L1_HcalHfMmOrPpOrPm_BptxAND", &L1_HcalHfMmOrPpOrPm_BptxAND, &b_L1_HcalHfMmOrPpOrPm_BptxAND);
   fChain->SetBranchAddress("L1_HcalHfMmOrPpOrPm_BptxAND_Prescl", &L1_HcalHfMmOrPpOrPm_BptxAND_Prescl, &b_L1_HcalHfMmOrPpOrPm_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_HcalHfMmpOrMpp", &L1_HcalHfMmpOrMpp, &b_L1_HcalHfMmpOrMpp);
   fChain->SetBranchAddress("L1_HcalHfMmpOrMpp_Prescl", &L1_HcalHfMmpOrMpp_Prescl, &b_L1_HcalHfMmpOrMpp_Prescl);
   fChain->SetBranchAddress("L1_HcalHfMmpOrMpp_BptxAND", &L1_HcalHfMmpOrMpp_BptxAND, &b_L1_HcalHfMmpOrMpp_BptxAND);
   fChain->SetBranchAddress("L1_HcalHfMmpOrMpp_BptxAND_Prescl", &L1_HcalHfMmpOrMpp_BptxAND_Prescl, &b_L1_HcalHfMmpOrMpp_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_HcalHoTotalOR", &L1_HcalHoTotalOR, &b_L1_HcalHoTotalOR);
   fChain->SetBranchAddress("L1_HcalHoTotalOR_Prescl", &L1_HcalHoTotalOR_Prescl, &b_L1_HcalHoTotalOR_Prescl);
   fChain->SetBranchAddress("L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND", &L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND, &b_L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND);
   fChain->SetBranchAddress("L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl", &L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl, &b_L1_MuOpen_ZdcCalo_NotBscMinBiasThresh2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND", &L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND, &b_L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND);
   fChain->SetBranchAddress("L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl", &L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl, &b_L1_MuOpen_ZdcCalo_NotHcalHfCoincidencePm_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1", &L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1, &b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl", &L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl, &b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold1_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2", &L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2, &b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl", &L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl, &b_L1_NotBsc2_BptxAND_BscMinBiasInnerThreshold2_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasOR", &L1_NotBsc2_BptxAND_BscMinBiasOR, &b_L1_NotBsc2_BptxAND_BscMinBiasOR);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl", &L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl, &b_L1_NotBsc2_BptxAND_BscMinBiasOR_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasThreshold1", &L1_NotBsc2_BptxAND_BscMinBiasThreshold1, &b_L1_NotBsc2_BptxAND_BscMinBiasThreshold1);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl", &L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl, &b_L1_NotBsc2_BptxAND_BscMinBiasThreshold1_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasThreshold2", &L1_NotBsc2_BptxAND_BscMinBiasThreshold2, &b_L1_NotBsc2_BptxAND_BscMinBiasThreshold2);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl", &L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl, &b_L1_NotBsc2_BptxAND_BscMinBiasThreshold2_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1", &L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1, &b_L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl", &L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl, &b_L1_NotBsc2_BptxAND_HcalHfCoincPmORBscMinBiasThresh1_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfCoincidencePm", &L1_NotBsc2_BptxAND_HcalHfCoincidencePm, &b_L1_NotBsc2_BptxAND_HcalHfCoincidencePm);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl", &L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl, &b_L1_NotBsc2_BptxAND_HcalHfCoincidencePm_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm", &L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm, &b_L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl", &L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl, &b_L1_NotBsc2_BptxAND_HcalHfMmOrPpOrPm_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfMmpOrMpp", &L1_NotBsc2_BptxAND_HcalHfMmpOrMpp, &b_L1_NotBsc2_BptxAND_HcalHfMmpOrMpp);
   fChain->SetBranchAddress("L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl", &L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl, &b_L1_NotBsc2_BptxAND_HcalHfMmpOrMpp_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasInnerThreshold1", &L1_NotBsc2_BscMinBiasInnerThreshold1, &b_L1_NotBsc2_BscMinBiasInnerThreshold1);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl", &L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl, &b_L1_NotBsc2_BscMinBiasInnerThreshold1_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasInnerThreshold2", &L1_NotBsc2_BscMinBiasInnerThreshold2, &b_L1_NotBsc2_BscMinBiasInnerThreshold2);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl", &L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl, &b_L1_NotBsc2_BscMinBiasInnerThreshold2_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasOR", &L1_NotBsc2_BscMinBiasOR, &b_L1_NotBsc2_BscMinBiasOR);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasOR_Prescl", &L1_NotBsc2_BscMinBiasOR_Prescl, &b_L1_NotBsc2_BscMinBiasOR_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasThreshold1", &L1_NotBsc2_BscMinBiasThreshold1, &b_L1_NotBsc2_BscMinBiasThreshold1);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasThreshold1_Prescl", &L1_NotBsc2_BscMinBiasThreshold1_Prescl, &b_L1_NotBsc2_BscMinBiasThreshold1_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasThreshold2", &L1_NotBsc2_BscMinBiasThreshold2, &b_L1_NotBsc2_BscMinBiasThreshold2);
   fChain->SetBranchAddress("L1_NotBsc2_BscMinBiasThreshold2_Prescl", &L1_NotBsc2_BscMinBiasThreshold2_Prescl, &b_L1_NotBsc2_BscMinBiasThreshold2_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_HcalHfCoincidencePm", &L1_NotBsc2_HcalHfCoincidencePm, &b_L1_NotBsc2_HcalHfCoincidencePm);
   fChain->SetBranchAddress("L1_NotBsc2_HcalHfCoincidencePm_Prescl", &L1_NotBsc2_HcalHfCoincidencePm_Prescl, &b_L1_NotBsc2_HcalHfCoincidencePm_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_HcalHfMmOrPpOrPm", &L1_NotBsc2_HcalHfMmOrPpOrPm, &b_L1_NotBsc2_HcalHfMmOrPpOrPm);
   fChain->SetBranchAddress("L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl", &L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl, &b_L1_NotBsc2_HcalHfMmOrPpOrPm_Prescl);
   fChain->SetBranchAddress("L1_NotBsc2_HcalHfMmpOrMpp", &L1_NotBsc2_HcalHfMmpOrMpp, &b_L1_NotBsc2_HcalHfMmpOrMpp);
   fChain->SetBranchAddress("L1_NotBsc2_HcalHfMmpOrMpp_Prescl", &L1_NotBsc2_HcalHfMmpOrMpp_Prescl, &b_L1_NotBsc2_HcalHfMmpOrMpp_Prescl);
   fChain->SetBranchAddress("L1_SingleEG12", &L1_SingleEG12, &b_L1_SingleEG12);
   fChain->SetBranchAddress("L1_SingleEG12_Prescl", &L1_SingleEG12_Prescl, &b_L1_SingleEG12_Prescl);
   fChain->SetBranchAddress("L1_SingleEG15", &L1_SingleEG15, &b_L1_SingleEG15);
   fChain->SetBranchAddress("L1_SingleEG15_Prescl", &L1_SingleEG15_Prescl, &b_L1_SingleEG15_Prescl);
   fChain->SetBranchAddress("L1_SingleEG18", &L1_SingleEG18, &b_L1_SingleEG18);
   fChain->SetBranchAddress("L1_SingleEG18_Prescl", &L1_SingleEG18_Prescl, &b_L1_SingleEG18_Prescl);
   fChain->SetBranchAddress("L1_SingleEG20", &L1_SingleEG20, &b_L1_SingleEG20);
   fChain->SetBranchAddress("L1_SingleEG20_Prescl", &L1_SingleEG20_Prescl, &b_L1_SingleEG20_Prescl);
   fChain->SetBranchAddress("L1_SingleEG22", &L1_SingleEG22, &b_L1_SingleEG22);
   fChain->SetBranchAddress("L1_SingleEG22_Prescl", &L1_SingleEG22_Prescl, &b_L1_SingleEG22_Prescl);
   fChain->SetBranchAddress("L1_SingleEG2_BptxAND", &L1_SingleEG2_BptxAND, &b_L1_SingleEG2_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG2_BptxAND_Prescl", &L1_SingleEG2_BptxAND_Prescl, &b_L1_SingleEG2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG30", &L1_SingleEG30, &b_L1_SingleEG30);
   fChain->SetBranchAddress("L1_SingleEG30_Prescl", &L1_SingleEG30_Prescl, &b_L1_SingleEG30_Prescl);
   fChain->SetBranchAddress("L1_SingleEG5_BptxAND", &L1_SingleEG5_BptxAND, &b_L1_SingleEG5_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG5_BptxAND_Prescl", &L1_SingleEG5_BptxAND_Prescl, &b_L1_SingleEG5_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG8_BptxAND", &L1_SingleEG8_BptxAND, &b_L1_SingleEG8_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG8_BptxAND_Prescl", &L1_SingleEG8_BptxAND_Prescl, &b_L1_SingleEG8_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet128", &L1_SingleJet128, &b_L1_SingleJet128);
   fChain->SetBranchAddress("L1_SingleJet128_Prescl", &L1_SingleJet128_Prescl, &b_L1_SingleJet128_Prescl);
   fChain->SetBranchAddress("L1_SingleJet128_BptxAND", &L1_SingleJet128_BptxAND, &b_L1_SingleJet128_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet128_BptxAND_Prescl", &L1_SingleJet128_BptxAND_Prescl, &b_L1_SingleJet128_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet128_Central_BptxAND", &L1_SingleJet128_Central_BptxAND, &b_L1_SingleJet128_Central_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet128_Central_BptxAND_Prescl", &L1_SingleJet128_Central_BptxAND_Prescl, &b_L1_SingleJet128_Central_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16", &L1_SingleJet16, &b_L1_SingleJet16);
   fChain->SetBranchAddress("L1_SingleJet16_Prescl", &L1_SingleJet16_Prescl, &b_L1_SingleJet16_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16_Central_NotETT140_BptxAND", &L1_SingleJet16_Central_NotETT140_BptxAND, &b_L1_SingleJet16_Central_NotETT140_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet16_Central_NotETT140_BptxAND_Prescl", &L1_SingleJet16_Central_NotETT140_BptxAND_Prescl, &b_L1_SingleJet16_Central_NotETT140_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo", &L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo, &b_L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl", &L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl, &b_L1_SingleJet20_Central_NotBptxOR_NotMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo", &L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo, &b_L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl", &L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl, &b_L1_SingleJet32_Central_NotBptxOR_NotMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36_BptxAND", &L1_SingleJet36_BptxAND, &b_L1_SingleJet36_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl", &L1_SingleJet36_BptxAND_Prescl, &b_L1_SingleJet36_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet52_BptxAND", &L1_SingleJet52_BptxAND, &b_L1_SingleJet52_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl", &L1_SingleJet52_BptxAND_Prescl, &b_L1_SingleJet52_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet52_Central_BptxAND", &L1_SingleJet52_Central_BptxAND, &b_L1_SingleJet52_Central_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet52_Central_BptxAND_Prescl", &L1_SingleJet52_Central_BptxAND_Prescl, &b_L1_SingleJet52_Central_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet68_BptxAND", &L1_SingleJet68_BptxAND, &b_L1_SingleJet68_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet68_BptxAND_Prescl", &L1_SingleJet68_BptxAND_Prescl, &b_L1_SingleJet68_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet80_BptxAND", &L1_SingleJet80_BptxAND, &b_L1_SingleJet80_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet80_BptxAND_Prescl", &L1_SingleJet80_BptxAND_Prescl, &b_L1_SingleJet80_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet80_Central", &L1_SingleJet80_Central, &b_L1_SingleJet80_Central);
   fChain->SetBranchAddress("L1_SingleJet80_Central_Prescl", &L1_SingleJet80_Central_Prescl, &b_L1_SingleJet80_Central_Prescl);
   fChain->SetBranchAddress("L1_SingleJet80_Central_BptxAND", &L1_SingleJet80_Central_BptxAND, &b_L1_SingleJet80_Central_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet80_Central_BptxAND_Prescl", &L1_SingleJet80_Central_BptxAND_Prescl, &b_L1_SingleJet80_Central_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet92", &L1_SingleJet92, &b_L1_SingleJet92);
   fChain->SetBranchAddress("L1_SingleJet92_Prescl", &L1_SingleJet92_Prescl, &b_L1_SingleJet92_Prescl);
   fChain->SetBranchAddress("L1_SingleJet92_BptxAND", &L1_SingleJet92_BptxAND, &b_L1_SingleJet92_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet92_BptxAND_Prescl", &L1_SingleJet92_BptxAND_Prescl, &b_L1_SingleJet92_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet92_Central", &L1_SingleJet92_Central, &b_L1_SingleJet92_Central);
   fChain->SetBranchAddress("L1_SingleJet92_Central_Prescl", &L1_SingleJet92_Central_Prescl, &b_L1_SingleJet92_Central_Prescl);
   fChain->SetBranchAddress("L1_SingleJet92_Central_BptxAND", &L1_SingleJet92_Central_BptxAND, &b_L1_SingleJet92_Central_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet92_Central_BptxAND_Prescl", &L1_SingleJet92_Central_BptxAND_Prescl, &b_L1_SingleJet92_Central_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu0_BptxAND", &L1_SingleMu0_BptxAND, &b_L1_SingleMu0_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu0_BptxAND_Prescl", &L1_SingleMu0_BptxAND_Prescl, &b_L1_SingleMu0_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu10", &L1_SingleMu10, &b_L1_SingleMu10);
   fChain->SetBranchAddress("L1_SingleMu10_Prescl", &L1_SingleMu10_Prescl, &b_L1_SingleMu10_Prescl);
   fChain->SetBranchAddress("L1_SingleMu14", &L1_SingleMu14, &b_L1_SingleMu14);
   fChain->SetBranchAddress("L1_SingleMu14_Prescl", &L1_SingleMu14_Prescl, &b_L1_SingleMu14_Prescl);
   fChain->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   fChain->SetBranchAddress("L1_SingleMu20_Prescl", &L1_SingleMu20_Prescl, &b_L1_SingleMu20_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3_BptxAND", &L1_SingleMu3_BptxAND, &b_L1_SingleMu3_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu3_BptxAND_Prescl", &L1_SingleMu3_BptxAND_Prescl, &b_L1_SingleMu3_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu5_BptxAND", &L1_SingleMu5_BptxAND, &b_L1_SingleMu5_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu5_BptxAND_Prescl", &L1_SingleMu5_BptxAND_Prescl, &b_L1_SingleMu5_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   fChain->SetBranchAddress("L1_SingleMu7_Prescl", &L1_SingleMu7_Prescl, &b_L1_SingleMu7_Prescl);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo", &L1_SingleMuBeamHalo, &b_L1_SingleMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo_Prescl", &L1_SingleMuBeamHalo_Prescl, &b_L1_SingleMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen_BptxAND", &L1_SingleMuOpen_BptxAND, &b_L1_SingleMuOpen_BptxAND);
   fChain->SetBranchAddress("L1_SingleMuOpen_BptxAND_Prescl", &L1_SingleMuOpen_BptxAND_Prescl, &b_L1_SingleMuOpen_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ZdcCaloMinus", &L1_ZdcCaloMinus, &b_L1_ZdcCaloMinus);
   fChain->SetBranchAddress("L1_ZdcCaloMinus_Prescl", &L1_ZdcCaloMinus_Prescl, &b_L1_ZdcCaloMinus_Prescl);
   fChain->SetBranchAddress("L1_ZdcCaloMinus_BptxAND", &L1_ZdcCaloMinus_BptxAND, &b_L1_ZdcCaloMinus_BptxAND);
   fChain->SetBranchAddress("L1_ZdcCaloMinus_BptxAND_Prescl", &L1_ZdcCaloMinus_BptxAND_Prescl, &b_L1_ZdcCaloMinus_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ZdcCaloPlus", &L1_ZdcCaloPlus, &b_L1_ZdcCaloPlus);
   fChain->SetBranchAddress("L1_ZdcCaloPlus_Prescl", &L1_ZdcCaloPlus_Prescl, &b_L1_ZdcCaloPlus_Prescl);
   fChain->SetBranchAddress("L1_ZdcCaloPlus_BptxAND", &L1_ZdcCaloPlus_BptxAND, &b_L1_ZdcCaloPlus_BptxAND);
   fChain->SetBranchAddress("L1_ZdcCaloPlus_BptxAND_Prescl", &L1_ZdcCaloPlus_BptxAND_Prescl, &b_L1_ZdcCaloPlus_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ZdcCaloPlus_ZdcCaloMinus", &L1_ZdcCaloPlus_ZdcCaloMinus, &b_L1_ZdcCaloPlus_ZdcCaloMinus);
   fChain->SetBranchAddress("L1_ZdcCaloPlus_ZdcCaloMinus_Prescl", &L1_ZdcCaloPlus_ZdcCaloMinus_Prescl, &b_L1_ZdcCaloPlus_ZdcCaloMinus_Prescl);
   fChain->SetBranchAddress("L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND", &L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND, &b_L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND);
   fChain->SetBranchAddress("L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl", &L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl, &b_L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   fChain->SetBranchAddress("L1_ZeroBias_Prescl", &L1_ZeroBias_Prescl, &b_L1_ZeroBias_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0", &L1Tech_BPTX_PreBPTX_v0, &b_L1Tech_BPTX_PreBPTX_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_PreBPTX.v0_Prescl", &L1Tech_BPTX_PreBPTX_v0_Prescl, &b_L1Tech_BPTX_PreBPTX_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_minus.v0", &L1Tech_BPTX_minus_v0, &b_L1Tech_BPTX_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_minus.v0_Prescl", &L1Tech_BPTX_minus_v0_Prescl, &b_L1Tech_BPTX_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0", &L1Tech_BPTX_minus_AND_not_plus_v0, &b_L1Tech_BPTX_minus_AND_not_plus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_minus_AND_not_plus.v0_Prescl", &L1Tech_BPTX_minus_AND_not_plus_v0_Prescl, &b_L1Tech_BPTX_minus_AND_not_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus.v0", &L1Tech_BPTX_plus_v0, &b_L1Tech_BPTX_plus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus.v0_Prescl", &L1Tech_BPTX_plus_v0_Prescl, &b_L1Tech_BPTX_plus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0", &L1Tech_BPTX_plus_AND_NOT_minus_v0, &b_L1Tech_BPTX_plus_AND_NOT_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_NOT_minus.v0_Prescl", &L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl, &b_L1Tech_BPTX_plus_AND_NOT_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0", &L1Tech_BPTX_plus_AND_minus_v0, &b_L1Tech_BPTX_plus_AND_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus.v0_Prescl", &L1Tech_BPTX_plus_AND_minus_v0_Prescl, &b_L1Tech_BPTX_plus_AND_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0", &L1Tech_BPTX_plus_AND_minus_instance1_v0, &b_L1Tech_BPTX_plus_AND_minus_instance1_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_AND_minus_instance1.v0_Prescl", &L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl, &b_L1Tech_BPTX_plus_AND_minus_instance1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0", &L1Tech_BPTX_plus_OR_minus_v0, &b_L1Tech_BPTX_plus_OR_minus_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_plus_OR_minus.v0_Prescl", &L1Tech_BPTX_plus_OR_minus_v0_Prescl, &b_L1Tech_BPTX_plus_OR_minus_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BPTX_quiet.v0", &L1Tech_BPTX_quiet_v0, &b_L1Tech_BPTX_quiet_v0);
   fChain->SetBranchAddress("L1Tech_BPTX_quiet.v0_Prescl", &L1Tech_BPTX_quiet_v0_Prescl, &b_L1Tech_BPTX_quiet_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_HighMultiplicity.v0", &L1Tech_BSC_HighMultiplicity_v0, &b_L1Tech_BSC_HighMultiplicity_v0);
   fChain->SetBranchAddress("L1Tech_BSC_HighMultiplicity.v0_Prescl", &L1Tech_BSC_HighMultiplicity_v0_Prescl, &b_L1Tech_BSC_HighMultiplicity_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_inner.v0", &L1Tech_BSC_halo_beam1_inner_v0, &b_L1Tech_BSC_halo_beam1_inner_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_inner.v0_Prescl", &L1Tech_BSC_halo_beam1_inner_v0_Prescl, &b_L1Tech_BSC_halo_beam1_inner_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_outer.v0", &L1Tech_BSC_halo_beam1_outer_v0, &b_L1Tech_BSC_halo_beam1_outer_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam1_outer.v0_Prescl", &L1Tech_BSC_halo_beam1_outer_v0_Prescl, &b_L1Tech_BSC_halo_beam1_outer_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_inner.v0", &L1Tech_BSC_halo_beam2_inner_v0, &b_L1Tech_BSC_halo_beam2_inner_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_inner.v0_Prescl", &L1Tech_BSC_halo_beam2_inner_v0_Prescl, &b_L1Tech_BSC_halo_beam2_inner_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_outer.v0", &L1Tech_BSC_halo_beam2_outer_v0, &b_L1Tech_BSC_halo_beam2_outer_v0);
   fChain->SetBranchAddress("L1Tech_BSC_halo_beam2_outer.v0_Prescl", &L1Tech_BSC_halo_beam2_outer_v0_Prescl, &b_L1Tech_BSC_halo_beam2_outer_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_OR.v0", &L1Tech_BSC_minBias_OR_v0, &b_L1Tech_BSC_minBias_OR_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_OR.v0_Prescl", &L1Tech_BSC_minBias_OR_v0_Prescl, &b_L1Tech_BSC_minBias_OR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold1.v0", &L1Tech_BSC_minBias_inner_threshold1_v0, &b_L1Tech_BSC_minBias_inner_threshold1_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold1.v0_Prescl", &L1Tech_BSC_minBias_inner_threshold1_v0_Prescl, &b_L1Tech_BSC_minBias_inner_threshold1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold2.v0", &L1Tech_BSC_minBias_inner_threshold2_v0, &b_L1Tech_BSC_minBias_inner_threshold2_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_inner_threshold2.v0_Prescl", &L1Tech_BSC_minBias_inner_threshold2_v0_Prescl, &b_L1Tech_BSC_minBias_inner_threshold2_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold1.v0", &L1Tech_BSC_minBias_threshold1_v0, &b_L1Tech_BSC_minBias_threshold1_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold1.v0_Prescl", &L1Tech_BSC_minBias_threshold1_v0_Prescl, &b_L1Tech_BSC_minBias_threshold1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold2.v0", &L1Tech_BSC_minBias_threshold2_v0, &b_L1Tech_BSC_minBias_threshold2_v0);
   fChain->SetBranchAddress("L1Tech_BSC_minBias_threshold2.v0_Prescl", &L1Tech_BSC_minBias_threshold2_v0_Prescl, &b_L1Tech_BSC_minBias_threshold2_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam1.v0", &L1Tech_BSC_splash_beam1_v0, &b_L1Tech_BSC_splash_beam1_v0);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam1.v0_Prescl", &L1Tech_BSC_splash_beam1_v0_Prescl, &b_L1Tech_BSC_splash_beam1_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam2.v0", &L1Tech_BSC_splash_beam2_v0, &b_L1Tech_BSC_splash_beam2_v0);
   fChain->SetBranchAddress("L1Tech_BSC_splash_beam2.v0_Prescl", &L1Tech_BSC_splash_beam2_v0_Prescl, &b_L1Tech_BSC_splash_beam2_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0", &L1Tech_CASTOR_HaloMuon_v0, &b_L1Tech_CASTOR_HaloMuon_v0);
   fChain->SetBranchAddress("L1Tech_CASTOR_HaloMuon.v0_Prescl", &L1Tech_CASTOR_HaloMuon_v0_Prescl, &b_L1Tech_CASTOR_HaloMuon_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_DT_GlobalOR.v0", &L1Tech_DT_GlobalOR_v0, &b_L1Tech_DT_GlobalOR_v0);
   fChain->SetBranchAddress("L1Tech_DT_GlobalOR.v0_Prescl", &L1Tech_DT_GlobalOR_v0_Prescl, &b_L1Tech_DT_GlobalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St1Sect45_down.v0", &L1Tech_FSC_St1Sect45_down_v0, &b_L1Tech_FSC_St1Sect45_down_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St1Sect45_down.v0_Prescl", &L1Tech_FSC_St1Sect45_down_v0_Prescl, &b_L1Tech_FSC_St1Sect45_down_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St1Sect45_upp.v0", &L1Tech_FSC_St1Sect45_upp_v0, &b_L1Tech_FSC_St1Sect45_upp_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St1Sect45_upp.v0_Prescl", &L1Tech_FSC_St1Sect45_upp_v0_Prescl, &b_L1Tech_FSC_St1Sect45_upp_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St1Sect56_down.v0", &L1Tech_FSC_St1Sect56_down_v0, &b_L1Tech_FSC_St1Sect56_down_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St1Sect56_down.v0_Prescl", &L1Tech_FSC_St1Sect56_down_v0_Prescl, &b_L1Tech_FSC_St1Sect56_down_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St1Sect56_upp.v0", &L1Tech_FSC_St1Sect56_upp_v0, &b_L1Tech_FSC_St1Sect56_upp_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St1Sect56_upp.v0_Prescl", &L1Tech_FSC_St1Sect56_upp_v0_Prescl, &b_L1Tech_FSC_St1Sect56_upp_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St2Sect45_down.v0", &L1Tech_FSC_St2Sect45_down_v0, &b_L1Tech_FSC_St2Sect45_down_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St2Sect45_down.v0_Prescl", &L1Tech_FSC_St2Sect45_down_v0_Prescl, &b_L1Tech_FSC_St2Sect45_down_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St2Sect45_upp.v0", &L1Tech_FSC_St2Sect45_upp_v0, &b_L1Tech_FSC_St2Sect45_upp_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St2Sect45_upp.v0_Prescl", &L1Tech_FSC_St2Sect45_upp_v0_Prescl, &b_L1Tech_FSC_St2Sect45_upp_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St2Sect56_down.v0", &L1Tech_FSC_St2Sect56_down_v0, &b_L1Tech_FSC_St2Sect56_down_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St2Sect56_down.v0_Prescl", &L1Tech_FSC_St2Sect56_down_v0_Prescl, &b_L1Tech_FSC_St2Sect56_down_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St2Sect56_upp.v0", &L1Tech_FSC_St2Sect56_upp_v0, &b_L1Tech_FSC_St2Sect56_upp_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St2Sect56_upp.v0_Prescl", &L1Tech_FSC_St2Sect56_upp_v0_Prescl, &b_L1Tech_FSC_St2Sect56_upp_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_downLeft.v0", &L1Tech_FSC_St3Sect45_downLeft_v0, &b_L1Tech_FSC_St3Sect45_downLeft_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_downLeft.v0_Prescl", &L1Tech_FSC_St3Sect45_downLeft_v0_Prescl, &b_L1Tech_FSC_St3Sect45_downLeft_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_downRight.v0", &L1Tech_FSC_St3Sect45_downRight_v0, &b_L1Tech_FSC_St3Sect45_downRight_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_downRight.v0_Prescl", &L1Tech_FSC_St3Sect45_downRight_v0_Prescl, &b_L1Tech_FSC_St3Sect45_downRight_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_uppLeft.v0", &L1Tech_FSC_St3Sect45_uppLeft_v0, &b_L1Tech_FSC_St3Sect45_uppLeft_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_uppLeft.v0_Prescl", &L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl, &b_L1Tech_FSC_St3Sect45_uppLeft_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_uppRight.v0", &L1Tech_FSC_St3Sect45_uppRight_v0, &b_L1Tech_FSC_St3Sect45_uppRight_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect45_uppRight.v0_Prescl", &L1Tech_FSC_St3Sect45_uppRight_v0_Prescl, &b_L1Tech_FSC_St3Sect45_uppRight_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect56_downLeft.v0", &L1Tech_FSC_St3Sect56_downLeft_v0, &b_L1Tech_FSC_St3Sect56_downLeft_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect56_downLeft.v0_Prescl", &L1Tech_FSC_St3Sect56_downLeft_v0_Prescl, &b_L1Tech_FSC_St3Sect56_downLeft_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect56_downRight.v0", &L1Tech_FSC_St3Sect56_downRight_v0, &b_L1Tech_FSC_St3Sect56_downRight_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect56_downRight.v0_Prescl", &L1Tech_FSC_St3Sect56_downRight_v0_Prescl, &b_L1Tech_FSC_St3Sect56_downRight_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect56_uppLeft.v0", &L1Tech_FSC_St3Sect56_uppLeft_v0, &b_L1Tech_FSC_St3Sect56_uppLeft_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect56_uppLeft.v0_Prescl", &L1Tech_FSC_St3Sect56_uppLeft_v0_Prescl, &b_L1Tech_FSC_St3Sect56_uppLeft_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect56_uppRight.v0", &L1Tech_FSC_St3Sect56_uppRight_v0, &b_L1Tech_FSC_St3Sect56_uppRight_v0);
   fChain->SetBranchAddress("L1Tech_FSC_St3Sect56_uppRight.v0_Prescl", &L1Tech_FSC_St3Sect56_uppRight_v0_Prescl, &b_L1Tech_FSC_St3Sect56_uppRight_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0", &L1Tech_HCAL_HBHE_totalOR_v0, &b_L1Tech_HCAL_HBHE_totalOR_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HBHE_totalOR.v0_Prescl", &L1Tech_HCAL_HBHE_totalOR_v0_Prescl, &b_L1Tech_HCAL_HBHE_totalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1", &L1Tech_HCAL_HF_MMP_or_MPP_v1, &b_L1Tech_HCAL_HF_MMP_or_MPP_v1);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_MMP_or_MPP.v1_Prescl", &L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl, &b_L1Tech_HCAL_HF_MMP_or_MPP_v1_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2", &L1Tech_HCAL_HF_coincidence_PM_v2, &b_L1Tech_HCAL_HF_coincidence_PM_v2);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_coincidence_PM.v2_Prescl", &L1Tech_HCAL_HF_coincidence_PM_v2_Prescl, &b_L1Tech_HCAL_HF_coincidence_PM_v2_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0", &L1Tech_HCAL_HF_single_channel_v0, &b_L1Tech_HCAL_HF_single_channel_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HF_single_channel.v0_Prescl", &L1Tech_HCAL_HF_single_channel_v0_Prescl, &b_L1Tech_HCAL_HF_single_channel_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0", &L1Tech_HCAL_HO_totalOR_v0, &b_L1Tech_HCAL_HO_totalOR_v0);
   fChain->SetBranchAddress("L1Tech_HCAL_HO_totalOR.v0_Prescl", &L1Tech_HCAL_HO_totalOR_v0_Prescl, &b_L1Tech_HCAL_HO_totalOR_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RB0_Cosmics.v0", &L1Tech_RPC_TTU_RB0_Cosmics_v0, &b_L1Tech_RPC_TTU_RB0_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RB0_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RB0_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus1_Cosmics.v0", &L1Tech_RPC_TTU_RBminus1_Cosmics_v0, &b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus1_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBminus1_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus2_Cosmics.v0", &L1Tech_RPC_TTU_RBminus2_Cosmics_v0, &b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBminus2_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBminus2_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0", &L1Tech_RPC_TTU_RBplus1_Cosmics_v0, &b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus1_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBplus1_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0", &L1Tech_RPC_TTU_RBplus2_Cosmics_v0, &b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBplus2_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_RBplus2_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBst1_collisions.v0", &L1Tech_RPC_TTU_RBst1_collisions_v0, &b_L1Tech_RPC_TTU_RBst1_collisions_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_RBst1_collisions.v0_Prescl", &L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl, &b_L1Tech_RPC_TTU_RBst1_collisions_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0", &L1Tech_RPC_TTU_barrel_Cosmics_v0, &b_L1Tech_RPC_TTU_barrel_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_barrel_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_barrel_Cosmics_v0_Prescl);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0", &L1Tech_RPC_TTU_pointing_Cosmics_v0, &b_L1Tech_RPC_TTU_pointing_Cosmics_v0);
   fChain->SetBranchAddress("L1Tech_RPC_TTU_pointing_Cosmics.v0_Prescl", &L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl, &b_L1Tech_RPC_TTU_pointing_Cosmics_v0_Prescl);
   Notify();
}

Bool_t HLT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HLT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HLT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HLT_cxx
