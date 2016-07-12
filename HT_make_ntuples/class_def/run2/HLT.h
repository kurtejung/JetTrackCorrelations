//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  7 22:33:03 2016 by ROOT version 6.02/13
// from TTree HltTree/
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root
//////////////////////////////////////////////////////////

#ifndef HLT_h
#define HLT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class HLT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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
   ULong64_t       Event;
   Int_t           LumiBlock;
   Int_t           Run;
   Int_t           Bx;
   Int_t           Orbit;
   Double_t        AvgInstDelLumi;
   Int_t           HLTriggerFirstPath;
   Int_t           HLTriggerFirstPath_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity60_v3;
   Int_t           HLT_PixelTracks_Multiplicity60_v3_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity85_v2;
   Int_t           HLT_PixelTracks_Multiplicity85_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity110_v2;
   Int_t           HLT_PixelTracks_Multiplicity110_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity135_v2;
   Int_t           HLT_PixelTracks_Multiplicity135_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity160_v2;
   Int_t           HLT_PixelTracks_Multiplicity160_v2_Prescl;
   Int_t           HLT_Physics_v2;
   Int_t           HLT_Physics_v2_Prescl;
   Int_t           DST_Physics_v1;
   Int_t           DST_Physics_v1_Prescl;
   Int_t           HLT_Random_v1;
   Int_t           HLT_Random_v1_Prescl;
   Int_t           HLT_EcalCalibration_v1;
   Int_t           HLT_EcalCalibration_v1_Prescl;
   Int_t           HLT_HcalCalibration_v1;
   Int_t           HLT_HcalCalibration_v1_Prescl;
   Int_t           AlCa_EcalPhiSym_v3;
   Int_t           AlCa_EcalPhiSym_v3_Prescl;
   Int_t           HLT_L1Tech6_BPTX_MinusOnly_v1;
   Int_t           HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl;
   Int_t           HLT_L1Tech5_BPTX_PlusOnly_v2;
   Int_t           HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl;
   Int_t           HLT_L1Tech7_NoBPTX_v1;
   Int_t           HLT_L1Tech7_NoBPTX_v1_Prescl;
   Int_t           HLT_L1TOTEM1_MinBias_v1;
   Int_t           HLT_L1TOTEM1_MinBias_v1_Prescl;
   Int_t           HLT_L1TOTEM2_ZeroBias_v1;
   Int_t           HLT_L1TOTEM2_ZeroBias_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1OR_v1;
   Int_t           HLT_L1MinimumBiasHF1OR_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF2OR_v1;
   Int_t           HLT_L1MinimumBiasHF2OR_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF1AND_v1;
   Int_t           HLT_L1MinimumBiasHF1AND_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF2AND_v1;
   Int_t           HLT_L1MinimumBiasHF2AND_v1_Prescl;
   Int_t           HLT_L1MinimumBiasHF2ORNoBptxGating_v1;
   Int_t           HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl;
   Int_t           AlCa_RPCMuonNoTriggers_v2;
   Int_t           AlCa_RPCMuonNoTriggers_v2_Prescl;
   Int_t           AlCa_RPCMuonNoHits_v2;
   Int_t           AlCa_RPCMuonNoHits_v2_Prescl;
   Int_t           AlCa_RPCMuonNormalisation_v2;
   Int_t           AlCa_RPCMuonNormalisation_v2_Prescl;
   Int_t           AlCa_LumiPixels_Random_v1;
   Int_t           AlCa_LumiPixels_Random_v1_Prescl;
   Int_t           AlCa_LumiPixels_ZeroBias_v2;
   Int_t           AlCa_LumiPixels_ZeroBias_v2_Prescl;

   Int_t           HLT_AK4CaloJet40_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet40_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet60_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet60_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet80_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet80_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet100_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet100_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet110_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet110_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet120_Eta5p1_v1;
   Int_t           HLT_AK4CaloJet120_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet150_v1;
   Int_t           HLT_AK4CaloJet150_v1_Prescl;
   Int_t           HLT_AK4PFJet40_Eta5p1_v1;
   Int_t           HLT_AK4PFJet40_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet60_Eta5p1_v1;
   Int_t           HLT_AK4PFJet60_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet80_Eta5p1_v1;
   Int_t           HLT_AK4PFJet80_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet100_Eta5p1_v1;
   Int_t           HLT_AK4PFJet100_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet110_Eta5p1_v1;
   Int_t           HLT_AK4PFJet110_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4PFJet120_Eta5p1_v1;
   Int_t           HLT_AK4PFJet120_Eta5p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet80_Jet35_Eta1p1_v1;
   Int_t           HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet80_Jet35_Eta0p7_v1;
   Int_t           HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl;
   Int_t           HLT_AK4CaloJet100_Jet35_Eta1p1_v1;
   Int_t           HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl;
   Int_t           HLT_AK4CaloJet100_Jet35_Eta0p7_v1;
   Int_t           HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl;
   Int_t           HLT_AK4CaloJet80_45_45_Eta2p1_v1;
   Int_t           HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl;

   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1;
   Int_t           HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl;
   Int_t           HLT_HIL1DoubleMu0_v1;
   Int_t           HLT_HIL1DoubleMu0_v1_Prescl;
   Int_t           HLT_HIL1DoubleMu10_v1;
   Int_t           HLT_HIL1DoubleMu10_v1_Prescl;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v1;
   Int_t           HLT_HIL2DoubleMu0_NHitQ_v1_Prescl;
   Int_t           HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1;
   Int_t           HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl;
   Int_t           HLT_HIL3DoubleMu0_OS_m7to14_v1;
   Int_t           HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl;
   Int_t           HLT_HIL2Mu3_NHitQ10_v1;
   Int_t           HLT_HIL2Mu3_NHitQ10_v1_Prescl;
   Int_t           HLT_HIL3Mu3_NHitQ15_v1;
   Int_t           HLT_HIL3Mu3_NHitQ15_v1_Prescl;
   Int_t           HLT_HIL2Mu5_NHitQ10_v1;
   Int_t           HLT_HIL2Mu5_NHitQ10_v1_Prescl;
   Int_t           HLT_HIL3Mu5_NHitQ15_v1;
   Int_t           HLT_HIL3Mu5_NHitQ15_v1_Prescl;
   Int_t           HLT_HIL2Mu7_NHitQ10_v1;
   Int_t           HLT_HIL2Mu7_NHitQ10_v1_Prescl;
   Int_t           HLT_HIL3Mu7_NHitQ15_v1;
   Int_t           HLT_HIL3Mu7_NHitQ15_v1_Prescl;
   Int_t           HLT_HIL2Mu15_v1;
   Int_t           HLT_HIL2Mu15_v1_Prescl;
   Int_t           HLT_HIL3Mu15_v1;
   Int_t           HLT_HIL3Mu15_v1_Prescl;
   Int_t           HLT_HIL2Mu20_v1;
   Int_t           HLT_HIL2Mu20_v1_Prescl;
   Int_t           HLT_HIL3Mu20_v1;
   Int_t           HLT_HIL3Mu20_v1_Prescl;
   Int_t           HLT_FullTrack18ForPPRef_v2;
   Int_t           HLT_FullTrack18ForPPRef_v2_Prescl;
   Int_t           HLT_FullTrack24ForPPRef_v2;
   Int_t           HLT_FullTrack24ForPPRef_v2_Prescl;
   Int_t           HLT_FullTrack34ForPPRef_v3;
   Int_t           HLT_FullTrack34ForPPRef_v3_Prescl;
   Int_t           HLT_FullTrack45ForPPRef_v2;
   Int_t           HLT_FullTrack45ForPPRef_v2_Prescl;
   Int_t           HLT_FullTrack53ForPPRef_v2;
   Int_t           HLT_FullTrack53ForPPRef_v2_Prescl;

   Int_t           HLT_HIL1CastorMediumJet_v1;
   Int_t           HLT_HIL1CastorMediumJet_v1_Prescl;
   Int_t           HLT_HICastorMediumJetPixel_SingleTrack_v1;
   Int_t           HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt8_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt15_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt20_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt30_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt40_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt50_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt60_v1;
   Int_t           HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl;

   Int_t           HLT_AK4PFBJetBCSV60_Eta2p1_v1;
   Int_t           HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFBJetBCSV80_Eta2p1_v1;
   Int_t           HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFDJet60_Eta2p1_v1;
   Int_t           HLT_AK4PFDJet60_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFDJet80_Eta2p1_v1;
   Int_t           HLT_AK4PFDJet80_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFBJetBSSV60_Eta2p1_v1;
   Int_t           HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl;
   Int_t           HLT_AK4PFBJetBSSV80_Eta2p1_v1;
   Int_t           HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl;
   Int_t           HLTriggerFinalPath;
   Int_t           HLTriggerFinalPath_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity60_v1;
   Int_t           HLT_PixelTracks_Multiplicity60_v1_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity60_v2;
   Int_t           HLT_PixelTracks_Multiplicity60_v2_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity85_v1;
   Int_t           HLT_PixelTracks_Multiplicity85_v1_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity110_v1;
   Int_t           HLT_PixelTracks_Multiplicity110_v1_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity135_v1;
   Int_t           HLT_PixelTracks_Multiplicity135_v1_Prescl;
   Int_t           HLT_PixelTracks_Multiplicity160_v1;
   Int_t           HLT_PixelTracks_Multiplicity160_v1_Prescl;
   Int_t           HLT_Physics_v1;
   Int_t           HLT_Physics_v1_Prescl;
   Int_t           HLT_L1Tech5_BPTX_PlusOnly_v1;
   Int_t           HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl;

   Int_t           HLT_FullTrack18ForPPRef_v1;
   Int_t           HLT_FullTrack18ForPPRef_v1_Prescl;
   Int_t           HLT_FullTrack18ForPPRef_v3;
   Int_t           HLT_FullTrack18ForPPRef_v3_Prescl;
   Int_t           HLT_FullTrack24ForPPRef_v1;
   Int_t           HLT_FullTrack24ForPPRef_v1_Prescl;
   Int_t           HLT_FullTrack24ForPPRef_v3;
   Int_t           HLT_FullTrack24ForPPRef_v3_Prescl;
   Int_t           HLT_FullTrack34ForPPRef_v1;
   Int_t           HLT_FullTrack34ForPPRef_v1_Prescl;
   Int_t           HLT_FullTrack34ForPPRef_v2;
   Int_t           HLT_FullTrack34ForPPRef_v2_Prescl;
   Int_t           HLT_FullTrack34ForPPRef_v4;
   Int_t           HLT_FullTrack34ForPPRef_v4_Prescl;
   Int_t           HLT_FullTrack45ForPPRef_v1;
   Int_t           HLT_FullTrack45ForPPRef_v1_Prescl;
   Int_t           HLT_FullTrack45ForPPRef_v3;
   Int_t           HLT_FullTrack45ForPPRef_v3_Prescl;
   Int_t           HLT_FullTrack53ForPPRef_v1;
   Int_t           HLT_FullTrack53ForPPRef_v1_Prescl;
   Int_t           HLT_FullTrack53ForPPRef_v3;
   Int_t           HLT_FullTrack53ForPPRef_v3_Prescl;
   Int_t           L1_AlwaysTrue;
   Int_t           L1_AlwaysTrue_Prescl;
   Int_t           L1_CastorHighJet_BptxAND;
   Int_t           L1_CastorHighJet_BptxAND_Prescl;
   Int_t           L1_CastorMediumJet_BptxAND;
   Int_t           L1_CastorMediumJet_BptxAND_Prescl;
   Int_t           L1_DoubleMu0_BptxAND;
   Int_t           L1_DoubleMu0_BptxAND_Prescl;
   Int_t           L1_DoubleMu0_MinimumBiasHF1_AND;
   Int_t           L1_DoubleMu0_MinimumBiasHF1_AND_Prescl;
   Int_t           L1_DoubleMu3_BptxAND;
   Int_t           L1_DoubleMu3_BptxAND_Prescl;
   Int_t           L1_DoubleMuOpen_BptxAND;
   Int_t           L1_DoubleMuOpen_BptxAND_Prescl;
   Int_t           L1_DoubleMuOpen_NotMinimumBiasHF2_AND;
   Int_t           L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl;
   Int_t           L1_DoubleJet20;
   Int_t           L1_DoubleJet20_Prescl;
   Int_t           L1_DoubleJet28_BptxAND;
   Int_t           L1_DoubleJet28_BptxAND_Prescl;
   Int_t           L1_DoubleJet32;
   Int_t           L1_DoubleJet32_Prescl;
   Int_t           L1_ETT130;
   Int_t           L1_ETT130_Prescl;
   Int_t           L1_ETT15_BptxAND;
   Int_t           L1_ETT15_BptxAND_Prescl;
   Int_t           L1_ETT20_BptxAND;
   Int_t           L1_ETT20_BptxAND_Prescl;
   Int_t           L1_ETT30_BptxAND;
   Int_t           L1_ETT30_BptxAND_Prescl;
   Int_t           L1_ETT40;
   Int_t           L1_ETT40_Prescl;
   Int_t           L1_ETT40_BptxAND;
   Int_t           L1_ETT40_BptxAND_Prescl;
   Int_t           L1_ETT45_BptxAND;
   Int_t           L1_ETT45_BptxAND_Prescl;
   Int_t           L1_ETT50_BptxAND;
   Int_t           L1_ETT50_BptxAND_Prescl;
   Int_t           L1_ETT55_BptxAND;
   Int_t           L1_ETT55_BptxAND_Prescl;
   Int_t           L1_ETT60_BptxAND;
   Int_t           L1_ETT60_BptxAND_Prescl;
   Int_t           L1_ETT70_BptxAND;
   Int_t           L1_ETT70_BptxAND_Prescl;
   Int_t           L1_ETT90_BptxAND;
   Int_t           L1_ETT90_BptxAND_Prescl;
   Int_t           L1_MinimumBiasHF1_AND;
   Int_t           L1_MinimumBiasHF1_AND_Prescl;
   Int_t           L1_MinimumBiasHF1_OR;
   Int_t           L1_MinimumBiasHF1_OR_Prescl;
   Int_t           L1_MinimumBiasHF1_XOR;
   Int_t           L1_MinimumBiasHF1_XOR_Prescl;
   Int_t           L1_MinimumBiasHF2_AND;
   Int_t           L1_MinimumBiasHF2_AND_Prescl;
   Int_t           L1_MinimumBiasHF2_OR;
   Int_t           L1_MinimumBiasHF2_OR_Prescl;
   Int_t           L1_MinimumBiasHF2_OR_NoBptxGating;
   Int_t           L1_MinimumBiasHF2_OR_NoBptxGating_Prescl;
   Int_t           L1_MuOpen_NotMinimumBiasHF2_AND;
   Int_t           L1_MuOpen_NotMinimumBiasHF2_AND_Prescl;
   Int_t           L1_NotMinimumBiasHF1_OR;
   Int_t           L1_NotMinimumBiasHF1_OR_Prescl;
   Int_t           L1_NotMinimumBiasHF2_AND;
   Int_t           L1_NotMinimumBiasHF2_AND_Prescl;
   Int_t           L1_NotZdcOR_BptxAND;
   Int_t           L1_NotZdcOR_BptxAND_Prescl;
   Int_t           L1_SingleEG12_BptxAND;
   Int_t           L1_SingleEG12_BptxAND_Prescl;
   Int_t           L1_SingleEG15_BptxAND;
   Int_t           L1_SingleEG15_BptxAND_Prescl;
   Int_t           L1_SingleEG20;
   Int_t           L1_SingleEG20_Prescl;
   Int_t           L1_SingleEG20_BptxAND;
   Int_t           L1_SingleEG20_BptxAND_Prescl;
   Int_t           L1_SingleEG2_BptxAND;
   Int_t           L1_SingleEG2_BptxAND_Prescl;
   Int_t           L1_SingleEG30_BptxAND;
   Int_t           L1_SingleEG30_BptxAND_Prescl;
   Int_t           L1_SingleEG5;
   Int_t           L1_SingleEG5_Prescl;
   Int_t           L1_SingleEG5_BptxAND;
   Int_t           L1_SingleEG5_BptxAND_Prescl;
   Int_t           L1_SingleEG7_BptxAND;
   Int_t           L1_SingleEG7_BptxAND_Prescl;
   Int_t           L1_SingleMu12_BptxAND;
   Int_t           L1_SingleMu12_BptxAND_Prescl;
   Int_t           L1_SingleMu16_BptxAND;
   Int_t           L1_SingleMu16_BptxAND_Prescl;
   Int_t           L1_SingleMu16_MinimumBiasHF1_AND;
   Int_t           L1_SingleMu16_MinimumBiasHF1_AND_Prescl;
   Int_t           L1_SingleMu3_BptxAND;
   Int_t           L1_SingleMu3_BptxAND_Prescl;
   Int_t           L1_SingleMu3_MinimumBiasHF1_AND;
   Int_t           L1_SingleMu3_MinimumBiasHF1_AND_Prescl;
   Int_t           L1_SingleMu3p5;
   Int_t           L1_SingleMu3p5_Prescl;
   Int_t           L1_SingleMu5_BptxAND;
   Int_t           L1_SingleMu5_BptxAND_Prescl;
   Int_t           L1_SingleMu7_BptxAND;
   Int_t           L1_SingleMu7_BptxAND_Prescl;
   Int_t           L1_SingleMuBeamHalo;
   Int_t           L1_SingleMuBeamHalo_Prescl;
   Int_t           L1_SingleMuOpen;
   Int_t           L1_SingleMuOpen_Prescl;
   Int_t           L1_SingleMuOpen_BptxAND;
   Int_t           L1_SingleMuOpen_BptxAND_Prescl;
   Int_t           L1_SingleMuOpen_NotBptxOR;
   Int_t           L1_SingleMuOpen_NotBptxOR_Prescl;
   Int_t           L1_SingleJet12_BptxAND;
   Int_t           L1_SingleJet12_BptxAND_Prescl;
   Int_t           L1_SingleJet16;
   Int_t           L1_SingleJet16_Prescl;
   Int_t           L1_SingleJet16_BptxAND;
   Int_t           L1_SingleJet16_BptxAND_Prescl;
   Int_t           L1_SingleJet200;
   Int_t           L1_SingleJet200_Prescl;
   Int_t           L1_SingleJet20_BptxAND;
   Int_t           L1_SingleJet20_BptxAND_Prescl;
   Int_t           L1_SingleJet24_BptxAND;
   Int_t           L1_SingleJet24_BptxAND_Prescl;
   Int_t           L1_SingleJet28_BptxAND;
   Int_t           L1_SingleJet28_BptxAND_Prescl;
   Int_t           L1_SingleJet32_BptxAND;
   Int_t           L1_SingleJet32_BptxAND_Prescl;
   Int_t           L1_SingleJet36;
   Int_t           L1_SingleJet36_Prescl;
   Int_t           L1_SingleJet36_BptxAND;
   Int_t           L1_SingleJet36_BptxAND_Prescl;
   Int_t           L1_SingleJet40_BptxAND;
   Int_t           L1_SingleJet40_BptxAND_Prescl;
   Int_t           L1_SingleJet44_BptxAND;
   Int_t           L1_SingleJet44_BptxAND_Prescl;
   Int_t           L1_SingleJet48_BptxAND;
   Int_t           L1_SingleJet48_BptxAND_Prescl;
   Int_t           L1_SingleJet52_BptxAND;
   Int_t           L1_SingleJet52_BptxAND_Prescl;
   Int_t           L1_SingleJet60_BptxAND;
   Int_t           L1_SingleJet60_BptxAND_Prescl;
   Int_t           L1_SingleJet68;
   Int_t           L1_SingleJet68_Prescl;
   Int_t           L1_SingleJet68_BptxAND;
   Int_t           L1_SingleJet68_BptxAND_Prescl;
   Int_t           L1_SingleJet8_BptxAND;
   Int_t           L1_SingleJet8_BptxAND_Prescl;
   Int_t           L1_SingleJetC20_NotBptxOR;
   Int_t           L1_SingleJetC20_NotBptxOR_Prescl;
   Int_t           L1_SingleJetC32_NotBptxOR;
   Int_t           L1_SingleJetC32_NotBptxOR_Prescl;
   Int_t           L1_TOTEM_0;
   Int_t           L1_TOTEM_0_Prescl;
   Int_t           L1_TOTEM_1;
   Int_t           L1_TOTEM_1_Prescl;
   Int_t           L1_TOTEM_2;
   Int_t           L1_TOTEM_2_Prescl;
   Int_t           L1_TOTEM_3;
   Int_t           L1_TOTEM_3_Prescl;
   Int_t           L1_ZdcOR;
   Int_t           L1_ZdcOR_Prescl;
   Int_t           L1_ZdcOR_BptxAND;
   Int_t           L1_ZdcOR_BptxAND_Prescl;
   Int_t           L1_ZdcXOR;
   Int_t           L1_ZdcXOR_Prescl;
   Int_t           L1_ZdcXOR_BptxAND;
   Int_t           L1_ZdcXOR_BptxAND_Prescl;
   Int_t           L1_ZeroBias;
   Int_t           L1_ZeroBias_Prescl;


   ///// PbPb Trigger Decisions ////

   Int_t           HLT_HIPuAK4CaloJet40_Eta5p1_v1;
   Int_t           HLT_HIPuAK4CaloJet40_Eta5p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet60_Eta5p1_v1;
   Int_t           HLT_HIPuAK4CaloJet60_Eta5p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet80_Eta5p1_v1;
   Int_t           HLT_HIPuAK4CaloJet80_Eta5p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1;
   Int_t           HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet100_Eta5p1_v1;
   Int_t           HLT_HIPuAK4CaloJet100_Eta5p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet110_Eta5p1_v1;
   Int_t           HLT_HIPuAK4CaloJet110_Eta5p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet120_Eta5p1_v1;
   Int_t           HLT_HIPuAK4CaloJet120_Eta5p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet150_Eta5p1_v1;
   Int_t           HLT_HIPuAK4CaloJet150_Eta5p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1;
   Int_t           HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1;
   Int_t           HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1;
   Int_t           HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1;
   Int_t           HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1;
   Int_t           HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1;
   Int_t           HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1;
   Int_t           HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1;
   Int_t           HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1;
   Int_t           HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1;
   Int_t           HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1;
   Int_t           HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1;
   Int_t           HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1;
   Int_t           HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloDJet60_Eta2p1_v1;
   Int_t           HLT_HIPuAK4CaloDJet60_Eta2p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloDJet80_Eta2p1_v1;
   Int_t           HLT_HIPuAK4CaloDJet80_Eta2p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1;
   Int_t           HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1;
   Int_t           HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1;
   Int_t           HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1_Prescl;
   Int_t           HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1;
   Int_t           HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1_Prescl;
   

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
   TBranch        *b_Event;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Bx;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_AvgInstDelLumi;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLTriggerFirstPath_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v3;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v3_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity110_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity110_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity135_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity135_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity160_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity160_v2_Prescl;   //!
   TBranch        *b_HLT_Physics_v2;   //!
   TBranch        *b_HLT_Physics_v2_Prescl;   //!
   TBranch        *b_DST_Physics_v1;   //!
   TBranch        *b_DST_Physics_v1_Prescl;   //!
   TBranch        *b_HLT_Random_v1;   //!
   TBranch        *b_HLT_Random_v1_Prescl;   //!
   TBranch        *b_HLT_EcalCalibration_v1;   //!
   TBranch        *b_HLT_EcalCalibration_v1_Prescl;   //!
   TBranch        *b_HLT_HcalCalibration_v1;   //!
   TBranch        *b_HLT_HcalCalibration_v1_Prescl;   //!
   TBranch        *b_AlCa_EcalPhiSym_v3;   //!
   TBranch        *b_AlCa_EcalPhiSym_v3_Prescl;   //!
   TBranch        *b_HLT_L1Tech6_BPTX_MinusOnly_v1;   //!
   TBranch        *b_HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl;   //!
   TBranch        *b_HLT_L1Tech5_BPTX_PlusOnly_v2;   //!
   TBranch        *b_HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl;   //!
   TBranch        *b_HLT_L1Tech7_NoBPTX_v1;   //!
   TBranch        *b_HLT_L1Tech7_NoBPTX_v1_Prescl;   //!
   TBranch        *b_HLT_L1TOTEM1_MinBias_v1;   //!
   TBranch        *b_HLT_L1TOTEM1_MinBias_v1_Prescl;   //!
   TBranch        *b_HLT_L1TOTEM2_ZeroBias_v1;   //!
   TBranch        *b_HLT_L1TOTEM2_ZeroBias_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1OR_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2OR_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2OR_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1AND_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF1AND_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2AND_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2AND_v1_Prescl;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2ORNoBptxGating_v1;   //!
   TBranch        *b_HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers_v2;   //!
   TBranch        *b_AlCa_RPCMuonNoTriggers_v2_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNoHits_v2;   //!
   TBranch        *b_AlCa_RPCMuonNoHits_v2_Prescl;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation_v2;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation_v2_Prescl;   //!
   TBranch        *b_AlCa_LumiPixels_Random_v1;   //!
   TBranch        *b_AlCa_LumiPixels_Random_v1_Prescl;   //!
   TBranch        *b_AlCa_LumiPixels_ZeroBias_v2;   //!
   TBranch        *b_AlCa_LumiPixels_ZeroBias_v2_Prescl;   //!

   TBranch        *b_HLT_AK4CaloJet40_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet40_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet60_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet60_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet80_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet80_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet100_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet100_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet110_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet110_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet120_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet120_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet150_v1;   //!
   TBranch        *b_HLT_AK4CaloJet150_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet40_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet40_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet60_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet60_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet80_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet80_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet100_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet100_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet110_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet110_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFJet120_Eta5p1_v1;   //!
   TBranch        *b_HLT_AK4PFJet120_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet80_Jet35_Eta1p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet80_Jet35_Eta0p7_v1;   //!
   TBranch        *b_HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet100_Jet35_Eta1p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet100_Jet35_Eta0p7_v1;   //!
   TBranch        *b_HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl;   //!
   TBranch        *b_HLT_AK4CaloJet80_45_45_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl;   //!

   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1;   //!
   TBranch        *b_HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_v1;   //!
   TBranch        *b_HLT_HIL1DoubleMu0_v1_Prescl;   //!
   TBranch        *b_HLT_HIL1DoubleMu10_v1;   //!
   TBranch        *b_HLT_HIL1DoubleMu10_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v1;   //!
   TBranch        *b_HLT_HIL2DoubleMu0_NHitQ_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3DoubleMu0_OS_m7to14_v1;   //!
   TBranch        *b_HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ10_v1;   //!
   TBranch        *b_HLT_HIL2Mu3_NHitQ10_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu3_NHitQ15_v1;   //!
   TBranch        *b_HLT_HIL3Mu3_NHitQ15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu5_NHitQ10_v1;   //!
   TBranch        *b_HLT_HIL2Mu5_NHitQ10_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu5_NHitQ15_v1;   //!
   TBranch        *b_HLT_HIL3Mu5_NHitQ15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu7_NHitQ10_v1;   //!
   TBranch        *b_HLT_HIL2Mu7_NHitQ10_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu7_NHitQ15_v1;   //!
   TBranch        *b_HLT_HIL3Mu7_NHitQ15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu15_v1;   //!
   TBranch        *b_HLT_HIL2Mu15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu15_v1;   //!
   TBranch        *b_HLT_HIL3Mu15_v1_Prescl;   //!
   TBranch        *b_HLT_HIL2Mu20_v1;   //!
   TBranch        *b_HLT_HIL2Mu20_v1_Prescl;   //!
   TBranch        *b_HLT_HIL3Mu20_v1;   //!
   TBranch        *b_HLT_HIL3Mu20_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v2_Prescl;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v2_Prescl;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v3_Prescl;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v2_Prescl;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v2_Prescl;   //!

   TBranch        *b_HLT_HIL1CastorMediumJet_v1;   //!
   TBranch        *b_HLT_HIL1CastorMediumJet_v1_Prescl;   //!
   TBranch        *b_HLT_HICastorMediumJetPixel_SingleTrack_v1;   //!
   TBranch        *b_HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt8_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt15_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt20_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt30_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt40_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt50_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt60_v1;   //!
   TBranch        *b_HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl;   //!

   TBranch        *b_HLT_AK4PFBJetBCSV60_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFBJetBCSV80_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFDJet60_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFDJet60_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFDJet80_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFDJet80_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFBJetBSSV60_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_AK4PFBJetBSSV80_Eta2p1_v1;   //!
   TBranch        *b_HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFinalPath_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v1_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v2;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity60_v2_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity85_v1_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity110_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity110_v1_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity135_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity135_v1_Prescl;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity160_v1;   //!
   TBranch        *b_HLT_PixelTracks_Multiplicity160_v1_Prescl;   //!
   TBranch        *b_HLT_Physics_v1;   //!
   TBranch        *b_HLT_Physics_v1_Prescl;   //!
   TBranch        *b_HLT_L1Tech5_BPTX_PlusOnly_v1;   //!
   TBranch        *b_HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl;   //!

   TBranch        *b_HLT_FullTrack18ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack18ForPPRef_v3_Prescl;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack24ForPPRef_v3_Prescl;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v2;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v2_Prescl;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v4;   //!
   TBranch        *b_HLT_FullTrack34ForPPRef_v4_Prescl;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack45ForPPRef_v3_Prescl;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v1;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v1_Prescl;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v3;   //!
   TBranch        *b_HLT_FullTrack53ForPPRef_v3_Prescl;   //!
   TBranch        *b_L1_AlwaysTrue;   //!
   TBranch        *b_L1_AlwaysTrue_Prescl;   //!
   TBranch        *b_L1_CastorHighJet_BptxAND;   //!
   TBranch        *b_L1_CastorHighJet_BptxAND_Prescl;   //!
   TBranch        *b_L1_CastorMediumJet_BptxAND;   //!
   TBranch        *b_L1_CastorMediumJet_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_BptxAND;   //!
   TBranch        *b_L1_DoubleMu0_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMu0_MinimumBiasHF1_AND;   //!
   TBranch        *b_L1_DoubleMu0_MinimumBiasHF1_AND_Prescl;   //!
   TBranch        *b_L1_DoubleMu3_BptxAND;   //!
   TBranch        *b_L1_DoubleMu3_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND;   //!
   TBranch        *b_L1_DoubleMuOpen_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleMuOpen_NotMinimumBiasHF2_AND;   //!
   TBranch        *b_L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl;   //!
   TBranch        *b_L1_DoubleJet20;   //!
   TBranch        *b_L1_DoubleJet20_Prescl;   //!
   TBranch        *b_L1_DoubleJet28_BptxAND;   //!
   TBranch        *b_L1_DoubleJet28_BptxAND_Prescl;   //!
   TBranch        *b_L1_DoubleJet32;   //!
   TBranch        *b_L1_DoubleJet32_Prescl;   //!
   TBranch        *b_L1_ETT130;   //!
   TBranch        *b_L1_ETT130_Prescl;   //!
   TBranch        *b_L1_ETT15_BptxAND;   //!
   TBranch        *b_L1_ETT15_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT20_BptxAND;   //!
   TBranch        *b_L1_ETT20_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT30_BptxAND;   //!
   TBranch        *b_L1_ETT30_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT40;   //!
   TBranch        *b_L1_ETT40_Prescl;   //!
   TBranch        *b_L1_ETT40_BptxAND;   //!
   TBranch        *b_L1_ETT40_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT45_BptxAND;   //!
   TBranch        *b_L1_ETT45_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT50_BptxAND;   //!
   TBranch        *b_L1_ETT50_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT55_BptxAND;   //!
   TBranch        *b_L1_ETT55_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT60_BptxAND;   //!
   TBranch        *b_L1_ETT60_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT70_BptxAND;   //!
   TBranch        *b_L1_ETT70_BptxAND_Prescl;   //!
   TBranch        *b_L1_ETT90_BptxAND;   //!
   TBranch        *b_L1_ETT90_BptxAND_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF1_AND;   //!
   TBranch        *b_L1_MinimumBiasHF1_AND_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF1_OR;   //!
   TBranch        *b_L1_MinimumBiasHF1_OR_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF1_XOR;   //!
   TBranch        *b_L1_MinimumBiasHF1_XOR_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF2_AND;   //!
   TBranch        *b_L1_MinimumBiasHF2_AND_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF2_OR;   //!
   TBranch        *b_L1_MinimumBiasHF2_OR_Prescl;   //!
   TBranch        *b_L1_MinimumBiasHF2_OR_NoBptxGating;   //!
   TBranch        *b_L1_MinimumBiasHF2_OR_NoBptxGating_Prescl;   //!
   TBranch        *b_L1_MuOpen_NotMinimumBiasHF2_AND;   //!
   TBranch        *b_L1_MuOpen_NotMinimumBiasHF2_AND_Prescl;   //!
   TBranch        *b_L1_NotMinimumBiasHF1_OR;   //!
   TBranch        *b_L1_NotMinimumBiasHF1_OR_Prescl;   //!
   TBranch        *b_L1_NotMinimumBiasHF2_AND;   //!
   TBranch        *b_L1_NotMinimumBiasHF2_AND_Prescl;   //!
   TBranch        *b_L1_NotZdcOR_BptxAND;   //!
   TBranch        *b_L1_NotZdcOR_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG12_BptxAND;   //!
   TBranch        *b_L1_SingleEG12_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG15_BptxAND;   //!
   TBranch        *b_L1_SingleEG15_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG20;   //!
   TBranch        *b_L1_SingleEG20_Prescl;   //!
   TBranch        *b_L1_SingleEG20_BptxAND;   //!
   TBranch        *b_L1_SingleEG20_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG2_BptxAND;   //!
   TBranch        *b_L1_SingleEG2_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG30_BptxAND;   //!
   TBranch        *b_L1_SingleEG30_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG5;   //!
   TBranch        *b_L1_SingleEG5_Prescl;   //!
   TBranch        *b_L1_SingleEG5_BptxAND;   //!
   TBranch        *b_L1_SingleEG5_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleEG7_BptxAND;   //!
   TBranch        *b_L1_SingleEG7_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu12_BptxAND;   //!
   TBranch        *b_L1_SingleMu12_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu16_BptxAND;   //!
   TBranch        *b_L1_SingleMu16_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu16_MinimumBiasHF1_AND;   //!
   TBranch        *b_L1_SingleMu16_MinimumBiasHF1_AND_Prescl;   //!
   TBranch        *b_L1_SingleMu3_BptxAND;   //!
   TBranch        *b_L1_SingleMu3_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu3_MinimumBiasHF1_AND;   //!
   TBranch        *b_L1_SingleMu3_MinimumBiasHF1_AND_Prescl;   //!
   TBranch        *b_L1_SingleMu3p5;   //!
   TBranch        *b_L1_SingleMu3p5_Prescl;   //!
   TBranch        *b_L1_SingleMu5_BptxAND;   //!
   TBranch        *b_L1_SingleMu5_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMu7_BptxAND;   //!
   TBranch        *b_L1_SingleMu7_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMuBeamHalo;   //!
   TBranch        *b_L1_SingleMuBeamHalo_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen_BptxAND;   //!
   TBranch        *b_L1_SingleMuOpen_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleJet12_BptxAND;   //!
   TBranch        *b_L1_SingleJet12_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet16;   //!
   TBranch        *b_L1_SingleJet16_Prescl;   //!
   TBranch        *b_L1_SingleJet16_BptxAND;   //!
   TBranch        *b_L1_SingleJet16_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet200;   //!
   TBranch        *b_L1_SingleJet200_Prescl;   //!
   TBranch        *b_L1_SingleJet20_BptxAND;   //!
   TBranch        *b_L1_SingleJet20_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet24_BptxAND;   //!
   TBranch        *b_L1_SingleJet24_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet28_BptxAND;   //!
   TBranch        *b_L1_SingleJet28_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet32_BptxAND;   //!
   TBranch        *b_L1_SingleJet32_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet36;   //!
   TBranch        *b_L1_SingleJet36_Prescl;   //!
   TBranch        *b_L1_SingleJet36_BptxAND;   //!
   TBranch        *b_L1_SingleJet36_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet40_BptxAND;   //!
   TBranch        *b_L1_SingleJet40_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet44_BptxAND;   //!
   TBranch        *b_L1_SingleJet44_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet48_BptxAND;   //!
   TBranch        *b_L1_SingleJet48_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet52_BptxAND;   //!
   TBranch        *b_L1_SingleJet52_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet60_BptxAND;   //!
   TBranch        *b_L1_SingleJet60_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet68;   //!
   TBranch        *b_L1_SingleJet68_Prescl;   //!
   TBranch        *b_L1_SingleJet68_BptxAND;   //!
   TBranch        *b_L1_SingleJet68_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet8_BptxAND;   //!
   TBranch        *b_L1_SingleJet8_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJetC20_NotBptxOR;   //!
   TBranch        *b_L1_SingleJetC20_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_SingleJetC32_NotBptxOR;   //!
   TBranch        *b_L1_SingleJetC32_NotBptxOR_Prescl;   //!
   TBranch        *b_L1_TOTEM_0;   //!
   TBranch        *b_L1_TOTEM_0_Prescl;   //!
   TBranch        *b_L1_TOTEM_1;   //!
   TBranch        *b_L1_TOTEM_1_Prescl;   //!
   TBranch        *b_L1_TOTEM_2;   //!
   TBranch        *b_L1_TOTEM_2_Prescl;   //!
   TBranch        *b_L1_TOTEM_3;   //!
   TBranch        *b_L1_TOTEM_3_Prescl;   //!
   TBranch        *b_L1_ZdcOR;   //!
   TBranch        *b_L1_ZdcOR_Prescl;   //!
   TBranch        *b_L1_ZdcOR_BptxAND;   //!
   TBranch        *b_L1_ZdcOR_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZdcXOR;   //!
   TBranch        *b_L1_ZdcXOR_Prescl;   //!
   TBranch        *b_L1_ZdcXOR_BptxAND;   //!
   TBranch        *b_L1_ZdcXOR_BptxAND_Prescl;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_Prescl;   //!


   //// PbPb TRIGGER BRANCHES //////
      TBranch        *b_HLT_HIPuAK4CaloJet40_Eta5p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet40_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet60_Eta5p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet60_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Eta5p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Eta5p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet110_Eta5p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet110_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet120_Eta5p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet120_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet150_Eta5p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet150_Eta5p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloDJet60_Eta2p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloDJet60_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloDJet80_Eta2p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloDJet80_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1_Prescl;   //!
   TBranch        *b_HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1;   //!
   TBranch        *b_HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1_Prescl;   //!

   HLT(TTree *tree=0, bool ispp=0);
   virtual ~HLT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, bool ispp);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

HLT::HLT(TTree *tree, bool ispp) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://eoscms//eos/cms/store/group/phys_heavyions/kjung/pp5TeV_HighPtJet80AOD_May2016_IVFVtx/HighPtJet80/crab_pp5Tev_HighPtJet80_AOD_recalibJP_IVFVtx/160507_222047/0000/HiForestAOD_data_89.root:/hltanalysis");
      dir->GetObject("HltTree",tree);

   }
   Init(tree, ispp);
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
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HLT::Init(TTree *tree, bool ispp)
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
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Bx", &Bx, &b_Bx);
   fChain->SetBranchAddress("Orbit", &Orbit, &b_Orbit);
   fChain->SetBranchAddress("AvgInstDelLumi", &AvgInstDelLumi, &b_AvgInstDelLumi);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("HLTriggerFirstPath_Prescl", &HLTriggerFirstPath_Prescl, &b_HLTriggerFirstPath_Prescl);

   if(ispp){
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v3", &HLT_PixelTracks_Multiplicity60_v3, &b_HLT_PixelTracks_Multiplicity60_v3);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v3_Prescl", &HLT_PixelTracks_Multiplicity60_v3_Prescl, &b_HLT_PixelTracks_Multiplicity60_v3_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85_v2", &HLT_PixelTracks_Multiplicity85_v2, &b_HLT_PixelTracks_Multiplicity85_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85_v2_Prescl", &HLT_PixelTracks_Multiplicity85_v2_Prescl, &b_HLT_PixelTracks_Multiplicity85_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110_v2", &HLT_PixelTracks_Multiplicity110_v2, &b_HLT_PixelTracks_Multiplicity110_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110_v2_Prescl", &HLT_PixelTracks_Multiplicity110_v2_Prescl, &b_HLT_PixelTracks_Multiplicity110_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135_v2", &HLT_PixelTracks_Multiplicity135_v2, &b_HLT_PixelTracks_Multiplicity135_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135_v2_Prescl", &HLT_PixelTracks_Multiplicity135_v2_Prescl, &b_HLT_PixelTracks_Multiplicity135_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160_v2", &HLT_PixelTracks_Multiplicity160_v2, &b_HLT_PixelTracks_Multiplicity160_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160_v2_Prescl", &HLT_PixelTracks_Multiplicity160_v2_Prescl, &b_HLT_PixelTracks_Multiplicity160_v2_Prescl);
   fChain->SetBranchAddress("HLT_Physics_v2", &HLT_Physics_v2, &b_HLT_Physics_v2);
   fChain->SetBranchAddress("HLT_Physics_v2_Prescl", &HLT_Physics_v2_Prescl, &b_HLT_Physics_v2_Prescl);
   fChain->SetBranchAddress("DST_Physics_v1", &DST_Physics_v1, &b_DST_Physics_v1);
   fChain->SetBranchAddress("DST_Physics_v1_Prescl", &DST_Physics_v1_Prescl, &b_DST_Physics_v1_Prescl);
   fChain->SetBranchAddress("HLT_Random_v1", &HLT_Random_v1, &b_HLT_Random_v1);
   fChain->SetBranchAddress("HLT_Random_v1_Prescl", &HLT_Random_v1_Prescl, &b_HLT_Random_v1_Prescl);
   fChain->SetBranchAddress("HLT_EcalCalibration_v1", &HLT_EcalCalibration_v1, &b_HLT_EcalCalibration_v1);
   fChain->SetBranchAddress("HLT_EcalCalibration_v1_Prescl", &HLT_EcalCalibration_v1_Prescl, &b_HLT_EcalCalibration_v1_Prescl);
   fChain->SetBranchAddress("HLT_HcalCalibration_v1", &HLT_HcalCalibration_v1, &b_HLT_HcalCalibration_v1);
   fChain->SetBranchAddress("HLT_HcalCalibration_v1_Prescl", &HLT_HcalCalibration_v1_Prescl, &b_HLT_HcalCalibration_v1_Prescl);
   fChain->SetBranchAddress("AlCa_EcalPhiSym_v3", &AlCa_EcalPhiSym_v3, &b_AlCa_EcalPhiSym_v3);
   fChain->SetBranchAddress("AlCa_EcalPhiSym_v3_Prescl", &AlCa_EcalPhiSym_v3_Prescl, &b_AlCa_EcalPhiSym_v3_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech6_BPTX_MinusOnly_v1", &HLT_L1Tech6_BPTX_MinusOnly_v1, &b_HLT_L1Tech6_BPTX_MinusOnly_v1);
   fChain->SetBranchAddress("HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl", &HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl, &b_HLT_L1Tech6_BPTX_MinusOnly_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech5_BPTX_PlusOnly_v2", &HLT_L1Tech5_BPTX_PlusOnly_v2, &b_HLT_L1Tech5_BPTX_PlusOnly_v2);
   fChain->SetBranchAddress("HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl", &HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl, &b_HLT_L1Tech5_BPTX_PlusOnly_v2_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech7_NoBPTX_v1", &HLT_L1Tech7_NoBPTX_v1, &b_HLT_L1Tech7_NoBPTX_v1);
   fChain->SetBranchAddress("HLT_L1Tech7_NoBPTX_v1_Prescl", &HLT_L1Tech7_NoBPTX_v1_Prescl, &b_HLT_L1Tech7_NoBPTX_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1TOTEM1_MinBias_v1", &HLT_L1TOTEM1_MinBias_v1, &b_HLT_L1TOTEM1_MinBias_v1);
   fChain->SetBranchAddress("HLT_L1TOTEM1_MinBias_v1_Prescl", &HLT_L1TOTEM1_MinBias_v1_Prescl, &b_HLT_L1TOTEM1_MinBias_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1TOTEM2_ZeroBias_v1", &HLT_L1TOTEM2_ZeroBias_v1, &b_HLT_L1TOTEM2_ZeroBias_v1);
   fChain->SetBranchAddress("HLT_L1TOTEM2_ZeroBias_v1_Prescl", &HLT_L1TOTEM2_ZeroBias_v1_Prescl, &b_HLT_L1TOTEM2_ZeroBias_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_v1", &HLT_L1MinimumBiasHF1OR_v1, &b_HLT_L1MinimumBiasHF1OR_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1OR_v1_Prescl", &HLT_L1MinimumBiasHF1OR_v1_Prescl, &b_HLT_L1MinimumBiasHF1OR_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2OR_v1", &HLT_L1MinimumBiasHF2OR_v1, &b_HLT_L1MinimumBiasHF2OR_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2OR_v1_Prescl", &HLT_L1MinimumBiasHF2OR_v1_Prescl, &b_HLT_L1MinimumBiasHF2OR_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1AND_v1", &HLT_L1MinimumBiasHF1AND_v1, &b_HLT_L1MinimumBiasHF1AND_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF1AND_v1_Prescl", &HLT_L1MinimumBiasHF1AND_v1_Prescl, &b_HLT_L1MinimumBiasHF1AND_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2AND_v1", &HLT_L1MinimumBiasHF2AND_v1, &b_HLT_L1MinimumBiasHF2AND_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2AND_v1_Prescl", &HLT_L1MinimumBiasHF2AND_v1_Prescl, &b_HLT_L1MinimumBiasHF2AND_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2ORNoBptxGating_v1", &HLT_L1MinimumBiasHF2ORNoBptxGating_v1, &b_HLT_L1MinimumBiasHF2ORNoBptxGating_v1);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl", &HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl, &b_HLT_L1MinimumBiasHF2ORNoBptxGating_v1_Prescl);
  
   fChain->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1", &HLT_AK4CaloJet40_Eta5p1_v1, &b_HLT_AK4CaloJet40_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1_Prescl", &HLT_AK4CaloJet40_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet40_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1", &HLT_AK4CaloJet60_Eta5p1_v1, &b_HLT_AK4CaloJet60_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1_Prescl", &HLT_AK4CaloJet60_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet60_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1", &HLT_AK4CaloJet80_Eta5p1_v1, &b_HLT_AK4CaloJet80_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1_Prescl", &HLT_AK4CaloJet80_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet80_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1", &HLT_AK4CaloJet100_Eta5p1_v1, &b_HLT_AK4CaloJet100_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1_Prescl", &HLT_AK4CaloJet100_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet100_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet110_Eta5p1_v1", &HLT_AK4CaloJet110_Eta5p1_v1, &b_HLT_AK4CaloJet110_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet110_Eta5p1_v1_Prescl", &HLT_AK4CaloJet110_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet110_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet120_Eta5p1_v1", &HLT_AK4CaloJet120_Eta5p1_v1, &b_HLT_AK4CaloJet120_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet120_Eta5p1_v1_Prescl", &HLT_AK4CaloJet120_Eta5p1_v1_Prescl, &b_HLT_AK4CaloJet120_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet150_v1", &HLT_AK4CaloJet150_v1, &b_HLT_AK4CaloJet150_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet150_v1_Prescl", &HLT_AK4CaloJet150_v1_Prescl, &b_HLT_AK4CaloJet150_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet40_Eta5p1_v1", &HLT_AK4PFJet40_Eta5p1_v1, &b_HLT_AK4PFJet40_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet40_Eta5p1_v1_Prescl", &HLT_AK4PFJet40_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet40_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet60_Eta5p1_v1", &HLT_AK4PFJet60_Eta5p1_v1, &b_HLT_AK4PFJet60_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet60_Eta5p1_v1_Prescl", &HLT_AK4PFJet60_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet60_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1", &HLT_AK4PFJet80_Eta5p1_v1, &b_HLT_AK4PFJet80_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1_Prescl", &HLT_AK4PFJet80_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet80_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet100_Eta5p1_v1", &HLT_AK4PFJet100_Eta5p1_v1, &b_HLT_AK4PFJet100_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet100_Eta5p1_v1_Prescl", &HLT_AK4PFJet100_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet100_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet110_Eta5p1_v1", &HLT_AK4PFJet110_Eta5p1_v1, &b_HLT_AK4PFJet110_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet110_Eta5p1_v1_Prescl", &HLT_AK4PFJet110_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet110_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFJet120_Eta5p1_v1", &HLT_AK4PFJet120_Eta5p1_v1, &b_HLT_AK4PFJet120_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFJet120_Eta5p1_v1_Prescl", &HLT_AK4PFJet120_Eta5p1_v1_Prescl, &b_HLT_AK4PFJet120_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta1p1_v1", &HLT_AK4CaloJet80_Jet35_Eta1p1_v1, &b_HLT_AK4CaloJet80_Jet35_Eta1p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl", &HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl, &b_HLT_AK4CaloJet80_Jet35_Eta1p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta0p7_v1", &HLT_AK4CaloJet80_Jet35_Eta0p7_v1, &b_HLT_AK4CaloJet80_Jet35_Eta0p7_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl", &HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl, &b_HLT_AK4CaloJet80_Jet35_Eta0p7_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta1p1_v1", &HLT_AK4CaloJet100_Jet35_Eta1p1_v1, &b_HLT_AK4CaloJet100_Jet35_Eta1p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl", &HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl, &b_HLT_AK4CaloJet100_Jet35_Eta1p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta0p7_v1", &HLT_AK4CaloJet100_Jet35_Eta0p7_v1, &b_HLT_AK4CaloJet100_Jet35_Eta0p7_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl", &HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl, &b_HLT_AK4CaloJet100_Jet35_Eta0p7_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_45_45_Eta2p1_v1", &HLT_AK4CaloJet80_45_45_Eta2p1_v1, &b_HLT_AK4CaloJet80_45_45_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl", &HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl, &b_HLT_AK4CaloJet80_45_45_Eta2p1_v1_Prescl);

   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1", &HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl", &HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet40Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1", &HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl", &HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet60Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1", &HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl", &HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet80Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1", &HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl", &HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_AK4CaloJet100Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton10Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton15Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton20Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton30Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1", &HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1, &b_HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl", &HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl, &b_HLT_HIL2Mu3Eta2p5_HIPhoton40Eta1p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_v1", &HLT_HIL1DoubleMu0_v1, &b_HLT_HIL1DoubleMu0_v1);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu0_v1_Prescl", &HLT_HIL1DoubleMu0_v1_Prescl, &b_HLT_HIL1DoubleMu0_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu10_v1", &HLT_HIL1DoubleMu10_v1, &b_HLT_HIL1DoubleMu10_v1);
   fChain->SetBranchAddress("HLT_HIL1DoubleMu10_v1_Prescl", &HLT_HIL1DoubleMu10_v1_Prescl, &b_HLT_HIL1DoubleMu10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v1", &HLT_HIL2DoubleMu0_NHitQ_v1, &b_HLT_HIL2DoubleMu0_NHitQ_v1);
   fChain->SetBranchAddress("HLT_HIL2DoubleMu0_NHitQ_v1_Prescl", &HLT_HIL2DoubleMu0_NHitQ_v1_Prescl, &b_HLT_HIL2DoubleMu0_NHitQ_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1", &HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1, &b_HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1);
   fChain->SetBranchAddress("HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl", &HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl, &b_HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3DoubleMu0_OS_m7to14_v1", &HLT_HIL3DoubleMu0_OS_m7to14_v1, &b_HLT_HIL3DoubleMu0_OS_m7to14_v1);
   fChain->SetBranchAddress("HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl", &HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl, &b_HLT_HIL3DoubleMu0_OS_m7to14_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ10_v1", &HLT_HIL2Mu3_NHitQ10_v1, &b_HLT_HIL2Mu3_NHitQ10_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu3_NHitQ10_v1_Prescl", &HLT_HIL2Mu3_NHitQ10_v1_Prescl, &b_HLT_HIL2Mu3_NHitQ10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu3_NHitQ15_v1", &HLT_HIL3Mu3_NHitQ15_v1, &b_HLT_HIL3Mu3_NHitQ15_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu3_NHitQ15_v1_Prescl", &HLT_HIL3Mu3_NHitQ15_v1_Prescl, &b_HLT_HIL3Mu3_NHitQ15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu5_NHitQ10_v1", &HLT_HIL2Mu5_NHitQ10_v1, &b_HLT_HIL2Mu5_NHitQ10_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu5_NHitQ10_v1_Prescl", &HLT_HIL2Mu5_NHitQ10_v1_Prescl, &b_HLT_HIL2Mu5_NHitQ10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu5_NHitQ15_v1", &HLT_HIL3Mu5_NHitQ15_v1, &b_HLT_HIL3Mu5_NHitQ15_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu5_NHitQ15_v1_Prescl", &HLT_HIL3Mu5_NHitQ15_v1_Prescl, &b_HLT_HIL3Mu5_NHitQ15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu7_NHitQ10_v1", &HLT_HIL2Mu7_NHitQ10_v1, &b_HLT_HIL2Mu7_NHitQ10_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu7_NHitQ10_v1_Prescl", &HLT_HIL2Mu7_NHitQ10_v1_Prescl, &b_HLT_HIL2Mu7_NHitQ10_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu7_NHitQ15_v1", &HLT_HIL3Mu7_NHitQ15_v1, &b_HLT_HIL3Mu7_NHitQ15_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu7_NHitQ15_v1_Prescl", &HLT_HIL3Mu7_NHitQ15_v1_Prescl, &b_HLT_HIL3Mu7_NHitQ15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v1", &HLT_HIL2Mu15_v1, &b_HLT_HIL2Mu15_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu15_v1_Prescl", &HLT_HIL2Mu15_v1_Prescl, &b_HLT_HIL2Mu15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu15_v1", &HLT_HIL3Mu15_v1, &b_HLT_HIL3Mu15_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu15_v1_Prescl", &HLT_HIL3Mu15_v1_Prescl, &b_HLT_HIL3Mu15_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL2Mu20_v1", &HLT_HIL2Mu20_v1, &b_HLT_HIL2Mu20_v1);
   fChain->SetBranchAddress("HLT_HIL2Mu20_v1_Prescl", &HLT_HIL2Mu20_v1_Prescl, &b_HLT_HIL2Mu20_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIL3Mu20_v1", &HLT_HIL3Mu20_v1, &b_HLT_HIL3Mu20_v1);
   fChain->SetBranchAddress("HLT_HIL3Mu20_v1_Prescl", &HLT_HIL3Mu20_v1_Prescl, &b_HLT_HIL3Mu20_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v2", &HLT_FullTrack18ForPPRef_v2, &b_HLT_FullTrack18ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v2_Prescl", &HLT_FullTrack18ForPPRef_v2_Prescl, &b_HLT_FullTrack18ForPPRef_v2_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v2", &HLT_FullTrack24ForPPRef_v2, &b_HLT_FullTrack24ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v2_Prescl", &HLT_FullTrack24ForPPRef_v2_Prescl, &b_HLT_FullTrack24ForPPRef_v2_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v3", &HLT_FullTrack34ForPPRef_v3, &b_HLT_FullTrack34ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v3_Prescl", &HLT_FullTrack34ForPPRef_v3_Prescl, &b_HLT_FullTrack34ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v2", &HLT_FullTrack45ForPPRef_v2, &b_HLT_FullTrack45ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v2_Prescl", &HLT_FullTrack45ForPPRef_v2_Prescl, &b_HLT_FullTrack45ForPPRef_v2_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v2", &HLT_FullTrack53ForPPRef_v2, &b_HLT_FullTrack53ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v2_Prescl", &HLT_FullTrack53ForPPRef_v2_Prescl, &b_HLT_FullTrack53ForPPRef_v2_Prescl);

   fChain->SetBranchAddress("HLT_HIL1CastorMediumJet_v1", &HLT_HIL1CastorMediumJet_v1, &b_HLT_HIL1CastorMediumJet_v1);
   fChain->SetBranchAddress("HLT_HIL1CastorMediumJet_v1_Prescl", &HLT_HIL1CastorMediumJet_v1_Prescl, &b_HLT_HIL1CastorMediumJet_v1_Prescl);
   fChain->SetBranchAddress("HLT_HICastorMediumJetPixel_SingleTrack_v1", &HLT_HICastorMediumJetPixel_SingleTrack_v1, &b_HLT_HICastorMediumJetPixel_SingleTrack_v1);
   fChain->SetBranchAddress("HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl", &HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl, &b_HLT_HICastorMediumJetPixel_SingleTrack_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt8_v1", &HLT_DmesonPPTrackingGlobal_Dpt8_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt8_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt8_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt15_v1", &HLT_DmesonPPTrackingGlobal_Dpt15_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt15_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt15_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt20_v1", &HLT_DmesonPPTrackingGlobal_Dpt20_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt20_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt20_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt30_v1", &HLT_DmesonPPTrackingGlobal_Dpt30_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt30_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt30_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt40_v1", &HLT_DmesonPPTrackingGlobal_Dpt40_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt40_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt40_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt50_v1", &HLT_DmesonPPTrackingGlobal_Dpt50_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt50_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt50_v1_Prescl);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt60_v1", &HLT_DmesonPPTrackingGlobal_Dpt60_v1, &b_HLT_DmesonPPTrackingGlobal_Dpt60_v1);
   fChain->SetBranchAddress("HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl", &HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl, &b_HLT_DmesonPPTrackingGlobal_Dpt60_v1_Prescl);

   fChain->SetBranchAddress("HLT_AK4PFBJetBCSV60_Eta2p1_v1", &HLT_AK4PFBJetBCSV60_Eta2p1_v1, &b_HLT_AK4PFBJetBCSV60_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl", &HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl, &b_HLT_AK4PFBJetBCSV60_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFBJetBCSV80_Eta2p1_v1", &HLT_AK4PFBJetBCSV80_Eta2p1_v1, &b_HLT_AK4PFBJetBCSV80_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl", &HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl, &b_HLT_AK4PFBJetBCSV80_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFDJet60_Eta2p1_v1", &HLT_AK4PFDJet60_Eta2p1_v1, &b_HLT_AK4PFDJet60_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFDJet60_Eta2p1_v1_Prescl", &HLT_AK4PFDJet60_Eta2p1_v1_Prescl, &b_HLT_AK4PFDJet60_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFDJet80_Eta2p1_v1", &HLT_AK4PFDJet80_Eta2p1_v1, &b_HLT_AK4PFDJet80_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFDJet80_Eta2p1_v1_Prescl", &HLT_AK4PFDJet80_Eta2p1_v1_Prescl, &b_HLT_AK4PFDJet80_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFBJetBSSV60_Eta2p1_v1", &HLT_AK4PFBJetBSSV60_Eta2p1_v1, &b_HLT_AK4PFBJetBSSV60_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl", &HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl, &b_HLT_AK4PFBJetBSSV60_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_AK4PFBJetBSSV80_Eta2p1_v1", &HLT_AK4PFBJetBSSV80_Eta2p1_v1, &b_HLT_AK4PFBJetBSSV80_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl", &HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl, &b_HLT_AK4PFBJetBSSV80_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("HLTriggerFinalPath_Prescl", &HLTriggerFinalPath_Prescl, &b_HLTriggerFinalPath_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v1", &HLT_PixelTracks_Multiplicity60_v1, &b_HLT_PixelTracks_Multiplicity60_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v1_Prescl", &HLT_PixelTracks_Multiplicity60_v1_Prescl, &b_HLT_PixelTracks_Multiplicity60_v1_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v2", &HLT_PixelTracks_Multiplicity60_v2, &b_HLT_PixelTracks_Multiplicity60_v2);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60_v2_Prescl", &HLT_PixelTracks_Multiplicity60_v2_Prescl, &b_HLT_PixelTracks_Multiplicity60_v2_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85_v1", &HLT_PixelTracks_Multiplicity85_v1, &b_HLT_PixelTracks_Multiplicity85_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85_v1_Prescl", &HLT_PixelTracks_Multiplicity85_v1_Prescl, &b_HLT_PixelTracks_Multiplicity85_v1_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110_v1", &HLT_PixelTracks_Multiplicity110_v1, &b_HLT_PixelTracks_Multiplicity110_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110_v1_Prescl", &HLT_PixelTracks_Multiplicity110_v1_Prescl, &b_HLT_PixelTracks_Multiplicity110_v1_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135_v1", &HLT_PixelTracks_Multiplicity135_v1, &b_HLT_PixelTracks_Multiplicity135_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135_v1_Prescl", &HLT_PixelTracks_Multiplicity135_v1_Prescl, &b_HLT_PixelTracks_Multiplicity135_v1_Prescl);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160_v1", &HLT_PixelTracks_Multiplicity160_v1, &b_HLT_PixelTracks_Multiplicity160_v1);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160_v1_Prescl", &HLT_PixelTracks_Multiplicity160_v1_Prescl, &b_HLT_PixelTracks_Multiplicity160_v1_Prescl);
   fChain->SetBranchAddress("HLT_Physics_v1", &HLT_Physics_v1, &b_HLT_Physics_v1);
   fChain->SetBranchAddress("HLT_Physics_v1_Prescl", &HLT_Physics_v1_Prescl, &b_HLT_Physics_v1_Prescl);
   fChain->SetBranchAddress("HLT_L1Tech5_BPTX_PlusOnly_v1", &HLT_L1Tech5_BPTX_PlusOnly_v1, &b_HLT_L1Tech5_BPTX_PlusOnly_v1);
   fChain->SetBranchAddress("HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl", &HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl, &b_HLT_L1Tech5_BPTX_PlusOnly_v1_Prescl);

   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v1", &HLT_FullTrack18ForPPRef_v1, &b_HLT_FullTrack18ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v1_Prescl", &HLT_FullTrack18ForPPRef_v1_Prescl, &b_HLT_FullTrack18ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v3", &HLT_FullTrack18ForPPRef_v3, &b_HLT_FullTrack18ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack18ForPPRef_v3_Prescl", &HLT_FullTrack18ForPPRef_v3_Prescl, &b_HLT_FullTrack18ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v1", &HLT_FullTrack24ForPPRef_v1, &b_HLT_FullTrack24ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v1_Prescl", &HLT_FullTrack24ForPPRef_v1_Prescl, &b_HLT_FullTrack24ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v3", &HLT_FullTrack24ForPPRef_v3, &b_HLT_FullTrack24ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack24ForPPRef_v3_Prescl", &HLT_FullTrack24ForPPRef_v3_Prescl, &b_HLT_FullTrack24ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v1", &HLT_FullTrack34ForPPRef_v1, &b_HLT_FullTrack34ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v1_Prescl", &HLT_FullTrack34ForPPRef_v1_Prescl, &b_HLT_FullTrack34ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v2", &HLT_FullTrack34ForPPRef_v2, &b_HLT_FullTrack34ForPPRef_v2);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v2_Prescl", &HLT_FullTrack34ForPPRef_v2_Prescl, &b_HLT_FullTrack34ForPPRef_v2_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v4", &HLT_FullTrack34ForPPRef_v4, &b_HLT_FullTrack34ForPPRef_v4);
   fChain->SetBranchAddress("HLT_FullTrack34ForPPRef_v4_Prescl", &HLT_FullTrack34ForPPRef_v4_Prescl, &b_HLT_FullTrack34ForPPRef_v4_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v1", &HLT_FullTrack45ForPPRef_v1, &b_HLT_FullTrack45ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v1_Prescl", &HLT_FullTrack45ForPPRef_v1_Prescl, &b_HLT_FullTrack45ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v3", &HLT_FullTrack45ForPPRef_v3, &b_HLT_FullTrack45ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack45ForPPRef_v3_Prescl", &HLT_FullTrack45ForPPRef_v3_Prescl, &b_HLT_FullTrack45ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v1", &HLT_FullTrack53ForPPRef_v1, &b_HLT_FullTrack53ForPPRef_v1);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v1_Prescl", &HLT_FullTrack53ForPPRef_v1_Prescl, &b_HLT_FullTrack53ForPPRef_v1_Prescl);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v3", &HLT_FullTrack53ForPPRef_v3, &b_HLT_FullTrack53ForPPRef_v3);
   fChain->SetBranchAddress("HLT_FullTrack53ForPPRef_v3_Prescl", &HLT_FullTrack53ForPPRef_v3_Prescl, &b_HLT_FullTrack53ForPPRef_v3_Prescl);
   fChain->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue, &b_L1_AlwaysTrue);
   fChain->SetBranchAddress("L1_AlwaysTrue_Prescl", &L1_AlwaysTrue_Prescl, &b_L1_AlwaysTrue_Prescl);
   fChain->SetBranchAddress("L1_CastorHighJet_BptxAND", &L1_CastorHighJet_BptxAND, &b_L1_CastorHighJet_BptxAND);
   fChain->SetBranchAddress("L1_CastorHighJet_BptxAND_Prescl", &L1_CastorHighJet_BptxAND_Prescl, &b_L1_CastorHighJet_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_CastorMediumJet_BptxAND", &L1_CastorMediumJet_BptxAND, &b_L1_CastorMediumJet_BptxAND);
   fChain->SetBranchAddress("L1_CastorMediumJet_BptxAND_Prescl", &L1_CastorMediumJet_BptxAND_Prescl, &b_L1_CastorMediumJet_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu0_BptxAND", &L1_DoubleMu0_BptxAND, &b_L1_DoubleMu0_BptxAND);
   fChain->SetBranchAddress("L1_DoubleMu0_BptxAND_Prescl", &L1_DoubleMu0_BptxAND_Prescl, &b_L1_DoubleMu0_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu0_MinimumBiasHF1_AND", &L1_DoubleMu0_MinimumBiasHF1_AND, &b_L1_DoubleMu0_MinimumBiasHF1_AND);
   fChain->SetBranchAddress("L1_DoubleMu0_MinimumBiasHF1_AND_Prescl", &L1_DoubleMu0_MinimumBiasHF1_AND_Prescl, &b_L1_DoubleMu0_MinimumBiasHF1_AND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMu3_BptxAND", &L1_DoubleMu3_BptxAND, &b_L1_DoubleMu3_BptxAND);
   fChain->SetBranchAddress("L1_DoubleMu3_BptxAND_Prescl", &L1_DoubleMu3_BptxAND_Prescl, &b_L1_DoubleMu3_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMuOpen_BptxAND", &L1_DoubleMuOpen_BptxAND, &b_L1_DoubleMuOpen_BptxAND);
   fChain->SetBranchAddress("L1_DoubleMuOpen_BptxAND_Prescl", &L1_DoubleMuOpen_BptxAND_Prescl, &b_L1_DoubleMuOpen_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleMuOpen_NotMinimumBiasHF2_AND", &L1_DoubleMuOpen_NotMinimumBiasHF2_AND, &b_L1_DoubleMuOpen_NotMinimumBiasHF2_AND);
   fChain->SetBranchAddress("L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl", &L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl, &b_L1_DoubleMuOpen_NotMinimumBiasHF2_AND_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet20", &L1_DoubleJet20, &b_L1_DoubleJet20);
   fChain->SetBranchAddress("L1_DoubleJet20_Prescl", &L1_DoubleJet20_Prescl, &b_L1_DoubleJet20_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet28_BptxAND", &L1_DoubleJet28_BptxAND, &b_L1_DoubleJet28_BptxAND);
   fChain->SetBranchAddress("L1_DoubleJet28_BptxAND_Prescl", &L1_DoubleJet28_BptxAND_Prescl, &b_L1_DoubleJet28_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_DoubleJet32", &L1_DoubleJet32, &b_L1_DoubleJet32);
   fChain->SetBranchAddress("L1_DoubleJet32_Prescl", &L1_DoubleJet32_Prescl, &b_L1_DoubleJet32_Prescl);
   fChain->SetBranchAddress("L1_ETT130", &L1_ETT130, &b_L1_ETT130);
   fChain->SetBranchAddress("L1_ETT130_Prescl", &L1_ETT130_Prescl, &b_L1_ETT130_Prescl);
   fChain->SetBranchAddress("L1_ETT15_BptxAND", &L1_ETT15_BptxAND, &b_L1_ETT15_BptxAND);
   fChain->SetBranchAddress("L1_ETT15_BptxAND_Prescl", &L1_ETT15_BptxAND_Prescl, &b_L1_ETT15_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT20_BptxAND", &L1_ETT20_BptxAND, &b_L1_ETT20_BptxAND);
   fChain->SetBranchAddress("L1_ETT20_BptxAND_Prescl", &L1_ETT20_BptxAND_Prescl, &b_L1_ETT20_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT30_BptxAND", &L1_ETT30_BptxAND, &b_L1_ETT30_BptxAND);
   fChain->SetBranchAddress("L1_ETT30_BptxAND_Prescl", &L1_ETT30_BptxAND_Prescl, &b_L1_ETT30_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT40", &L1_ETT40, &b_L1_ETT40);
   fChain->SetBranchAddress("L1_ETT40_Prescl", &L1_ETT40_Prescl, &b_L1_ETT40_Prescl);
   fChain->SetBranchAddress("L1_ETT40_BptxAND", &L1_ETT40_BptxAND, &b_L1_ETT40_BptxAND);
   fChain->SetBranchAddress("L1_ETT40_BptxAND_Prescl", &L1_ETT40_BptxAND_Prescl, &b_L1_ETT40_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT45_BptxAND", &L1_ETT45_BptxAND, &b_L1_ETT45_BptxAND);
   fChain->SetBranchAddress("L1_ETT45_BptxAND_Prescl", &L1_ETT45_BptxAND_Prescl, &b_L1_ETT45_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT50_BptxAND", &L1_ETT50_BptxAND, &b_L1_ETT50_BptxAND);
   fChain->SetBranchAddress("L1_ETT50_BptxAND_Prescl", &L1_ETT50_BptxAND_Prescl, &b_L1_ETT50_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT55_BptxAND", &L1_ETT55_BptxAND, &b_L1_ETT55_BptxAND);
   fChain->SetBranchAddress("L1_ETT55_BptxAND_Prescl", &L1_ETT55_BptxAND_Prescl, &b_L1_ETT55_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT60_BptxAND", &L1_ETT60_BptxAND, &b_L1_ETT60_BptxAND);
   fChain->SetBranchAddress("L1_ETT60_BptxAND_Prescl", &L1_ETT60_BptxAND_Prescl, &b_L1_ETT60_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT70_BptxAND", &L1_ETT70_BptxAND, &b_L1_ETT70_BptxAND);
   fChain->SetBranchAddress("L1_ETT70_BptxAND_Prescl", &L1_ETT70_BptxAND_Prescl, &b_L1_ETT70_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ETT90_BptxAND", &L1_ETT90_BptxAND, &b_L1_ETT90_BptxAND);
   fChain->SetBranchAddress("L1_ETT90_BptxAND_Prescl", &L1_ETT90_BptxAND_Prescl, &b_L1_ETT90_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_AND", &L1_MinimumBiasHF1_AND, &b_L1_MinimumBiasHF1_AND);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_AND_Prescl", &L1_MinimumBiasHF1_AND_Prescl, &b_L1_MinimumBiasHF1_AND_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_OR", &L1_MinimumBiasHF1_OR, &b_L1_MinimumBiasHF1_OR);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_OR_Prescl", &L1_MinimumBiasHF1_OR_Prescl, &b_L1_MinimumBiasHF1_OR_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_XOR", &L1_MinimumBiasHF1_XOR, &b_L1_MinimumBiasHF1_XOR);
   fChain->SetBranchAddress("L1_MinimumBiasHF1_XOR_Prescl", &L1_MinimumBiasHF1_XOR_Prescl, &b_L1_MinimumBiasHF1_XOR_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_AND", &L1_MinimumBiasHF2_AND, &b_L1_MinimumBiasHF2_AND);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_AND_Prescl", &L1_MinimumBiasHF2_AND_Prescl, &b_L1_MinimumBiasHF2_AND_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_OR", &L1_MinimumBiasHF2_OR, &b_L1_MinimumBiasHF2_OR);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_OR_Prescl", &L1_MinimumBiasHF2_OR_Prescl, &b_L1_MinimumBiasHF2_OR_Prescl);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_OR_NoBptxGating", &L1_MinimumBiasHF2_OR_NoBptxGating, &b_L1_MinimumBiasHF2_OR_NoBptxGating);
   fChain->SetBranchAddress("L1_MinimumBiasHF2_OR_NoBptxGating_Prescl", &L1_MinimumBiasHF2_OR_NoBptxGating_Prescl, &b_L1_MinimumBiasHF2_OR_NoBptxGating_Prescl);
   fChain->SetBranchAddress("L1_MuOpen_NotMinimumBiasHF2_AND", &L1_MuOpen_NotMinimumBiasHF2_AND, &b_L1_MuOpen_NotMinimumBiasHF2_AND);
   fChain->SetBranchAddress("L1_MuOpen_NotMinimumBiasHF2_AND_Prescl", &L1_MuOpen_NotMinimumBiasHF2_AND_Prescl, &b_L1_MuOpen_NotMinimumBiasHF2_AND_Prescl);
   fChain->SetBranchAddress("L1_NotMinimumBiasHF1_OR", &L1_NotMinimumBiasHF1_OR, &b_L1_NotMinimumBiasHF1_OR);
   fChain->SetBranchAddress("L1_NotMinimumBiasHF1_OR_Prescl", &L1_NotMinimumBiasHF1_OR_Prescl, &b_L1_NotMinimumBiasHF1_OR_Prescl);
   fChain->SetBranchAddress("L1_NotMinimumBiasHF2_AND", &L1_NotMinimumBiasHF2_AND, &b_L1_NotMinimumBiasHF2_AND);
   fChain->SetBranchAddress("L1_NotMinimumBiasHF2_AND_Prescl", &L1_NotMinimumBiasHF2_AND_Prescl, &b_L1_NotMinimumBiasHF2_AND_Prescl);
   fChain->SetBranchAddress("L1_NotZdcOR_BptxAND", &L1_NotZdcOR_BptxAND, &b_L1_NotZdcOR_BptxAND);
   fChain->SetBranchAddress("L1_NotZdcOR_BptxAND_Prescl", &L1_NotZdcOR_BptxAND_Prescl, &b_L1_NotZdcOR_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG12_BptxAND", &L1_SingleEG12_BptxAND, &b_L1_SingleEG12_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG12_BptxAND_Prescl", &L1_SingleEG12_BptxAND_Prescl, &b_L1_SingleEG12_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG15_BptxAND", &L1_SingleEG15_BptxAND, &b_L1_SingleEG15_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG15_BptxAND_Prescl", &L1_SingleEG15_BptxAND_Prescl, &b_L1_SingleEG15_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG20", &L1_SingleEG20, &b_L1_SingleEG20);
   fChain->SetBranchAddress("L1_SingleEG20_Prescl", &L1_SingleEG20_Prescl, &b_L1_SingleEG20_Prescl);
   fChain->SetBranchAddress("L1_SingleEG20_BptxAND", &L1_SingleEG20_BptxAND, &b_L1_SingleEG20_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG20_BptxAND_Prescl", &L1_SingleEG20_BptxAND_Prescl, &b_L1_SingleEG20_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG2_BptxAND", &L1_SingleEG2_BptxAND, &b_L1_SingleEG2_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG2_BptxAND_Prescl", &L1_SingleEG2_BptxAND_Prescl, &b_L1_SingleEG2_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG30_BptxAND", &L1_SingleEG30_BptxAND, &b_L1_SingleEG30_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG30_BptxAND_Prescl", &L1_SingleEG30_BptxAND_Prescl, &b_L1_SingleEG30_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG5", &L1_SingleEG5, &b_L1_SingleEG5);
   fChain->SetBranchAddress("L1_SingleEG5_Prescl", &L1_SingleEG5_Prescl, &b_L1_SingleEG5_Prescl);
   fChain->SetBranchAddress("L1_SingleEG5_BptxAND", &L1_SingleEG5_BptxAND, &b_L1_SingleEG5_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG5_BptxAND_Prescl", &L1_SingleEG5_BptxAND_Prescl, &b_L1_SingleEG5_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleEG7_BptxAND", &L1_SingleEG7_BptxAND, &b_L1_SingleEG7_BptxAND);
   fChain->SetBranchAddress("L1_SingleEG7_BptxAND_Prescl", &L1_SingleEG7_BptxAND_Prescl, &b_L1_SingleEG7_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu12_BptxAND", &L1_SingleMu12_BptxAND, &b_L1_SingleMu12_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu12_BptxAND_Prescl", &L1_SingleMu12_BptxAND_Prescl, &b_L1_SingleMu12_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu16_BptxAND", &L1_SingleMu16_BptxAND, &b_L1_SingleMu16_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu16_BptxAND_Prescl", &L1_SingleMu16_BptxAND_Prescl, &b_L1_SingleMu16_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu16_MinimumBiasHF1_AND", &L1_SingleMu16_MinimumBiasHF1_AND, &b_L1_SingleMu16_MinimumBiasHF1_AND);
   fChain->SetBranchAddress("L1_SingleMu16_MinimumBiasHF1_AND_Prescl", &L1_SingleMu16_MinimumBiasHF1_AND_Prescl, &b_L1_SingleMu16_MinimumBiasHF1_AND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3_BptxAND", &L1_SingleMu3_BptxAND, &b_L1_SingleMu3_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu3_BptxAND_Prescl", &L1_SingleMu3_BptxAND_Prescl, &b_L1_SingleMu3_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3_MinimumBiasHF1_AND", &L1_SingleMu3_MinimumBiasHF1_AND, &b_L1_SingleMu3_MinimumBiasHF1_AND);
   fChain->SetBranchAddress("L1_SingleMu3_MinimumBiasHF1_AND_Prescl", &L1_SingleMu3_MinimumBiasHF1_AND_Prescl, &b_L1_SingleMu3_MinimumBiasHF1_AND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu3p5", &L1_SingleMu3p5, &b_L1_SingleMu3p5);
   fChain->SetBranchAddress("L1_SingleMu3p5_Prescl", &L1_SingleMu3p5_Prescl, &b_L1_SingleMu3p5_Prescl);
   fChain->SetBranchAddress("L1_SingleMu5_BptxAND", &L1_SingleMu5_BptxAND, &b_L1_SingleMu5_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu5_BptxAND_Prescl", &L1_SingleMu5_BptxAND_Prescl, &b_L1_SingleMu5_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMu7_BptxAND", &L1_SingleMu7_BptxAND, &b_L1_SingleMu7_BptxAND);
   fChain->SetBranchAddress("L1_SingleMu7_BptxAND_Prescl", &L1_SingleMu7_BptxAND_Prescl, &b_L1_SingleMu7_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo", &L1_SingleMuBeamHalo, &b_L1_SingleMuBeamHalo);
   fChain->SetBranchAddress("L1_SingleMuBeamHalo_Prescl", &L1_SingleMuBeamHalo_Prescl, &b_L1_SingleMuBeamHalo_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   fChain->SetBranchAddress("L1_SingleMuOpen_Prescl", &L1_SingleMuOpen_Prescl, &b_L1_SingleMuOpen_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen_BptxAND", &L1_SingleMuOpen_BptxAND, &b_L1_SingleMuOpen_BptxAND);
   fChain->SetBranchAddress("L1_SingleMuOpen_BptxAND_Prescl", &L1_SingleMuOpen_BptxAND_Prescl, &b_L1_SingleMuOpen_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR, &b_L1_SingleMuOpen_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR_Prescl", &L1_SingleMuOpen_NotBptxOR_Prescl, &b_L1_SingleMuOpen_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_SingleJet12_BptxAND", &L1_SingleJet12_BptxAND, &b_L1_SingleJet12_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet12_BptxAND_Prescl", &L1_SingleJet12_BptxAND_Prescl, &b_L1_SingleJet12_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16", &L1_SingleJet16, &b_L1_SingleJet16);
   fChain->SetBranchAddress("L1_SingleJet16_Prescl", &L1_SingleJet16_Prescl, &b_L1_SingleJet16_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16_BptxAND", &L1_SingleJet16_BptxAND, &b_L1_SingleJet16_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet16_BptxAND_Prescl", &L1_SingleJet16_BptxAND_Prescl, &b_L1_SingleJet16_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet200", &L1_SingleJet200, &b_L1_SingleJet200);
   fChain->SetBranchAddress("L1_SingleJet200_Prescl", &L1_SingleJet200_Prescl, &b_L1_SingleJet200_Prescl);
   fChain->SetBranchAddress("L1_SingleJet20_BptxAND", &L1_SingleJet20_BptxAND, &b_L1_SingleJet20_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet20_BptxAND_Prescl", &L1_SingleJet20_BptxAND_Prescl, &b_L1_SingleJet20_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet24_BptxAND", &L1_SingleJet24_BptxAND, &b_L1_SingleJet24_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet24_BptxAND_Prescl", &L1_SingleJet24_BptxAND_Prescl, &b_L1_SingleJet24_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet28_BptxAND", &L1_SingleJet28_BptxAND, &b_L1_SingleJet28_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet28_BptxAND_Prescl", &L1_SingleJet28_BptxAND_Prescl, &b_L1_SingleJet28_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet32_BptxAND", &L1_SingleJet32_BptxAND, &b_L1_SingleJet32_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet32_BptxAND_Prescl", &L1_SingleJet32_BptxAND_Prescl, &b_L1_SingleJet32_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36", &L1_SingleJet36, &b_L1_SingleJet36);
   fChain->SetBranchAddress("L1_SingleJet36_Prescl", &L1_SingleJet36_Prescl, &b_L1_SingleJet36_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36_BptxAND", &L1_SingleJet36_BptxAND, &b_L1_SingleJet36_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl", &L1_SingleJet36_BptxAND_Prescl, &b_L1_SingleJet36_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet40_BptxAND", &L1_SingleJet40_BptxAND, &b_L1_SingleJet40_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet40_BptxAND_Prescl", &L1_SingleJet40_BptxAND_Prescl, &b_L1_SingleJet40_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet44_BptxAND", &L1_SingleJet44_BptxAND, &b_L1_SingleJet44_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet44_BptxAND_Prescl", &L1_SingleJet44_BptxAND_Prescl, &b_L1_SingleJet44_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet48_BptxAND", &L1_SingleJet48_BptxAND, &b_L1_SingleJet48_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet48_BptxAND_Prescl", &L1_SingleJet48_BptxAND_Prescl, &b_L1_SingleJet48_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet52_BptxAND", &L1_SingleJet52_BptxAND, &b_L1_SingleJet52_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl", &L1_SingleJet52_BptxAND_Prescl, &b_L1_SingleJet52_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet60_BptxAND", &L1_SingleJet60_BptxAND, &b_L1_SingleJet60_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet60_BptxAND_Prescl", &L1_SingleJet60_BptxAND_Prescl, &b_L1_SingleJet60_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet68", &L1_SingleJet68, &b_L1_SingleJet68);
   fChain->SetBranchAddress("L1_SingleJet68_Prescl", &L1_SingleJet68_Prescl, &b_L1_SingleJet68_Prescl);
   fChain->SetBranchAddress("L1_SingleJet68_BptxAND", &L1_SingleJet68_BptxAND, &b_L1_SingleJet68_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet68_BptxAND_Prescl", &L1_SingleJet68_BptxAND_Prescl, &b_L1_SingleJet68_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet8_BptxAND", &L1_SingleJet8_BptxAND, &b_L1_SingleJet8_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet8_BptxAND_Prescl", &L1_SingleJet8_BptxAND_Prescl, &b_L1_SingleJet8_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJetC20_NotBptxOR", &L1_SingleJetC20_NotBptxOR, &b_L1_SingleJetC20_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJetC20_NotBptxOR_Prescl", &L1_SingleJetC20_NotBptxOR_Prescl, &b_L1_SingleJetC20_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_SingleJetC32_NotBptxOR", &L1_SingleJetC32_NotBptxOR, &b_L1_SingleJetC32_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJetC32_NotBptxOR_Prescl", &L1_SingleJetC32_NotBptxOR_Prescl, &b_L1_SingleJetC32_NotBptxOR_Prescl);
   fChain->SetBranchAddress("L1_TOTEM_0", &L1_TOTEM_0, &b_L1_TOTEM_0);
   fChain->SetBranchAddress("L1_TOTEM_0_Prescl", &L1_TOTEM_0_Prescl, &b_L1_TOTEM_0_Prescl);
   fChain->SetBranchAddress("L1_TOTEM_1", &L1_TOTEM_1, &b_L1_TOTEM_1);
   fChain->SetBranchAddress("L1_TOTEM_1_Prescl", &L1_TOTEM_1_Prescl, &b_L1_TOTEM_1_Prescl);
   fChain->SetBranchAddress("L1_TOTEM_2", &L1_TOTEM_2, &b_L1_TOTEM_2);
   fChain->SetBranchAddress("L1_TOTEM_2_Prescl", &L1_TOTEM_2_Prescl, &b_L1_TOTEM_2_Prescl);
   fChain->SetBranchAddress("L1_TOTEM_3", &L1_TOTEM_3, &b_L1_TOTEM_3);
   fChain->SetBranchAddress("L1_TOTEM_3_Prescl", &L1_TOTEM_3_Prescl, &b_L1_TOTEM_3_Prescl);
   fChain->SetBranchAddress("L1_ZdcOR", &L1_ZdcOR, &b_L1_ZdcOR);
   fChain->SetBranchAddress("L1_ZdcOR_Prescl", &L1_ZdcOR_Prescl, &b_L1_ZdcOR_Prescl);
   fChain->SetBranchAddress("L1_ZdcOR_BptxAND", &L1_ZdcOR_BptxAND, &b_L1_ZdcOR_BptxAND);
   fChain->SetBranchAddress("L1_ZdcOR_BptxAND_Prescl", &L1_ZdcOR_BptxAND_Prescl, &b_L1_ZdcOR_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ZdcXOR", &L1_ZdcXOR, &b_L1_ZdcXOR);
   fChain->SetBranchAddress("L1_ZdcXOR_Prescl", &L1_ZdcXOR_Prescl, &b_L1_ZdcXOR_Prescl);
   fChain->SetBranchAddress("L1_ZdcXOR_BptxAND", &L1_ZdcXOR_BptxAND, &b_L1_ZdcXOR_BptxAND);
   fChain->SetBranchAddress("L1_ZdcXOR_BptxAND_Prescl", &L1_ZdcXOR_BptxAND_Prescl, &b_L1_ZdcXOR_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   fChain->SetBranchAddress("L1_ZeroBias_Prescl", &L1_ZeroBias_Prescl, &b_L1_ZeroBias_Prescl);
}


/// START PBPB TRIGGERS...

else{
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_v1", &HLT_HIPuAK4CaloJet40_Eta5p1_v1, &b_HLT_HIPuAK4CaloJet40_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_v1_Prescl", &HLT_HIPuAK4CaloJet40_Eta5p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet40_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_v1", &HLT_HIPuAK4CaloJet60_Eta5p1_v1, &b_HLT_HIPuAK4CaloJet60_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_v1_Prescl", &HLT_HIPuAK4CaloJet60_Eta5p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet60_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1", &HLT_HIPuAK4CaloJet80_Eta5p1_v1, &b_HLT_HIPuAK4CaloJet80_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1_Prescl", &HLT_HIPuAK4CaloJet80_Eta5p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet80_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1", &HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1, &b_HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1_Prescl", &HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1_Prescl, &b_HLT_HIPuAK4CaloJet80_Eta5p1ForZS_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1", &HLT_HIPuAK4CaloJet100_Eta5p1_v1, &b_HLT_HIPuAK4CaloJet100_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1_Prescl", &HLT_HIPuAK4CaloJet100_Eta5p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet100_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet110_Eta5p1_v1", &HLT_HIPuAK4CaloJet110_Eta5p1_v1, &b_HLT_HIPuAK4CaloJet110_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet110_Eta5p1_v1_Prescl", &HLT_HIPuAK4CaloJet110_Eta5p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet110_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet120_Eta5p1_v1", &HLT_HIPuAK4CaloJet120_Eta5p1_v1, &b_HLT_HIPuAK4CaloJet120_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet120_Eta5p1_v1_Prescl", &HLT_HIPuAK4CaloJet120_Eta5p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet120_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet150_Eta5p1_v1", &HLT_HIPuAK4CaloJet150_Eta5p1_v1, &b_HLT_HIPuAK4CaloJet150_Eta5p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet150_Eta5p1_v1_Prescl", &HLT_HIPuAK4CaloJet150_Eta5p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet150_Eta5p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1", &HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1, &b_HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1_Prescl", &HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1_Prescl, &b_HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1", &HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1, &b_HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1_Prescl", &HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1_Prescl, &b_HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1", &HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1, &b_HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1_Prescl", &HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1_Prescl, &b_HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1", &HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1, &b_HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1_Prescl", &HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1_Prescl, &b_HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1", &HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1, &b_HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1_Prescl", &HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1_Prescl, &b_HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1", &HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1, &b_HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1_Prescl", &HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1_Prescl, &b_HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1", &HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1, &b_HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1_Prescl", &HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1_Prescl, &b_HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1", &HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1, &b_HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1_Prescl", &HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1_Prescl, &b_HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1", &HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1, &b_HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1_Prescl", &HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet80_Jet35_Eta1p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1", &HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1, &b_HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1_Prescl", &HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1_Prescl, &b_HLT_HIPuAK4CaloJet80_Jet35_Eta0p7_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1", &HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1, &b_HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1_Prescl", &HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet100_Jet35_Eta1p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1", &HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1, &b_HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1_Prescl", &HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1_Prescl, &b_HLT_HIPuAK4CaloJet100_Jet35_Eta0p7_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1", &HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1, &b_HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1_Prescl", &HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1_Prescl, &b_HLT_HIPuAK4CaloJet80_45_45_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloDJet60_Eta2p1_v1", &HLT_HIPuAK4CaloDJet60_Eta2p1_v1, &b_HLT_HIPuAK4CaloDJet60_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloDJet60_Eta2p1_v1_Prescl", &HLT_HIPuAK4CaloDJet60_Eta2p1_v1_Prescl, &b_HLT_HIPuAK4CaloDJet60_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloDJet80_Eta2p1_v1", &HLT_HIPuAK4CaloDJet80_Eta2p1_v1, &b_HLT_HIPuAK4CaloDJet80_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloDJet80_Eta2p1_v1_Prescl", &HLT_HIPuAK4CaloDJet80_Eta2p1_v1_Prescl, &b_HLT_HIPuAK4CaloDJet80_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1", &HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1, &b_HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1_Prescl", &HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1_Prescl, &b_HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1", &HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1, &b_HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1_Prescl", &HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1_Prescl, &b_HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1", &HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1, &b_HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1_Prescl", &HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1_Prescl, &b_HLT_HIPuAK4CaloBJetSSV60_Eta2p1_v1_Prescl);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1", &HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1, &b_HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1);
   fChain->SetBranchAddress("HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1_Prescl", &HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1_Prescl, &b_HLT_HIPuAK4CaloBJetSSV80_Eta2p1_v1_Prescl);
}

   Notify();
}

Bool_t HLT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
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
