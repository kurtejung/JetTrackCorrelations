
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLine.h"

#include <iostream>

using namespace std;

double ReturnDivError(double num, double numErr, double den, double denErr){
	if(den==0 || num==0) return 0;
	double err2 = pow(num/den,2)*(pow(numErr/num,2)+pow(denErr/den,2));
	return sqrt(err2);
}

double ReturnSubError(double first, double firstErr, double sec, double secErr){
	double err2 = pow(first*firstErr,2)+pow(sec*secErr,2);
	return sqrt(err2);
}


void drawJetTrackProj(){

	bool doMixEvt = true; //If false, do the intermediate range subtraction instead
	
	const int nTrkPtBins = 4;
	const int jetPtBins = 1;
	const int centralityBins = 4;
	const int collType = 2;

	string histoTitle = "Data";
	string PbPbhistoTitle = "RecoJet_RecoTrack";
	string PbPbhistoToDivide = "GenJet_RecoTrack";

	const string trkCorr = "";
	const string trkCorr2 = "";

	//pp  bins are filled duplicate in every cent bin so just picking 1 is fine
	//const int centBins[centralityBins+1] = {0,10,30,50,100};
	const int centBins[centralityBins+1] = {100,50,30,10,0};
	const int jptBins[jetPtBins+1] = {100,300};
	const int trkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8};

	double nearSideLow = -0.785; //(-pi/4)
	double nearSideHigh = 0.785; //(+pi/4)

	if(!doMixEvt){
		nearSideLow = -1.57;
		nearSideHigh = 4.71;
	}

	double midRangeLow = 1.378; // 3pi/8
	double midRangeHigh = 1.768; // 5pi/8

	double awaySideLow = 2.356; //(3pi/4)
	double awaySideHigh = 3.93; //(5pi/4);
	
	TH2D *inclCorrPlots[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH2D *bjetCorrPlots[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclFG_nearPeak[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bFG_nearPeak[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *toMix_nearPeak[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclFG_awayPeak[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bFG_awayPeak[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclBG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bBG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *BGToMix[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH2D *mixPlots[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH2D *bmixPlots[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *mixProjEta[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bmixProjEta[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *mixProjEtaToDivide[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclPhiBG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bPhiBG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *ToMixPhiBG[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclPhiFG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bPhiFG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *ToMixPhiFG[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH2D *sidebandPlots[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH2D *bsidebandPlots[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH2D *inclCorrPlotsToDivide[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH2D *mixPlotsToDivide[centralityBins][jetPtBins][nTrkPtBins][collType];

	TFile *inIncl = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/Data2015_pp_fullMerge_withBjet.root");
	TFile *inB = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/Data2015_pp_fullMerge_withBjet.root");

	//PbPb Data
	//TFile *inInclPbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/Data2015_PbPb_fullMerge_withBjet.root");
	//TFile *inBPbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/Data2015_PbPb_fullMerge_withBjet.root");

	//PbPb MC
	TFile *inInclPbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_RecoRecoReduced.root");
	TFile *inBPbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_RecoRecoReduced.root");
	TFile *halliePbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_GenRecoReduced.root");
	//TFile *halliePbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/HydJet_RecoJet_GenTrack_Aug23.root");
	
	cout << "loaded files" << endl;

	TH1D *nEvents = (TH1D*)inIncl->Get(Form("%s_Nevents_dijets",histoTitle.c_str()));
	int InclEvents = nEvents->GetEntries();

	TH1D *bnEvents = (TH1D*)inB->Get(Form("%s_Nevents_dijets",histoTitle.c_str()));
	int bEvents = bnEvents->GetEntries();

	TH1D *nEventsPbPb = (TH1D*)inInclPbPb->Get(Form("%s_Nevents_dijets",PbPbhistoTitle.c_str()));
	int InclEvents_PbPb = nEvents->GetEntries();

	TH1D *bnEventsPbPb = (TH1D*)inBPbPb->Get(Form("%s_Nevents_dijets",PbPbhistoTitle.c_str()));
	int bEvents_PbPb = bnEvents->GetEntries();
	
	cout << "evt level content loaded" << endl;

	TH1D *nJets_PbPb[centralityBins];
	TH1D *nBjets_PbPb[centralityBins];
	TH1D *nJetsToDiv_PbPb[centralityBins];

	TH1D *nJets_pp = (TH1D*)inIncl->Get(Form("%s_all_jets_corrpTCent0_Cent10_Pt100_Pt300",histoTitle.c_str()))->Clone("nJets_pp");
	TH1D *nBjets_pp = (TH1D*)inB->Get(Form("%s_all_bjets_corrpTCent0_Cent10_Pt100_Pt300",histoTitle.c_str()))->Clone("nBjets_pp");

	for(int i=0; i<centralityBins; i++){
		
		cout << "at step1" << endl;

		cout << "getting " << Form("%s_all_jets_corrpTCent%d_Cent%d_Pt100_Pt300",PbPbhistoTitle.c_str(),centBins[i+1],centBins[i]) << endl;
		nJets_PbPb[i] = (TH1D*)inInclPbPb->Get(Form("%s_all_jets_corrpTCent%d_Cent%d_Pt100_Pt300",PbPbhistoTitle.c_str(),centBins[i+1],centBins[i]))->Clone(Form("nJets_PbPb_cent%d_%d",centBins[i+1],centBins[i]));
		cout << "getting " << Form("%s_all_bjets_corrpTCent%d_Cent%d_Pt100_Pt300",PbPbhistoTitle.c_str(),centBins[i+1],centBins[i]) << endl;
		nBjets_PbPb[i] = (TH1D*)inBPbPb->Get(Form("%s_all_bjets_corrpTCent%d_Cent%d_Pt100_Pt300",PbPbhistoTitle.c_str(),centBins[i+1],centBins[i]))->Clone(Form("nBjets_PbPb_cent%d_%d",centBins[i+1],centBins[i]));
		
		cout << "at step2" << endl;
		nJetsToDiv_PbPb[i] = (TH1D*)halliePbPb->Get(Form("%s_all_jets_corrpTCent%d_Cent%d_Pt100_Pt300",PbPbhistoToDivide.c_str(),centBins[i+1],centBins[i]))->Clone(Form("nJetsToDiv_PbPb_cent%d_%d",centBins[i+1],centBins[i]));


		for(int j=0; j<jetPtBins; j++){
			for(int k=0; k<nTrkPtBins; k++){
				cout << "start!" << endl;
				inclCorrPlots[i][j][k][0] = (TH2D*)inIncl->Get(Form("%s_hJetTrackSignalBackgroundCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",histoTitle.c_str(),centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_pp",centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				bjetCorrPlots[i][j][k][0] = (TH2D*)inB->Get(Form("%s_hbJetTrackSignalBackgroundCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",histoTitle.c_str(),centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("bJetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_pp",centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));

				cout << "loading " << Form("%s_hJetTrackSignalBackground%sCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",PbPbhistoTitle.c_str(),trkCorr.c_str(),centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]) << endl;
				inclCorrPlots[i][j][k][1] = (TH2D*)inInclPbPb->Get(Form("%s_hJetTrackSignalBackground%sCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",PbPbhistoTitle.c_str(),trkCorr.c_str(),centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				cout << "loading b" <<endl;
				bjetCorrPlots[i][j][k][1] = (TH2D*)inBPbPb->Get(Form("%s_hbJetTrackSignalBackground%sCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",PbPbhistoTitle.c_str(),trkCorr.c_str(),centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("bJetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				
				cout << "loading " << Form("%s_hJetTrackME%sCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",PbPbhistoTitle.c_str(),trkCorr.c_str(),centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]) << endl;
				mixPlots[i][j][k][1] = (TH2D*)inInclPbPb->Get(Form("%s_hJetTrackME%sCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",PbPbhistoTitle.c_str(),trkCorr.c_str(),centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorrMixed_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				cout << "loading b" << endl;
				bmixPlots[i][j][k][1] = (TH2D*)inBPbPb->Get(Form("%s_hbJetTrackME%sCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",PbPbhistoTitle.c_str(),trkCorr.c_str(),centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("bJetTrackCorrMixed_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));


				inclCorrPlotsToDivide[i][j][k][1] = (TH2D*)halliePbPb->Get(Form("%s_hJetTrackSignalBackground%sCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",PbPbhistoToDivide.c_str(),trkCorr2.c_str(),centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorrToDivide_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				mixPlotsToDivide[i][j][k][1] = (TH2D*)halliePbPb->Get(Form("%s_hJetTrackME%sCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",PbPbhistoToDivide.c_str(),trkCorr2.c_str(), centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorrMixedToDivide_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i+1],centBins[i],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));

				//if(i!=0) inclCorrPlotsToDivide[0][j][k][1]->Add(inclCorrPlotsToDivide[i][j][k][1]);
				//if(i!=0) mixPlotsToDivide[0][j][k][1]->Add(mixPlotsToDivide[i][j][k][1]);

				inclCorrPlots[i][j][k][0]->Rebin2D(2,2);
				bjetCorrPlots[i][j][k][0]->Rebin2D(2,2);

				inclCorrPlots[i][j][k][1]->Rebin2D(2,2);
				bjetCorrPlots[i][j][k][1]->Rebin2D(2,2);

				mixPlots[i][j][k][1]->Rebin2D(2,2);
				bmixPlots[i][j][k][1]->Rebin2D(2,2);

				inclCorrPlotsToDivide[i][j][k][1]->Rebin2D(2,2);
				mixPlotsToDivide[i][j][k][1]->Rebin2D(2,2);

				//inclCorrPlots[i][j][k][1]->Scale(2e12);
				//bjetCorrPlots[i][j][k][1]->Scale(2e12);

				//mixPlots[i][j][k][1]->Scale(2e12);
				//bmixPlots[i][j][k][1]->Scale(2e12);

				/*inclCorrPlotsToDivide[i][j][k][1]->Rebin2D(5,2);
				mixPlotsToDivide[i][j][k][1]->Rebin2D(5,2);*/

				//inclCorrPlotsToDivide[i][j][k][1]->Scale(5);
				//mixPlotsToDivide[i][j][k][1]->Scale(2e12);

				cout << "centrality " << centBins[i+1] << " " << centBins[i] << " trkPt " << trkPtBins[k] << " " << trkPtBins[k+1] << " entries kurt: "<< inclCorrPlots[i][j][k][1]->GetEntries() << " entries hallie: "<< inclCorrPlotsToDivide[i][j][k][1]->GetEntries() << endl;

				mixProjEta[i][j][k][1] = (TH1D*)mixPlots[i][j][k][1]->ProjectionX(Form("_pxMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,1));
				bmixProjEta[i][j][k][1] = (TH1D*)bmixPlots[i][j][k][1]->ProjectionX(Form("_bpxMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,1));
				mixProjEtaToDivide[i][j][k][1] = (TH1D*)mixPlotsToDivide[i][j][k][1]->ProjectionX(Form("_pxMixToDivide_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,1));

				//Project all the mixed histograms back to the profile in dEta
				for(int ixbin=1; ixbin<=mixPlots[i][j][k][1]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=mixPlots[i][j][k][1]->GetNbinsY(); iybin++){
						mixPlots[i][j][k][1]->SetBinContent(ixbin,iybin, mixProjEta[i][j][k][1]->GetBinContent(ixbin));
						mixPlots[i][j][k][1]->SetBinError(ixbin,iybin, mixProjEta[i][j][k][1]->GetBinError(ixbin)*sqrt(mixProjEta[i][j][k][1]->GetNbinsX()));

						bmixPlots[i][j][k][1]->SetBinContent(ixbin,iybin, bmixProjEta[i][j][k][1]->GetBinContent(ixbin));
						bmixPlots[i][j][k][1]->SetBinError(ixbin,iybin, bmixProjEta[i][j][k][1]->GetBinError(ixbin)*sqrt(bmixProjEta[i][j][k][1]->GetNbinsX()));
					}
				}
				for(int ixbin=1; ixbin<=mixPlotsToDivide[i][j][k][1]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=mixPlotsToDivide[i][j][k][1]->GetNbinsY(); iybin++){
						mixPlotsToDivide[i][j][k][1]->SetBinContent(ixbin,iybin, mixProjEtaToDivide[i][j][k][1]->GetBinContent(ixbin));
						mixPlotsToDivide[i][j][k][1]->SetBinError(ixbin,iybin, mixProjEtaToDivide[i][j][k][1]->GetBinError(ixbin)*sqrt(mixProjEtaToDivide[i][j][k][1]->GetNbinsX()));

					}
				}
				mixPlots[i][j][k][1]->Scale(1./mixProjEta[i][j][k][1]->GetBinContent(mixProjEta[i][j][k][1]->FindBin(0.)));
				bmixPlots[i][j][k][1]->Scale(1./bmixProjEta[i][j][k][1]->GetBinContent(bmixProjEta[i][j][k][1]->FindBin(0.)));
				mixPlotsToDivide[i][j][k][1]->Scale(1./mixProjEtaToDivide[i][j][k][1]->GetBinContent(mixProjEtaToDivide[i][j][k][1]->FindBin(0.)));

				//if(doMixEvt){
				//	inclCorrPlots[i][j][k][1]->Divide(mixPlots[i][j][k][1]);
				//	bjetCorrPlots[i][j][k][1]->Divide(bmixPlots[i][j][k][1]);
				//}
				
				for(int m=0; m<collType; m++){			

				    //now do the middle stuff as BG (or mix evts)
					int lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(midRangeLow);
					int highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(midRangeHigh);
					int lowBin2, highBin2;
					if(m) lowBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetYaxis()->FindBin(midRangeLow);
					if(m) highBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetYaxis()->FindBin(midRangeHigh);
					if(!doMixEvt || m==0){
						inclBG[i][j][k][m] = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionX(Form("_pxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
						bBG[i][j][k][m] = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionX(Form("_bpxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
						if(m){
							lowBin = inclCorrPlotsToDivide[i][j][k][m]->GetYaxis()->FindBin(midRangeLow);
							highBin = inclCorrPlotsToDivide[i][j][k][m]->GetYaxis()->FindBin(midRangeHigh);
							BGToMix[i][j][k][m] = (TH1D*)inclCorrPlotsToDivide[i][j][k][m]->ProjectionX(Form("_pxBGToMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin2, highBin2,"e");
						}
					}
					else{
						inclBG[i][j][k][m] = (TH1D*)mixPlots[i][j][k][m]->ProjectionX(Form("_pxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),0, -1,"e");
						bBG[i][j][k][m] = (TH1D*)bmixPlots[i][j][k][m]->ProjectionX(Form("_bpxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),0,-1,"e");
						BGToMix[i][j][k][m] = (TH1D*)mixPlotsToDivide[i][j][k][m]->ProjectionX(Form("_pxBGToMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),0, -1,"e");
					}


					if(!doMixEvt || m==0){
						inclBG[i][j][k][m]->Scale(1./inclBG[i][j][k][m]->GetBinContent(inclBG[i][j][k][m]->FindBin(0)));
						bBG[i][j][k][m]->Scale(1./bBG[i][j][k][m]->GetBinContent(bBG[i][j][k][m]->FindBin(0)));
						if(m) BGToMix[i][j][k][m]->Scale(1./BGToMix[i][j][k][m]->GetBinContent(BGToMix[i][j][k][m]->FindBin(0)));
					}

					//Divide out the acceptance via mixed events or by mid-region

					for(int ixbin=1; ixbin<=inclCorrPlots[i][j][k][m]->GetNbinsX(); ixbin++){
						for(int iybin=1; iybin<=inclCorrPlots[i][j][k][m]->GetNbinsY(); iybin++){
							if(inclBG[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin, inclCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin)/inclBG[i][j][k][m]->GetBinContent(ixbin));
							else inclCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin,0);
							if(inclBG[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlots[i][j][k][m]->SetBinError(ixbin,iybin, ReturnDivError(inclCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin), inclCorrPlots[i][j][k][m]->GetBinError(ixbin,iybin), inclBG[i][j][k][m]->GetBinContent(ixbin),inclBG[i][j][k][m]->GetBinError(ixbin)));

							if(bBG[i][j][k][m]->GetBinContent(ixbin)) bjetCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin, bjetCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin)/bBG[i][j][k][m]->GetBinContent(ixbin));
							else bjetCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin,0);
							if(bBG[i][j][k][m]->GetBinContent(ixbin)) bjetCorrPlots[i][j][k][m]->SetBinError(ixbin,iybin, ReturnDivError(bjetCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin), bjetCorrPlots[i][j][k][m]->GetBinError(ixbin,iybin), bBG[i][j][k][m]->GetBinContent(ixbin), bBG[i][j][k][m]->GetBinError(ixbin)));
						}
					}

					if(m==1){
						for(int ixbin=1; ixbin<=inclCorrPlotsToDivide[i][j][k][m]->GetNbinsX(); ixbin++){
							for(int iybin=1; iybin<=inclCorrPlotsToDivide[i][j][k][m]->GetNbinsY(); iybin++){
								if(BGToMix[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlotsToDivide[i][j][k][m]->SetBinContent(ixbin,iybin, inclCorrPlotsToDivide[i][j][k][m]->GetBinContent(ixbin,iybin)/BGToMix[i][j][k][m]->GetBinContent(ixbin));
								else inclCorrPlotsToDivide[i][j][k][m]->SetBinContent(ixbin,iybin,0);
								if(BGToMix[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlotsToDivide[i][j][k][m]->SetBinError(ixbin,iybin, ReturnDivError(inclCorrPlotsToDivide[i][j][k][m]->GetBinContent(ixbin,iybin), inclCorrPlotsToDivide[i][j][k][m]->GetBinError(ixbin,iybin), BGToMix[i][j][k][m]->GetBinContent(ixbin),BGToMix[i][j][k][m]->GetBinError(ixbin)));
							}
						}
					}

					if(m==0)inclCorrPlots[i][j][k][0]->Scale(1./nJets_pp->Integral());
					if(m==0)bjetCorrPlots[i][j][k][0]->Scale(1./nBjets_pp->Integral());

					if(m==1)inclCorrPlots[i][j][k][1]->Scale(1./nJets_PbPb[i]->Integral());
					if(m==1)bjetCorrPlots[i][j][k][1]->Scale(1./nBjets_PbPb[i]->Integral());
					if(m==1)inclCorrPlotsToDivide[i][j][k][1]->Scale(1./nJetsToDiv_PbPb[i]->Integral());
					if(m==1 && doMixEvt){
					
					    //Scale by the number of bins from eta (-3,3)
						//******* THIS PROCEDURE HERE MIGHT BE FUCKED AND IMPACTS YIELDS DRAMATICALLY!! ********** //

						inclCorrPlots[i][j][k][1]->Scale(inclCorrPlots[i][j][k][1]->GetNbinsX()*3./5.);
						bjetCorrPlots[i][j][k][1]->Scale(bjetCorrPlots[i][j][k][1]->GetNbinsX()*3./5.);
						inclCorrPlotsToDivide[i][j][k][1]->Scale(inclCorrPlotsToDivide[i][j][k][1]->GetNbinsX()*3./5.);
					}

					//Get the dphi projection to remove v2, etc
					lowBin = inclCorrPlots[i][j][k][m]->GetXaxis()->FindBin(-2.5);
					highBin = inclCorrPlots[i][j][k][m]->GetXaxis()->FindBin(-1.5);
					if(m) lowBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetXaxis()->FindBin(-2.5);
					if(m) highBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetXaxis()->FindBin(-1.5);
					int totalBins = highBin-lowBin;
					int totalBins2 = highBin2-lowBin2;
					inclPhiBG[i][j][k][m] = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionY(Form("_pyPhi_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					bPhiBG[i][j][k][m] = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionY(Form("_bpyPhi_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					if(m==1) ToMixPhiBG[i][j][k][m] = (TH1D*)inclCorrPlotsToDivide[i][j][k][m]->ProjectionY(Form("_pyPhiToMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin2, highBin2,"e");
					lowBin = inclCorrPlots[i][j][k][m]->GetXaxis()->FindBin(1.5);
					highBin = inclCorrPlots[i][j][k][m]->GetXaxis()->FindBin(2.5);
					if(m) lowBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetXaxis()->FindBin(1.5);
					if(m) highBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetXaxis()->FindBin(2.5);
					totalBins+=(highBin-lowBin);
					if(m) totalBins2+=(highBin2-lowBin2);
					TH1D *tmp1 = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionY(Form("_pyPhiTmp1_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					TH1D *tmp2 = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionY(Form("_bpyPhiTmp2_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					TH1D *tmp3;
					if(m==1) tmp3 = (TH1D*)inclCorrPlotsToDivide[i][j][k][m]->ProjectionY(Form("_pyPhiTmp3ToMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin2, highBin2,"e");
					inclPhiBG[i][j][k][m]->Add(tmp1);
					bPhiBG[i][j][k][m]->Add(tmp2);
					if(m==1) ToMixPhiBG[i][j][k][m]->Add(tmp3);
					
					inclPhiBG[i][j][k][m]->Scale(1./(double)(totalBins+2));
					bPhiBG[i][j][k][m]->Scale(1./(double)(totalBins+2));
					if(m==1) ToMixPhiBG[i][j][k][m]->Scale(1./(double)(totalBins2+2));

					lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(-5);
					highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(5);
					totalBins = highBin-lowBin;
					if(m) lowBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetXaxis()->FindBin(-5);
					if(m) highBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetXaxis()->FindBin(5);
					if(m) totalBins2 = highBin2-lowBin2;
					inclPhiFG[i][j][k][m] = inclCorrPlots[i][j][k][m]->ProjectionY(Form("_pyPhiFG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					bPhiFG[i][j][k][m] = bjetCorrPlots[i][j][k][m]->ProjectionY(Form("_bpyPhiFG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					if(m==1) ToMixPhiFG[i][j][k][m] = inclCorrPlotsToDivide[i][j][k][m]->ProjectionY(Form("_pyPhiFGToMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin2, highBin2,"e");
					inclPhiFG[i][j][k][m]->Scale(1./inclPhiFG[i][j][k][m]->GetNbinsX());
					bPhiFG[i][j][k][m]->Scale(1./bPhiFG[i][j][k][m]->GetNbinsX());
					if(m==1) ToMixPhiFG[i][j][k][m]->Scale(1./ToMixPhiFG[i][j][k][m]->GetNbinsX());

					// do background subtraction via dphi sidebands
					assert(inclCorrPlots[i][j][k][m]->GetNbinsX() == inclPhiBG[i][j][k][m]->GetNbinsX());
					if(m==1) assert(inclCorrPlotsToDivide[i][j][k][m]->GetNbinsX() == ToMixPhiBG[i][j][k][m]->GetNbinsX());


					for(int ixbin=1; ixbin<=inclCorrPlots[i][j][k][m]->GetNbinsX(); ixbin++){
						for(int iybin=1; iybin<=inclCorrPlots[i][j][k][m]->GetNbinsY(); iybin++){
							if(m==1 && inclCorrPlots[i][j][k][m]->GetXaxis()->GetBinCenter(ixbin) > -0.1 && inclCorrPlots[i][j][k][m]->GetXaxis()->GetBinCenter(ixbin) < 0.1 && inclCorrPlots[i][j][k][m]->GetYaxis()->GetBinCenter(iybin) > -0.1 && inclCorrPlots[i][j][k][m]->GetYaxis()->GetBinCenter(iybin) < 0.1){
								cout << "(x,y): (" << inclCorrPlots[i][j][k][m]->GetXaxis()->GetBinLowEdge(ixbin) << "-" << inclCorrPlots[i][j][k][m]->GetXaxis()->GetBinLowEdge(ixbin)+inclCorrPlots[i][j][k][m]->GetXaxis()->GetBinWidth(ixbin) << ", " << inclCorrPlots[i][j][k][m]->GetYaxis()->GetBinLowEdge(iybin) << "-" << inclCorrPlots[i][j][k][m]->GetYaxis()->GetBinLowEdge(iybin)+inclCorrPlots[i][j][k][m]->GetYaxis()->GetBinWidth(iybin) << ")" << endl;
								cout << "incl Point " << ixbin << " " << iybin << " content" << inclCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin) << " to be subtracted by " << inclPhiBG[i][j][k][m]->GetBinContent(iybin) << endl;
							}
							if(inclPhiBG[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin, inclCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin) - inclPhiBG[i][j][k][m]->GetBinContent(iybin));
							else inclCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin,0);
							if(inclPhiBG[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlots[i][j][k][m]->SetBinError(ixbin,iybin, sqrt(inclPhiBG[i][j][k][m]->GetNbinsX())*ReturnSubError(inclCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin), inclCorrPlots[i][j][k][m]->GetBinError(ixbin,iybin), inclPhiBG[i][j][k][m]->GetBinContent(iybin), inclPhiBG[i][j][k][m]->GetBinError(iybin)));

							if(bPhiBG[i][j][k][m]->GetBinContent(ixbin)) bjetCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin, bjetCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin) - bPhiBG[i][j][k][m]->GetBinContent(iybin));
							else bjetCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin,0);
							if(bPhiBG[i][j][k][m]->GetBinContent(ixbin)) bjetCorrPlots[i][j][k][m]->SetBinError(ixbin,iybin, sqrt(bPhiBG[i][j][k][m]->GetNbinsX())*ReturnSubError(bjetCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin), bjetCorrPlots[i][j][k][m]->GetBinError(ixbin,iybin), bPhiBG[i][j][k][m]->GetBinContent(iybin), bPhiBG[i][j][k][m]->GetBinError(iybin)));
						}
					}
					if(m==1){
						for(int ixbin=1; ixbin<=inclCorrPlotsToDivide[i][j][k][m]->GetNbinsX(); ixbin++){
							for(int iybin=1; iybin<=inclCorrPlotsToDivide[i][j][k][m]->GetNbinsY(); iybin++){
								if(inclCorrPlotsToDivide[i][j][k][m]->GetXaxis()->GetBinCenter(ixbin) > -0.1 && inclCorrPlotsToDivide[i][j][k][m]->GetXaxis()->GetBinCenter(ixbin) < 0.1 && inclCorrPlotsToDivide[i][j][k][m]->GetYaxis()->GetBinCenter(iybin) > -0.1 && inclCorrPlotsToDivide[i][j][k][m]->GetYaxis()->GetBinCenter(iybin) < 0.1){
									cout << "hallie Point " << ixbin << " " << iybin << " content" << inclCorrPlotsToDivide[i][j][k][m]->GetBinContent(ixbin,iybin) << " to be subtracted by " << ToMixPhiBG[i][j][k][m]->GetBinContent(iybin) << endl;
								}
								if(ToMixPhiBG[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlotsToDivide[i][j][k][m]->SetBinContent(ixbin,iybin, inclCorrPlotsToDivide[i][j][k][m]->GetBinContent(ixbin,iybin) - ToMixPhiBG[i][j][k][m]->GetBinContent(iybin));
								else inclCorrPlotsToDivide[i][j][k][m]->SetBinContent(ixbin,iybin,0);
								if(ToMixPhiBG[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlotsToDivide[i][j][k][m]->SetBinError(ixbin,iybin, sqrt(ToMixPhiBG[i][j][k][m]->GetNbinsX())*ReturnSubError(inclCorrPlotsToDivide[i][j][k][m]->GetBinContent(ixbin,iybin), inclCorrPlotsToDivide[i][j][k][m]->GetBinError(ixbin,iybin), ToMixPhiBG[i][j][k][m]->GetBinContent(iybin), ToMixPhiBG[i][j][k][m]->GetBinError(iybin)));
							}
						}
					}

					//Finally, project near-side peak from the fg/mix - bg correlations
					lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(nearSideLow);
					highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(nearSideHigh);

					if(m)lowBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetYaxis()->FindBin(nearSideLow);
					if(m) highBin2 = inclCorrPlotsToDivide[i][j][k][m]->GetYaxis()->FindBin(nearSideHigh);
					inclFG_nearPeak[i][j][k][m] = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionX(Form("_pxNear_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					bFG_nearPeak[i][j][k][m] = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionX(Form("_bpxNear_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					if(m==1) toMix_nearPeak[i][j][k][m] = (TH1D*)inclCorrPlotsToDivide[i][j][k][m]->ProjectionX(Form("_pxNearToMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin2, highBin2,"e");

					if(i==0 && j==0 && k==3 && m==1){
						for(int ibin=1; ibin<=inclFG_nearPeak[i][j][k][m]->GetNbinsX(); ibin++){
							cout << "inclJet content, bin" << ibin << " (center: " << inclFG_nearPeak[i][j][k][m]->GetBinCenter(ibin) << " ) " << inclFG_nearPeak[i][j][k][m]->GetBinContent(ibin) << endl;
							cout << "bJet content, bin" << ibin << " " << bFG_nearPeak[i][j][k][m]->GetBinContent(ibin) << endl;
						}
					} 

					//inclFG_nearPeak[i][j][k][m]->Rebin(2);
					//bFG_nearPeak[i][j][k][m]->Rebin(2);

					for(int ibin=0; ibin<inclFG_nearPeak[i][j][k][m]->GetNbinsX(); ibin++){
						inclFG_nearPeak[i][j][k][m]->SetBinContent(ibin, inclFG_nearPeak[i][j][k][m]->GetBinContent(ibin)/inclFG_nearPeak[i][j][k][m]->GetBinWidth(ibin));
						inclFG_nearPeak[i][j][k][m]->SetBinError(ibin, inclFG_nearPeak[i][j][k][m]->GetBinError(ibin)/inclFG_nearPeak[i][j][k][m]->GetBinWidth(ibin));
						bFG_nearPeak[i][j][k][m]->SetBinContent(ibin, bFG_nearPeak[i][j][k][m]->GetBinContent(ibin)/bFG_nearPeak[i][j][k][m]->GetBinWidth(ibin));
						bFG_nearPeak[i][j][k][m]->SetBinError(ibin, bFG_nearPeak[i][j][k][m]->GetBinError(ibin)/bFG_nearPeak[i][j][k][m]->GetBinWidth(ibin));
					}
					if(m==1){
						for(int ibin=0; ibin<toMix_nearPeak[i][j][k][m]->GetNbinsX(); ibin++){
							toMix_nearPeak[i][j][k][m]->SetBinContent(ibin, toMix_nearPeak[i][j][k][m]->GetBinContent(ibin)/toMix_nearPeak[i][j][k][m]->GetBinWidth(ibin));
							toMix_nearPeak[i][j][k][m]->SetBinError(ibin, toMix_nearPeak[i][j][k][m]->GetBinError(ibin)/toMix_nearPeak[i][j][k][m]->GetBinWidth(ibin));
						}
					}
				}
			}
		}
	}

	////***************************************************************************/////
	/////                       START PLOTTING STUFF NOW                            ////
	////***************************************************************************/////

	//Check the ratios of the dEta Yields
	TH1D *dEtaProjForRatio[centralityBins][nTrkPtBins];
	TH1D *dEtaProjForRatioDen[centralityBins][nTrkPtBins];
	TLatex *etaProjLabel[centralityBins][nTrkPtBins];
	TLegend *leg1 = new TLegend(0.5,0.5,0.9,0.9);
	TLine *line1 = new TLine(-3.5,0,3.5,0);
	line1->SetLineStyle(2);
	TCanvas *cRatios = new TCanvas("cRatios","",1200,1200);
	cRatios->Divide(centralityBins, nTrkPtBins);
	for(int i=0; i<centralityBins; i++){
		for(int j=0; j<nTrkPtBins; j++){
			cRatios->cd(j*centralityBins+i+1);
			dEtaProjForRatio[i][j] = (TH1D*)inclFG_nearPeak[i][0][j][1]->Clone(Form("dEtaProjForRatio_%d_%d",i,j));//->ProjectionX();
			dEtaProjForRatioDen[i][j] = (TH1D*)toMix_nearPeak[i][0][j][1]->Clone(Form("dEtaProjForRatioDen_%d_%d",i,j));//->ProjectionX();
			//dEtaProjForRatio[i][j]->Divide(dEtaProjForRatioDen[i][j]);
			dEtaProjForRatio[i][j]->GetXaxis()->SetLabelSize(0.06);
			dEtaProjForRatio[i][j]->GetYaxis()->SetLabelSize(0.06);
			dEtaProjForRatio[i][j]->Rebin(2);
			dEtaProjForRatio[i][j]->GetXaxis()->SetRangeUser(-3.5,3.5);
			dEtaProjForRatio[i][j]->GetXaxis()->SetTitleSize(0.06);
			dEtaProjForRatio[i][j]->GetXaxis()->SetTitleOffset(1.3);
			dEtaProjForRatio[i][j]->SetXTitle("Track-Jet #Delta#eta");
			//dEtaProjForRatio[i][j]->SetMaximum(3);
			dEtaProjForRatio[i][j]->Draw("");
			cout << "kurt peak bin err: "<< dEtaProjForRatio[i][j]->GetBinError(dEtaProjForRatio[i][j]->FindBin(0)) << endl;
			dEtaProjForRatioDen[i][j]->SetLineColor(2);
			dEtaProjForRatioDen[i][j]->Rebin(2);
			//dEtaProjForRatioDen[i][j]->Scale(4*1.3);
			dEtaProjForRatioDen[i][j]->Draw("same");
			//ToMixPhiBG[i][0][j][1]->Rebin(4);
			//ToMixPhiBG[i][0][j][1]->Draw("Same");
			line1->Draw("same");
			cout << "hallie peak bin err: "<< dEtaProjForRatioDen[i][j]->GetBinError(dEtaProjForRatioDen[i][j]->FindBin(0)) << endl;
			etaProjLabel[i][j] = new TLatex(-3.3, dEtaProjForRatio[i][j]->GetMaximum()*0.8, Form("Cent %d to %d, pT %d to %d",centBins[i],centBins[i+1],trkPtBins[j],trkPtBins[j+1]));
			etaProjLabel[i][j]->SetTextSize(0.08);
			etaProjLabel[i][j]->Draw();
			if(i==0 && j==0) leg1->AddEntry(dEtaProjForRatio[i][j],"Kurt PbPb Data","l");
			if(i==0 && j==0) leg1->AddEntry(dEtaProjForRatioDen[i][j],"Hallie PbPb Data","l");
			if(i==centralityBins-1 && j==nTrkPtBins-1) leg1->Draw("same");
		}
	}


	//check the dphi projections
	TCanvas *phi = new TCanvas("phi","",1200,600);
	phi->Divide(4,1);
	phi->cd(1);
	inclPhiFG[0][0][0][0]->Draw();
	inclPhiBG[0][0][0][0]->SetLineColor(2);
	inclPhiBG[0][0][0][0]->Draw("same");
	for(int i=0; i<3; i++){
		phi->cd(i+2);
		inclPhiFG[0][0][i][1]->Draw();
		inclPhiBG[0][0][i][1]->SetLineColor(2);
		inclPhiBG[0][0][i][1]->Draw("same");
	}

	//check the deta projections
	TCanvas *eta = new TCanvas("eta","",1200,600);
	eta->Divide(4,1);
	eta->cd(1);
	bFG_nearPeak[0][0][0][0]->Draw();
	int lowBin = inclFG_nearPeak[0][0][0][0]->FindBin(-1.5);
	int highBin = inclFG_nearPeak[0][0][0][0]->FindBin(1.5);
	//inclBG[0][0][0][0]->Scale((inclFG_nearPeak[0][0][0][0]->Integral(5,lowBin)+inclFG_nearPeak[0][0][0][0]->Integral(highBin,inclFG_nearPeak[0][0][0][0]->GetNbinsX()-5)) / (inclBG[0][0][0][0]->Integral(5,lowBin)+inclBG[0][0][0][0]->Integral(highBin,inclBG[0][0][0][0]->GetNbinsX()-5)));
	inclBG[0][0][0][0]->SetLineColor(2);
	bBG[0][0][0][0]->Draw("same");
	for(int i=0; i<3; i++){
		eta->cd(i+2);
		bFG_nearPeak[0][0][i][1]->Draw();
		int lowBin = inclFG_nearPeak[0][0][i][1]->FindBin(-3);
		int highBin = inclFG_nearPeak[0][0][i][1]->FindBin(-2);
		//inclBG[0][0][i][1]->Scale((inclFG_nearPeak[0][0][i][1]->Integral(5,lowBin)+inclFG_nearPeak[0][0][i][1]->Integral(highBin,inclFG_nearPeak[0][0][i][1]->GetNbinsX()-5)) / (inclBG[0][0][i][1]->Integral(5,lowBin)+inclBG[0][0][i][1]->Integral(highBin,inclBG[0][0][i][1]->GetNbinsX()-5)));
		inclBG[0][0][i][1]->SetLineColor(2);
		inclBG[0][0][i][1]->Scale(inclFG_nearPeak[0][0][i][1]->Integral(lowBin,highBin)/inclBG[0][0][i][1]->Integral(lowBin,highBin));
		bBG[0][0][i][1]->Draw("same");
	}
	/*for(int i=1; i<centralityBins; i++){
		for(int j=0; j<jetPtBins; j++){
			for(int k=0; k<nTrkPtBins; k++){
				inclFG_nearPeak[0][j][k][1]->Add(inclFG_nearPeak[i][j][k][1]);
				bFG_nearPeak[0][j][k][1]->Add(bFG_nearPeak[i][j][k][1]);
				inclBG[0][j][k][1]->Add(inclBG[i][j][k][1]);
				bBG[0][j][k][1]->Add(bBG[i][j][k][1]);
				inclFG_awayPeak[0][j][k][1]->Add(inclFG_awayPeak[i][j][k][1]);
				bFG_awayPeak[0][j][k][1]->Add(bFG_awayPeak[i][j][k][1]);
			}
		}
	}*/

	int tptBin = 0;
	double bFrac = 1./0.035/0.90; //  1/bfrac/taggingPurity

	TLatex *l1[3];

	TCanvas *cc = new TCanvas("cc","",1200,600);
	cc->Divide(5,1);
	int dummy=1;
	for(int tptBin=0; tptBin<nTrkPtBins; tptBin++){
		cc->cd(dummy);
		cc->GetPad(dummy++)->SetLogy();
		l1[dummy-1] = new TLatex(0,1.1,Form("%d < Track p_{T} < %d",trkPtBins[tptBin],trkPtBins[tptBin+1]));
	//inclFG_nearPeak[0][0][tptBin]->Scale(1./(double)nEvents);
	//inclBG[0][0][tptBin]->Scale(1./(double)nEvents);
		for(int collision=0; collision<2; collision++){
			for(int icent=0; icent<centralityBins; icent++){
				inclFG_nearPeak[icent][0][tptBin][collision]->SetXTitle("Jet #eta");
				inclFG_nearPeak[icent][0][tptBin][collision]->SetYTitle("Jet Signal (mix-evt corrected)");
				//if(!icent && collision) inclFG_nearPeak[icent][0][tptBin][collision]->Draw();

				inclBG[icent][0][tptBin][collision]->SetMarkerColor(6);
				inclBG[icent][0][tptBin][collision]->SetLineColor(6);

				inclBG[icent][0][tptBin][collision]->Scale(0.9*inclFG_nearPeak[icent][0][tptBin][collision]->Integral(20,30)/inclBG[icent][0][tptBin][collision]->Integral(20,30));
				//if(!icent && collision) inclFG_nearPeak[icent][0][tptBin][collision]->Divide(inclBG[icent][0][tptBin][collision]);
				if(!icent && collision) inclFG_nearPeak[icent][0][tptBin][collision]->Draw();

				//if(!icent && collision) inclBG[icent][0][tptBin][collision]->Draw("same");

				bFG_nearPeak[icent][0][tptBin][collision]->SetMarkerStyle(24);
				//bFG_nearPeak[icent][0][tptBin][collision]->SetMarkerColor(kRed);
				//bFG_nearPeak[icent][0][tptBin][collision]->SetLineColor(kRed);
	//bFG_nearPeak[0][0][tptBin][collision]->Scale(bFrac);
				//if(!icent && collision) bFG_nearPeak[icent][0][tptBin][collision]->Draw("same");

			}
		}
		l1[dummy-1]->Draw("Same");
	}


	///*************************************************////
	///     Check the mixing output
	///*************************************************////


	TCanvas *checkMix = new TCanvas("checkMix","",800,800);
	checkMix->Divide(2,1);
	checkMix->cd(1);
	//inclPhiFG[0][0][0][1]->Draw("");
	//inclPhiBG[0][0][0][1]->Draw("same");
	//mixPlots[3][0][0][1]->Draw("surf1");
	TH1D *tt_proj = inclCorrPlots[3][0][0][1]->ProjectionX();
	TH1D *tt_proj2 = inclCorrPlotsToDivide[3][0][0][1]->ProjectionX();
	tt_proj->Divide(tt_proj2);
	tt_proj->Draw();
	tt_proj2->SetLineColor(2);
	//tt_proj2->Draw("same");
	checkMix->cd(2);
	//inclCorrPlots[2][0][0][1]->Draw("colz");

	TH1D *tmix_proj = mixPlots[3][0][0][1]->ProjectionX();
	TH1D *tmix_proj2 = mixPlotsToDivide[3][0][0][1]->ProjectionX();
	tmix_proj->Divide(tmix_proj2);
	tmix_proj->Draw();
	//mixPlotsToDivide[3][0][0][1]->Divide(mixPlots[3][0][0][1]);
	//mixPlotsToDivide[3][0][0][1]->Draw("surf1");


	///*************************************************////
	///     Draw the delta-Eta (inclusive) plots
	///*************************************************////

	int fillColors[nTrkPtBins] = {kRed-3, kGreen-3, kBlue-3, kOrange-3};//, kViolet-5};
	TLegend *leg = new TLegend(0.1,0.1,0.9,0.9);
	TLatex *centBinLabels[3];

	THStack *pbpb_eta[centralityBins];
	THStack *pbpb_eta_ratio[centralityBins];
	TH1D *inclClone[centralityBins];
	TH1D *bClone[centralityBins];
	THStack *pp_eta = new THStack("pp_eta","");
	for(int itrkpt=nTrkPtBins-1; itrkpt>=0; itrkpt--){
		inclFG_nearPeak[0][0][itrkpt][0]->SetFillColor(fillColors[itrkpt]);
		//inclFG_nearPeak[0][0][itrkpt][0]->GetXaxis()->SetRangeUser(-1.6,1.6);
		pp_eta->Add(inclFG_nearPeak[0][0][itrkpt][0]);
		leg->AddEntry(inclFG_nearPeak[0][0][itrkpt][0],Form("%d < p_{T}^{assoc} < %d",trkPtBins[itrkpt],trkPtBins[itrkpt+1]));
	}

	for(int ibin=0; ibin<centralityBins; ibin++){
		pbpb_eta[ibin] = new THStack(Form("pbpb_eta_centBin%d",ibin),"");
		pbpb_eta_ratio[ibin] = new THStack(Form("pbpb_eta_ratio_centBin%d",ibin),"");
		for(int itrkpt=nTrkPtBins-1; itrkpt>=0; itrkpt--){
			inclFG_nearPeak[ibin][0][itrkpt][1]->SetLineColor(1);
			//inclFG_nearPeak[ibin][0][itrkpt][1]->GetXaxis()->SetRangeUser(-1.6,1.6);
			inclFG_nearPeak[ibin][0][itrkpt][1]->SetFillColor(fillColors[itrkpt]);
			pbpb_eta[ibin]->Add(inclFG_nearPeak[ibin][0][itrkpt][1]);

			string newname(inclFG_nearPeak[ibin][0][itrkpt][1]->GetName());
			newname.append("_clone");
			inclClone[ibin] = (TH1D*)(inclFG_nearPeak[ibin][0][itrkpt][1])->Clone(newname.c_str());
			inclClone[ibin]->Add(inclFG_nearPeak[ibin][0][itrkpt][0],-1);
			//inclClone[ibin]->Divide(bFG_nearPeak[ibin][0][itrkpt][1]);
			pbpb_eta_ratio[ibin]->Add(inclClone[ibin]);
		}
	}

	TCanvas *dEtaCorr = new TCanvas("dEtaCorr","",1200,800);
	dEtaCorr->Divide(centralityBins,2);
	dEtaCorr->cd(1);
	pp_eta->Draw("hist");
	pp_eta->SetMaximum(20);//pp_eta->GetMaximum()*1.2);
	pp_eta->SetMinimum(-1);
	pp_eta->GetXaxis()->SetLimits(-1.6,1.6);
	pp_eta->GetXaxis()->SetTitle("#Delta#eta");
	pp_eta->Draw("hist");
	TLatex *pplabel = new TLatex(-1.3,pp_eta->GetMaximum()*1.1,"pp Data, inclusive jets");
	pplabel->Draw("same");
	for(int ibin=0; ibin<centralityBins-1; ibin++){
		dEtaCorr->cd(ibin+2);
		pbpb_eta[ibin]->Draw("hist");
		pbpb_eta[ibin]->SetMaximum(20.);//pbpb_eta[ibin]->GetMaximum()*1.2);
		pbpb_eta[ibin]->SetMinimum(-1);
		pbpb_eta[ibin]->GetXaxis()->SetLimits(-1.6,1.6);
		pbpb_eta[ibin]->GetXaxis()->SetTitle("Track-Jet #Delta#eta");
		pbpb_eta[ibin]->Draw("hist");
		centBinLabels[ibin] = new TLatex(-1.3,pbpb_eta[ibin]->GetMaximum()*1.1,Form("PbPb, %d-%d%%",centBins[ibin],centBins[ibin+1]));
		centBinLabels[ibin]->Draw("same");

		dEtaCorr->cd(ibin+2+centralityBins);
		pbpb_eta_ratio[ibin]->Draw("hist");
		pbpb_eta_ratio[ibin]->GetXaxis()->SetLimits(-1.6,1.6);
		pbpb_eta_ratio[ibin]->SetMaximum(10);
		pbpb_eta_ratio[ibin]->SetMinimum(-4);
		pbpb_eta_ratio[ibin]->GetXaxis()->SetTitle("Track-Jet #Delta#eta");
		pbpb_eta_ratio[ibin]->Draw("hist");
	}
	dEtaCorr->cd(5);
	leg->Draw();

	///*************************************************////
	///     Draw the delta-eta (b-jet) plots
	///*************************************************////

	TLatex *bcentBinLabels[3];
	THStack *bjet_pbpb_eta[centralityBins];
	THStack *bjet_pbpb_eta_ratio[centralityBins];
	THStack *bjet_pp_eta = new THStack("bjet_pp_eta","");
	for(int itrkpt=nTrkPtBins-1; itrkpt>=0; itrkpt--){
		bFG_nearPeak[0][0][itrkpt][0]->SetFillColor(fillColors[itrkpt]);
		bFG_nearPeak[0][0][itrkpt][0]->GetXaxis()->SetRangeUser(-1.6,1.6);
		bjet_pp_eta->Add(bFG_nearPeak[0][0][itrkpt][0]);
	}

	for(int ibin=0; ibin<centralityBins; ibin++){
		bjet_pbpb_eta[ibin] = new THStack(Form("b_pbpb_eta_centBin%d",ibin),"");
		bjet_pbpb_eta_ratio[ibin] = new THStack(Form("b_pbpb_eta_ratio_centBin%d",ibin),"");
		for(int itrkpt=nTrkPtBins-1; itrkpt>=0; itrkpt--){
			bFG_nearPeak[ibin][0][itrkpt][1]->SetLineColor(1);
			bFG_nearPeak[ibin][0][itrkpt][1]->GetXaxis()->SetRangeUser(-1.6,1.6);
			bFG_nearPeak[ibin][0][itrkpt][1]->SetFillColor(fillColors[itrkpt]);
			bjet_pbpb_eta[ibin]->Add(bFG_nearPeak[ibin][0][itrkpt][1]);

			string newname(bFG_nearPeak[ibin][0][itrkpt][1]->GetName());
			newname.append("_bjet_clone");
			bClone[ibin] = (TH1D*)(bFG_nearPeak[ibin][0][itrkpt][1])->Clone(newname.c_str());
			bClone[ibin]->Add(bFG_nearPeak[ibin][0][itrkpt][0],-1);
			bjet_pbpb_eta_ratio[ibin]->Add(bClone[ibin]);
		}
	}

	TCanvas *bdEtaCorr = new TCanvas("bdEtaCorr","",1200,800);
	bdEtaCorr->Divide(centralityBins,2);
	bdEtaCorr->cd(1);
	bjet_pp_eta->Draw("hist");
	bjet_pp_eta->SetMaximum(20.);//pp_eta->GetMaximum()*1.2);
	bjet_pp_eta->SetMinimum(-1.);
	bjet_pp_eta->GetXaxis()->SetLimits(-1.6,1.6);
	bjet_pp_eta->Draw("hist");
	TLatex *bpplabel = new TLatex(-1.3,bjet_pp_eta->GetMaximum()*1.1,"pp Data, b-jets");
	bpplabel->Draw("same");
	for(int ibin=0; ibin<centralityBins-1; ibin++){
		bdEtaCorr->cd(ibin+2);
		bjet_pbpb_eta[ibin]->Draw("hist");
		bjet_pbpb_eta[ibin]->SetMaximum(20.);//bjet_pbpb_eta[ibin]->GetMaximum()*1.2);
		bjet_pbpb_eta[ibin]->SetMinimum(-1);
		bjet_pbpb_eta[ibin]->GetXaxis()->SetLimits(-1.6,1.6);
		bjet_pbpb_eta[ibin]->Draw("hist");
		bcentBinLabels[ibin] = new TLatex(-1.3,bjet_pbpb_eta[ibin]->GetMaximum()*1.1,Form("PbPb, %d-%d%%",centBins[ibin],centBins[ibin+1]));
		bcentBinLabels[ibin]->Draw("same");

		bdEtaCorr->cd(ibin+2+centralityBins);
		bjet_pbpb_eta_ratio[ibin]->Draw("hist");
		bjet_pbpb_eta_ratio[ibin]->GetXaxis()->SetLimits(-1.6,1.6);
		bjet_pbpb_eta_ratio[ibin]->SetMaximum(10);
		bjet_pbpb_eta_ratio[ibin]->SetMinimum(-4);
		bjet_pbpb_eta_ratio[ibin]->GetXaxis()->SetTitle("Track-Jet #Delta#eta");
		bjet_pbpb_eta_ratio[ibin]->Draw("hist");
	}


	bdEtaCorr->cd(5);
	leg->Draw();


}
