
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
	
	const int nTrkPtBins = 5;
	const int jetPtBins = 1;
	const int centralityBins = 4;
	const int collType = 2;

	//pp  bins are filled duplicate in every cent bin so just picking 1 is fine
	const int centBins[centralityBins+1] = {0,10,30,50,100};
	const int jptBins[jetPtBins+1] = {100,300};
	const int trkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};

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

	TH1D *inclFG_awayPeak[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bFG_awayPeak[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclBG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bBG[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH2D *mixPlots[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH2D *bmixPlots[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *mixProjEta[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bmixProjEta[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclPhiBG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bPhiBG[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclPhiFG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bPhiFG[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH2D *sidebandPlots[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH2D *bsidebandPlots[centralityBins][jetPtBins][nTrkPtBins][collType];

	TFile *inIncl = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/pp_Data_5TeV_InclusiveJet.root");
	TFile *inB = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/pp_Data_5TeV_csv0p9Filter.root");

	TFile *inInclPbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/Data2015_PbPb_fullMerge_withBjet.root");
	TFile *inBPbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/Data2015_PbPb_fullMerge_withBjet.root");

	TH1D *nEvents = (TH1D*)inIncl->Get("Data_Nevents_dijets");
	int InclEvents = nEvents->GetEntries();

	TH1D *bnEvents = (TH1D*)inB->Get("Data_Nevents_dijets");
	int bEvents = bnEvents->GetEntries();

	TH1D *nEventsPbPb = (TH1D*)inInclPbPb->Get("Data_Nevents_dijets");
	int InclEvents_PbPb = nEvents->GetEntries();

	TH1D *bnEventsPbPb = (TH1D*)inBPbPb->Get("Data_Nevents_dijets");
	int bEvents_PbPb = bnEvents->GetEntries();

	TH1D *nJets_PbPb[centralityBins];
	TH1D *nBjets_PbPb[centralityBins];

	TH1D *nJets_pp = (TH1D*)inIncl->Get("Data_all_jets_corrpTCent0_Cent10_Pt100_Pt300")->Clone("nJets_pp");
	TH1D *nBjets_pp = (TH1D*)inB->Get("Data_all_jets_corrpTCent0_Cent10_Pt100_Pt300")->Clone("nBjets_pp");

	for(int i=0; i<centralityBins; i++){

		nJets_PbPb[i] = (TH1D*)inInclPbPb->Get(Form("Data_all_jets_corrpTCent%d_Cent%d_Pt100_Pt300",centBins[i],centBins[i+1]))->Clone(Form("nJets_PbPb_cent%d_%d",centBins[i],centBins[i+1]));
		nBjets_PbPb[i] = (TH1D*)inBPbPb->Get(Form("Data_all_bjets_corrpTCent%d_Cent%d_Pt100_Pt300",centBins[i],centBins[i+1]))->Clone(Form("nBjets_PbPb_cent%d_%d",centBins[i],centBins[i+1]));

		for(int j=0; j<jetPtBins; j++){
			for(int k=0; k<nTrkPtBins; k++){
				inclCorrPlots[i][j][k][0] = (TH2D*)inIncl->Get(Form("Data_hJetTrackSignalBackgroundCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_pp",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				bjetCorrPlots[i][j][k][0] = (TH2D*)inB->Get(Form("Data_hJetTrackSignalBackgroundCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("bJetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_pp",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));

				inclCorrPlots[i][j][k][1] = (TH2D*)inInclPbPb->Get(Form("Data_hJetTrackSignalBackgroundCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				bjetCorrPlots[i][j][k][1] = (TH2D*)inBPbPb->Get(Form("Data_hbJetTrackSignalBackgroundCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("bJetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));

				mixPlots[i][j][k][1] = (TH2D*)inInclPbPb->Get(Form("Data_hJetTrackMECent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorrMixed_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				bmixPlots[i][j][k][1] = (TH2D*)inBPbPb->Get(Form("Data_hbJetTrackMECent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("bJetTrackCorrMixed_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));



				//inclCorrPlots[i][j][k][0]->Scale(1./nJets_pp->GetEntries());
				//bjetCorrPlots[i][j][k][0]->Scale(1./nBjets_pp->GetEntries());

				//inclCorrPlots[i][j][k][1]->Scale(1./nJets_PbPb[i]->GetEntries());
				//bjetCorrPlots[i][j][k][1]->Scale(1./nBjets_PbPb[i]->GetEntries());

				inclCorrPlots[i][j][k][0]->Rebin2D(2,2);
				bjetCorrPlots[i][j][k][0]->Rebin2D(2,2);

				inclCorrPlots[i][j][k][1]->Rebin2D(2,2);
				bjetCorrPlots[i][j][k][1]->Rebin2D(2,2);

				mixPlots[i][j][k][1]->Rebin2D(2,2);
				bmixPlots[i][j][k][1]->Rebin2D(2,2);

				mixProjEta[i][j][k][1] = (TH1D*)mixPlots[i][j][k][1]->ProjectionX(Form("_pxMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,1));
				bmixProjEta[i][j][k][1] = (TH1D*)bmixPlots[i][j][k][1]->ProjectionX(Form("_bpxMix_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,1));

				//Project all the mixed histograms back to the profile in dEta
				for(int ixbin=1; ixbin<=mixPlots[i][j][k][1]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=mixPlots[i][j][k][1]->GetNbinsY(); iybin++){
						mixPlots[i][j][k][1]->SetBinContent(ixbin,iybin, mixProjEta[i][j][k][1]->GetBinContent(ixbin));
						mixPlots[i][j][k][1]->SetBinError(ixbin,iybin, mixProjEta[i][j][k][1]->GetBinError(ixbin));

						bmixPlots[i][j][k][1]->SetBinContent(ixbin,iybin, bmixProjEta[i][j][k][1]->GetBinContent(ixbin));
						bmixPlots[i][j][k][1]->SetBinError(ixbin,iybin, bmixProjEta[i][j][k][1]->GetBinError(ixbin));
					}
				}
				mixPlots[i][j][k][1]->Scale(1./mixPlots[i][j][k][1]->GetBinContent(mixPlots[i][j][k][1]->FindBin(0.)));
				bmixPlots[i][j][k][1]->Scale(1./bmixPlots[i][j][k][1]->GetBinContent(bmixPlots[i][j][k][1]->FindBin(0.)));

				//if(doMixEvt){
				//	inclCorrPlots[i][j][k][1]->Divide(mixPlots[i][j][k][1]);
				//	bjetCorrPlots[i][j][k][1]->Divide(bmixPlots[i][j][k][1]);
				//}
				
				for(int m=0; m<collType; m++){			

				    //now do the middle stuff as BG (or mix evts)
					int lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(midRangeLow);
					int highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(midRangeHigh);
					if(!doMixEvt || m==0){
						inclBG[i][j][k][m] = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionX(Form("_pxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
						bBG[i][j][k][m] = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionX(Form("_bpxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					}
					else{
						inclBG[i][j][k][m] = (TH1D*)mixPlots[i][j][k][m]->ProjectionX(Form("_pxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),0, -1,"e");
						bBG[i][j][k][m] = (TH1D*)bmixPlots[i][j][k][m]->ProjectionX(Form("_bpxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),0,-1,"e");
					}


					if(!doMixEvt || m==0){
						inclBG[i][j][k][m]->Scale(1./inclBG[i][j][k][m]->GetBinContent(inclBG[i][j][k][m]->FindBin(0)));
						bBG[i][j][k][m]->Scale(1./bBG[i][j][k][m]->GetBinContent(bBG[i][j][k][m]->FindBin(0)));
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

					if(m==0)inclCorrPlots[i][j][k][0]->Scale(1./nJets_pp->GetEntries());
					if(m==0)bjetCorrPlots[i][j][k][0]->Scale(1./nBjets_pp->GetEntries());

					if(m==1)inclCorrPlots[i][j][k][1]->Scale(1./nJets_PbPb[i]->GetEntries());
					if(m==1)bjetCorrPlots[i][j][k][1]->Scale(1./nBjets_PbPb[i]->GetEntries());
					if(m==1 && doMixEvt){
						inclCorrPlots[i][j][k][1]->Scale(40.);
						bjetCorrPlots[i][j][k][1]->Scale(40.);
					}


					//Get the dphi projection to remove v2, etc
					lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(-3.5);
					highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(-2);
					int totalBins = highBin-lowBin;
					inclPhiBG[i][j][k][m] = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionY(Form("_pyPhi_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					bPhiBG[i][j][k][m] = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionY(Form("_bpyPhi_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(2);
					highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(3.5);
					totalBins+=(highBin-lowBin);
					TH1D *tmp1 = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionY(Form("_pyPhiTmp1_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					TH1D *tmp2 = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionY(Form("_bpyPhiTmp2_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					inclPhiBG[i][j][k][m]->Add(tmp1);
					bPhiBG[i][j][k][m]->Add(tmp2);
					
					inclPhiBG[i][j][k][m]->Scale(0.9999/(double)(totalBins+1));
					bPhiBG[i][j][k][m]->Scale(0.9999/(double)(totalBins+1));

					lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(-5);
					highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(5);
					totalBins = highBin-lowBin;
					inclPhiFG[i][j][k][m] = inclCorrPlots[i][j][k][m]->ProjectionY(Form("_pyPhiFG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					bPhiFG[i][j][k][m] = bjetCorrPlots[i][j][k][m]->ProjectionY(Form("_bpyPhiFG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					inclPhiFG[i][j][k][m]->Scale(1./totalBins);
					bPhiFG[i][j][k][m]->Scale(1./totalBins);

					// do background subtraction via dphi sidebands
					for(int ixbin=1; ixbin<=inclCorrPlots[i][j][k][m]->GetNbinsX(); ixbin++){
						for(int iybin=1; iybin<=inclCorrPlots[i][j][k][m]->GetNbinsY(); iybin++){
							if(inclPhiBG[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin, inclCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin) - inclPhiBG[i][j][k][m]->GetBinContent(iybin));
							else inclCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin,0);
							if(inclPhiBG[i][j][k][m]->GetBinContent(ixbin)) inclCorrPlots[i][j][k][m]->SetBinError(ixbin,iybin, ReturnSubError(inclCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin), inclCorrPlots[i][j][k][m]->GetBinError(ixbin,iybin), inclPhiBG[i][j][k][m]->GetBinContent(iybin), inclPhiBG[i][j][k][m]->GetBinError(iybin)));

							if(bPhiBG[i][j][k][m]->GetBinContent(ixbin)) bjetCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin, bjetCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin) - bPhiBG[i][j][k][m]->GetBinContent(iybin));
							else bjetCorrPlots[i][j][k][m]->SetBinContent(ixbin,iybin,0);
							if(bPhiBG[i][j][k][m]->GetBinContent(ixbin)) bjetCorrPlots[i][j][k][m]->SetBinError(ixbin,iybin, ReturnSubError(bjetCorrPlots[i][j][k][m]->GetBinContent(ixbin,iybin), bjetCorrPlots[i][j][k][m]->GetBinError(ixbin,iybin), bPhiBG[i][j][k][m]->GetBinContent(iybin), bPhiBG[i][j][k][m]->GetBinError(iybin)));
						}
					}

					//Finally, project near-side peak from the fg/mix - bg correlations
					lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(nearSideLow);
					highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(nearSideHigh);
					inclFG_nearPeak[i][j][k][m] = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionX(Form("_pxNear_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					bFG_nearPeak[i][j][k][m] = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionX(Form("_bpxNear_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");

					//inclFG_nearPeak[i][j][k][m]->Rebin(2);
					//bFG_nearPeak[i][j][k][m]->Rebin(2);

					for(int ibin=0; ibin<inclFG_nearPeak[i][j][k][m]->GetNbinsX(); ibin++){
						inclFG_nearPeak[i][j][k][m]->SetBinContent(ibin, inclFG_nearPeak[i][j][k][m]->GetBinContent(ibin)/inclFG_nearPeak[i][j][k][m]->GetBinWidth(ibin));
						inclFG_nearPeak[i][j][k][m]->SetBinError(ibin, inclFG_nearPeak[i][j][k][m]->GetBinError(ibin)/inclFG_nearPeak[i][j][k][m]->GetBinWidth(ibin));
						bFG_nearPeak[i][j][k][m]->SetBinContent(ibin, bFG_nearPeak[i][j][k][m]->GetBinContent(ibin)/bFG_nearPeak[i][j][k][m]->GetBinWidth(ibin));
						bFG_nearPeak[i][j][k][m]->SetBinError(ibin, bFG_nearPeak[i][j][k][m]->GetBinError(ibin)/bFG_nearPeak[i][j][k][m]->GetBinWidth(ibin));
					}
				}
			}
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
	inclCorrPlots[0][0][0][1]->Draw("surf1");
	checkMix->cd(2);
	//inclCorrPlots[2][0][0][1]->Draw("colz");
	bjetCorrPlots[0][0][0][1]->Draw("surf1");


	///*************************************************////
	///     Draw the delta-Eta (inclusive) plots
	///*************************************************////

	int fillColors[nTrkPtBins] = {kRed-3, kGreen-3, kBlue-3, kOrange-3, kViolet-5};
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
	pp_eta->SetMaximum(30);//pp_eta->GetMaximum()*1.2);
	pp_eta->SetMinimum(-1);
	pp_eta->GetXaxis()->SetLimits(-1.6,1.6);
	pp_eta->GetXaxis()->SetTitle("#Delta#eta");
	pp_eta->Draw("hist");
	TLatex *pplabel = new TLatex(-1.3,pp_eta->GetMaximum()*1.1,"pp Data, inclusive jets");
	pplabel->Draw("same");
	for(int ibin=0; ibin<centralityBins-1; ibin++){
		dEtaCorr->cd(ibin+2);
		pbpb_eta[ibin]->Draw("hist");
		pbpb_eta[ibin]->SetMaximum(30.);//pbpb_eta[ibin]->GetMaximum()*1.2);
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
	bjet_pp_eta->SetMaximum(30.);//pp_eta->GetMaximum()*1.2);
	bjet_pp_eta->SetMinimum(-1.);
	bjet_pp_eta->GetXaxis()->SetLimits(-1.6,1.6);
	bjet_pp_eta->Draw("hist");
	TLatex *bpplabel = new TLatex(-1.3,bjet_pp_eta->GetMaximum()*1.1,"pp Data, b-jets");
	bpplabel->Draw("same");
	for(int ibin=0; ibin<centralityBins-1; ibin++){
		bdEtaCorr->cd(ibin+2);
		bjet_pbpb_eta[ibin]->Draw("hist");
		bjet_pbpb_eta[ibin]->SetMaximum(30.);//bjet_pbpb_eta[ibin]->GetMaximum()*1.2);
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
