

#include "Rtypes.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH1D.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include <iostream>

using namespace std;

void drawJetShapes(){
	
	const bool shapeNormalize = false;
	
	const int nCentBins=2;
	const int trkPtBins=7;
	
	const string xTrkBins[trkPtBins+1] = {"TrkPt0p5", "TrkPt0p7", "TrkPt1", "TrkPt2", "TrkPt4", "TrkPt8", "TrkPt16", "TrkPt999" };
	const string xCentBins[nCentBins+1] = {"Cent0", "Cent30", "Cent100"};
	
	const double mean_pts[trkPtBins] = {0.55,0.844,1.35,2.35,5.07,9.72,41.};
	const double ptWidth[trkPtBins] = {0.2,0.3,1.,2.,4.,8.,100.};
	
	const int colors[trkPtBins] = {kBlack, kBlue-9, kYellow-9, kOrange+1, kGreen+3, kRed+1, kRed+3};
	const double newXbins[15] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};
	const double newEtaBins[18] = {-2,-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,2.0};
	
	TH1D *PbPb_data[nCentBins][trkPtBins];
	TH1D *PbPb_incl_data[nCentBins][trkPtBins];
	TH1D *PbPbHisto_offTrk[nCentBins][trkPtBins];
	TH1D *pp_data[trkPtBins];
	TH1D *pp_MC[trkPtBins];
	TH1D *pp_incl[trkPtBins];
	TH1D *pp_MC_incl[trkPtBins];
	
	TH1D *pp_eta[trkPtBins];
	TH1D *pp_eta_incl[trkPtBins];
	
	THStack *PbPbStack[nCentBins];
	THStack *PbPbInclStack[nCentBins];
	THStack *ppStack;
	THStack *ppMCStack;
	THStack *ppInclStack;
	THStack *ppInclMCStack;
	
	THStack *ppEtaStack;
	THStack *ppInclEtaStack;
	
	TH1D *ppSumHisto;
	TH1D *ppInclSumHisto;
	
	TH1D *ppEtaSumHisto;
	TH1D *ppInclEtaSumHisto;
	
	TH1D *ppMCSumHisto;
	TH1D *ppMCInclSumHisto;
	TH1D *PbPbMCSumHisto[nCentBins];
	TH1D *PbPbMCInclSumHisto[nCentBins];
	TH1D *PbPbSumHisto[nCentBins];
	TH1D *PbPbppRatios[nCentBins];
	TH1D *PbPbbInclRatios[nCentBins];
	TH1D *PbPbInclSumHisto[nCentBins];
	
	TH1D *PbPbSumHisto_offTrk[nCentBins];
	
	TH1D *ppBInclRatio;
	
	TH1D *ppBInclEtaRatio;
	
	TH1D *ppMC_qm17Shape;
	TH1D *ppMC_qm17Shape_data;
	TH1D *PbPb_qm17Shape_data;
	
	TFile *finQM17 = new TFile("ppMC_QM17test_bgSubtractedCorrs_dR_fullHistos.root");
	//TFile *finQM17 = new TFile("Data_bgSubtractedCorrs_dR_fullHistos.root");
	double xbin[11] = {0.5,0.7,1,2,3,4,8,12,16,20,999};
	ppMC_qm17Shape = (TH1D*)finQM17->Get("JetShape2_Yield_jet_BkgSub_pTweightedInclusive_pp_TrkPt1_TrkPt2")->Clone("ppMC_qm17Shape");
	for(int i=3; i<10; i++){
		cout << "getting " << Form("JetShape2_Yield_jet_BkgSub_pTweightedInclusive_pp_TrkPt%g_TrkPt%g",xbin[i],xbin[i+1]) << endl;
		ppMC_qm17Shape->Add((TH1D*)finQM17->Get(Form("JetShape2_Yield_jet_BkgSub_pTweightedInclusive_pp_TrkPt%g_TrkPt%g",xbin[i],xbin[i+1])));	
	}
	for(int ibin=0; ibin<ppMC_qm17Shape->GetNbinsX(); ibin++){
		ppMC_qm17Shape->SetBinContent(ibin, ppMC_qm17Shape->GetBinContent(ibin)/ppMC_qm17Shape->GetBinWidth(ibin));
		ppMC_qm17Shape->SetBinError(ibin, ppMC_qm17Shape->GetBinError(ibin)/ppMC_qm17Shape->GetBinWidth(ibin));
	}
	ppMC_qm17Shape = (TH1D*)ppMC_qm17Shape->Rebin(14,ppMC_qm17Shape->GetName(), newXbins);
	
	TFile *finQM17data = new TFile("/Users/kjung/JetTrack2016/results_plotting/Jet_Shapes_pTweighted.root");
	double xbin2[11] = {0.5,0.7,1,2,3,4,8,12,16,20,300};
	ppMC_qm17Shape_data = (TH1D*)finQM17data->Get("JetShape2_Yield_BkgSub_pTweightedInclusive_pp_TrkPt1_TrkPt2")->Clone("ppMC_qm17Shape_data");
	PbPb_qm17Shape_data = (TH1D*)finQM17data->Get("JetShape2_Yield_BkgSub_pTweightedInclusive_Cent10_Cent30_TrkPt1_TrkPt2");
	for(int i=3; i<10; i++){
		cout << "getting " << Form("JetShape2_Yield_BkgSub_pTweightedInclusive_pp_TrkPt%g_TrkPt%g",xbin2[i],xbin2[i+1]) << endl;
		ppMC_qm17Shape_data->Add((TH1D*)finQM17data->Get(Form("JetShape2_Yield_BkgSub_pTweightedInclusive_pp_TrkPt%g_TrkPt%g",xbin2[i],xbin2[i+1])));
		PbPb_qm17Shape_data->Add((TH1D*)finQM17data->Get(Form("JetShape2_Yield_BkgSub_pTweightedInclusive_Cent10_Cent30_TrkPt%g_TrkPt%g",xbin2[i],xbin2[i+1])));
	}
	for(int ibin=0; ibin<ppMC_qm17Shape_data->GetNbinsX(); ibin++){
		//ppMC_qm17Shape_data->SetBinContent(ibin, ppMC_qm17Shape_data->GetBinContent(ibin)/ppMC_qm17Shape_data->GetBinWidth(ibin));
		//ppMC_qm17Shape_data->SetBinError(ibin, ppMC_qm17Shape_data->GetBinError(ibin)/ppMC_qm17Shape_data->GetBinWidth(ibin));
	}
	ppMC_qm17Shape_data = (TH1D*)ppMC_qm17Shape_data->Rebin(14,ppMC_qm17Shape_data->GetName(), newXbins);
	PbPb_qm17Shape_data = (TH1D*)PbPb_qm17Shape_data->Rebin(14,PbPb_qm17Shape_data->GetName(), newXbins);
	
	TFile *finb = new TFile("BJetData_bgSubtractedCorrs_dR_fullHistos_looseTrkCuts_wideBinning.root");
	TFile *finb_badTrk = new TFile("BJetData_bgSubtractedCorrs_dR_fullHistos_looseTrkCuts_officialJTCs_wideBinning.root");
	TFile *fin = new TFile("InclJetData_bgSubtractedCorrs_dR_fullHistos_looseTrkCuts_officialJTCs.root");
	TFile *finHMC = new TFile("HallieppMC_bgSubtractedCorrs_dR_fullHistos.root");
	TFile *finMC = new TFile("ppMC_bjet_bgSubtractedCorrs_dR_fullHistos.root");
	TFile *finInclMC = new TFile("ppMC_PbPbMCMixing_bgSubtractedCorrs_dR_fullHistos.root");
	
	TFile *finPbPbMC = new TFile("PbPbMC_bjet_bgSubtractedCorrs_dR_fullHistos.root");
	TFile *finInclPbPbMC = new TFile("PbPbMC_PbPbMCMixing_bgSubtractedCorrs_dR_fullHistos.root");
		
	cout << "cp1" << endl;
	for(int i=0; i<nCentBins; i++){
		for(int j=1; j<trkPtBins; j++){
			string toGet = Form("JetShape2_Yield_bjet_BkgSub_pTweightedInclusive_%s_%s_%s_%s",xCentBins[i].c_str(),xCentBins[i+1].c_str(), xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
			cout << toGet << endl;
			PbPb_data[i][j] = (TH1D*)finb->Get(toGet.c_str())->Clone(toGet.c_str());
			PbPb_data[i][j] = (TH1D*)PbPb_data[i][j]->Rebin(14,PbPb_data[i][j]->GetName(), newXbins);
			///PbPb_data[i][j] = (TH1D*)PbPb_data[i][j]->Rebin(17, PbPb_data[i][j]->GetName(), newEtaBins);
			PbPb_data[i][j]->SetFillColor(colors[j]);
			PbPb_data[i][j]->SetLineColor(1);
			//PbPb_data[i][j]->GetXaxis()->SetRangeUser(-2,2);
			PbPb_data[i][j]->GetXaxis()->SetRangeUser(0,1);
			//PbPb_data[i][j]->Scale(1./ptWidth[j]);
			
			toGet = Form("JetShape2_Yield_bjet_BkgSub_pTweightedInclusive_%s_%s_%s_%s",xCentBins[i].c_str(),xCentBins[i+1].c_str(), xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
			PbPbHisto_offTrk[i][j] = (TH1D*)finb_badTrk->Get(toGet.c_str())->Clone(toGet.c_str());
			PbPbHisto_offTrk[i][j] = (TH1D*)PbPbHisto_offTrk[i][j]->Rebin(14,PbPbHisto_offTrk[i][j]->GetName(), newXbins);
			
			toGet = Form("JetShape2_Yield_jet_BkgSub_pTweightedInclusive_%s_%s_%s_%s",xCentBins[i].c_str(),xCentBins[i+1].c_str(), xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
			PbPb_incl_data[i][j] = (TH1D*)fin->Get(toGet.c_str())->Clone(toGet.c_str());
			PbPb_incl_data[i][j] = (TH1D*)PbPb_incl_data[i][j]->Rebin(14,PbPb_incl_data[i][j]->GetName(), newXbins);
			PbPb_incl_data[i][j]->SetFillColor(colors[j]);
			PbPb_incl_data[i][j]->SetLineColor(1);
			PbPb_incl_data[i][j]->GetXaxis()->SetRangeUser(0,1);
			
			if(j==1){
				PbPbMCInclSumHisto[i] = (TH1D*)finInclPbPbMC->Get(toGet.c_str())->Clone(toGet.c_str());
				toGet = Form("JetShape2_Yield_bjet_BkgSub_pTweightedInclusive_%s_%s_%s_%s",xCentBins[i].c_str(),xCentBins[i+1].c_str(), xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
				PbPbMCSumHisto[i] = (TH1D*)finPbPbMC->Get(toGet.c_str())->Clone(toGet.c_str());

			}
			else{
				PbPbMCInclSumHisto[i]->Add((TH1D*)finInclPbPbMC->Get(toGet.c_str()));
				toGet = Form("JetShape2_Yield_bjet_BkgSub_pTweightedInclusive_%s_%s_%s_%s",xCentBins[i].c_str(),xCentBins[i+1].c_str(), xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
				PbPbMCSumHisto[i]->Add((TH1D*)finPbPbMC->Get(toGet.c_str()));
			}
		}
	}
	cout << "cp2" << endl;
	for(int j=0; j<trkPtBins; j++){
		string toGet = Form("JetShape2_Yield_bjet_BkgSub_pTweightedInclusive_pp_%s_%s",xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
		cout << toGet << endl;
		pp_data[j] = (TH1D*)finb->Get(toGet.c_str())->Clone(toGet.c_str());
		pp_data[j] = (TH1D*)pp_data[j]->Rebin(14,pp_data[j]->GetName(), newXbins);
		//pp_data[j] = (TH1D*)pp_data[j]->Rebin(17, pp_data[j]->GetName(), newEtaBins);
		pp_data[j]->SetFillColor(colors[j]);
		//pp_data[j]->GetXaxis()->SetRangeUser(-2,2);
		pp_data[j]->GetXaxis()->SetRangeUser(0,1);
		pp_data[j]->SetLineColor(1);
		//pp_data[j]->Scale(1./ptWidth[j]);
		
		toGet = Form("JetEta_Yield_bjet_BkgSub_pTweightedInclusive_pp_%s_%s",xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
		pp_eta[j] = (TH1D*)finb->Get(toGet.c_str())->Clone(toGet.c_str());
		pp_eta[j] = (TH1D*)pp_eta[j]->Rebin(17, pp_eta[j]->GetName(), newEtaBins);
		pp_eta[j]->SetFillColor(colors[j]);
		pp_eta[j]->GetXaxis()->SetRangeUser(-2,2);
		pp_eta[j]->SetLineColor(1);
		//pp_eta[j]->Scale(1./ptWidth[j]);
	}
	
	cout << "cp3" << endl;
	for(int j=0; j<trkPtBins; j++){
		string toGet = Form("JetShape2_Yield_jet_BkgSub_pTweightedInclusive_pp_%s_%s",xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
		pp_incl[j] = (TH1D*)fin->Get(toGet.c_str())->Clone("pp1");
		pp_incl[j] = (TH1D*)pp_incl[j]->Rebin(14,pp_incl[j]->GetName(), newXbins);
		//pp_incl[j] = (TH1D*)pp_incl[j]->Rebin(17, pp_incl[j]->GetName(), newEtaBins);
		pp_incl[j]->SetFillColor(colors[j]);
		//pp_incl[j]->GetXaxis()->SetRangeUser(-2,2);
		pp_incl[j]->GetXaxis()->SetRangeUser(0,1);
		pp_incl[j]->SetLineColor(1);
		//pp_incl[j]->Scale(1./ptWidth[j]);
		
		cout << "bin " << j << " data int: "<< pp_incl[j]->Integral() << endl;
		
		toGet = Form("JetEta_Yield_jet_BkgSub_pTweightedInclusive_pp_%s_%s",xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
		pp_eta_incl[j] = (TH1D*)fin->Get(toGet.c_str())->Clone(toGet.c_str());
		pp_eta_incl[j] = (TH1D*)pp_eta_incl[j]->Rebin(17, pp_eta_incl[j]->GetName(), newEtaBins);
		pp_eta_incl[j]->SetFillColor(colors[j]);
		pp_eta_incl[j]->GetXaxis()->SetRangeUser(-2,2);
		pp_eta_incl[j]->SetLineColor(1);
		//pp_eta_incl[j]->Scale(1./ptWidth[j]);
	}
	
	cout << "cp4" << endl;
	for(int j=0; j<trkPtBins; j++){
		string toGet = Form("JetShape2_Yield_bjet_BkgSub_pTweightedInclusive_pp_%s_%s",xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
		string renameGet = toGet+"_bMC";
		pp_MC[j] = (TH1D*)finMC->Get(toGet.c_str())->Clone(renameGet.c_str());
		pp_MC[j] = (TH1D*)pp_MC[j]->Rebin(14,pp_MC[j]->GetName(), newXbins);
		pp_MC[j]->SetFillColor(colors[j]);
		pp_MC[j]->GetXaxis()->SetRangeUser(0,1);
		//pp_data[j]->Scale(1./ptWidth[j]);
		
		toGet = Form("JetShape2_Yield_jet_BkgSub_pTweightedInclusive_pp_%s_%s",xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
		renameGet = toGet+"_inclMC";
		pp_MC_incl[j] = (TH1D*)((TH1D*)finInclMC->Get(toGet.c_str()))->Clone(renameGet.c_str());
		cout << "bin " << j << "MC integral " << pp_MC_incl[j]->Integral() << " data int: "<< pp_incl[j]->Integral() << endl;
		
		pp_MC_incl[j] = (TH1D*)pp_MC_incl[j]->Rebin(14,pp_MC_incl[j]->GetName(), newXbins);
		pp_MC_incl[j]->SetFillColor(colors[j]);
		pp_MC_incl[j]->GetXaxis()->SetRangeUser(0,1);
		
		
	}

	cout << "cp5" << endl;
	//Stack and sum them...
	for(int i=0; i<nCentBins; i++){
		PbPbStack[i] = new THStack(Form("PbPbStack_Cent%d",i),"");
		PbPbInclStack[i] = new THStack(Form("PbPbInclStack_Cent%d",i),"");
		for(int j=2; j<trkPtBins; j++){
			for(int ibin=1; ibin<=PbPb_data[i][j]->GetNbinsX(); ibin++){
				if(PbPb_data[i][j]->GetBinContent(ibin)>0){
					PbPb_data[i][j]->SetBinContent(ibin, PbPb_data[i][j]->GetBinContent(ibin)/PbPb_data[i][j]->GetBinWidth(ibin));
					PbPb_data[i][j]->SetBinError(ibin, PbPb_data[i][j]->GetBinError(ibin)/PbPb_data[i][j]->GetBinWidth(ibin));
					PbPbHisto_offTrk[i][j]->SetBinContent(ibin, PbPbHisto_offTrk[i][j]->GetBinContent(ibin)/PbPbHisto_offTrk[i][j]->GetBinWidth(ibin));
					PbPbHisto_offTrk[i][j]->SetBinError(ibin, PbPbHisto_offTrk[i][j]->GetBinError(ibin)/PbPbHisto_offTrk[i][j]->GetBinWidth(ibin));
				}
				else{
					PbPb_data[i][j]->SetBinContent(ibin, 0);
					PbPb_data[i][j]->SetBinError(ibin, 0);
				}
				
				if(PbPb_incl_data[i][j]->GetBinContent(ibin)>0){
					PbPb_incl_data[i][j]->SetBinContent(ibin, PbPb_incl_data[i][j]->GetBinContent(ibin)/PbPb_incl_data[i][j]->GetBinWidth(ibin));
					PbPb_incl_data[i][j]->SetBinError(ibin, PbPb_incl_data[i][j]->GetBinError(ibin)/PbPb_incl_data[i][j]->GetBinWidth(ibin));
				}
				else{
					PbPb_incl_data[i][j]->SetBinContent(ibin,0);
					PbPb_incl_data[i][j]->SetBinError(ibin,0);
				}
			}
			PbPbStack[i]->Add(PbPb_data[i][j],"hist");
			PbPbInclStack[i]->Add(PbPb_incl_data[i][j],"hist");
		}
		for(int ibin=1; ibin<=PbPbMCSumHisto[i]->GetNbinsX(); ibin++){
			if(PbPbMCSumHisto[i]->GetBinContent(ibin)>0){
				PbPbMCSumHisto[i]->SetBinContent(ibin, PbPbMCSumHisto[i]->GetBinContent(ibin)/PbPbMCSumHisto[i]->GetBinWidth(ibin));
				PbPbMCSumHisto[i]->SetBinError(ibin, PbPbMCSumHisto[i]->GetBinError(ibin)/PbPbMCSumHisto[i]->GetBinWidth(ibin));
			}
			else{
				PbPbMCSumHisto[i]->SetBinContent(ibin, 0);
				PbPbMCSumHisto[i]->SetBinError(ibin, 0);
			}
			if(PbPbMCInclSumHisto[i]->GetBinContent(ibin)>0){
				PbPbMCInclSumHisto[i]->SetBinContent(ibin, PbPbMCInclSumHisto[i]->GetBinContent(ibin)/PbPbMCInclSumHisto[i]->GetBinWidth(ibin));
				PbPbMCInclSumHisto[i]->SetBinError(ibin, PbPbMCInclSumHisto[i]->GetBinError(ibin)/PbPbMCInclSumHisto[i]->GetBinWidth(ibin));
			}
			else{
				PbPbMCInclSumHisto[i]->SetBinContent(ibin, 0);
				PbPbMCInclSumHisto[i]->SetBinError(ibin, 0);
			}
		}
		PbPbSumHisto[i] = (TH1D*)PbPb_data[i][2]->Clone(Form("PbPbSumHisto_cent%d",i));
		PbPbInclSumHisto[i] = (TH1D*)PbPb_incl_data[i][2]->Clone(Form("PbPbSumInclHisto_cent%d",i));
		PbPbSumHisto_offTrk[i] = (TH1D*)PbPbHisto_offTrk[i][2]->Clone(Form("PbPbSumHisto_offTrk_cent%d",i));
		for(int j=3; j<trkPtBins; j++){
			PbPbSumHisto[i]->Add(PbPb_data[i][j]);
			PbPbInclSumHisto[i]->Add(PbPb_incl_data[i][j]);
			PbPbSumHisto_offTrk[i]->Add(PbPbHisto_offTrk[i][j]);
		}
		
	}
	
	cout << "cp6 "<< endl;
	ppStack = new THStack("ppStack","");
	ppInclStack = new THStack("ppInclStack","");
	ppEtaStack = new THStack("ppEtaStack","");
	ppInclEtaStack = new THStack("ppInclEtaStack","");
	
	for(int j=2; j<trkPtBins; j++){
		for(int ibin=1; ibin<=pp_data[j]->GetNbinsX(); ibin++){
			if(pp_data[j]->GetBinContent(ibin)>0){
				pp_data[j]->SetBinContent(ibin, pp_data[j]->GetBinContent(ibin)/pp_data[j]->GetBinWidth(ibin));
				pp_data[j]->SetBinError(ibin, pp_data[j]->GetBinError(ibin)/pp_data[j]->GetBinWidth(ibin));
			}
			else{
				pp_data[j]->SetBinContent(ibin, 0);
				pp_data[j]->SetBinError(ibin,0);
			}
			
			if(pp_MC[j]->GetBinContent(ibin)>0){
				pp_MC[j]->SetBinContent(ibin, pp_MC[j]->GetBinContent(ibin)/pp_MC[j]->GetBinWidth(ibin));
				pp_MC[j]->SetBinError(ibin, pp_MC[j]->GetBinError(ibin)/pp_MC[j]->GetBinWidth(ibin));
			}
			else{
				pp_MC[j]->SetBinContent(ibin, 0);
				pp_MC[j]->SetBinError(ibin, 0);
			}
			
			if(pp_incl[j]->GetBinContent(ibin)>0){
				pp_incl[j]->SetBinContent(ibin, pp_incl[j]->GetBinContent(ibin)/pp_incl[j]->GetBinWidth(ibin));
				pp_incl[j]->SetBinError(ibin, pp_incl[j]->GetBinError(ibin)/pp_incl[j]->GetBinWidth(ibin));	
			}
			else{
				pp_incl[j]->SetBinContent(ibin,0);
				pp_incl[j]->SetBinError(ibin,0);
			}
			
			if(pp_MC_incl[j]->GetBinContent(ibin)>0){
				pp_MC_incl[j]->SetBinContent(ibin, pp_MC_incl[j]->GetBinContent(ibin)/pp_MC_incl[j]->GetBinWidth(ibin));
				pp_MC_incl[j]->SetBinError(ibin, pp_MC_incl[j]->GetBinError(ibin)/pp_MC_incl[j]->GetBinWidth(ibin));
			}
			else{
				pp_MC_incl[j]->SetBinContent(ibin, 0);
				pp_MC_incl[j]->SetBinError(ibin, 0);
			}
		}
		
		for(int ibin=1; ibin<=pp_eta[j]->GetNbinsX(); ibin++){
			pp_eta[j]->SetBinContent(ibin, pp_eta[j]->GetBinContent(ibin)/pp_eta[j]->GetBinWidth(ibin));
			pp_eta[j]->SetBinError(ibin, pp_eta[j]->GetBinError(ibin)/pp_eta[j]->GetBinWidth(ibin));
			
			pp_eta_incl[j]->SetBinContent(ibin, pp_eta_incl[j]->GetBinContent(ibin)/pp_eta_incl[j]->GetBinWidth(ibin));
			pp_eta_incl[j]->SetBinError(ibin, pp_eta_incl[j]->GetBinError(ibin)/pp_eta_incl[j]->GetBinWidth(ibin));
		}
	}
	cout << "cp7" << endl;
	
	ppSumHisto = (TH1D*)pp_data[2]->Clone("ppSumHisto");
	ppMCSumHisto = (TH1D*)pp_MC[2]->Clone("ppMCSumHisto");
	ppMCInclSumHisto = (TH1D*)pp_MC_incl[2]->Clone("ppMCInclSumHisto");
	ppInclSumHisto = (TH1D*)pp_incl[2]->Clone("ppInclSumHisto");
	ppEtaSumHisto = (TH1D*)pp_eta[2]->Clone("ppEtaSumHisto");
	ppInclEtaSumHisto = (TH1D*)pp_eta_incl[2]->Clone("ppInclEtaSumHisto");
	for(int j=3; j<trkPtBins; j++){
		ppSumHisto->Add(pp_data[j]);
		ppMCSumHisto->Add(pp_MC[j]);
		ppMCInclSumHisto->Add(pp_MC_incl[j]);
		ppInclSumHisto->Add(pp_incl[j]);
		
		ppEtaSumHisto->Add(pp_eta[j]);
		ppInclEtaSumHisto->Add(pp_eta_incl[j]);
	}
	
	for(int j=2; j<trkPtBins; j++){
		if(shapeNormalize) pp_data[j]->Scale(1./ppSumHisto->Integral("width"));
		ppStack->Add(pp_data[j],"hist");
		if(shapeNormalize) pp_incl[j]->Scale(1./ppInclSumHisto->Integral("width"));
		ppInclStack->Add(pp_incl[j],"hist");
		
		ppEtaStack->Add(pp_eta[j],"hist");
		ppInclEtaStack->Add(pp_eta_incl[j],"hist");
	}
	if(shapeNormalize){
		ppSumHisto->Scale(1./ppSumHisto->Integral("width"));
		ppInclSumHisto->Scale(1./ppInclSumHisto->Integral("width"));
		ppMCSumHisto->Scale(1./ppMCSumHisto->Integral("width"));
		ppMCInclSumHisto->Scale(1./ppMCInclSumHisto->Integral("width"));
	}
	
	cout << "cp8" << endl;
	TH1D *PbPbMCbInclRatios[nCentBins];
	for(int i=0; i<nCentBins; i++){
		PbPbppRatios[i] = (TH1D*)PbPbSumHisto[i]->Clone(Form("PbPbppRatio_cent%d",i));
		PbPbbInclRatios[i] = (TH1D*)PbPbSumHisto[i]->Clone(Form("PbPbbInclRatios_cent%d",i));
		cout << " cent bin " << i << " " << PbPbppRatios[i]->GetBinContent(1) << endl;
		PbPbppRatios[i]->Divide(ppSumHisto);
		PbPbSumHisto_offTrk[i]->Divide(ppSumHisto);
		//PbPbppRatios[i]->Divide(ppMC_qm17Shape_data);
		PbPbbInclRatios[i]->Divide(PbPbInclSumHisto[i]);
		
		PbPbMCbInclRatios[i] = (TH1D*)PbPbMCSumHisto[i]->Clone(Form("PbPbMCbInclRatio_cent%d",i));
		PbPbMCbInclRatios[i]->Divide(PbPbMCInclSumHisto[i]);
	}
	
	ppBInclRatio = (TH1D*)ppSumHisto->Clone("ppBInclRatio");
	ppBInclRatio->Divide(ppInclSumHisto);
	TH1D *ppMCbInclRatio = (TH1D*)ppMCSumHisto->Clone("ppMCbInclRatio");
	ppMCbInclRatio->Divide(ppMCInclSumHisto);
	ppMCbInclRatio->SetLineColor(kBlue+3);
	ppMCbInclRatio->SetLineStyle(2);
	ppMCbInclRatio->SetFillColor(0);
	
	//ppBInclEtaRatio = (TH1D*)ppEtaSumHisto->Clone("ppBInclEtaRatio");
	//ppBInclEtaRatio->Divide(ppInclEtaSumHisto);
	
	cout << "cp9" << endl;
	
	TLatex *titles[4];
	for(int i=0; i<nCentBins; i++){
		if(i==0) titles[i] = new TLatex(0.2,5e2,"PbPb, 0-30% Centrality");
		if(i==1) titles[i] = new TLatex(0.2,5e2,"PbPb, 30-100% Centrality");
	}
	TLegend *leg1 = new TLegend(0.1,0.1,0.9,0.9);
	leg1->AddEntry(pp_data[2],"1 < p_{T} < 2 GeV","f");
	leg1->AddEntry(pp_data[3],"2 < p_{T} < 4 GeV","f");
	leg1->AddEntry(pp_data[4],"4 < p_{T} < 8 GeV","f");
	leg1->AddEntry(pp_data[5],"8 < p_{T} < 16 GeV","f");
	leg1->AddEntry(pp_data[6],"p_{T} > 16 GeV","f");
	
	//*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*//
	
	TH1D *qm17PbPbRatio = (TH1D*)PbPb_qm17Shape_data->Clone("qm17PbPbRatio");
	qm17PbPbRatio->Divide(ppMC_qm17Shape_data);
	
	TCanvas *cc = new TCanvas("cc","",1000,600);
	cc->Divide(nCentBins+1,2);
	for(int i=0; i<nCentBins+1; i++){
		cc->cd(i+1);
		cc->GetPad(i+1)->SetLogy();
		if(i==0){
			ppStack->SetMaximum(3000);
			ppStack->SetMinimum(1);
			ppStack->Draw("");
			if(shapeNormalize) ppStack->GetHistogram()->GetYaxis()->SetTitle("b-jet pp #rho(r)");
			else ppStack->GetHistogram()->GetYaxis()->SetTitle("b-jet pp #Rho (r)");
			ppStack->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			ppStack->GetHistogram()->GetXaxis()->SetTitle("Jet-Track #DeltaR");
			ppSumHisto->Draw("same");
			ppMCSumHisto->SetMarkerColor(kBlue-5);
			ppMCSumHisto->SetMarkerStyle(25);
			ppMCSumHisto->Draw("Same");
			ppMCInclSumHisto->SetMarkerStyle(25);
			ppMCInclSumHisto->SetMarkerColor(kBlue-3);
			ppMCInclSumHisto->Draw("Same");
			ppMC_qm17Shape->SetMarkerColor(kGreen-4);
			ppMC_qm17Shape->Draw("same");
			ppMC_qm17Shape_data->SetMarkerStyle(20);
			ppMC_qm17Shape_data->SetMarkerColor(kOrange+2);
			ppMC_qm17Shape_data->Draw("same");
		} 
		else{
			//PbPbStack[5-i]->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			PbPbStack[i-1]->SetMaximum(3000);
			//PbPbStack[i-1]->SetMinimum(ppStack->GetMinimum());
			PbPbStack[i-1]->SetMinimum(1);
			PbPbStack[i-1]->Draw("");
			PbPbStack[i-1]->GetHistogram()->GetXaxis()->SetTitle("Jet-Track #DeltaR");
			PbPbStack[i-1]->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			PbPbSumHisto[i-1]->Draw("same");
			titles[i-1]->Draw("Same");
			PbPb_qm17Shape_data->SetMarkerColor(kOrange+2);
			PbPb_qm17Shape_data->SetMarkerStyle(20);
			if(i==1) PbPb_qm17Shape_data->Draw("same");
			if(i==1) PbPbInclSumHisto[i-1]->Draw("same");
		}
		//cc->GetPad(i)->Update();
	}
	cc->cd(4);
	leg1->Draw();
	//cc->cd(6);
	//TH1D *ppDataMCRatio = (TH1D*)ppSumHisto->Clone("ppDataMCRatio");
	//ppDataMCRatio->Divide(ppMCSumHisto);
	//ppDataMCRatio->GetYaxis()->SetTitle("Data / Kurt MC");
	//ppDataMCRatio->Draw();
	
	for(int i=5; i<=6; i++){
		cc->cd(i);
		PbPbppRatios[(i-5)]->SetMaximum(3.5);
		PbPbppRatios[(i-5)]->SetMinimum(0);
		PbPbppRatios[(i-5)]->SetYTitle("b-jet PbPb / b-jet pp");
		PbPbppRatios[(i-5)]->SetXTitle("Jet-Track #DeltaR");
		PbPbppRatios[(i-5)]->Draw();
		qm17PbPbRatio->SetMarkerColor(kOrange+2);
		qm17PbPbRatio->SetMarkerStyle(20);
		if(i==5) qm17PbPbRatio->Draw("Same");
		PbPbSumHisto_offTrk[(i-5)]->SetMarkerColor(kGreen+2);
		if(i==5) PbPbSumHisto_offTrk[(i-5)]->Draw("same");
		
		
	}
	
	//*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*//
	
	TCanvas *cc4 = new TCanvas("cc4","",1000,600);
	cc4->Divide(nCentBins+1,2);
	for(int i=0; i<nCentBins+1; i++){
		cc4->cd(i+1);
		cc4->GetPad(i+1)->SetLogy();
		if(i==0){}
		else{
			//PbPbStack[5-i]->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			PbPbStack[i-1]->SetMaximum(3000);
			//PbPbStack[i-1]->SetMinimum(ppStack->GetMinimum());
			PbPbStack[i-1]->SetMinimum(1);
			PbPbStack[i-1]->Draw("");
			PbPbStack[i-1]->GetHistogram()->GetXaxis()->SetTitle("Jet-Track #DeltaR");
			PbPbStack[i-1]->GetHistogram()->GetYaxis()->SetTitle("b-jet #rho(r)");
			PbPbStack[i-1]->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			PbPbSumHisto[i-1]->Draw("same");
			PbPbMCSumHisto[i-1]->SetMarkerColor(4);
			PbPbMCSumHisto[i-1]->Draw("same");
			titles[i-1]->Draw("Same");
		}
		//cc->GetPad(i)->Update();
	}
	cc4->cd(4);
	leg1->Draw();
	
	for(int i=5; i<=6; i++){
		cc4->cd(i);
		PbPbbInclRatios[(i-5)]->SetMaximum(3.5);
		PbPbbInclRatios[(i-5)]->SetMinimum(0);
		PbPbbInclRatios[(i-5)]->SetYTitle("b-jet PbPb / Incl-Jet PbPb");
		PbPbbInclRatios[(i-5)]->SetXTitle("Jet-Track #DeltaR");
		PbPbbInclRatios[(i-5)]->Draw();
		PbPbMCbInclRatios[(i-5)]->SetMarkerColor(4);
		PbPbMCbInclRatios[(i-5)]->Draw("same");
		
	}
	
	//*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*//
	
	TCanvas *cc100 = new TCanvas("cc100","",1000,600);
	cc100->Divide(nCentBins+1,2);
	for(int i=0; i<nCentBins+1; i++){
		cc100->cd(i+1);
		cc100->GetPad(i+1)->SetLogy();
		if(i==0){}
		else{
			//PbPbStack[5-i]->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			PbPbInclStack[i-1]->SetMaximum(3000);
			//PbPbStack[i-1]->SetMinimum(ppStack->GetMinimum());
			PbPbInclStack[i-1]->SetMinimum(1);
			PbPbInclStack[i-1]->Draw("");
			PbPbInclStack[i-1]->GetHistogram()->GetXaxis()->SetTitle("Jet-Track #DeltaR");
			PbPbInclStack[i-1]->GetHistogram()->GetYaxis()->SetTitle("inclusive jet #rho(r)");
			PbPbInclStack[i-1]->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			PbPbInclSumHisto[i-1]->Draw("same");
			PbPbMCInclSumHisto[i-1]->SetMarkerColor(4);
			PbPbMCInclSumHisto[i-1]->Draw("same");
			titles[i-1]->Draw("Same");
		}
		//cc->GetPad(i)->Update();
	}
	cc100->cd(4);
	leg1->Draw();
	
	for(int i=5; i<=6; i++){
		cc4->cd(i);
		
	}
	
	//*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*//
	
	TCanvas *cc2 = new TCanvas("cc2","",800,600);
	cc2->Divide(2,2);
	cc2->cd(1);
	ppStack->SetMinimum(0.5);
	//ppStack->Draw("");
	//ppSumHisto->Draw("same");
	ppMCSumHisto->SetMarkerStyle(25);
	ppMCSumHisto->SetMarkerColor(kBlue-3);
	//ppMCSumHisto->Draw("same");
	
	cc2->cd(2);
	ppInclStack->SetMaximum(3000);
	ppInclStack->SetMinimum(0.5);
	ppInclStack->Draw("");
	if(shapeNormalize) ppInclStack->GetHistogram()->GetYaxis()->SetTitle("Incl Jet pp #rho(r)");
	else ppInclStack->GetHistogram()->GetYaxis()->SetTitle("Incl Jet pp #Rho (r)");
	ppInclStack->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
	ppInclStack->GetHistogram()->GetXaxis()->SetTitle("Jet-Track #DeltaR");
	ppInclSumHisto->Draw("same");
	ppMCInclSumHisto->SetMarkerStyle(25);
	ppMCInclSumHisto->SetMarkerColor(kBlue-3);
	ppMCInclSumHisto->Draw("Same");
	ppMC_qm17Shape->SetMarkerColor(kGreen-4);
	ppMC_qm17Shape->Draw("same");
	ppMC_qm17Shape_data->SetMarkerStyle(20);
	ppMC_qm17Shape_data->SetMarkerColor(kOrange+2);
	ppMC_qm17Shape_data->Draw("same");
	
	cc2->cd(3);
	leg1->Draw();
	
	cc2->cd(4);
	ppBInclRatio->SetMaximum(3);
	ppBInclRatio->SetMinimum(0);
	ppBInclRatio->GetXaxis()->SetLabelSize(0.06);
	ppBInclRatio->GetYaxis()->SetLabelSize(0.1);
	ppBInclRatio->SetYTitle("b-jet / Inclusive-Jet");
	//ppBInclRatio->Draw("");
	//ppMCbInclRatio->Draw("hist,same");
	TH1D *dataRatio = (TH1D*)ppMC_qm17Shape_data->Clone("dataRatio");
	dataRatio->Divide(ppInclSumHisto);
	dataRatio->SetYTitle("Ratio to Data");
	dataRatio->GetXaxis()->SetLabelSize(0.06);
	dataRatio->GetYaxis()->SetLabelSize(0.06);
	dataRatio->Draw();
	TH1D *qm17MCratio = (TH1D*)ppMC_qm17Shape->Clone("qm17MCratio");
	qm17MCratio->Divide(ppInclSumHisto);
	qm17MCratio->Draw("Same");
	TH1D *stdMCratio = (TH1D*)ppMCInclSumHisto->Clone("stdMCratio");
	stdMCratio->Divide(ppInclSumHisto);
	stdMCratio->Draw("same");
	
	
	
	//*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*//
	
	TCanvas *cc3 = new TCanvas("cc3","",800,600);
	cc3->Divide(2,2);
	cc3->cd(1);
	ppEtaStack->SetMinimum(0.5);
	ppEtaStack->SetMaximum(3000);
	ppEtaStack->Draw("");
	ppEtaStack->GetHistogram()->GetYaxis()->SetTitle("b-jet pp #Delta#eta");
	ppEtaStack->GetHistogram()->GetXaxis()->SetRangeUser(-2,2);
	ppEtaStack->GetHistogram()->GetXaxis()->SetTitle("Jet-Track #DeltaR");
	ppEtaStack->Draw("");
	ppEtaSumHisto->Draw("same");
	
	/*cc3->cd(2);
	ppInclEtaStack->SetMaximum(3000);
	ppInclEtaStack->SetMinimum(0.5);
	ppInclEtaStack->Draw("");
	ppInclEtaStack->GetHistogram()->GetYaxis()->SetTitle("Incl Jet pp #Delta#eta");
	ppInclEtaStack->GetHistogram()->GetXaxis()->SetRangeUser(-2,2);
	ppInclEtaStack->GetHistogram()->GetXaxis()->SetTitle("Jet-Track #DeltaR");
	ppInclEtaSumHisto->Draw("same");
	
	cc3->cd(3);
	leg1->Draw();
	
	cc3->cd(4);
	ppBInclEtaRatio->SetMaximum(3);
	ppBInclEtaRatio->SetMinimum(0);
	ppBInclEtaRatio->GetXaxis()->SetLabelSize(0.06);
	ppBInclEtaRatio->GetYaxis()->SetLabelSize(0.1);
	ppBInclEtaRatio->SetYTitle("b-jet / Inclusive-Jet");
	ppBInclEtaRatio->Draw("");*/
	
	
}